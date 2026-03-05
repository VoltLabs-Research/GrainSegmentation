#include <volt/grain_segmentation_service.h>
#include <volt/core/frame_adapter.h>
#include <volt/core/analysis_result.h>
#include <volt/utilities/concurrence/parallel_system.h>
#include <volt/utilities/json_utils.h>
#include <spdlog/spdlog.h>
#include <map>

namespace Volt{

using namespace Volt::Particles;

GrainSegmentationService::GrainSegmentationService()
    : _rmsd(0.10f),
      _identificationMode(StructureAnalysis::Mode::PTM),
      _adoptOrphanAtoms(true),
      _minGrainAtomCount(100),
      _handleCoherentInterfaces(true),
      _outputBonds(false){}

void GrainSegmentationService::setIdentificationMode(StructureAnalysis::Mode mode){
    _identificationMode = mode;
}

void GrainSegmentationService::setRMSD(float rmsd){
    _rmsd = rmsd;
}

void GrainSegmentationService::setParameters(
    bool adoptOrphanAtoms,
    int minGrainAtomCount,
    bool handleCoherentInterfaces,
    bool outputBonds
){
    _adoptOrphanAtoms = adoptOrphanAtoms;
    _minGrainAtomCount = minGrainAtomCount;
    _handleCoherentInterfaces = handleCoherentInterfaces;
    _outputBonds = outputBonds;
}

json GrainSegmentationService::compute(const LammpsParser::Frame &frame, const std::string &outputFilename){
    if(frame.natoms <= 0)
        return AnalysisResult::failure("Invalid number of atoms");

    auto positions = FrameAdapter::createPositionPropertyShared(frame);
    if(!positions)
        return AnalysisResult::failure("Failed to create position property");

    std::vector<Matrix3> preferredOrientations;
    preferredOrientations.push_back(Matrix3::Identity());

    auto structuretypes = std::make_unique<ParticleProperty>(frame.natoms, DataType::Int, 1, 0, true);
    AnalysisContext context(
        positions.get(),
        frame.simulationCell,
        LATTICE_BCC,
        nullptr,
        structuretypes.get(),
        std::move(preferredOrientations)
    );

    auto structureAnalysis = std::make_unique<StructureAnalysis>(
        context,
        false,
        _identificationMode,
        _rmsd
    );

    {
        PROFILE("Identify Structures");
        structureAnalysis->identifyStructures();
    }

    std::vector<int> extractedStructureTypes;
    extractedStructureTypes.reserve(frame.natoms);
    for(int i = 0; i < frame.natoms; i++){
        extractedStructureTypes.push_back(structureAnalysis->context().structureTypes->getInt(i));
    }

    if(!outputFilename.empty()){
        if(_identificationMode == StructureAnalysis::Mode::PTM){
            spdlog::warn("PTM data export to file not available in standalone package");
        }
        return performGrainSegmentation(frame, *structureAnalysis, extractedStructureTypes, outputFilename);
    }

    return AnalysisResult::failure("No output filename specified");
}

json GrainSegmentationService::performGrainSegmentation(
    const LammpsParser::Frame &frame,
    const StructureAnalysis &structureAnalysis,
    const std::vector<int>& structureTypes,
    const std::string &outputFile
){
    spdlog::info("Starting grain segmentation analysis...");

    try{
        const auto& ctx = structureAnalysis.context();
        if(!ctx.ptmOrientation || !ctx.correspondencesCode){
            spdlog::error("PTM orientation data not available. Grain segmentation requires PTM mode.");
            return AnalysisResult::failure("Grain segmentation requires PTM mode with orientation output enabled.");
        }

        auto positions = std::make_shared<ParticleProperty>(frame.natoms, ParticleProperty::PositionProperty, 0, true);
        for(size_t i = 0; i < frame.positions.size() && i < static_cast<size_t>(frame.natoms); i++){
            positions->setPoint3(i, frame.positions[i]);
        }

        auto structures = std::make_shared<ParticleProperty>(frame.natoms, DataType::Int, 1, 0, false);
        for(size_t i = 0; i < structureTypes.size(); i++){
            structures->setInt(i, structureTypes[i]);
        }

        auto orientations = std::make_shared<ParticleProperty>(frame.natoms, DataType::Double, 4, 0, false);
        for(size_t i = 0; i < static_cast<size_t>(frame.natoms); i++){
            for(int c = 0; c < 4; c++){
                orientations->setDoubleComponent(i, c, ctx.ptmOrientation->getDoubleComponent(i, c));
            }
        }

        auto correspondences = std::make_shared<ParticleProperty>(frame.natoms, DataType::Int64, 1, 0, false);
        {
            auto* src = reinterpret_cast<const uint64_t*>(ctx.correspondencesCode->data());
            auto* dst = reinterpret_cast<uint64_t*>(correspondences->data());
            std::copy(src, src + frame.natoms, dst);
        }

        spdlog::info("Running GrainSegmentationEngine1...");
        auto engine1 = std::make_shared<GrainSegmentationEngine1>(
            positions,
            structures,
            orientations,
            correspondences,
            &frame.simulationCell,
            _handleCoherentInterfaces,
            _outputBonds
        );

        engine1->perform();

        spdlog::info("GrainSegmentationEngine1 complete. Suggested merging threshold: {:.4f}", engine1->suggestedMergingThreshold());
        spdlog::info("Running GrainSegmentationEngine2...");

        GrainSegmentationEngine2 engine2(
            engine1,
            _adoptOrphanAtoms,
            static_cast<size_t>(_minGrainAtomCount),
            true
        );

        engine2.perform();
        spdlog::info("Found {} grains", engine2.grainCount());

        auto atomClusters = engine2.atomClusters();
        std::vector<int> grainIds(frame.natoms, 0);
        for(size_t i = 0; i < static_cast<size_t>(frame.natoms); i++){
            grainIds[i] = atomClusters->getInt(i);
        }

        // Build grain center-of-mass map
        std::map<int, Point3> grainCenters;
        std::map<int, int> grainAtomCount;
        for(int i = 0; i < frame.natoms; i++){
            int gid = grainIds[i];
            if(i < static_cast<int>(frame.positions.size())){
                const auto& p = frame.positions[i];
                grainCenters[gid] = grainCenters[gid] + Vector3(p.x(), p.y(), p.z());
                grainAtomCount[gid]++;
            }
        }

        // Build grains sub_listing
        json grainsArray = json::array();
        for(const auto &grain : engine2.grains()){
            json grainInfo;
            grainInfo["id"] = grain.id;
            grainInfo["size"] = grain.size;
            grainInfo["orientation"] = {
                grain.orientation.x(),
                grain.orientation.y(),
                grain.orientation.z(),
                grain.orientation.w()
            };
            // Center of mass
            if(grainAtomCount.count(grain.id) && grainAtomCount[grain.id] > 0){
                int cnt = grainAtomCount[grain.id];
                const auto& c = grainCenters[grain.id];
                grainInfo["pos"] = { c.x() / cnt, c.y() / cnt, c.z() / cnt };
            } else {
                grainInfo["pos"] = {0.0, 0.0, 0.0};
            }
            grainsArray.push_back(grainInfo);
        }

        json result;
        result["main_listing"] = {
            { "total_grains", static_cast<int>(engine2.grainCount()) },
            { "merging_threshold", engine1->suggestedMergingThreshold() }
        };
        result["sub_listings"] = { { "grains", grainsArray } };

        const std::string msgpackPath = outputFile + "_grains.msgpack";
        if(JsonUtils::writeJsonMsgpackToFile(result, msgpackPath, false)){
            spdlog::info("Exported grain data to: {}", msgpackPath);
        }else{
            spdlog::warn("Could not write grains msgpack: {}", msgpackPath);
        }

        return result;
    }catch(const std::exception& e){
        spdlog::error("Grain segmentation error: {}", e.what());
        return AnalysisResult::failure(std::string("Grain segmentation failed: ") + e.what());
    }
}

}

