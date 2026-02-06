#pragma once

#include <volt/core/volt.h>
#include <volt/core/lammps_parser.h>
#include <nlohmann/json.hpp>
#include <volt/core/particle_property.h>
#include <volt/structures/crystal_structure_types.h>
#include <volt/analysis/structure_analysis.h>
#include <volt/analysis/analysis_context.h>
#include <volt/grain_segmentation_engine.h>
#include <string>

namespace Volt{
using json = nlohmann::json;

class GrainSegmentationService{
public:
    GrainSegmentationService();
    
    void setIdentificationMode(StructureAnalysis::Mode mode);
    void setRMSD(float rmsd);

    void setParameters(
        bool adoptOrphanAtoms,
        int minGrainAtomCount,
        bool handleCoherentInterfaces,
        bool outputBonds
    );

    json compute(
        const LammpsParser::Frame &frame,
        const std::string &outputFilename = ""
    );

private:
    float _rmsd;
    StructureAnalysis::Mode _identificationMode;

    bool _adoptOrphanAtoms;
    int _minGrainAtomCount;
    bool _handleCoherentInterfaces;
    bool _outputBonds;

    json performGrainSegmentation(
        const LammpsParser::Frame &frame,
        const StructureAnalysis& structureAnalysis,
        const std::vector<int>& structureTypes,
        const std::string& outputFile
    );
};

}
