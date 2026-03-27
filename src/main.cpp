#include <volt/cli/common.h>
#include <volt/grain_segmentation_service.h>
#include <oneapi/tbb/global_control.h>
#include <tbb/info.h>

#include <algorithm>
#include <fstream>
#include <set>

using namespace Volt;
using namespace Volt::CLI;

void showUsage(const std::string& name) {
    printUsageHeader(name, "Volt - Grain Segmentation");
    std::cerr
        << "  --rmsd <float>                        RMSD threshold for PTM. [default: 0.1]\n"
        << "  --minGrainAtomCount <int>             Minimum atoms per grain. [default: 100]\n"
        << "  --adoptOrphanAtoms <true|false>       Adopt orphan atoms. [default: true]\n"
        << "  --handleCoherentInterfaces <true|false> Handle coherent interfaces. [default: true]\n"
        << "  --outputBonds                         Output neighbor bonds. [default: false]\n"
        << "  --threads <int>                       Max worker threads (TBB/OMP). [default: auto]\n";
    printHelpOption();
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        showUsage(argv[0]);
        return 1;
    }
    
    std::string filename, outputBase;
    auto opts = parseArgs(argc, argv, filename, outputBase);
    
    if (hasOption(opts, "--help") || filename.empty()) {
        showUsage(argv[0]);
        return filename.empty() ? 1 : 0;
    }
    
    if (!hasOption(opts, "--threads")) {
        const int maxAvailableThreads = static_cast<int>(oneapi::tbb::info::default_concurrency());
        int physicalCores = 0;
        std::ifstream cpuinfo("/proc/cpuinfo");
        if (cpuinfo.is_open()) {
            std::set<std::pair<int, int>> physicalCoreIds;
            int fallbackCpuCores = 0;
            int physicalId = -1;
            int coreId = -1;
            std::string line;
            while (std::getline(cpuinfo, line)) {
                if (line.empty()) {
                    if (physicalId >= 0 && coreId >= 0) {
                        physicalCoreIds.emplace(physicalId, coreId);
                    }
                    physicalId = -1;
                    coreId = -1;
                    continue;
                }
                if (line.rfind("physical id", 0) == 0) {
                    physicalId = std::stoi(line.substr(line.find(':') + 1));
                } else if (line.rfind("core id", 0) == 0) {
                    coreId = std::stoi(line.substr(line.find(':') + 1));
                } else if (line.rfind("cpu cores", 0) == 0) {
                    fallbackCpuCores = std::max(fallbackCpuCores, std::stoi(line.substr(line.find(':') + 1)));
                }
            }
            if (physicalId >= 0 && coreId >= 0) {
                physicalCoreIds.emplace(physicalId, coreId);
            }
            physicalCores = !physicalCoreIds.empty()
                ? static_cast<int>(physicalCoreIds.size())
                : fallbackCpuCores;
        }
        int defaultThreads = maxAvailableThreads;
        if (physicalCores > 0) {
            defaultThreads = std::min(maxAvailableThreads, physicalCores);
        }
        opts["--threads"] = std::to_string(std::max(1, defaultThreads));
    }

    const int requestedThreads = getInt(opts, "--threads");
    oneapi::tbb::global_control parallelControl(
        oneapi::tbb::global_control::max_allowed_parallelism,
        static_cast<std::size_t>(std::max(1, requestedThreads))
    );
    initLogging("grain-segmentation");
    spdlog::info("Using {} threads (OneTBB)", requestedThreads);
    
    LammpsParser::Frame frame;
    if (!parseFrame(filename, frame)) return 1;
    
    outputBase = deriveOutputBase(filename, outputBase);
    spdlog::info("Output base: {}", outputBase);
    
    bool adoptOrphanAtoms = getString(opts, "--adoptOrphanAtoms", "true") == "true";
    int minGrainAtomCount = getInt(opts, "--minGrainAtomCount", 100);
    bool handleCoherentInterfaces = getString(opts, "--handleCoherentInterfaces", "true") == "true";
    bool outputBonds = hasOption(opts, "--outputBonds");
    
    spdlog::info("Grain segmentation parameters:");
    spdlog::info("  - adoptOrphanAtoms: {}", adoptOrphanAtoms);
    spdlog::info("  - minGrainAtomCount: {}", minGrainAtomCount);
    spdlog::info("  - handleCoherentInterfaces: {}", handleCoherentInterfaces);
    spdlog::info("  - outputBonds: {}", outputBonds);
    
    GrainSegmentationService analyzer;
    analyzer.setRMSD(getDouble(opts, "--rmsd", 0.1f));
    analyzer.setParameters(
        adoptOrphanAtoms,
        minGrainAtomCount,
        handleCoherentInterfaces,
        outputBonds
    );
    
    spdlog::info("Starting grain segmentation...");
    json result = analyzer.compute(frame, outputBase);
    
    if (result.value("is_failed", false)) {
        spdlog::error("Analysis failed: {}", result.value("error", "Unknown error"));
        return 1;
    }
    return 0;
}
