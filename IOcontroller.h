#pragma once

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/stat.h>
#endif

class IoController
{
public:
    static void mkdir_absence(const char* outFolder)
    {
#if defined(_WIN32)
        CreateDirectoryA(outFolder, nullptr); // can be used on Windows
#else
        mkdir(outFolder, 0733); // can be used on non-Windows
#endif
    }

    /// Save a serialized file
    template <class T>
    static void SaveFile(const std::string filename, const T& output)
    {
        std::ofstream outfile(filename, std::ios::binary);

        if (!outfile.eof() && !outfile.fail())
        {
            StreamType res;
            serialize(output, res);
            outfile.write(reinterpret_cast<char*>(&res[0]), res.size());
            outfile.close();
            res.clear();
            std::cout << "Save file successfully: " << filename << '\n';
        }
        else
        {
            std::cout << "Save file failed: " + filename << '\n';
            exit(1);
        }
    }

    /// Load a serialized file
    template <class T>
    static void LoadFile(const std::string filename, T& input)
    {
        std::ifstream infile(filename, std::ios::binary);

        if (!infile.eof() && !infile.fail())
        {
            infile.seekg(0, std::ios_base::end);
            const std::streampos fileSize = infile.tellg();
            infile.seekg(0, std::ios_base::beg);
            std::vector<uint8_t> res(fileSize);
            infile.read(reinterpret_cast<char*>(&res[0]), fileSize);
            infile.close();
            input.clear();
            auto it = res.cbegin();
            input = deserialize<T>(it, res.cend());
            res.clear();
        }
        else
        {
            std::cout << "Cannot open file: " + filename << '\n';
            exit(1);
        }
    }

    /// Save graph structure to a file
    static void SaveGraphStruct(const std::string graphName, const Graph& vecGraph, const bool isReverse)
    {
        std::string postfix = ".vec.graph";

        if (isReverse) postfix = ".vec.rvs.graph";

        const std::string filename = graphName + postfix;
        SaveFile(filename, vecGraph);
    }

    /// Load graph structure from a file
    static void LoadGraphStruct(const std::string graphName, Graph& vecGraph, const bool isReverse)
    {
        std::string postfix = ".vec.graph";

        if (isReverse) postfix = ".vec.rvs.graph";

        const std::string filename = graphName + postfix;
        LoadFile(filename, vecGraph);
    }

    static void SaveGraphProbDist(const std::string graphName, int dist)
    {
        std::ofstream outFile(graphName + ".probdist");
        outFile<< dist;
    }


    static int LoadGraphProbDist(const std::string graphName)
    {
        std::string filename = graphName + ".probdist";
        std::ifstream infile(filename);
        int probDist = WEIGHTS;

        if (!infile.is_open())
        {
            std::cout << "The file \"" + filename + "\" can NOT be opened\n";
            return probDist;
        }

        infile >> probDist;
        infile.close();
        std::cout << "probability distribution: " << probDist << std::endl;
        return probDist;
    }

    /// New member for RFIM
    /// Save the community structure to a file
    static void SaveCommStruct(const std::string graphName, const NodeComms& vecComms,
        const CommLists& vecCommLists, const CommFairList& commFairList)
    {
        SaveFile(graphName + ".veccomm", vecComms);
        SaveFile(graphName + ".commvec", vecCommLists);
        SaveFile(graphName + ".commfair", commFairList);
    }

    /// Load the community structure from a file
    static void LoadCommStruct(const std::string& graphName,
        NodeComms& vecComms,
        CommLists& vecCommLists,
        CommFairList& fairList)
    {
        LoadFile(graphName + ".veccomm", vecComms);
        LoadFile(graphName + ".commvec", vecCommLists);
        LoadFile(graphName + ".commfair", fairList);
    }

    /// Get out-file name
    // static std::string BuildOutFileName(const std::string graphName, const std::string algName, const int seedsize,
    //                                     const std::string probDist, const float probEdge, const double saturateRatio = -1.0)
    // {
    //     std::string outFileName;

    //     if (probDist == "uniform")
    //     {
    //         outFileName = graphName + "_" + algName + "_k" + std::to_string(seedsize) + "_" + probDist + "_" + std::to_string(probEdge);
    //     }

    //     outFileName = graphName + "_" + algName + "_k" + std::to_string(seedsize) + "_" + probDist;

    //     if (saturateRatio != -1.0) outFileName += "_sr" + std::to_string(saturateRatio);

    //     return outFileName;
    // }

    /// Get out-file name
    static std::string BuildOutFileName(const std::string& graphName, const std::string& algName, const int seedsize,
        const std::string& probDist, const float probEdge, const double saturateRatio = -1.0, const double eps = -1.0,
        const double delta = -1.0)
    {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(3);

        std::string outFileName = graphName + "_" + algName + "_k" + std::to_string(seedsize) + "_" + probDist;

        if (probDist == "uniform")
        {
            oss.str(""); oss.clear();
            oss << probEdge;
            outFileName += oss.str();
        }

        if (saturateRatio != -1.0)
        {
            oss.str(""); oss.clear();
            oss << saturateRatio;
            outFileName += "_sr" + oss.str();
        }

        if (eps != -1.0)
        {
            oss.str(""); oss.clear();
            oss << eps;
            outFileName += "_eps" + oss.str();
        }

        if (delta != -1.0)
        {
            oss.str(""); oss.clear();
            oss << delta;
            outFileName += "_dlt" + oss.str();
        }

        return outFileName;
    }

    /// Print the results
    static void WriteResult(const std::string& outFileName, const TResult& resultObj, const std::string& outFolder, const FuncType funcType = IM)
    {
        if (funcType == IM)
        {
            const auto approx = resultObj.get_approximation();
            const auto runTime = resultObj.get_running_time();
            const auto influence = resultObj.get_influence();
            const auto influenceOriginal = resultObj.get_influence_original();
            const auto seedSize = resultObj.get_seed_size();
            const auto RRsetsSize = resultObj.get_RRsets_size();
            std::cout << "   --------------------" << std::endl;
            std::cout << "  |Approx.: " << approx << std::endl;
            std::cout << "  |Time (sec): " << runTime << std::endl;
            std::cout << "  |Influence: " << influence << std::endl;
            std::cout << "  |Self-estimated influence: " << influenceOriginal << std::endl;
            std::cout << "  |#Seeds: " << seedSize << std::endl;
            std::cout << "  |#RR sets: " << RRsetsSize << std::endl;
            std::cout << "   --------------------" << std::endl;
            mkdir_absence(outFolder.c_str());
            std::ofstream outFileNew(outFolder + "/" + outFileName);

            if (outFileNew.is_open())
            {
                outFileNew << "Approx.: " << approx << std::endl;
                outFileNew << "Time (sec): " << runTime << std::endl;
                outFileNew << "Influence: " << influence << std::endl;
                outFileNew << "Self-estimated influence: " << influenceOriginal << std::endl;
                outFileNew << "#Seeds: " << seedSize << std::endl;
                outFileNew << "#RR sets: " << RRsetsSize << std::endl;
                outFileNew.close();
            }
        }

        else if (funcType == RFIM)
        {
            const auto approx = resultObj.get_approximation();
            const auto runTime = resultObj.get_running_time();
            const auto influence = resultObj.get_influence();
            const auto influenceOriginal = resultObj.get_influence_original();
            const auto seedSize = resultObj.get_seed_size();
            const auto RRsetsSize = resultObj.get_RRsets_size();
            const auto fairInfluence = resultObj.get_fair_inf();
            const auto fairInfluenceOriginal = resultObj.get_fair_inf_ori();
            const auto rScore = resultObj.get_robustness_score();
            const auto rScoreOri = resultObj.get_robustness_score_ori();
            const auto objMC = resultObj.get_resMC();

            std::cout << "   --------------------" << std::endl;
            std::cout << "  |Approx.: " << approx << std::endl;
            std::cout << "  |Time (sec): " << runTime << std::endl;
            // std::cout << "  |Influence: " << influence << std::endl;
            // std::cout << "  |Self-estimated influence: " << influenceOriginal << std::endl;
            std::cout << "  |Objective rho: " << fairInfluence << std::endl;
            // std::cout << "  |Fair influence original: " << fairInfluenceOriginal << std::endl;
            std::cout << "  |Robustness score (Saturate Greedy): " << rScore << std::endl;
            std::cout << "  |Self-estimated Robustness score (Saturate Greedy): " << rScoreOri << std::endl;
            std::cout << "  |Objective rho (MC simulation): " << objMC << std::endl;
            std::cout << "  |#Seeds: " << seedSize << std::endl;
            std::cout << "  |#RR sets: " << RRsetsSize << std::endl;
            std::cout << "   --------------------" << std::endl;
            mkdir_absence(outFolder.c_str());
            std::ofstream outFileNew(outFolder + "/" + outFileName);

            if (outFileNew.is_open())
            {
                outFileNew << "Approx.: " << approx << std::endl;
                outFileNew << "Time (sec): " << runTime << std::endl;
                outFileNew << "Influence: " << influence << std::endl;
                outFileNew << "Self-estimated influence: " << influenceOriginal << std::endl;
                // outFileNew << "Self-estimated fair influence: " << fairInfluenceOriginal << std::endl;
                outFileNew << "Robustness score (Saturate Greedy): " << rScore << std::endl;
                outFileNew << "Self-estimated robustness score (Saturate Greedy): " << rScoreOri << std::endl;
                outFileNew << "Objective rho: " << fairInfluence << std::endl;
                outFileNew << "Objective rho (MC simulation): " << objMC << std::endl;
                outFileNew << "#Seeds: " << seedSize << std::endl;
                outFileNew << "#RR sets: " << RRsetsSize << std::endl;
                outFileNew.close();
            }
        }
    }

    /// Print the seeds
    static void WriteOrderSeeds(const std::string& outFileName, const TResult& resultObj, const std::string& outFolder)
    {
        auto vecSeed = resultObj.get_seed_vec();
        mkdir_absence(outFolder.c_str());
        const auto outpath = outFolder + "/seed";
        mkdir_absence(outpath.c_str());
        std::ofstream outFile(outpath + "/seed_" + outFileName);

        for (auto i = 0; i < vecSeed.size(); i++)
        {
            outFile << vecSeed[i] << '\n';
        }

        outFile.close();
    }

    /// Return a filename that does not yet exist by appending _2, _3, etc.
    static std::string GetNonExistFileName(const std::string& basePath, const std::string& filename)
    {
        std::string fullPath = basePath + "/" + filename;

        std::ifstream test(fullPath);
        if (!test.fail()) {
            test.close();

            int count = 2;
            while (true) {
                std::string newName = filename + "_" + std::to_string(count);
                std::string newFullPath = basePath + "/" + newName;

                std::ifstream testNew(newFullPath);
                if (testNew.fail()) {
                    return newName;
                }
                ++count;
            }
        }

        return filename;
    }
};

using TIO = IoController;
using PIO = std::shared_ptr<IoController>;
