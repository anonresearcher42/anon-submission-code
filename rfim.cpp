#include "stdafx.h"
#include "SFMT/dSFMT/dSFMT.c"
#include "alg.cpp"

void init_random_seed()
{
    // Randomize the seed for generating random numbers
    dsfmt_gv_init_gen_rand(static_cast<uint32_t>(time(nullptr)));
}

int main(int argc, char* argv[])
{

    TArgument Arg(argc, argv);

    if (Arg._probDist == PROB_DIST_ERROR)
    {
        LogInfo("The input probability distribution is not supported:", Arg._probDistStr);
        LogInfo("The supported probability distribution: weights, wc, uniform, skewed");
        return 0;
    }

    if (Arg._func == FUNC_ERROR)
    {
        LogInfo("The input func is not supported: ", Arg._funcStr);
        LogInfo("The supported func: format, im, fim, rfim");
        return 0;
    }

    init_random_seed();

    const std::string infilename = Arg._dir + "/" + Arg._graphname;
    const std::string commFilenames = Arg._dir + "/" + Arg._graphname + ".comm";
    if (Arg._func == FORMAT)
    {
        // Format the graph and the community
        GraphBase::FormatGraph(infilename, Arg._probDist, Arg._wcVar, Arg._probEdge, Arg._skewType, commFilenames);
        return 0;
    }

    bool verbose = true;
    if (!verbose) {
        std::cout.setstate(std::ios_base::failbit);
    }

    Timer mainTimer("main");
    // Load the reverse graph
    Graph graph;
    std::cout << "Loading graph: " << infilename << std::endl;
    GraphBase::LoadGraph(graph, infilename);
    int probDist = GraphBase::LoadGraphProbDist(infilename);

    // Load the original graph
    Graph graphOri;
    std::cout << "Loading original graph: " << infilename << std::endl;
    GraphBase::LoadGraph(graphOri, infilename, false);

    // Initialize a result object to record the results
    TResult tRes;

    // pAlg: pointer to the algorithm instance (TAlg)
    PAlg pAlg;

    NodeComms vecComm;
    CommLists commVec;
    CommFairList commFair;

    // Load the community structures
    if (Arg._func == FIM || Arg._func == RFIM)
    {
        std::cout << "Loading community: " << commFilenames << std::endl;
        GraphBase::LoadComm(vecComm, commVec, commFair, infilename);
        pAlg = std::make_unique<TAlg>(graph, vecComm, commVec, commFair, tRes, graphOri);
    }
    else // Arg._func == IM
    {
        pAlg = std::make_unique<TAlg>(graph, tRes);
    }
    pAlg->set_vanilla_sample(Arg._vanilla); // true
    pAlg->set_prob_dist((ProbDist)probDist); // Set propagation model
    auto delta = Arg._delta;
    if (delta < 0) delta = 1.0 / graph.size();

    FairType _fairType = FAIR_ERROR;
    if (Arg._isLinearFair) _fairType = LINEAR;
    else _fairType = SCALE_POWER;
    pAlg->set_fair_type(_fairType);

    int seedSize = Arg._seedsize;
    std::cout << "seedSize k=" << seedSize << std::endl;
    Arg.build_outfilename(seedSize, (ProbDist)probDist, graph);
    std::cout << "---The Begin of " << Arg._outFileName << "---\n";

    switch (Arg._func)
    {
    case IM:
        pAlg->subsimWithHIST(seedSize, Arg._eps, delta);
        break;
    case FIM:
        std::cout << " Not support for FIM" << std::endl;
        break;
    case RFIM:
        switch (Arg._methodType)
        {
            case SATURATE_GREEDY_HIST:
                pAlg->saturateGreedyWithHIST(seedSize, Arg._eps, delta, Arg._gamma, Arg._saturateRatio);
                break;
            case SINGLE_GREEDY:
                pAlg->singleGreedyWithHIST(seedSize, Arg._eps, delta, Arg._gamma, Arg._saturateRatio);
                break;
            case ALL_GREEDY:
                pAlg->allGreedyWithHIST(seedSize, Arg._eps, delta, Arg._gamma, Arg._saturateRatio);
                break;
            case SATURATE_GREEDY_CELF:
                pAlg->saturateGreedyWithCELF(seedSize, Arg._eps, delta, Arg._gamma, Arg._saturateRatio);
                break;
            default:
                std::cerr<< "Error: Unknown method type" << std::endl;
                break;
        }
        break;
    default:
        std::cerr<< "Error: Unknown function type" << std::endl;
        break;
    }

    Arg._outFileName = TIO::GetNonExistFileName(Arg._resultFolder, Arg._outFileName);
    TIO::WriteResult(Arg._outFileName, tRes, Arg._resultFolder, Arg._func);
    TIO::WriteOrderSeeds(Arg._outFileName, tRes, Arg._resultFolder);
    std::cout << "---The End of " << Arg._outFileName << "---\n";
    pAlg->RefreshHypergraph();
    return 0;
}