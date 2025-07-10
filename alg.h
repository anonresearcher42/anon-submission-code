#pragma once

class Alg
{
private:
    /// _numV: number of nodes in the graph.
    uint32_t _numV;
    /// _numE: number of edges in the graph.
    size_t _numE;

    /// _numMetric: number of metrics.
    uint32_t _numMetric = 1;
    /// _numComm: number of communities under each metric.
    std::vector<uint32_t> _numComm;

    /// _numRRsets: number of RR/G-RR sets.
    size_t _numRRsets = 0;
    /// Upper bound in the last round for __mode=1.
    double _boundLast = DBL_MAX;
    /// The minimum upper bound among all rounds for __model=2.
    double _boundMin = DBL_MAX;

    // number of MC simulations
    int _numMC = 10000;

    // original seed budget k
    int _k = 0;
    // saturation ratio
    double _beta = 0.0;
    // sigma_min
    double _fInfMin = 0.0;
    // sigma_max
    double _fInfMax = 0.0;
    // H_min
    double _rFairInfMin = 0.0;
    // _alpha_min
    double _alphaMin = DBL_MAX;

    /// Two hyper-graphs, one is used for selecting seeds and the other is used for validating influence.
    THyperGraph _hyperGraph, _hyperGraphVldt;
    THyperGraph _hyperGraphOri;
    /// Result object.
    TResult &_res;
    /// Seed set.
    Nodelist _vecSeed;
    ProbDist _probDist = WC;
    /// Best seed set S^{*} in saturateGreedyWithHIST
    // Nodelist _vecBestSeed;

    /// The greedy (instead of optimal) fair influence \sigma(S^{g}_{q}) for each metric for each size of seed |S|.
    /// _vecOptFairInf[q][k]: (for FIM) the greedy fair influence for metric q and seed size k,
    /// which is used for Saturate Greedy in RFIM.
    std::vector<double> _vecOptFairInf;
    /// The current fair influence under the constructed G-RR sets and selected seeds _vecSeed.
    double _fairInf = 0.0;

    // the parameter c in Saturate Greedy
    double _c = 0.5;

    double _baseNumRRsets = 0.0;

    std::vector<uint32_t> _vecOutDegree;
    std::vector<uint32_t> _vecVldtInf;

    double _lowerDeg = 0.0;

    double _threshold = DBL_MAX;

    /// Fair influence for each metric, for calculating the roughly lower / upper bound of H
    std::vector<double> _vecFairInf;
    /// Fair influence for each metric, for calculating the upper bound of H
    std::vector<double> _vecFairInfVldt;
    /// Fair influence for differen size of seeds for each metric, for the roughly lower bound of H
    std::vector<std::vector<double>> _vecFairInfDiffSize;
    /// Tighter upper bound of Fair influence for each metric
    std::vector<double> _vecBoundMinFairInf;

    /// The general coverage of the G-RR sets for different size of seeds
    std::vector<double> _vecSelfVldtGFCoverage;

    /// Maximum coverage by lazy updating.
    double MaxCoverVanilla(const int targetSize);
    double MaxCoverOutDegPrority(const int targetSize);
    double MaxCoverIMSentinel(std::vector<uint32_t> &seedSet, const int targetSize);
    double MaxCoverSentinelSet(const int targetSize, const int totalTargetSize);

    /// Maximum coverage by maintaining the top-k marginal coverage.
    double MaxCoverTopK(const int targetSize);
    /// Maximum coverage.
    double MaxCover(const int targetSize);

    /// Maximum Coverage Greedy selection for the top-k marginal (robust) fair influence.
    double MaxCoverGreedyMarginalFairSentinelSet(const int targetSize, const int totalTargetSize);
    double MaxCoverGreedyMarginalFairIMSentinel(std::vector<uint32_t> &seedSet, const int targetSize);
    /// Calculate the unbiased estimate fair influence under the constructed G-RR sets + node u.
    double ComputeFairMarginalGain(uint32_t u, bool usingVldtGraph = false);
    /// Update the G-RR sets coverage status (including _vecFairInf) after adding node u.
    void UpdateCoverageStatus(uint32_t u, std::vector<bool>& needRecompute);

    /// Adjust the general coverage bound.
    double AdjustGeneralBound(const double bound, const size_t numR);
    /// Adjust the general influence bound.
    double AdjustGeneralInfluenceBound(const double gInfBound, const size_t numR);


public:
    Alg(const Graph &graph, TResult &tRes): _hyperGraph(graph), _hyperGraphVldt(graph), _res(tRes), _hyperGraphOri(graph)
    {
        _numV = _hyperGraph.get_nodes();
        _numE = _hyperGraph.get_edges();
        _vecOutDegree = std::vector<uint32_t>(_numV);

        for (auto &nbrs : graph)
        {
            for (auto &node : nbrs)
            {
                _vecOutDegree[node.first]++;
            }
        }
    }

    Alg(const Graph &graph, const NodeComms &vecComm, const CommLists &commVec,
        const CommFairList &commFair, TResult &tRes, const Graph &graphOri, const FairType fairType = LINEAR)
        : _hyperGraph(graph, vecComm, commVec, commFair, fairType),
        _hyperGraphVldt(graph, vecComm, commVec, commFair, fairType), _res(tRes),
        _hyperGraphOri(graphOri, vecComm, commVec, commFair, fairType)
    {
        _numMetric = vecComm.size();

        _numV = _hyperGraph.get_nodes();
        _numE = _hyperGraph.get_edges();
        _vecOutDegree = std::vector<uint32_t>(_numV);

        for (auto &nbrs : graph)
        {
            for (auto &node : nbrs)
            {
                _vecOutDegree[node.first]++;
            }
        }

        _vecOptFairInf.resize(_numMetric);
        _numComm.resize(_numMetric);
        for (auto q = 0; q < _numMetric; ++q)
        {
            _numComm[q] = commVec[q].size();
        }

        _vecFairInf.resize(_numMetric);
        _vecFairInfVldt.resize(_numMetric);
        _vecFairInfDiffSize.resize(_numMetric);
        _vecBoundMinFairInf.resize(_numMetric);

        _hyperGraph.setPtr(&_c, &_fairInf, &_vecOptFairInf, &_vecFairInf, &_vecFairInfVldt, &_vecFairInfDiffSize);
        _hyperGraphVldt.setPtr(&_c, &_fairInf, &_vecOptFairInf, &_vecFairInf, &_vecFairInfVldt, &_vecFairInfDiffSize);
        _hyperGraphOri.setPtr(&_c, &_fairInf, &_vecOptFairInf, &_vecFairInf, &_vecFairInfVldt, &_vecFairInfDiffSize);

        _alphaMin = _hyperGraph.get_alpha_min();
    }

    ~Alg()
    {
    }

    /// Set cascade model.
    void set_prob_dist(const ProbDist weight);
    void set_vanilla_sample(const bool isVanilla);

    void RefreshHypergraph()
    {
        _hyperGraph.RefreshHypergraph();
        _hyperGraphVldt.RefreshHypergraph();
    }

    void set_fair_type(const FairType fairType)
    {
        _hyperGraph.set_fair_type(fairType);
        _hyperGraphVldt.set_fair_type(fairType);
        _hyperGraphOri.set_fair_type(fairType);
    }

    /// Evaluate influence spread for the seed set constructed
    double EfficInfVldtAlg();
    /// Evaluate influence spread for a given seed set
    double EfficInfVldtAlg(const Nodelist vecSeed);

    /// Evaluate the general fair influence of the seed set constructed
    double EfficFairInfVldtAlg();
    /// Evaluate the general fair influence of a given seed set
    double EfficFairInfVldtAlg(const Nodelist vecSeed);

    double estimateRRSize();

    double subsimOnly(const int targetSize, const double epsilon, const double delta);
    double IncreaseR2(std::unordered_set<uint32_t> &connSet, double a, double upperOPT, double targetAppr);
    // double IncreaseR2(std::unordered_set<uint32_t> &connSet, double a, double upperOPT, double targetAppr, const int totalTargetSize);

    double FindRemSet(const int targetSize, const double epsilon, const double targeEpsilon, const double delta);
    double FindDynamSub(const int totalTargetSize, const double epsilon, const double delta);

    double subsimWithHIST(const int targetSize, const double epsilon, const double delta);

    // return the upper bound of H
    double calculateUpperBound(const double a1);
    // return the lower bound of H
    double calculateLowerBound(const double a2);
    // return the estimated lower bound (by R1) of H
    double calculateLowerBoundR1(const double a2, const int k);

    // member functions for RFIM

    // Sentinal set selection phase for G-HIST
    double FindGeneralDynamSub(const int totalTargetSize, const double epsilon, const double delta);
    // Remaining set selection phase for G-HIST
    double FindGeneralRemSet(const int targetSize, const double epsilon, const double targetEpsilon, const double delta);
    // SG-HIST
    double SaturateGeneralizedHIST(const int targetSize, const double epsilon, const double delta);
    // S-HIST
    double saturateGreedyWithHIST(const int targetSize = 100, const double epsilon = 0.1, const double delta = 0.05, const double gamma = 0.1, const double saturateRatio = -1.0);
    // baseline: single greedy
    double singleGreedyWithHIST(const int targetSize = 100, const double epsilon = 0.1, const double delta = 0.05, const double gamma = 0.1, const double saturateRatio = -1.0);
    // baseline: all greedy
    double allGreedyWithHIST(const int targetSize = 100, const double epsilon = 0.1, const double delta = 0.05, const double gamma = 0.1, const double saturateRatio = -1.0);

    // baseline: saturate greedy with CELF, this baseline is not included in our work, since they are very inefficient
    double generalCELF(const int targetSize, const double epsilon, const double delta);
    double generalCELFMintss(const double targetThreshold, const double epsilon, const double delta);
    double saturateGreedyWithCELF(const int targetSize = 100, const double epsilon = 0.1, const double delta = 0.05, const double gamma = 0.1, const double saturateRatio = -1.0);
};

using TAlg = Alg;
using PAlg = std::unique_ptr<TAlg>;
