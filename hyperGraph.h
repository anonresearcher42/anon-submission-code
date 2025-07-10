#pragma once

class HyperGraph
{
private:
    /// _numV: number of nodes in the graph.
    uint32_t _numV;
    /// _numE: number of edges in the graph.
    size_t _numE = 0;
    /// _numRRsets: number of RR sets.
    size_t _numRRsets = 0;
    std::vector<bool> _vecVisitBool;
    std::vector<bool> _vecVisitBoolTmp;
    Nodelist _vecVisitNode;

    size_t _hit = 0;
    double _numSamplesEval = 0;
    double _hyperedgeAvgEval = 0.0;

    bool _isVanilla = false;

    /// The type of function when sampling RR sets
    FuncType _funcType = FUNC_ERROR;
    /// The id of metric when sampling RR sets, only works for FIM
    int _metricID = -1;
    /// The number of metrics, only works for RFIM
    int _numMetric = 0;

    /// The method of calauting fairness
    FairType _fairType = LINEAR;

    /// if the MaxCover uses Saturate Greedy, only works for _funcType=RFIM 
    bool _isSaturateGreedy = true;

    /// The greedy (instead of optimal) fair influence \sigma(S^{g}_{q}) for each metric for each size of seed |S|.
    /// (*_pvecOptFairInf)[q][k]: (for FIM) the greedy fair influence for metric q and seed size k,
    /// which is used for Saturate Greedy in RFIM.
    std::vector<double>* _pVecOptFairInf = nullptr;

    /// the parameter c
    const double* _p_c = nullptr;
    /// the current fair influence
    const double *_p_fairInf = nullptr;

    /// Generalized version of hit for FIM/RFIM.
    std::vector<std::vector<size_t>> _gHit;

    GraphEdgeStatus _edgeStatus;

    /// Initialization
    void InitHypergraph()
    {
        _numV = (uint32_t)_graph.size();

        for (auto& nbrs : _graph) _numE += nbrs.size();

        if (_pcommVec) {
            _commSize.clear();
            _commSize.resize(_pcommVec->size());
            for (size_t q = 0; q < _pcommVec->size(); ++q) {
                _commSize[q].resize((*_pcommVec)[q].size());
                for (size_t j = 0; j < (*_pcommVec)[q].size(); ++j) {
                    _commSize[q][j] = (*_pcommVec)[q][j].size();
                }
            }
        }

        _FRsets = FRsets(_numV);
        _vecVisitBool = std::vector<bool>(_numV);
        _vecVisitNode = Nodelist(_numV);

        _vecVisitBoolTmp = std::vector<bool>(_numV);

        _GFRsets.resize(_numMetric);
        _GRRsets.resize(_numMetric);
        _gHit.resize(_numMetric);
        _numActComm.resize(_numMetric);

        for (auto q = 0; q < _numMetric; ++q)
        {
            _GFRsets[q].resize(_commSize[q].size());
            _GRRsets[q].resize(_commSize[q].size());
            _gHit[q].resize(_commSize[q].size());
            _numActComm[q].resize(_commSize[q].size());
            for (auto j = 0; j < _commSize[q].size(); ++j)
            {
                _GFRsets[q][j].resize(_numV, FRset());
            }
        }

        _edgeStatus.resize(_numV);
        for (uint32_t u = 0; u < _numV; ++u)
		{
			_edgeStatus[u].resize(_graph[u].size(), UNCHECK);
		}
    }

    void refresh_edgeStatus(bool usingVecVisitBool = false)
    {
        if (usingVecVisitBool)
        {
            for (uint32_t u = 0; u < _graph.size(); ++u)
            {
                if (!_vecVisitBool[u]) continue;
                std::fill(_edgeStatus[u].begin(), _edgeStatus[u].end(), UNCHECK);
            }
        }
        else
        {
            for (uint32_t u = 0; u < _graph.size(); ++u)
            {
                std::fill(_edgeStatus[u].begin(), _edgeStatus[u].end(), UNCHECK);
            }
        }
    }

public:
    /// _graph: reverse graph
    const Graph& _graph;
    // Graph _graph;
    /// _FRsets: forward cover sets, _FRsets[i] is the node sets that node i can reach
    FRsets _FRsets;
    /// _RRsets: reverse cover sets, _RRsets[i] is the node set that can reach node i
    RRsets _RRsets;

    /// _GFRsets: generalized forward cover sets, _GFRsets[q][j][i] is indexes of G-hyperedges (G-RR sets) that node i can reach in q-th metric, j-th community
    std::vector<std::vector<FRsets>> _GFRsets;
    /// _GRRsets: generalized reverse cover sets, _GRRsets[q][j][i] is the rr-set that can reach node i in q-th metric, j-th community
    std::vector<std::vector<RRsets>> _GRRsets;

    /// _GIsCovered: generalized coverage status, _GIsCovered[q][j][i] is the coverage status of hyperedge i in q-th metric, j-th community
    std::vector<std::vector<std::vector<bool>>> _GIsCovered;
    /// _numGCovered: the number of covered edges under current seeds, _numGCovered[q][j] is the number in q-th metric, j-th community
    std::vector<std::vector<size_t>> _numGCovered;
    /// _GCoverage: The number of marginal coverage of the hyperedge, _GCoverage[q][j][i] is the number of marginal coverage of node i in q-th metric, j-th community
    std::vector<std::vector<std::vector<size_t>>> _GCoverage;

    /// _vecComm && _commVec: the community structures
    // (*_pvecComm)[q][j], the community ID of q-th metric, j-th node
    const NodeComms* _pvecComm = nullptr;
    // (*_pcommVec)[q][j], the set of node id's in q-th metric, j-th community
    const CommLists* _pcommVec = nullptr;

    // (*_pcommFair)[q][j], the fair parameter of q-th metric, j-th community
    const CommFairList* _pcommFair = nullptr;

    // (*_pvecFairInf)[q]: the fair influence of q-th metric
    std::vector<double>* _pvecFairInf = nullptr;
    // (*_pvecFairInfVldt)[q]: the fair influence of q-th metric (vdlt)
    std::vector<double>* _pvecFairInfVldt=nullptr;
    // (*_pvecFairInfDiffSize)[q]: the fair influence of q-th metric (diff size)
    std::vector<std::vector<double>>* _pvecFairInfDiffSize=nullptr;

    // _commSize: the size of each community under each metric
    std::vector<std::vector<uint32_t>> _commSize;

    // _uStart: the start node of the G-RR set
    std::vector<std::vector<uint32_t>> _uStart;

    // the number of activated nodes in each community for each metric.
    std::vector<std::vector<uint32_t>> _numActComm;

    ProbDist _probDist = WEIGHTS;

    explicit HyperGraph(const Graph& graph)
        : _graph(graph), _pvecComm(nullptr), _pcommVec(nullptr), _pcommFair(nullptr)
    {
        InitHypergraph();
        set_func_type(IM);
    }

    explicit HyperGraph(const Graph& graph, const NodeComms& vecComm,
        const CommLists& commVec, const CommFairList& commFair, const FairType fairType = LINEAR)
        : _graph(graph), _pvecComm(&vecComm), _pcommVec(&commVec), _pcommFair(&commFair)
    {
        _numMetric = vecComm.size();
        InitHypergraph();
        set_func_type(RFIM);
        set_fair_type(fairType);
    }

    /// Set cascade model
    void set_prob_dist(const ProbDist dist)
    {
        _probDist = dist;
    }

    void set_vanilla_sample(const bool isVanilla)
    {
        _isVanilla = isVanilla;
    }

    /// Returns the number of nodes in the graph.
    uint32_t get_nodes() const
    {
        return _numV;
    }

    /// Returns the number of edges in the graph.
    size_t get_edges() const
    {
        return _numE;
    }

    /// Returns the number of RR sets (or G-RR sets) in the graph.
    size_t get_RR_sets_size() const
    {
        return _numRRsets;
    }

    /// Set the function type (IM, FIM or RFIM), metric id (only works for FIM), fair type(LINEAR or SCALE_POWER)
    void set_func_type(const FuncType funcType, const int metricID = -1, const FairType fairType = FAIR_ERROR)
    {
        _funcType = funcType;
        if (metricID != -1 && funcType == FIM) _metricID = metricID;
        if ((funcType == FIM || funcType == RFIM) && fairType != FAIR_ERROR) set_fair_type(fairType);
        return ;
    }

    FuncType get_func_type()
    {
        return _funcType;
    }

    int get_metricID()
    {
        return _metricID;
    }

    /// Set the fair type
    void set_fair_type(const FairType fairType)
    {
        _fairType = fairType;
        return ;
    }

    void set_opt_fair_inf_ptr(std::vector<double>* p)
    {
        _pVecOptFairInf = p;
    }

    void set_c_for_RFIM(const double* c)
    {
        _p_c = c;
    }

    void set_currFairInf_ptr(const double * p)
    {
        _p_fairInf = p;
    }

    double get_alpha_min()
    {
        double alphaMin = DBL_MAX;
        for (auto q = 0; q < _numMetric; ++q)
        {
            for (auto j = 0; j < _commSize[q].size(); ++j)
            {
                alphaMin = (alphaMin > (*_pcommFair)[q][j]) ? (*_pcommFair)[q][j] : alphaMin;
            }
        }
        return alphaMin;
    }

    void setPtr(const double* c, const double * currFairInf, std::vector<double>* vecOptFairInf,
        std::vector<double>* vecFairInf, std::vector<double>* vecFairInfVldt, std::vector<std::vector<double>>* vecFairInfDiffSize)
    {
        set_c_for_RFIM(c);
        set_currFairInf_ptr(currFairInf);
        set_opt_fair_inf_ptr(vecOptFairInf);
        _pvecFairInf = vecFairInf;
        _pvecFairInfVldt = vecFairInfVldt;
        _pvecFairInfDiffSize = vecFairInfDiffSize;
    }

    void set_isSaturateGreedy(const bool isSaturateGreedy)
    { // only works for RFIM
        _isSaturateGreedy = isSaturateGreedy;
    }

    bool get_isSaturateGreedy()
    {
        return _isSaturateGreedy;
    }

    /// Initialize the G-coverage status
    void InitGcoverage()
    {
        _GIsCovered.resize(_numMetric);
        _numGCovered.resize(_numMetric);
        _GCoverage.resize(_numMetric);
        for (auto q = 0; q < _numMetric; ++q)
        {
            _GIsCovered[q].resize(_commSize[q].size());
            _numGCovered[q].resize(_commSize[q].size());
            _GCoverage[q].resize(_commSize[q].size());
            std::fill(_numGCovered[q].begin(), _numGCovered[q].end(), 0);
            for (auto j = 0; j < _commSize[q].size(); ++j)
            {
                _GIsCovered[q][j].resize(_numRRsets);
                std::fill(_GIsCovered[q][j].begin(), _GIsCovered[q][j].end(), false);
                _GCoverage[q][j].resize(_numV);
                for (auto i = 0; i < _numV; i++)
                {
                    _GCoverage[q][j][i] = _GFRsets[q][j][i].size();
                }
                _gHit[q][j] = 0;
            }
        }
        _hit = 0;
    }

    void InitGHit()
    {
        for (auto q = 0; q < _numMetric; ++q)
        {
            for (auto j = 0; j < _commSize[q].size(); ++j)
            {
                _gHit[q][j] = 0;
            }
        }
    }

    /// Generate a set of n RR/G-RR sets
    void BuildRRsets(const size_t numSamples)
    {
        void (*func)(const uint32_t uStart, const size_t hyperIdx);

        if (numSamples > SIZE_MAX)
        {
            std::cout << "Error:R too large" << std::endl;
            exit(1);
        }

        const auto prevSize = _numRRsets;
        _numRRsets = _numRRsets > numSamples ? _numRRsets : numSamples;

        if (_isVanilla) // FIM/RFIM cannot use SUBSIM because they generate G-RR sets
        {
            // std::cout << "Sample RR set by vanilla method" << std::endl;
            if (_funcType == IM)
            {
                for (auto i = prevSize; i < numSamples; i++) BuildOneRRset(dsfmt_gv_genrand_uint32_range(_numV), i);
            }
            else if (_funcType == FIM || _funcType == RFIM)
            {
                for (auto i = prevSize; i < numSamples; i++) BuildOneRRsetGeneralized(i);
            }
            else
            {
                std::cout << "Error: unknown function type" << std::endl;
                exit(1);
            }
            return ;
        }
        if (_probDist == WC)
        {
            // std::cout << "Sample RR sets in WC model" << std::endl;
            for (auto i = prevSize; i < numSamples; i++)
            {
                BuildOneRRsetWeighted(dsfmt_gv_genrand_uint32_range(_numV), i);
            }
            return ;
        }
        else if (_probDist == UNIFORM)
        {
            // std::cout << "Sample RR sets in uniform model" << std::endl;

            for (auto i = prevSize; i < numSamples; i++)
            {
                BuildOneRRsetConstant(dsfmt_gv_genrand_uint32_range(_numV), i);
            }

            return ;
        }
        else
        {
            for (auto i = prevSize; i < numSamples; i++)
            {
                BuildOneRRset(dsfmt_gv_genrand_uint32_range(_numV), i);;
            }
        }
        return ;
    }

    /// Return the average size of RR (G-RR) sets
    double EvalHyperedgeAvg()
    {
        return _hyperedgeAvgEval;
    }

    void display_hyperedge_stat()
    {
        size_t total_hyperedges = 0;

        for (size_t i = 0; i < _numRRsets; i++)
        {
            total_hyperedges += _RRsets[i].size();
        }

        double ave_hyperedge_size = 1.0 * total_hyperedges / _numRRsets;
        double var = 0.0;

        for (size_t i = 0; i < _numRRsets; i++)
        {
            double diff = _RRsets[i].size() - ave_hyperedge_size;
            var += (diff * diff);
        }

        var = var / _numRRsets;
        std::cout << "final average RR set size: " << ave_hyperedge_size << ", final variance: " << var << std::endl;
        return ;
    }

    // Calculate the average size of hyperedges (RRset or G-RRset)
    double HyperedgeAvg()
    {
        size_t totalHyperedges = 0;
        double avgSize = 0.0;

        switch (_funcType)
        {
            case IM:
            {
                for (size_t i = 0; i < _numRRsets; i++)
                {
                    totalHyperedges += _RRsets[i].size();
                }
                avgSize = 1.0 * totalHyperedges / _numRRsets;
                break;
            }
            case FIM:
            {
                int q = _metricID;
                for (auto j = 0; j < _commSize[q].size(); ++j)
                {
                    for (auto& rrset : _GRRsets[q][j])
                    {
                        totalHyperedges += rrset.size();
                    }
                }
                avgSize = 1.0 * totalHyperedges / _numRRsets / _commSize[q].size();
                break;
            }
            case RFIM:
            {
                int totalCommNum = 0;
                for (auto q = 0; q < _numMetric; ++q)
                {
                    totalCommNum += _commSize[q].size();
                    for (auto j = 0; j < _commSize[q].size(); ++j)
                    {
                        for (auto& rrset : _GRRsets[q][j])
                        {
                            totalHyperedges += rrset.size();
                        }
                    }
                }
                avgSize = 1.0 * totalHyperedges / _numRRsets / totalCommNum;
                break;
            }
            default:
            {
                std::cout << "Error: unknown function type" << std::endl;
                exit(1);
            }
        }

        return avgSize;
    }

    // Generate one RR set
    void BuildOneRRset(const uint32_t uStart, const size_t hyperIdx)
    {
        size_t numVisitNode = 0, currIdx = 0;
        _FRsets[uStart].push_back(hyperIdx);
        _vecVisitNode[numVisitNode++] = uStart;
        _vecVisitBool[uStart] = true;

        while (currIdx < numVisitNode)
        {
            const auto expand = _vecVisitNode[currIdx++];

            for (auto& nbr : _graph[expand])
            {
                const auto nbrId = nbr.first;

                if (_vecVisitBool[nbrId])
                    continue;

                const auto randDouble = dsfmt_gv_genrand_open_close();

                if (randDouble > nbr.second)
                    continue;

                _vecVisitNode[numVisitNode++] = nbrId;
                _vecVisitBool[nbrId] = true;
                _FRsets[nbrId].push_back(hyperIdx);
            }
        }

        for (int i = 0; i < numVisitNode; i++) _vecVisitBool[_vecVisitNode[i]] = false;

        _RRsets.push_back(RRset(_vecVisitNode.begin(), _vecVisitNode.begin() + numVisitNode));
    }

    // Generate one Generalized RR set for FIM or RFIM
    void BuildOneRRsetGeneralized(const size_t hyperIdx)
    {
        uint32_t numVisitNode = 0, currIdx = 0, numVisitNodeTmp;
        switch (_funcType)
        {
            case FIM:
            {
                auto q = _metricID;
                for (auto j = 0; j < _commSize[q].size(); ++j) // select j starting nodes
                {
                    if (_commSize[q][j] == 0) continue;
                    const uint32_t uStart = (*_pcommVec)[q][j][dsfmt_gv_genrand_uint32_range(_commSize[q][j])];
                    _GFRsets[q][j][uStart].push_back(hyperIdx);
                    _GRRsets[q][j].push_back(RRset(1, uStart));
                    if (_vecVisitBool[uStart]) continue;
                    _vecVisitNode[numVisitNode++] = uStart;
                    _vecVisitBool[uStart] = true;
                }
                while (currIdx < numVisitNode) // generate a G-RR set from the selected starting nodes
                {
                    const auto expand = _vecVisitNode[currIdx++];
                    if (_graph[expand].size() == 0) continue;
                    uint32_t idx = 0;
                    for (auto& nbr : _graph[expand])
                    {
                        const auto nbrId = nbr.first;
                        const auto randDouble = dsfmt_gv_genrand_open_close();
                        if (randDouble > nbr.second) {_edgeStatus[expand][idx++]=FAIL; continue;}
                        else {_edgeStatus[expand][idx++]=PASS;}
                        if (_vecVisitBool[nbrId]) continue;
                        _vecVisitNode[numVisitNode++] = nbrId;
                        _vecVisitBool[nbrId] = true;
                    }
                }
                for (auto j = 0; j < _commSize[q].size(); ++j) // update _GRRsets and _GFRsets
                {
                    currIdx = 0, numVisitNodeTmp = 0;

                    if (_commSize[q][j] <= 0) continue;
                    const auto uStart = _GRRsets[q][j][hyperIdx][numVisitNodeTmp];
                    _vecVisitBoolTmp[uStart] = true;
                    _vecVisitNode[numVisitNodeTmp++] = uStart;

                    while(currIdx < numVisitNodeTmp)
                    {
                        const auto expand = _vecVisitNode[currIdx++];
                        uint32_t idx = 0;
                        for (auto& nbr : _graph[expand])
                        {
                            if (_vecVisitBool[nbr.first] && !_vecVisitBoolTmp[nbr.first] && _edgeStatus[expand][idx] == PASS)
                            {
                                _vecVisitNode[numVisitNodeTmp++] = nbr.first;
                                _GFRsets[q][j][nbr.first].push_back(hyperIdx);
                                _vecVisitBoolTmp[nbr.first] = true;
                            }
                            idx++;
                        }
                    }
                    for (int i = 0; i < numVisitNodeTmp; i++) _vecVisitBoolTmp[_vecVisitNode[i]] = false;
                    _GRRsets[q][j][hyperIdx].insert(_GRRsets[q][j][hyperIdx].end(), _vecVisitNode.begin() + 1, _vecVisitNode.begin() + numVisitNodeTmp);
                }
                refresh_edgeStatus(true);
                std::fill(_vecVisitBool.begin(), _vecVisitBool.end(), false);
                break;
            }
            case RFIM:
            {
                for (auto q = 0; q < _numMetric; ++q)
                {
                    numVisitNode = 0;
                    for (auto j = 0; j < _commSize[q].size(); ++j) // select j starting nodes
                    {
                        if (_commSize[q][j] == 0) continue;
                        const uint32_t uStart = (*_pcommVec)[q][j][dsfmt_gv_genrand_uint32_range(_commSize[q][j])];
                        _GFRsets[q][j][uStart].push_back(hyperIdx);
                        _GRRsets[q][j].push_back(RRset(1, uStart));
                        if (_vecVisitBool[uStart]) continue;
                        _vecVisitNode[numVisitNode++] = uStart;
                        _vecVisitBool[uStart] = true;
                    }
                    while (currIdx < numVisitNode) // generate a G-RR set from the selected starting nodes
                    {
                        const auto expand = _vecVisitNode[currIdx++];
                        if (_graph[expand].size() == 0) continue;
                        uint32_t idx = 0;
                        for (auto& nbr : _graph[expand])
                        {
                            const auto nbrId = nbr.first;
                            const auto randDouble = dsfmt_gv_genrand_open_close();
                            if (randDouble > nbr.second) {_edgeStatus[expand][idx++]=FAIL; continue;}
                            else {_edgeStatus[expand][idx++]=PASS;}
                            if (_vecVisitBool[nbrId]) continue;
                            _vecVisitNode[numVisitNode++] = nbrId;
                            _vecVisitBool[nbrId] = true;
                        }
                    }
                    for (auto j = 0; j < _commSize[q].size(); ++j) // update _GRRsets and _GFRsets
                    {
                        currIdx = 0, numVisitNodeTmp = 0;
                        if (_commSize[q][j] <= 0) continue;
                        const auto uStart = _GRRsets[q][j][hyperIdx][numVisitNodeTmp];
                        _vecVisitBoolTmp[uStart] = true;
                        _vecVisitNode[numVisitNodeTmp++] = uStart;

                        while(currIdx < numVisitNodeTmp)
                        {
                            uint32_t idx = 0;
                            const auto expand = _vecVisitNode[currIdx++];
                            for (auto& nbr : _graph[expand])
                            {
                                if (_vecVisitBool[nbr.first] && ! _vecVisitBoolTmp[nbr.first] && _edgeStatus[expand][idx] == PASS)
                                {
                                    _vecVisitNode[numVisitNodeTmp++] = nbr.first;
                                    _GFRsets[q][j][nbr.first].push_back(hyperIdx);
                                    _vecVisitBoolTmp[nbr.first] = true;
                                }
                                idx++;
                            }
                        }
                        for (int i = 0; i < numVisitNodeTmp; i++) _vecVisitBoolTmp[_vecVisitNode[i]] = false;
                        _GRRsets[q][j][hyperIdx].insert(_GRRsets[q][j][hyperIdx].end(), _vecVisitNode.begin() + 1, _vecVisitNode.begin() + numVisitNodeTmp);
                    }
                    refresh_edgeStatus(true);
                    std::fill(_vecVisitBool.begin(), _vecVisitBool.end(), false);
                }
                break;
            }
            default:
            {
                std::cout << "Error: Unknown function type" << std::endl;
                exit(1);
            }
        }
        return ;
    }

    // Independent cascade with weighted probability
    void BuildOneRRsetWeighted(const uint32_t uStart, const size_t hyperIdx)
    {
        size_t numVisitNode = 0, currIdx = 0;
        _FRsets[uStart].push_back(hyperIdx);
        _vecVisitNode[numVisitNode++] = uStart;
        _vecVisitBool[uStart] = true;

        while (currIdx < numVisitNode)
        {
            const auto expand = _vecVisitNode[currIdx++];
            {
                if (_graph[expand].size() == 0) continue;

                double p =  _graph[expand][0].second;
                double log2Prob = Logarithm(1 - p);

                if (p < 1)
                {
                    double prob = dsfmt_gv_genrand_open_close();
                    int startPos = Logarithm(prob) / log2Prob;
                    int endPos = _graph[expand].size();

                    while (startPos < endPos)
                    {
                        const auto nbrId = _graph[expand][startPos].first;

                        if (_vecVisitBool[nbrId])
                        {
                            int increment = Logarithm(dsfmt_gv_genrand_open_close()) / log2Prob;
                            startPos += (increment + 1);
                            continue;
                        }

                        _vecVisitNode[numVisitNode++] = nbrId;
                        _vecVisitBool[nbrId] = true;
                        _FRsets[nbrId].push_back(hyperIdx);
                        int increment = Logarithm(dsfmt_gv_genrand_open_close()) / log2Prob;
                        startPos += increment + 1;
                    }
                }
                else
                {
                    for (auto& nbr : _graph[expand])
                    {
                        const auto nbrId = nbr.first;

                        if (_vecVisitBool[nbrId])
                            continue;

                        _vecVisitNode[numVisitNode++] = nbrId;
                        _vecVisitBool[nbrId] = true;
                        _FRsets[nbrId].push_back(hyperIdx);
                    }
                }
            }
        }

        for (int i = 0; i < numVisitNode; i++) _vecVisitBool[_vecVisitNode[i]] = false;

        _RRsets.push_back(RRset(_vecVisitNode.begin(), _vecVisitNode.begin() + numVisitNode));
    }

    /* independent cascade with constant probability */
    void BuildOneRRsetConstant(const uint32_t uStart, const size_t hyperIdx)
    {
        size_t numVisitNode = 0, currIdx = 0;
        _FRsets[uStart].push_back(hyperIdx);
        _vecVisitNode[numVisitNode++] = uStart;
        _vecVisitBool[uStart] = true;
        const double p = _graph[0][0].second;
        const double const_prob = Logarithm(1 - p);

        if (p == 1)
        {
            while (currIdx < numVisitNode)
            {
                const auto expand = _vecVisitNode[currIdx++];

                //std::cout<<_graph[expand].size()<<std::endl;
                if (_graph[expand].size() == 0) continue;

                for (auto& nbr : _graph[expand])
                {
                    const auto nbrId = nbr.first;

                    //std::cout<<nbr.first<<" "<<nbr.second<<std::endl;
                    if (_vecVisitBool[nbrId])
                        continue;

                    _vecVisitNode[numVisitNode++] = nbrId;
                    _vecVisitBool[nbrId] = true;
                    _FRsets[nbrId].push_back(hyperIdx);
                }
            }
        }
        else
        {
            while (currIdx < numVisitNode)
            {
                const auto expand = _vecVisitNode[currIdx++];

                if (0 == _graph[expand].size())
                {
                    continue;
                }

                int startPos = Logarithm(dsfmt_gv_genrand_open_close()) / const_prob;
                int endPos = _graph[expand].size();

                while (startPos < endPos)
                {
                    //std::cout<<"enter loop"<<std::endl;
                    const auto nbrId = _graph[expand][startPos].first;

                    if (!_vecVisitBool[nbrId])
                    {
                        _vecVisitNode[numVisitNode++] = nbrId;
                        _vecVisitBool[nbrId] = true;
                        _FRsets[nbrId].push_back(hyperIdx);
                    }

                    int increment = Logarithm(dsfmt_gv_genrand_open_close()) / const_prob;
                    startPos += increment + 1;
                }
            }
        }

        for (int i = 0; i < numVisitNode; i++) _vecVisitBool[_vecVisitNode[i]] = false;

        _RRsets.push_back(RRset(_vecVisitNode.begin(), _vecVisitNode.begin() + numVisitNode));
    }

    // Evaluate the influence spread of a seed set on current generated RR sets
    double CalculateInf(const Nodelist& vecSeed)
    {
        std::vector<bool> vecBoolVst = std::vector<bool>(_numRRsets);
        std::vector<bool> vecBoolSeed(_numV);

        for (auto seed : vecSeed) vecBoolSeed[seed] = true;

        for (auto seed : vecSeed)
        {
            for (auto node : _FRsets[seed])
            {
                vecBoolVst[node] = true;
            }
        }

        size_t count = std::count(vecBoolVst.begin(), vecBoolVst.end(), true);
        return 1.0 * count * _numV / _numRRsets;
    }

    // Compute the general coverage of a given seed set
    double ComputeFairInf(const Nodelist& vecSeed)
    {
        std::vector<bool> vecBoolVst = std::vector<bool>(_numRRsets);
        double fInf = 0.0, alpha;

        for (auto q = 0; q < _numMetric; ++q)
        {
            for (auto j = 0; j < _commSize[q].size(); ++j)
            {
                if (_commSize[q][j] == 0) continue;
                alpha = (*_pcommFair)[q][j];
                for (auto seed : vecSeed)
                {
                    for (auto hyperID : _GFRsets[q][j][seed])
                    {
                        vecBoolVst[hyperID] = true;
                    }
                }
                size_t count = std::count(vecBoolVst.begin(), vecBoolVst.end(), true);
                // fInf += alpha * count * _numV / _numRRsets;
                fInf += alpha * count * _commSize[q][j] / _numRRsets;
                std::fill(vecBoolVst.begin(), vecBoolVst.end(), false);
            }
        }
        return fInf;
    }

    void RefreshActComm()
    {
        switch (_funcType)
        {
            case RFIM:
            {
                for (int q = 0; q < _numMetric; q++)
                {
                    for (int j = 0; j < _commSize[q].size(); ++j) _numActComm[q][j] = 0;
                }
                break;
            }
            case FIM:
            {
                for (int j = 0; j < _commSize[_metricID].size(); ++j) _numActComm[_metricID][j] = 0;
                break;
            }
            default:
            {
                std::cout << "Error: Unknown function type" << std::endl;
                exit(1);
            }
        }
        return ;
    }

    // Compute the conservative lower bound fair influence of a seed set
    double ComputeConsLB(const Nodelist& vecSeed)
    {
        RefreshActComm();
        double fInf = 0.0, alpha;
        switch (_funcType)
        {
            case RFIM:
            {
                double currFInf;
                fInf = _isSaturateGreedy ? 0.0 : (double)_numV;
                for (int q = 0; q < _numMetric; q++)
                {
                    currFInf = 0.0;
                    for (auto seed : vecSeed)
                    {
                        _numActComm[q][(*_pvecComm)[q][seed]]++;
                    }
                    for (int j = 0; j < _commSize[q].size(); ++j)
                    {
                        alpha = (*_pcommFair)[q][j];
                        switch (_fairType)
                        {
                            case LINEAR:
                                currFInf += alpha * _numActComm[q][j];
                                break;
                            case SCALE_POWER:
                                currFInf += _commSize[q][j] * pow((double)_numActComm[q][j] / _commSize[q][j], alpha);
                                break;
                            default:
                                std::cout << "Error: Unknown fair type" << std::endl;
                                exit(1);
                        }
                    }
                    if (_isSaturateGreedy)
                    {
                        fInf += std::min((*_p_c), currFInf / (*_pVecOptFairInf)[q]);
                    }
                    else // simple greedy (baseline)
                    {
                        if (currFInf / (*_pVecOptFairInf)[q] < fInf) fInf = currFInf / (*_pVecOptFairInf)[q];
                    }
                }
                break;
            }

            case FIM:
            {
                int q = _metricID;
                for (auto seed : vecSeed)
                {
                    _numActComm[q][(*_pvecComm)[q][seed]]++;
                }
                for (int j = 0; j < _commSize[q].size(); ++j)
                {
                    alpha = (*_pcommFair)[q][j];
                    switch (_fairType)
                    {
                        case LINEAR:
                            fInf += alpha * _numActComm[q][j];
                            break;
                        case SCALE_POWER:
                            fInf += _commSize[q][j] * pow((double)_numActComm[q][j] / _commSize[q][j], alpha);
                            break;
                        default:
                            std::cout << "Error: Unknown fair type" << std::endl;
                            exit(1);
                    }
                }
                break;
            }

            default:
            {
                std::cout << "Error: Unknown function type" << std::endl;
                exit(1);
            }
        }
        return fInf;
    }

    // Compute the conservative upper bound fair influence of when generating numR G-RR sets
    double ComputeConsUB(const size_t numR)
    {
        double fInf = 0.0, alpha;
        switch (_funcType)
        {
            case RFIM:
            {
                double currFInf;
                fInf = _isSaturateGreedy ? 0.0 : (double)_numV;
                for (int q = 0; q < _numMetric; q++)
                {
                    currFInf = 0.0;
                    for (int j = 0; j < _commSize[q].size(); ++j)
                    {
                        alpha = (*_pcommFair)[q][j];
                        switch (_fairType)
                        {
                            case LINEAR:
                                currFInf += alpha * numR;
                                break;
                            case SCALE_POWER:
                                currFInf += _commSize[q][j] * pow((double)numR / _commSize[q][j], alpha);
                                break;
                            default:
                                std::cout << "Error: Unknown fair type" << std::endl;
                                exit(1);
                        }
                    }
                    if (_isSaturateGreedy)
                    {
                        fInf += std::min((*_p_c), currFInf / (*_pVecOptFairInf)[q]);
                    }
                    else // simple greedy (baseline)
                    {
                        if (currFInf / (*_pVecOptFairInf)[q] < fInf) fInf = currFInf / (*_pVecOptFairInf)[q];
                    }
                }
                break;
            }

            case FIM:
            {
                int q = _metricID;
                for (int j = 0; j < _commSize[q].size(); ++j)
                {
                    alpha = (*_pcommFair)[q][j];
                    switch (_fairType)
                    {
                        case LINEAR:
                            fInf += alpha * numR;
                            break;
                        case SCALE_POWER:
                            fInf += _commSize[q][j] * pow((double)numR / _commSize[q][j], alpha);
                            break;
                        default:
                            std::cout << "Error: Unknown fair type" << std::endl;
                            exit(1);
                    }
                }
                break;
            }

            default:
            {
                std::cout << "Error: Unknown function type" << std::endl;
                exit(1);
            }
        }

        return fInf;
    }

    // Compute the marginal fair influence of a node v given the current seed set
    double ComputeMarginalFairInf(const uint32_t v)
    {
        double fInf = 0.0, alpha;

        switch (_funcType)
        {
            case RFIM:
            {
                double currMarginalFInf; // The marginal fair influence after adding node v in metric q
                double currFInf; // The fair influence before adding node v in mertic q
                fInf = _isSaturateGreedy ? 0.0 : (double)_numV;
                for (auto q = 0; q < _numMetric; ++q)
                {
                    currMarginalFInf = 0.0; currFInf = 0.0;
                    for (auto j = 0; j < _commSize[q].size(); ++j)
                    {
                        if (_commSize[q][j] == 0) continue;
                        alpha = (*_pcommFair)[q][j];
                        switch (_fairType)
                        {
                        case LINEAR:
                            currMarginalFInf += alpha * _GCoverage[q][j][v] * _commSize[q][j] / _numRRsets;
                            currFInf += alpha * _numGCovered[q][j] * _commSize[q][j] / _numRRsets;
                            break;
                        case SCALE_POWER:
                            currMarginalFInf += alpha * _commSize[q][j] * computeMarginalScalePower(_numGCovered[q][j], _GCoverage[q][j][v], alpha, _numRRsets);
                            currFInf += _commSize[q][j] * (1 - alpha * computeScalePower(_numGCovered[q][j], alpha, _numRRsets));
                            break;
                        default:
                            std::cout << "Error: Unknown fair type" << std::endl;
                            exit(1);
                        }
                    }
                    if (_isSaturateGreedy)
                    {
                        fInf += std::min((*_p_c), (currMarginalFInf + currFInf) / (*_pVecOptFairInf)[q])
                            - std::min((*_p_c), currFInf / (*_pVecOptFairInf)[q]);
                    }
                    else // single greedy (baseline)
                    {
                        if (currMarginalFInf / (*_pVecOptFairInf)[q] < fInf) fInf = currMarginalFInf / (*_pVecOptFairInf)[q];
                    }
                }
                break;
            }

            case FIM:
            {
                int q = _metricID;
                for (auto j = 0; j < _commSize[q].size(); ++j)
                {
                    if (_commSize[q][j] == 0) continue;
                    alpha = (*_pcommFair)[q][j];
                    switch (_fairType)
                    {
                    case LINEAR: // LINEAR fair do not need to use _numGCovered
                        fInf += alpha * _GCoverage[q][j][v] * _commSize[q][j] / _numRRsets;
                        break;
                    case SCALE_POWER:
                        fInf += alpha * _commSize[q][j] * computeMarginalScalePower(_numGCovered[q][j], _GCoverage[q][j][v], alpha, _numRRsets);
                        break;
                    default:
                        std::cout << "Error: Unknown fair type" << std::endl;
                        exit(1);
                    }
                }
                break;
            }

            default:
            {
                std::cout << "Error: Unknown function type" << std::endl;
                exit(1);
            }
        }

        return fInf;
    }

    // Update the G-coverage status
    void UpdateCoverageStatus(const uint32_t v, std::vector<bool>& needRecompute)
    {
        for (auto q = 0; q < _numMetric; ++q)
        {
            if (_funcType == RFIM) (*_pvecFairInf)[q] = 0.0;
            for (auto j = 0; j < _commSize[q].size(); ++j)
            {
                if (_commSize[q][j] == 0) continue;
                for (auto hyperID : _GFRsets[q][j][v])
                {
                    if (_GIsCovered[q][j][hyperID]) continue;
                    _GIsCovered[q][j][hyperID] = true;
                    _numGCovered[q][j]++;
                    for (auto node : _GRRsets[q][j][hyperID])
                    {
                        needRecompute[node] = true;
                        _GCoverage[q][j][node]--;
                    }
                }
                double alpha = (*_pcommFair)[q][j];
                if (_funcType == RFIM) (*_pvecFairInf)[q] += alpha * _numGCovered[q][j] * _commSize[q][j] / _numRRsets;
            }
        }
        return ;
    }

    // Efficiently estimate the influence spread with sampling error epsilon within probability 1-delta
    double EfficInfVldtAlg(const Nodelist& vecSeed, const double delta = 1e-3, const double eps = 0.01)
    {
        const double c = 2.0 * (exp(1.0) - 2.0);
        const double LambdaL = 1.0 + 2.0 * c * (1.0 + eps) * log(2.0 / delta) / (eps * eps);
        size_t numHyperEdge = 0;
        size_t numCoverd = 0;
        std::vector<bool> vecBoolSeed(_numV);

        for (auto seed : vecSeed) vecBoolSeed[seed] = true;

        while (numCoverd < LambdaL)
        {
            numHyperEdge++;
            size_t numVisitNode = 0, currIdx = 0;
            const auto uStart = dsfmt_gv_genrand_uint32_range(_numV);

            if (vecBoolSeed[uStart])
            {
                // Stop, this sample is covered
                numCoverd++;
                continue;
            }

            _vecVisitNode[numVisitNode++] = uStart;
            _vecVisitBool[uStart] = true;

            while (currIdx < numVisitNode)
            {
                const auto expand = _vecVisitNode[currIdx++];

                for (auto& nbr : _graph[expand])
                {
                    const auto nbrId = nbr.first;

                    if (_vecVisitBool[nbrId])
                        continue;

                    const auto randDouble = dsfmt_gv_genrand_open_close();

                    if (randDouble > nbr.second)
                        continue;

                    if (vecBoolSeed[nbrId])
                    {
                        // Stop, this sample is covered
                        numCoverd++;
                        goto postProcess;
                    }

                    _vecVisitNode[numVisitNode++] = nbrId;
                    _vecVisitBool[nbrId] = true;
                }
            }

postProcess:

            for (auto i = 0; i < numVisitNode; i++)
                _vecVisitBool[_vecVisitNode[i]] = false;
        }

        return 1.0 * numCoverd * _numV / numHyperEdge;
    }

    // Refresh the hypergraph (GRR, GFR and gHit)
    void RefreshHypergraph()
    {
        switch (_funcType)
        {
            case IM:
            {
                if (_RRsets.size() != 0)
                {
                    for (auto i = _numRRsets; i--;) RRset().swap(_RRsets[i]);
                    RRsets().swap(_RRsets);

                    for (auto i = _numV; i--;) FRset().swap(_FRsets[i]);
                }
                break;
            }

            case FIM:
            {
                int q = _metricID;
                if (_numRRsets > 0 || _GRRsets[q][0].size() != 0)
                {
                    for (auto j = _commSize[q].size(); j--;)
                    {
                        int currRRsetsSize = _GRRsets[q][j].size();
                        for (auto i = currRRsetsSize; i--;) RRset().swap(_GRRsets[q][j][i]);
                        RRsets().swap(_GRRsets[q][j]);

                        for (auto i = _numV; i--;) FRset().swap(_GFRsets[q][j][i]);

                        _gHit[q][j] = 0;
                    }
                }
                break;
            }

            case RFIM:
            {
                if (_numRRsets > 0 || _GRRsets[0][0].size() != 0)
                {
                    for (auto q = _numMetric; q--;) // refresh G-RR/G-FR sets and _gHit
                    {
                        for (auto j = _commSize[q].size(); j--;)
                        {
                            int currRRsetsSize = _GRRsets[q][j].size();
                            for (auto i = currRRsetsSize; i--;) RRset().swap(_GRRsets[q][j][i]);
                            RRsets().swap(_GRRsets[q][j]);

                            for (auto i = _numV; i--;) FRset().swap(_GFRsets[q][j][i]);

                            _gHit[q][j] = 0;
                        }
                    }
                }
                break;
            }

            default:
            {
                std::cout << "Error: Unknown function type" << std::endl;
                exit(1);
            }
        }
        _numRRsets = 0;
        _hit = 0;
    }

    // Release memory
    void ReleaseMemory()
    {
        RefreshHypergraph();
        std::vector<bool>().swap(_vecVisitBool);
        Nodelist().swap(_vecVisitNode);
        FRsets().swap(_FRsets);
    }

    void BuildOneRRsetEarlyStopByVanilla(std::unordered_set<uint32_t> &connSet, const uint32_t uStart, const size_t hyperIdx)
    {
        size_t numVisitNode = 0, currIdx = 0;
        _FRsets[uStart].push_back(hyperIdx);
        _vecVisitNode[numVisitNode++] = uStart;
        _vecVisitBool[uStart] = true;

        if (connSet.find(uStart) != connSet.end())
        {
            _hit++;
            goto finished;
        }

        while (currIdx < numVisitNode)
        {
            const auto expand = _vecVisitNode[currIdx++];

            if (0 == _graph[expand].size())
            {
                continue;
            }

            for (auto& nbr : _graph[expand])
            {
                const auto nbrId = nbr.first;

                if (_vecVisitBool[nbrId])
                    continue;

                const auto randDouble = dsfmt_gv_genrand_open_close();

                if (randDouble > nbr.second)
                    continue;

                _vecVisitNode[numVisitNode++] = nbrId;
                _vecVisitBool[nbrId] = true;
                _FRsets[nbrId].push_back(hyperIdx);

                if (connSet.find(nbrId) != connSet.end())
                {
                    _hit++;
                    goto finished;
                }
            }
        }

finished:

        for (int i = 0; i < numVisitNode; i++) _vecVisitBool[_vecVisitNode[i]] = false;

        _RRsets.push_back(RRset(_vecVisitNode.begin(), _vecVisitNode.begin() + numVisitNode));
    }

    void BuildOneRRsetEarlyStopBySubsim(std::unordered_set<uint32_t> &connSet, const uint32_t uStart, const size_t hyperIdx)
    {
        size_t numVisitNode = 0, currIdx = 0;
        _FRsets[uStart].push_back(hyperIdx);
        _vecVisitNode[numVisitNode++] = uStart;
        _vecVisitBool[uStart] = true;

        if (connSet.find(uStart) != connSet.end())
        {
            _hit++;
            goto finished;
        }

        while (currIdx < numVisitNode)
        {
            const auto expand = _vecVisitNode[currIdx++];
            const double p =  _graph[expand][0].second;

            if (0 == _graph[expand].size())
            {
                continue;
            }

            if (p >= 1.0)
            {
                for (auto& nbr : _graph[expand])
                {
                    const auto nbrId = nbr.first;

                    if (_vecVisitBool[nbrId])
                        continue;

                    _vecVisitNode[numVisitNode++] = nbrId;
                    _vecVisitBool[nbrId] = true;
                    _FRsets[nbrId].push_back(hyperIdx);
                }

                continue;
            }

            const double const_prob = Logarithm(1 - p);
            int startPos = Logarithm(dsfmt_gv_genrand_open_close()) / const_prob;
            int endPos = _graph[expand].size();

            while (startPos < endPos)
            {
                const auto nbrId = _graph[expand][startPos].first;

                if (!_vecVisitBool[nbrId])
                {
                    _vecVisitNode[numVisitNode++] = nbrId;
                    _vecVisitBool[nbrId] = true;
                    _FRsets[nbrId].push_back(hyperIdx);
                }

                if (connSet.find(nbrId) != connSet.end())
                {
                    _hit++;
                    goto finished;
                }

                int increment = Logarithm(dsfmt_gv_genrand_open_close()) / const_prob;
                startPos += increment + 1;
            }
        }

finished:

        for (int i = 0; i < numVisitNode; i++) _vecVisitBool[_vecVisitNode[i]] = false;

        _RRsets.push_back(RRset(_vecVisitNode.begin(), _vecVisitNode.begin() + numVisitNode));
    }

    void BuildOneRRsetGenerlizedEarlyStop(std::unordered_set<uint32_t> &connSet, const size_t hyperIdx)
    {
        uint32_t numVisitNode = 0, currIdx = 0, numVisitNodeTmp;
        bool hit;
        switch (_funcType)
        {
            case FIM:
            {
                auto q = _metricID;
                for (auto j = 0; j < _commSize[q].size(); ++j) // select j starting nodes
                {
                    if (_commSize[q][j] == 0) continue;
                    const uint32_t uStart = (*_pcommVec)[q][j][dsfmt_gv_genrand_uint32_range(_commSize[q][j])];
                    _GFRsets[q][j][uStart].push_back(hyperIdx);
                    _GRRsets[q][j].push_back(RRset(1, uStart));
                    if (_vecVisitBool[uStart]) continue;
                    _vecVisitNode[numVisitNode++] = uStart;
                    _vecVisitBool[uStart] = true;
                }
                while (currIdx < numVisitNode) // generate a G-RR set from the selected starting nodes
                {
                    const auto expand = _vecVisitNode[currIdx++];
                    if (connSet.find(expand) != connSet.end()) continue; // hit, stop detecting more neighbors
                    uint32_t idx = 0;
                    for (auto& nbr : _graph[expand])
                    {
                        const auto nbrId = nbr.first;
                        const auto randDouble = dsfmt_gv_genrand_open_close();
                        if (randDouble > nbr.second) {_edgeStatus[expand][idx++]=FAIL; continue;}
                        else {_edgeStatus[expand][idx++]=PASS;}
                        if (_vecVisitBool[nbrId]) continue;
                        _vecVisitNode[numVisitNode++] = nbrId;
                        _vecVisitBool[nbrId] = true;
                    }
                }
                for (auto j = 0; j < _commSize[q].size(); ++j) // update _GRRsets and _GFRsets
                {
                    currIdx = 0, numVisitNodeTmp = 0;
                    hit = false;

                    if (_commSize[q][j] <= 0) continue;
                    const auto uStart = _GRRsets[q][j][hyperIdx][numVisitNodeTmp];
                    if (connSet.find(uStart) != connSet.end()) // hit
                    {
                        _gHit[q][j]++;
                        hit = true;
                    }
                    _vecVisitBoolTmp[uStart] = true;
                    _vecVisitNode[numVisitNodeTmp++] = uStart;

                    while(currIdx < numVisitNodeTmp && ! hit)
                    {
                        const auto expand = _vecVisitNode[currIdx++];
                        uint32_t idx = 0;
                        for (auto& nbr : _graph[expand])
                        {
                            EdgeStatus es = _edgeStatus[expand][idx++];
                            if (_vecVisitBool[nbr.first] && !_vecVisitBoolTmp[nbr.first] && es==PASS)
                            {
                                _vecVisitNode[numVisitNodeTmp++] = nbr.first;
                                _GFRsets[q][j][nbr.first].push_back(hyperIdx);
                                _vecVisitBoolTmp[nbr.first] = true;

                                if (connSet.find(nbr.first) != connSet.end()) // hit
                                {
                                    _gHit[q][j]++;
                                    hit = true;
                                    break;
                                }
                            }
                        }
                    }
                    for (int i = 0; i < numVisitNodeTmp; i++) _vecVisitBoolTmp[_vecVisitNode[i]] = false;
                    _GRRsets[q][j][hyperIdx].insert(_GRRsets[q][j][hyperIdx].end(), _vecVisitNode.begin() + 1, _vecVisitNode.begin() + numVisitNodeTmp);
                }
                refresh_edgeStatus(true);
                std::fill(_vecVisitBool.begin(), _vecVisitBool.end(), false);
                break;
            }
            case RFIM:
            {
                for (auto q = 0; q < _numMetric; ++q) 
                {
                    numVisitNode = 0; currIdx = 0;
                    for (auto j = 0; j < _commSize[q].size(); ++j) // select j starting nodes
                    {
                        if (_commSize[q][j] == 0) continue;
                        const uint32_t uStart = (*_pcommVec)[q][j][dsfmt_gv_genrand_uint32_range(_commSize[q][j])];
                        _GFRsets[q][j][uStart].push_back(hyperIdx);
                        _GRRsets[q][j].push_back(RRset(1, uStart));
                        if (_vecVisitBool[uStart]) continue;
                        _vecVisitNode[numVisitNode++] = uStart;
                        _vecVisitBool[uStart] = true;
                    }
                    while (currIdx < numVisitNode) // generate a G-RR set from the selected starting nodes
                    {
                        const auto expand = _vecVisitNode[currIdx++];
                        if (connSet.find(expand) != connSet.end()) continue; // hit, stop detecting more neighbors
                        uint32_t idx = 0;
                        for (auto& nbr : _graph[expand])
                        {
                            const auto nbrId = nbr.first;
                            const auto randDouble = dsfmt_gv_genrand_open_close();
                            if (randDouble > nbr.second) {_edgeStatus[expand][idx++]=FAIL; continue;}
                            else {_edgeStatus[expand][idx++]=PASS;}
                            if (_vecVisitBool[nbrId]) continue;
                            _vecVisitNode[numVisitNode++] = nbrId;
                            _vecVisitBool[nbrId] = true;
                        }
                    }
                    for (auto j = 0; j < _commSize[q].size(); ++j)
                    {
                        currIdx = 0, numVisitNodeTmp = 0;
                        hit = false;

                        if (_commSize[q][j] <= 0) continue;
                        const auto uStart = _GRRsets[q][j][hyperIdx][numVisitNodeTmp];
                        if (connSet.find(uStart) != connSet.end()) // hit
                        {
                            _gHit[q][j]++;
                            hit = true;
                        }
                        _vecVisitBoolTmp[uStart] = true;
                        _vecVisitNode[numVisitNodeTmp++] = uStart;

                        while(currIdx < numVisitNodeTmp && ! hit)
                        {
                            const auto expand = _vecVisitNode[currIdx++];
                            uint32_t idx = 0;
                            for (auto& nbr : _graph[expand])
                            {
                                EdgeStatus es = _edgeStatus[expand][idx++];
                                if (_vecVisitBool[nbr.first] && ! _vecVisitBoolTmp[nbr.first] && es==PASS)
                                {
                                    _vecVisitNode[numVisitNodeTmp++] = nbr.first;
                                    _GFRsets[q][j][nbr.first].push_back(hyperIdx);
                                    _vecVisitBoolTmp[nbr.first] = true;

                                    if (connSet.find(nbr.first) != connSet.end()) // hit
                                    {
                                        _gHit[q][j]++;
                                        hit = true;
                                        break;
                                    }
                                }
                            }
                        }
                        for (int i = 0; i < numVisitNodeTmp; i++) _vecVisitBoolTmp[_vecVisitNode[i]] = false;
                        _GRRsets[q][j][hyperIdx].insert(_GRRsets[q][j][hyperIdx].end(), _vecVisitNode.begin() + 1, _vecVisitNode.begin() + numVisitNodeTmp);
                    }
                    refresh_edgeStatus(true);
                    std::fill(_vecVisitBool.begin(), _vecVisitBool.end(), false);
                }
                break;
            }
            default:
            {
                std::cout << "Error: Unknown function type" << std::endl;
                exit(1);
            }
        }
        return ;
    }

    void BuildRRsetsEarlyStop(std::unordered_set<uint32_t> &connSet, const int numSamples)
    {
        const auto prevSize = _numRRsets;
        _numRRsets = _numRRsets > numSamples ? _numRRsets : numSamples;

        if (_isVanilla)
        {
            // std::cout << "Sample RR sets with early stop by Vanilla method" << std::endl;
            if (_funcType == IM)
            {
                for (auto i = prevSize; i < numSamples; i++)
                BuildOneRRsetEarlyStopByVanilla(connSet, dsfmt_gv_genrand_uint32_range(_numV), i);
            }
            else if(_funcType == FIM || _funcType == RFIM) // G-RR set generation do not use SUBSIM
            {
                for (auto i = prevSize; i < numSamples; i++) BuildOneRRsetGenerlizedEarlyStop(connSet, i);
            }
            else
            {
                std::cout << "Error: unknown function type" << std::endl;
                exit(1);
            }
        }
        else
        {
            // std::cout << "Sample RR sets with early stop By SUBSIM" << std::endl;
            for (auto i = prevSize; i < numSamples; i++)
            {
                BuildOneRRsetEarlyStopBySubsim(connSet, dsfmt_gv_genrand_uint32_range(_numV), i);
            }
        }
    }

    double EvalSeedSetInfByVanilla(std::unordered_set<uint32_t> &connSet, const int numSamples)
    {
        uint32_t numCovered = 0;
        int64_t totalHyperedgeSize = 0;
        _numSamplesEval = numSamples;

        for (int i = 1; i < numSamples; i++)
        {
            uint32_t uStart = dsfmt_gv_genrand_uint32_range(_numV);
            size_t numVisitNode = 0, currIdx = 0;
            _vecVisitNode[numVisitNode++] = uStart;
            _vecVisitBool[uStart] = true;

            if (connSet.find(uStart) != connSet.end())
            {
                numCovered++;
                goto finished;
            }

            while (currIdx < numVisitNode)
            {
                const auto expand = _vecVisitNode[currIdx++];

                if (0 == _graph[expand].size())
                {
                    continue;
                }

                const double p =  _graph[expand][0].second;

                for (auto& nbr : _graph[expand])
                {
                    const auto nbrId = nbr.first;

                    if (_vecVisitBool[nbrId])
                        continue;

                    const auto randDouble = dsfmt_gv_genrand_open_close();

                    if (randDouble > nbr.second)
                        continue;

                    _vecVisitNode[numVisitNode++] = nbrId;
                    _vecVisitBool[nbrId] = true;

                    if (connSet.find(nbrId) != connSet.end())
                    {
                        numCovered++;
                        goto finished;
                    }
                }
            }

        finished:
            totalHyperedgeSize += numVisitNode;

            for (int i = 0; i < numVisitNode; i++) _vecVisitBool[_vecVisitNode[i]] = false;
        }

        _hyperedgeAvgEval = 1.0 * totalHyperedgeSize / numSamples;
        return 1.0 * numCovered * _numV / numSamples;
    }

    double EvalSeedSetInfBySubsim(std::unordered_set<uint32_t> &connSet, const int numSamples)
    {
        uint32_t numCovered = 0;
        int64_t totalHyperedgeSize = 0;
        _numSamplesEval = numSamples;

        for (int i = 1; i < numSamples; i++)
        {
            uint32_t uStart = dsfmt_gv_genrand_uint32_range(_numV);
            size_t numVisitNode = 0, currIdx = 0;
            _vecVisitNode[numVisitNode++] = uStart;
            _vecVisitBool[uStart] = true;

            if (connSet.find(uStart) != connSet.end())
            {
                numCovered++;
                goto finished;
            }

            while (currIdx < numVisitNode)
            {
                const auto expand = _vecVisitNode[currIdx++];

                if (0 == _graph[expand].size())
                {
                    continue;
                }

                const double p =  _graph[expand][0].second;

                if (p >= 1.0)
                {
                    for (auto& nbr : _graph[expand])
                    {
                        const auto nbrId = nbr.first;

                        if (_vecVisitBool[nbrId])
                            continue;

                        _vecVisitNode[numVisitNode++] = nbrId;
                        _vecVisitBool[nbrId] = true;

                        if (connSet.find(nbrId) != connSet.end())
                        {
                            numCovered++;
                            goto finished;
                        }
                    }

                    continue;
                }

                const double const_prob = Logarithm(1 - p);
                int startPos = Logarithm(dsfmt_gv_genrand_open_close()) / const_prob;
                int endPos = _graph[expand].size();

                while (startPos < endPos)
                {
                    const auto nbrId = _graph[expand][startPos].first;

                    if (!_vecVisitBool[nbrId])
                    {
                        _vecVisitNode[numVisitNode++] = nbrId;
                        _vecVisitBool[nbrId] = true;
                    }

                    if (connSet.find(nbrId) != connSet.end())
                    {
                        numCovered++;
                        goto finished;
                    }

                    int increment = Logarithm(dsfmt_gv_genrand_open_close()) / const_prob;
                    startPos += increment + 1;
                }
            }

finished:
            totalHyperedgeSize += numVisitNode;

            for (int i = 0; i < numVisitNode; i++) _vecVisitBool[_vecVisitNode[i]] = false;
        }

        _hyperedgeAvgEval = 1.0 * totalHyperedgeSize / numSamples;
        return 1.0 * numCovered * _numV / numSamples;
    }

    double EvalSeedSetGeneralizedFInf(std::unordered_set<uint32_t> &connSet, const size_t numSamples)
    {
        RefreshHypergraph();
        size_t numVisitNode, currIdx, numVisitNodeTmp, totalSizeHyperEdges = 0, numHyperEdges = 0;
        bool hit;
        for (size_t hyperIdx = 0; hyperIdx < numSamples; hyperIdx++) // sample numSamples G-RR sets and compute _gHit
        {
            numVisitNode = 0, currIdx = 0;
            switch (_funcType)
            {
                case FIM:
                {
                    auto q = _metricID;
                    for (auto j = 0; j < _commSize[q].size(); ++j) // select j starting nodes
                    {
                        if (_commSize[q][j] == 0) continue;
                        const uint32_t uStart = (*_pcommVec)[q][j][dsfmt_gv_genrand_uint32_range(_commSize[q][j])];
                        _GRRsets[q][j].push_back(RRset(1, uStart));
                        if (_vecVisitBool[uStart]) continue;
                        _vecVisitNode[numVisitNode++] = uStart;
                        _vecVisitBool[uStart] = true;
                    }
                    while (currIdx < numVisitNode) // generate a G-RR set from the selected starting nodes
                    {
                        const auto expand = _vecVisitNode[currIdx++];
                        if (connSet.find(expand) != connSet.end()) continue; // hit, stop detecting more neighbors
                        uint32_t idx = 0;
                        for (auto& nbr : _graph[expand])
                        {
                            const auto nbrId = nbr.first;
                            const auto randDouble = dsfmt_gv_genrand_open_close();
                            if (randDouble > nbr.second) {_edgeStatus[expand][idx++]=FAIL; continue;}
                            else {_edgeStatus[expand][idx++]=PASS;}
                            if (_vecVisitBool[nbrId]) continue;
                            _vecVisitNode[numVisitNode++] = nbrId;
                            _vecVisitBool[nbrId] = true;
                        }
                    }
                    for (auto j = 0; j < _commSize[q].size(); ++j) // reverse propagation from each source
                    {
                        currIdx = 0, numVisitNodeTmp = 0;
                        hit = false;
    
                        if (_commSize[q][j] <= 0) continue;
                        const auto uStart = _GRRsets[q][j][hyperIdx][numVisitNodeTmp];
                        if (connSet.find(uStart) != connSet.end()) // hit
                        {
                            _gHit[q][j]++;
                            hit = true;
                        }
                        _vecVisitBoolTmp[uStart] = true;
                        _vecVisitNode[numVisitNodeTmp++] = uStart;
    
                        while(currIdx < numVisitNodeTmp && ! hit)
                        {
                            const auto expand = _vecVisitNode[currIdx++];
                            uint32_t idx = 0;
                            for (auto& nbr : _graph[expand])
                            {
                                EdgeStatus es = _edgeStatus[expand][idx++];
                                if (_vecVisitBool[nbr.first] && ! _vecVisitBoolTmp[nbr.first] && es==PASS)
                                {
                                    _vecVisitNode[numVisitNodeTmp++] = nbr.first;
                                    _vecVisitBoolTmp[nbr.first] = true;
    
                                    if (connSet.find(nbr.first) != connSet.end()) // hit
                                    {
                                        _gHit[q][j]++;
                                        hit = true;
                                        break;
                                    }
                                }
                            }
                        }
                        for (int i = 0; i < numVisitNodeTmp; i++) _vecVisitBoolTmp[_vecVisitNode[i]] = false;
                        totalSizeHyperEdges += numVisitNodeTmp;
                        numHyperEdges++;
                    }
                    refresh_edgeStatus(true);
                    std::fill(_vecVisitBool.begin(), _vecVisitBool.end(), false);
                    break;
                }

                case RFIM:
                {
                    for (auto q = 0; q < _numMetric; ++q) // select q * j starting nodes
                    {
                        numVisitNode = 0;
                        currIdx = 0;
                        for (auto j = 0; j < _commSize[q].size(); ++j)
                        {
                            if (_commSize[q][j] == 0) continue;
                            const uint32_t uStart = (*_pcommVec)[q][j][dsfmt_gv_genrand_uint32_range(_commSize[q][j])];
                            // _GFRsets[q][j][uStart].push_back(hyperIdx);
                            _GRRsets[q][j].push_back(RRset(1, uStart));
                            if (_vecVisitBool[uStart]) continue;
                            _vecVisitNode[numVisitNode++] = uStart;
                            _vecVisitBool[uStart] = true;
                        }
                        while (currIdx < numVisitNode) // generate a G-RR set from the selected starting nodes
                        {
                            const auto expand = _vecVisitNode[currIdx++];
                            if (connSet.find(expand) != connSet.end()) continue; // hit, stop detecting more neighbors
                            uint32_t idx = 0;
                            for (auto& nbr : _graph[expand])
                            {
                                const auto nbrId = nbr.first;
                                const auto randDouble = dsfmt_gv_genrand_open_close();
                                if (randDouble > nbr.second) {_edgeStatus[expand][idx++]=FAIL; continue;}
                                else {_edgeStatus[expand][idx++]=PASS;}
                                if (_vecVisitBool[nbrId]) continue;
                                _vecVisitNode[numVisitNode++] = nbrId;
                                _vecVisitBool[nbrId] = true;
                            }
                        }
                        for (auto j = 0; j < _commSize[q].size(); ++j)
                        {
                            currIdx = 0, numVisitNodeTmp = 0;
                            hit = false;

                            if (_commSize[q][j] <= 0) continue;
                            const auto uStart = _GRRsets[q][j][hyperIdx][numVisitNodeTmp];
                            if (connSet.find(uStart) != connSet.end()) // hit
                            {
                                _gHit[q][j]++;
                                hit = true;
                            }
                            _vecVisitBoolTmp[uStart] = true;
                            _vecVisitNode[numVisitNodeTmp++] = uStart;

                            while(currIdx < numVisitNodeTmp && ! hit)
                            {
                                const auto expand = _vecVisitNode[currIdx++];
                                uint32_t idx = 0;
                                for (auto& nbr : _graph[expand])
                                {
                                    EdgeStatus es = _edgeStatus[expand][idx++];
                                    if (_vecVisitBool[nbr.first] && ! _vecVisitBoolTmp[nbr.first] && es==PASS)
                                    {
                                        _vecVisitNode[numVisitNodeTmp++] = nbr.first;
                                        _vecVisitBoolTmp[nbr.first] = true;

                                        if (connSet.find(nbr.first) != connSet.end()) // hit
                                        {
                                            _gHit[q][j]++;
                                            hit = true;
                                            break;
                                        }
                                    }
                                }
                            }
                            for (int i = 0; i < numVisitNodeTmp; i++) _vecVisitBoolTmp[_vecVisitNode[i]] = false;
                            totalSizeHyperEdges += numVisitNodeTmp;
                            numHyperEdges++;
                        }
                        refresh_edgeStatus(true);
                        std::fill(_vecVisitBool.begin(), _vecVisitBool.end(), false);
                    }
                    break;
                }

                default:
                {
                    std::cout << "Error: Unknown function type" << std::endl;
                    exit(1);
                }
            }
        }

        _hyperedgeAvgEval = (double)totalSizeHyperEdges / (double)numHyperEdges;

        switch (_funcType) // use _gHit to estimate the fair influence
        {
            case FIM:
            {
                int q = _metricID;
                double fInf = 0.0;
                for (int j = 0; j < _commSize[q].size(); ++j)
                {
                    if (_commSize[q][j] == 0) continue;
                    double alpha = (*_pcommFair)[q][j];
                    switch (_fairType)
                    {
                        case LINEAR:
                            fInf += alpha * _gHit[q][j] * _commSize[q][j] / numSamples;
                            break;
                        case SCALE_POWER:
                            fInf += _commSize[q][j] * (1 - alpha * computeScalePower(_gHit[q][j], alpha, numSamples));
                            break;
                        default:
                            std::cerr << "Error: Unknown fair type" << std::endl;
                            exit(1);
                    }
                }
                return fInf;
            }

            case RFIM:
            {
                double fInf = _isSaturateGreedy ? 0.0 : (double)_numV;
                for (int q = 0; q < _numMetric; ++q)
                {
                    double currMarginalFInf = 0.0;
                    for (int j = 0; j < _commSize[q].size(); ++j)
                    {
                        if (_commSize[q][j] == 0) continue;
                        double alpha = (*_pcommFair)[q][j];
                        switch (_fairType)
                        {
                            case LINEAR:
                                currMarginalFInf += alpha * _gHit[q][j] * _commSize[q][j] / numSamples;
                                break;
                            case SCALE_POWER:
                                currMarginalFInf += _commSize[q][j] * (1 - alpha * computeScalePower(_gHit[q][j], alpha, numSamples));
                                break;
                            default:
                                std::cerr << "Error: Unknown fair type" << std::endl;
                                exit(1);
                        }
                    }
                    if (_isSaturateGreedy)
                    {
                        fInf += std::min((*_p_c), currMarginalFInf / (*_pVecOptFairInf)[q]);
                        (*_pvecFairInfVldt)[q] = currMarginalFInf;
                    }
                    else
                    {
                        if (currMarginalFInf / (*_pVecOptFairInf)[q] < fInf) fInf = currMarginalFInf / (*_pVecOptFairInf)[q];
                        (*_pvecFairInfVldt)[q] = currMarginalFInf;
                    }
                }
                return fInf;
            }

            default:
            {
                std::cerr << "Error: Unknown function type" << std::endl;
                exit(1);
            }
        }
    }

    double gFInfEst(std::unordered_set<uint32_t> &connSet, const double epsilon, const double delta)
    {
        RefreshHypergraph();
        size_t numVisitNode, currIdx, hyperIdx = 0, numVisitNodeTmp, totalSizeHyperEdges = 0, numHyperEdges = 0, numSamples = 0;
        bool hit;
        double gFairInf = 0.0;
        double a = 0.0;
        double b = 1.0;
        double sumGFairInf = 0.0;

        double gCoverageThreshold = 2 * (b - a) * (1 + epsilon) * ((b - a) / b + epsilon / 3) * log(2 / delta) / pow2(epsilon);

        while (gFairInf * numSamples < gCoverageThreshold)
        {
            numVisitNode = 0; currIdx = 0; numSamples++;
            switch (_funcType) // sample G-RR sets and compute _gHit
            {
                case FIM:
                {
                    auto q = _metricID;
                    for (auto j = 0; j < _commSize[q].size(); ++j) // select j starting nodes
                    {
                        if (_commSize[q][j] == 0) continue;
                        const uint32_t uStart = (*_pcommVec)[q][j][dsfmt_gv_genrand_uint32_range(_commSize[q][j])];
                        _GRRsets[q][j].push_back(RRset(1, uStart));
                        if (_vecVisitBool[uStart]) continue;
                        _vecVisitNode[numVisitNode++] = uStart;
                        _vecVisitBool[uStart] = true;
                    }
                    while (currIdx < numVisitNode) // generate a G-RR set from the selected starting nodes
                    {
                        const auto expand = _vecVisitNode[currIdx++];
                        if (connSet.find(expand) != connSet.end()) continue; // hit, stop detecting more neighbors
                        for (auto& nbr : _graph[expand])
                        {
                            const auto nbrId = nbr.first;
                            if (_vecVisitBool[nbrId]) continue;
                            const auto randDouble = dsfmt_gv_genrand_open_close();
                            if (randDouble > nbr.second) continue;
                            _vecVisitNode[numVisitNode++] = nbrId;
                            _vecVisitBool[nbrId] = true;
                        }
                    }
                    for (auto j = 0; j < _commSize[q].size(); ++j) // reverse propagation from each source
                    {
                        currIdx = 0, numVisitNodeTmp = 0;
                        hit = false;
    
                        if (_commSize[q][j] <= 0) continue;
                        const auto uStart = _GRRsets[q][j][hyperIdx][numVisitNodeTmp];
                        if (connSet.find(uStart) != connSet.end()) // hit
                        {
                            _gHit[q][j]++;
                            hit = true;
                        }
                        _vecVisitBoolTmp[uStart] = true;
                        _vecVisitNode[numVisitNodeTmp++] = uStart;
    
                        while(currIdx < numVisitNodeTmp && ! hit)
                        {
                            const auto expand = _vecVisitNode[currIdx++];
                            for (auto& nbr : _graph[expand])
                            {
                                if (_vecVisitBool[nbr.first] && ! _vecVisitBoolTmp[nbr.first])
                                {
                                    _vecVisitNode[numVisitNodeTmp++] = nbr.first;
                                    _vecVisitBoolTmp[nbr.first] = true;
    
                                    if (connSet.find(nbr.first) != connSet.end()) // hit
                                    {
                                        _gHit[q][j]++;
                                        hit = true;
                                        break;
                                    }
                                }
                            }
                        }
                        for (int i = 0; i < numVisitNodeTmp; i++) _vecVisitBoolTmp[_vecVisitNode[i]] = false;
                    }
                    std::fill(_vecVisitBool.begin(), _vecVisitBool.end(), false);
                    break;
                }
                case RFIM:
                {
                    for (auto q = 0; q < _numMetric; ++q) // select q * j starting nodes
                    {
                        for (auto j = 0; j < _commSize[q].size(); ++j)
                        {
                            if (_commSize[q][j] == 0) continue;
                            const uint32_t uStart = (*_pcommVec)[q][j][dsfmt_gv_genrand_uint32_range(_commSize[q][j])];
                            _GRRsets[q][j].push_back(RRset(1, uStart));
                            if (_vecVisitBool[uStart]) continue;
                            _vecVisitNode[numVisitNode++] = uStart;
                            _vecVisitBool[uStart] = true;
                        }
                    }
                    while (currIdx < numVisitNode) // generate a G-RR set from the selected starting nodes
                    {
                        const auto expand = _vecVisitNode[currIdx++];
                        if (connSet.find(expand) != connSet.end()) continue; // hit, stop detecting more neighbors
                        for (auto& nbr : _graph[expand])
                        {
                            const auto nbrId = nbr.first;
                            if (_vecVisitBool[nbrId]) continue;
                            const auto randDouble = dsfmt_gv_genrand_open_close();
                            if (randDouble > nbr.second) continue;
                            _vecVisitNode[numVisitNode++] = nbrId;
                            _vecVisitBool[nbrId] = true;
                        }
                    }
                    for (auto q = 0; q < _numMetric; ++q) // update _gHit
                    {
                        for (auto j = 0; j < _commSize[q].size(); ++j)
                        {
                            currIdx = 0, numVisitNodeTmp = 0;
                            hit = false;

                            if (_commSize[q][j] <= 0) continue;
                            const auto uStart = _GRRsets[q][j][hyperIdx][numVisitNodeTmp];
                            if (connSet.find(uStart) != connSet.end()) // hit
                            {
                                _gHit[q][j]++;
                                hit = true;
                            }
                            _vecVisitBoolTmp[uStart] = true;
                            _vecVisitNode[numVisitNodeTmp++] = uStart;
    
                            while(currIdx < numVisitNodeTmp && ! hit)
                            {
                                const auto expand = _vecVisitNode[currIdx++];
                                for (auto& nbr : _graph[expand])
                                {
                                    if (_vecVisitBool[nbr.first] && ! _vecVisitBoolTmp[nbr.first])
                                    {
                                        _vecVisitNode[numVisitNodeTmp++] = nbr.first;
                                        _vecVisitBoolTmp[nbr.first] = true;

                                        if (connSet.find(nbr.first) != connSet.end()) // hit
                                        {
                                            _gHit[q][j]++;
                                            hit = true;
                                            break;
                                        }
                                    }
                                }
                            }
                            for (int i = 0; i < numVisitNodeTmp; i++) _vecVisitBoolTmp[_vecVisitNode[i]] = false;
                            totalSizeHyperEdges += numVisitNodeTmp;
                            numHyperEdges++;
                        }
                    }
                    std::fill(_vecVisitBool.begin(), _vecVisitBool.end(), false);

                    break;
                }
                default:
                {
                    std::cout << "Error: Unknown function type" << std::endl;
                    exit(1);
                }
            }

            switch (_funcType) // use _gHit to estimate the general influence
            {
                case FIM:
                {
                    int q = _metricID;
                    double fInf = 0.0;
                    for (int j = 0; j < _commSize[q].size(); ++j)
                    {
                        if (_commSize[q][j] == 0) continue;
                        double alpha = (*_pcommFair)[q][j];
                        switch (_fairType)
                        {
                            case LINEAR:
                                fInf += alpha * _gHit[q][j] * _commSize[q][j] / numSamples;
                                break;
                            case SCALE_POWER:
                                fInf += _commSize[q][j] * (1 - alpha * computeScalePower(_gHit[q][j], alpha, numSamples));
                                break;
                            default:
                                std::cerr << "Error: Unknown fair type" << std::endl;
                                exit(1);
                        }
                    }
                    gFairInf = fInf;
                    break;
                }

                case RFIM:
                {
                    double fInf = _isSaturateGreedy ? 0.0 : (double)_numV;
                    for (int q = 0; q < _numMetric; ++q)
                    {
                        double currMarginalFInf = 0.0;
                        for (int j = 0; j < _commSize[q].size(); ++j)
                        {
                            if (_commSize[q][j] == 0) continue;
                            double alpha = (*_pcommFair)[q][j];
                            switch (_fairType)
                            {
                                case LINEAR:
                                    currMarginalFInf += alpha * _gHit[q][j] * _commSize[q][j] / numSamples;
                                    break;
                                case SCALE_POWER:
                                    currMarginalFInf += _commSize[q][j] * (1 - alpha * computeScalePower(_gHit[q][j], alpha, numSamples));
                                    break;
                                default:
                                    std::cerr << "Error: Unknown fair type" << std::endl;
                                    exit(1);
                            }
                        }
                        if (_isSaturateGreedy)
                        {
                            fInf += std::min((*_p_c), currMarginalFInf / (*_pVecOptFairInf)[q]);
                        }
                        else
                        {
                            if (currMarginalFInf / (*_pVecOptFairInf)[q] < fInf) fInf = currMarginalFInf / (*_pVecOptFairInf)[q];
                        }
                    }
                    gFairInf = fInf;
                    break;
                }
                default:
                {
                    std::cerr << "Error: Unknown function type" << std::endl;
                    exit(1);
                }
            }
            hyperIdx++;
        }
        _numRRsets = numSamples;
        return gFairInf;
    }

    // Evaluate the (fair) influence of seed set by generating a given number of RR (G-RR) sets
    double EvalSeedSetInf(std::unordered_set<uint32_t> &connSet, const size_t numSamples)
    {
        double resInf = 0.0;
        if (_funcType == IM && _isVanilla)
        {
            resInf = EvalSeedSetInfByVanilla(connSet, numSamples);
        }
        else if (_funcType == IM && !_isVanilla)
        {
            resInf = EvalSeedSetInfBySubsim(connSet, numSamples);
        }
        else if (_funcType == FIM || _funcType == RFIM)
        {
            resInf = EvalSeedSetGeneralizedFInf(connSet, numSamples);
        }
        else
        {
            std::cout << "Error: unknown function type" << std::endl;
            exit(1);
        }
        this->_numRRsets = numSamples;
        return resInf;
    }

    // Calculate the fair influence spread of the current seed set on generated RR sets
    double CalculateInfEarlyStop(const bool isSaturate = true)
    {
        switch (_funcType)
        {
            case IM:
            {
                return 1.0 * _hit * _numV / _numRRsets;
            }
            case FIM:
            {
                int q = _metricID;
                double fInf = 0.0, alpha;
                for (int j = 0; j < _commSize[q].size(); ++j)
                {
                    if (_commSize[q][j] == 0) continue;
                    alpha = (*_pcommFair)[q][j];
                    switch (_fairType)
                    {
                        case LINEAR:
                            fInf += alpha * _gHit[q][j] * _commSize[q][j] / _numRRsets;
                            break;
                        case SCALE_POWER:
                            fInf += _commSize[q][j] * (1 - alpha * computeScalePower(_gHit[q][j], (*_pcommFair)[q][j], _numRRsets));
                            break;
                        default:
                            std::cout << "Error: Unknown fair type" << std::endl;
                            exit(1);
                    }
                }
                return fInf;
            }
            case RFIM:
            {
                double fInf = 0.0, alpha;
                double currMarginalFInf;
                fInf = (_isSaturateGreedy && isSaturate) ? 0.0 : (double)_numV;
                for (auto q = 0; q < _numMetric; ++q)
                {
                    currMarginalFInf = 0.0;
                    for (auto j = 0; j < _commSize[q].size(); ++j)
                    {
                        if (_commSize[q][j] == 0) continue;
                        alpha = (*_pcommFair)[q][j];
                        switch (_fairType)
                        {
                            case LINEAR:
                                currMarginalFInf += alpha * _gHit[q][j] * _commSize[q][j] / _numRRsets;
                                break;
                            case SCALE_POWER:
                                currMarginalFInf += _commSize[q][j] * (1 - alpha * computeScalePower(_gHit[q][j], (*_pcommFair)[q][j], _numRRsets));
                                break;
                            default:
                                std::cout << "Error: Unknown fair type" << std::endl;
                                exit(1);
                        }
                    }
                    if (_isSaturateGreedy && isSaturate)
                    {
                        fInf += std::min((*_p_c), currMarginalFInf / (*_pVecOptFairInf)[q]);
                        (*_pvecFairInfVldt)[q] = currMarginalFInf;
                    }
                    else // simple greedy (baseline)
                    {
                        if (currMarginalFInf / (*_pVecOptFairInf)[q] < fInf) fInf = currMarginalFInf / (*_pVecOptFairInf)[q];
                        (*_pvecFairInfVldt)[q] = currMarginalFInf;
                    }
                }
                return fInf;
            }
            default:
            {
                std::cout << "Error: Unknown function type" << std::endl;
                exit(1);
            }
        }
    }

    // MC simulation
    double MonteCarloSimulate(const std::unordered_set<uint32_t>& seedSet, int numSimulations)
    {
        // _numRRsets = 1; RefreshHypergraph();
        InitGHit();

        if (seedSet.empty()) return 0.0;

        double totalInfluence = 0.0;

        for (int sim = 0; sim < numSimulations; ++sim)
        {
            std::fill(_vecVisitBool.begin(), _vecVisitBool.end(), false);
            uint32_t visitIdx = 0;

            for (auto u : seedSet)
            {
                _vecVisitBool[u] = true;
                _vecVisitNode[visitIdx++] = u;
            }

            uint32_t curr = 0;
            while (curr < visitIdx)
            {
                uint32_t u = _vecVisitNode[curr++]; // u is the newly activated node

                // update _gHit for each newly activated node
                for (int q = 0; q < _numMetric; ++q)
                {
                    int comm = (*_pvecComm)[q][u];
                    if (comm != -1)
                    {
                        _gHit[q][comm]++;
                    }
                }

                for (const auto& [v, prob] : _graph[u])
                {
                    if (!_vecVisitBool[v] && dsfmt_gv_genrand_open_close() < prob)
                    {
                        _vecVisitBool[v] = true;
                        _vecVisitNode[visitIdx++] = v;
                    }
                }
            }
        }

        // use _gHit to calculate the G-Influence
        double totalGfairInf = (_funcType == FIM || _funcType == RFIM) ? 
            (_funcType == FIM ? 0.0 : (_isSaturateGreedy ? 0.0 : (double)_numV)) : 0.0;

        for (int q = 0; q < _numMetric; ++q)
        {
            if (_funcType == FIM && _metricID != q) continue;

            double fairInfQ = 0.0;
            for (int j = 0; j < _commSize[q].size(); ++j)
            {
                if (_commSize[q][j] == 0) continue;

                double alpha = (*_pcommFair)[q][j];
                double avgAct = _gHit[q][j] / (double)numSimulations;

                switch (_fairType)
                {
                    case LINEAR:
                        fairInfQ += alpha * avgAct;
                        break;
                    case SCALE_POWER:
                        fairInfQ += _commSize[q][j] * pow(avgAct / _commSize[q][j], alpha);
                        break;
                    default:
                        std::cerr << "Unknown fair type!" << std::endl;
                        exit(1);
                }
            }

            if (_funcType == FIM)
            {
                totalGfairInf = fairInfQ;
            }
            else if (_funcType == RFIM)
            {
                if (_isSaturateGreedy)
                {
                    totalGfairInf += std::min((*_p_c), fairInfQ / (*_pVecOptFairInf)[q]);
                }
                else
                {
                    totalGfairInf = std::min(totalGfairInf, fairInfQ / (*_pVecOptFairInf)[q]);
                }
            }
            else
            {
                std::cerr << "Unknown func type!" << std::endl;
                exit(1);
            }
        }

        // _numRRsets = 1; RefreshHypergraph();
        InitGHit();
        return totalGfairInf;
    }
};

using THyperGraph = HyperGraph;
using PHyperGraph = std::shared_ptr<THyperGraph>;
