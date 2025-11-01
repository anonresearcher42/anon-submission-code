#include "stdafx.h"

// code from OPIM
double Alg::MaxCoverVanilla(const int targetSize)
{
    // optimization with minimum upper bound among all rounds [Default].
    _boundLast = DBL_MAX, _boundMin = DBL_MAX;
    FRset coverage(_numV, 0);
    size_t maxDeg = 0;

    for (auto i = _numV; i--;)
    {
        const auto deg = _hyperGraph._FRsets[i].size();
        coverage[i] = deg;

        if (deg > maxDeg) maxDeg = deg;
    }

    // degMap: map degree to the nodes with this degree
    RRsets degMap(maxDeg + 1);

    for (auto i = _numV; i--;)
    {
        if (coverage[i] == 0) continue;

        degMap[coverage[i]].push_back(i);
    }

    size_t sumInf = 0;
    // check if an edge is removed
    std::vector<bool> edgeMark(_numRRsets, false);
    _vecSeed.clear();

    for (auto deg = maxDeg; deg > 0; deg--) // Enusre deg > 0
    {
        auto& vecNode = degMap[deg];

        for (auto idx = vecNode.size(); idx--;)
        {
            auto argmaxIdx = vecNode[idx];
            const auto currDeg = coverage[argmaxIdx];

            if (deg > currDeg)
            {
                degMap[currDeg].push_back(argmaxIdx);
                continue;
            }

            if (true)
            {
                // Find upper bound
                auto topk = targetSize;
                auto degBound = deg;
                FRset vecBound(targetSize);
                // Initialize vecBound
                auto idxBound = idx + 1;

                while (topk && idxBound--)
                {
                    vecBound[--topk] = coverage[degMap[degBound][idxBound]];
                }

                while (topk && --degBound)
                {
                    idxBound = degMap[degBound].size();

                    while (topk && idxBound--)
                    {
                        vecBound[--topk] = coverage[degMap[degBound][idxBound]];
                    }
                }

                MakeMinHeap(vecBound);
                // Find the top-k marginal coverage
                auto flag = topk == 0;

                while (flag && idxBound--)
                {
                    const auto currDegBound = coverage[degMap[degBound][idxBound]];

                    if (vecBound[0] >= degBound)
                    {
                        flag = false;
                    }
                    else if (vecBound[0] < currDegBound)
                    {
                        MinHeapReplaceMinValue(vecBound, currDegBound);
                    }
                }

                while (flag && --degBound)
                {
                    idxBound = degMap[degBound].size();

                    while (flag && idxBound--)
                    {
                        const auto currDegBound = coverage[degMap[degBound][idxBound]];

                        if (vecBound[0] >= degBound)
                        {
                            flag = false;
                        }
                        else if (vecBound[0] < currDegBound)
                        {
                            MinHeapReplaceMinValue(vecBound, currDegBound);
                        }
                    }
                }

                _boundLast = double(accumulate(vecBound.begin(), vecBound.end(), size_t(0)) + sumInf) * _numV / _numRRsets;

                if (_boundMin > _boundLast) _boundMin = _boundLast;
            }

            if (_vecSeed.size() >= targetSize)
            {
                // Top-k influential nodes constructed
                const auto finalInf = 1.0 * sumInf * _numV / _numRRsets;
                // std::cout << ">>>[greedy-lazy] influence: " << finalInf << ", min-bound: " << _boundMin <<
                //           ", last-bound: " << _boundLast << '\n';
                return finalInf;
            }

            sumInf += currDeg;
            _vecSeed.push_back(argmaxIdx);
            coverage[argmaxIdx] = 0;

            for (auto edgeIdx : _hyperGraph._FRsets[argmaxIdx])
            {
                if (edgeMark[edgeIdx]) continue;

                edgeMark[edgeIdx] = true;

                for (auto nodeIdx : _hyperGraph._RRsets[edgeIdx])
                {
                    if (coverage[nodeIdx] == 0) continue; // This node is seed, skip

                    coverage[nodeIdx]--;
                }
            }
        }

        degMap.pop_back();
    }

    return 1.0 * _numV; // All RR sets are covered.
}

// In the case of identical marginal，the node with largest out degree is chosen
double Alg::MaxCoverOutDegPrority(const int targetSize)
{
    _boundLast = DBL_MAX, _boundMin = DBL_MAX;
    FRset coverage(_numV, 0);
    size_t maxDeg = 0;

    for (auto i = _numV; i--;)
    {
        const auto deg = _hyperGraph._FRsets[i].size();
        coverage[i] = deg;

        if (deg > maxDeg) maxDeg = deg;
    }

    RRsets degMap(maxDeg + 1); // degMap: map degree to the nodes with this degree

    for (auto i = _numV; i--;)
    {
        if (coverage[i] == 0) continue;

        degMap[coverage[i]].push_back(i);
    }

    size_t sumInf = 0;
    // check if an edge is removed
    std::vector<bool> edgeMark(_numRRsets, false);
    _vecSeed.clear();

    for (auto deg = maxDeg; deg > 0; deg--) // Enusre deg > 0
    {
        auto& origVecNode = degMap[deg];
        std::vector<std::pair<uint32_t, uint32_t>> vecPair;

        for (auto idx = origVecNode.size(); idx--;)
        {
            auto node = origVecNode[idx];
            auto nodeCoverage = coverage[node];

            if (deg > nodeCoverage)
            {
                degMap[nodeCoverage].push_back(node);
                continue;
            }

            vecPair.push_back(std::make_pair(_vecOutDegree[node], node));
        }

        degMap.pop_back();

        if (vecPair.size() == 0)
        {
            continue;
        }

        /* sort nodes by their out-degre in ascending order */
        sort(vecPair.begin(), vecPair.end());
        std::vector<uint32_t> newVecNode;

        for (auto &nodePair : vecPair)
        {
            newVecNode.push_back(nodePair.second);
        }

        degMap.push_back(newVecNode);
        auto &vecNode = degMap[deg];

        for (auto idx = vecNode.size(); idx--;)
        {
            auto argmaxIdx = vecNode[idx];
            const auto currDeg = coverage[argmaxIdx];

            if (deg > currDeg)
            {
                degMap[currDeg].push_back(argmaxIdx);
                continue;
            }

            if (true)
            {
                // Find upper bound
                auto topk = targetSize;
                auto degBound = deg;
                FRset vecBound(targetSize);
                // Initialize vecBound
                auto idxBound = idx + 1;

                while (topk && idxBound--)
                {
                    vecBound[--topk] = coverage[degMap[degBound][idxBound]];
                }

                while (topk && --degBound)
                {
                    idxBound = degMap[degBound].size();

                    while (topk && idxBound--)
                    {
                        vecBound[--topk] = coverage[degMap[degBound][idxBound]];
                    }
                }

                MakeMinHeap(vecBound);
                auto flag = topk == 0;

                while (flag && idxBound--)
                {
                    const auto currDegBound = coverage[degMap[degBound][idxBound]];

                    if (vecBound[0] >= degBound)
                    {
                        flag = false;
                    }
                    else if (vecBound[0] < currDegBound)
                    {
                        MinHeapReplaceMinValue(vecBound, currDegBound);
                    }
                }

                while (flag && --degBound)
                {
                    idxBound = degMap[degBound].size();

                    while (flag && idxBound--)
                    {
                        const auto currDegBound = coverage[degMap[degBound][idxBound]];

                        if (vecBound[0] >= degBound)
                        {
                            flag = false;
                        }
                        else if (vecBound[0] < currDegBound)
                        {
                            MinHeapReplaceMinValue(vecBound, currDegBound);
                        }
                    }
                }

                _boundLast = double(accumulate(vecBound.begin(), vecBound.end(), size_t(0)) + sumInf) * _numV / _numRRsets;

                if (_boundMin > _boundLast) _boundMin = _boundLast;
            }

            if (_vecSeed.size() >= targetSize)
            {
                const auto finalInf = 1.0 * sumInf * _numV / _numRRsets;
                // std::cout << ">>>[greedy-lazy] influence: " << finalInf << ", min-bound: " << _boundMin <<
                //           ", last-bound: " << _boundLast << '\n';
                return finalInf;
            }

            sumInf += currDeg;
            _vecSeed.push_back(argmaxIdx);
            coverage[argmaxIdx] = 0;

            for (auto edgeIdx : _hyperGraph._FRsets[argmaxIdx])
            {
                if (edgeMark[edgeIdx]) continue;

                edgeMark[edgeIdx] = true;

                for (auto nodeIdx : _hyperGraph._RRsets[edgeIdx])
                {
                    if (coverage[nodeIdx] == 0) continue; // This node is seed, skip

                    coverage[nodeIdx]--;
                }
            }
        }

        degMap.pop_back();
    }

    return 1.0 * _numV; // All RR sets are covered.
}

// max cover used in the IM-sentinel Phase
double Alg::MaxCoverIMSentinel(std::vector<uint32_t> &seedSet, const int targetSize)
{
    // seedSet: the sentinel set obtained in the Sentinel Set Selection Phase
    std::unordered_set<uint32_t> subSeedSet(seedSet.begin(), seedSet.end());
    _boundLast = DBL_MAX, _boundMin = DBL_MAX;
    FRset coverage(_numV, 0);
    size_t maxDeg = 0;

    for (auto i = _numV; i--;)
    {
        const auto deg = _hyperGraph._FRsets[i].size();
        coverage[i] = deg;

        if (deg > maxDeg) maxDeg = deg;
    }

    RRsets degMap(maxDeg + 1); // degMap: map degree to the nodes with this degree

    for (auto i = _numV; i--;)
    {
        if (coverage[i] == 0) continue;

        degMap[coverage[i]].push_back(i);
    }

    size_t sumInf = 0;
    // check if an edge is removed
    std::vector<bool> edgeMark(_numRRsets, false);
    _vecSeed.clear();

    for (auto node : seedSet)
    {
        sumInf += coverage[node];
        _vecSeed.push_back(node);
        coverage[node] = 0;

        for (auto edgeIdx : _hyperGraph._FRsets[node])
        {
            if (edgeMark[edgeIdx]) continue;

            edgeMark[edgeIdx] = true;

            for (auto nodeIdx : _hyperGraph._RRsets[edgeIdx])
            {
                if (coverage[nodeIdx] == 0) continue; // This node is seed, skip

                coverage[nodeIdx]--;
            }
        }
    }

    for (auto deg = maxDeg; deg > 0; deg--) // Enusre deg > 0
    {
        auto& origVecNode = degMap[deg];
        std::vector<std::pair<uint32_t, uint32_t>> vecPair;

        for (auto idx = origVecNode.size(); idx--;)
        {
            auto node = origVecNode[idx];
            auto nodeCoverage = coverage[node];

            if (deg > nodeCoverage)
            {
                degMap[nodeCoverage].push_back(node);
                continue;
            }

            vecPair.push_back(std::make_pair(_vecOutDegree[node], node));
        }

        degMap.pop_back();

        if (vecPair.size() == 0)
        {
            continue;
        }

        /* sort nodes by their out-degre in ascending order */
        sort(vecPair.begin(), vecPair.end());
        std::vector<uint32_t> newVecNode;

        for (auto &nodePair : vecPair)
        {
            newVecNode.push_back(nodePair.second);
        }

        degMap.push_back(newVecNode);
        auto &vecNode = degMap[deg];

        for (auto idx = vecNode.size(); idx--;)
        {
            auto argmaxIdx = vecNode[idx];
            const auto currDeg = coverage[argmaxIdx];

            if (deg > currDeg)
            {
                degMap[currDeg].push_back(argmaxIdx);
                continue;
            }

            if (true)
            {
                // Find upper bound
                auto topk = targetSize;
                auto degBound = deg;
                FRset vecBound(targetSize);
                // Initialize vecBound
                auto idxBound = idx + 1;

                while (topk && idxBound--)
                {
                    vecBound[--topk] = coverage[degMap[degBound][idxBound]];
                }

                while (topk && --degBound)
                {
                    idxBound = degMap[degBound].size();

                    while (topk && idxBound--)
                    {
                        vecBound[--topk] = coverage[degMap[degBound][idxBound]];
                    }
                }

                MakeMinHeap(vecBound);
                auto flag = topk == 0;

                while (flag && idxBound--)
                {
                    const auto currDegBound = coverage[degMap[degBound][idxBound]];

                    if (vecBound[0] >= degBound)
                    {
                        flag = false;
                    }
                    else if (vecBound[0] < currDegBound)
                    {
                        MinHeapReplaceMinValue(vecBound, currDegBound);
                    }
                }

                while (flag && --degBound)
                {
                    idxBound = degMap[degBound].size();

                    while (flag && idxBound--)
                    {
                        const auto currDegBound = coverage[degMap[degBound][idxBound]];

                        if (vecBound[0] >= degBound)
                        {
                            flag = false;
                        }
                        else if (vecBound[0] < currDegBound)
                        {
                            MinHeapReplaceMinValue(vecBound, currDegBound);
                        }
                    }
                }

                _boundLast = double(accumulate(vecBound.begin(), vecBound.end(), size_t(0)) + sumInf) * _numV / _numRRsets;

                if (_boundMin > _boundLast) _boundMin = _boundLast;
            }

            if (_vecSeed.size() >= targetSize)
            {
                const auto finalInf = 1.0 * sumInf * _numV / _numRRsets;
                // std::cout << ">>>[greedy-lazy] influence: " << finalInf << ", min-bound: " << _boundMin <<
                //           ", last-bound: " << _boundLast << '\n';
                return finalInf;
            }

            sumInf += currDeg;
            _vecSeed.push_back(argmaxIdx);
            coverage[argmaxIdx] = 0;

            for (auto edgeIdx : _hyperGraph._FRsets[argmaxIdx])
            {
                if (edgeMark[edgeIdx]) continue;

                edgeMark[edgeIdx] = true;

                for (auto nodeIdx : _hyperGraph._RRsets[edgeIdx])
                {
                    if (coverage[nodeIdx] == 0) continue; // This node is seed, skip

                    coverage[nodeIdx]--;
                }
            }
        }

        degMap.pop_back();
    }

    return 1.0 * _numV; // All RR sets are covered.
}

double Alg::MaxCoverSentinelSet(const int targetSize, const int totalTargetSize)
{
    //targetSize: the size of the sentinel set
    //totalTargetSize: the total number of the seed set
    _boundLast = DBL_MAX, _boundMin = DBL_MAX;
    FRset coverage(_numV, 0);
    size_t maxDeg = 0;

    for (auto i = _numV; i--;)
    {
        const auto deg = _hyperGraph._FRsets[i].size();
        coverage[i] = deg;

        if (deg > maxDeg) maxDeg = deg;
    }

    RRsets degMap(maxDeg + 1); // degMap: map degree to the nodes with this degree

    for (auto i = _numV; i--;)
    {
        if (coverage[i] == 0) continue;

        degMap[coverage[i]].push_back(i);
    }

    size_t sumInf = 0;
    // check if an edge is removed
    std::vector<bool> edgeMark(_numRRsets, false);
    _vecSeed.clear();

    for (auto deg = maxDeg; deg > 0; deg--) // Enusre deg > 0
    {
        auto& origVecNode = degMap[deg];
        std::vector<std::pair<uint32_t, uint32_t>> vecPair;

        for (auto idx = origVecNode.size(); idx--;)
        {
            auto node = origVecNode[idx];
            auto nodeCoverage = coverage[node];

            if (deg > nodeCoverage)
            {
                degMap[nodeCoverage].push_back(node);
                continue;
            }

            vecPair.push_back(std::make_pair(_vecOutDegree[node], node));
        }

        degMap.pop_back();

        if (vecPair.size() == 0)
        {
            continue;
        }

        /* sort nodes by their out-degre in ascending order */
        sort(vecPair.begin(), vecPair.end());
        std::vector<uint32_t> newVecNode;

        for (auto &nodePair : vecPair)
        {
            newVecNode.push_back(nodePair.second);
        }

        degMap.push_back(newVecNode);
        auto &vecNode = degMap[deg];

        for (auto idx = vecNode.size(); idx--;)
        {
            auto argmaxIdx = vecNode[idx];
            const auto currDeg = coverage[argmaxIdx];

            if (deg > currDeg)
            {
                degMap[currDeg].push_back(argmaxIdx);
                continue;
            }

            if (true)
            {
                // Find upper bound
                auto topk = totalTargetSize;
                auto degBound = deg;
                FRset vecBound(totalTargetSize);
                // Initialize vecBound
                auto idxBound = idx + 1;

                while (topk && idxBound--)
                {
                    vecBound[--topk] = coverage[degMap[degBound][idxBound]];
                }

                while (topk && --degBound)
                {
                    idxBound = degMap[degBound].size();

                    while (topk && idxBound--)
                    {
                        vecBound[--topk] = coverage[degMap[degBound][idxBound]];
                    }
                }

                MakeMinHeap(vecBound);
                // Find the top-k marginal coverage
                auto flag = topk == 0;

                while (flag && idxBound--)
                {
                    const auto currDegBound = coverage[degMap[degBound][idxBound]];

                    if (vecBound[0] >= degBound)
                    {
                        flag = false;
                    }
                    else if (vecBound[0] < currDegBound)
                    {
                        MinHeapReplaceMinValue(vecBound, currDegBound);
                    }
                }

                while (flag && --degBound)
                {
                    idxBound = degMap[degBound].size();

                    while (flag && idxBound--)
                    {
                        const auto currDegBound = coverage[degMap[degBound][idxBound]];

                        if (vecBound[0] >= degBound)
                        {
                            flag = false;
                        }
                        else if (vecBound[0] < currDegBound)
                        {
                            MinHeapReplaceMinValue(vecBound, currDegBound);
                        }
                    }
                }

                _boundLast = double(accumulate(vecBound.begin(), vecBound.end(), size_t(0)) + sumInf) * _numV / _numRRsets;

                if (_boundMin > _boundLast) _boundMin = _boundLast;
            }

            if (_vecSeed.size() >= targetSize)
            {
                goto afterGreedy;
            }

            sumInf += currDeg;
            _vecSeed.push_back(argmaxIdx);
            _vecVldtInf.push_back(sumInf);
            coverage[argmaxIdx] = 0;

            for (auto edgeIdx : _hyperGraph._FRsets[argmaxIdx])
            {
                if (edgeMark[edgeIdx]) continue;

                edgeMark[edgeIdx] = true;

                for (auto nodeIdx : _hyperGraph._RRsets[edgeIdx])
                {
                    if (coverage[nodeIdx] == 0) continue; // This node is seed, skip

                    coverage[nodeIdx]--;
                }
            }
        }

        degMap.pop_back();
    }

afterGreedy:
    const auto finalInf = 1.0 * sumInf * _numV / _numRRsets;
    // std::cout << "  >>>[greedy-lazy] influence: " << finalInf << ", seed set: " << _vecSeed.size() << ", min-bound: " << _boundMin <<
    //           ", last-bound: " << _boundLast << std::endl;

    if (_vecSeed.size() == targetSize)
    {
        // if the sample size is sufficiently large, the result is reliable
        if (_numRRsets > 1000)
        {
            return finalInf;
        }
    }

    // if covering all the RR-sets, the sentinel set may include some nodes which cover only a small number of RR-sets.
    // such nodes should not be included.
    uint32_t threshold = 0.9 * sumInf;
    for (int i = _vecSeed.size() - 1; i > 0; i--)
    {
        if (_vecVldtInf[i - 1] >= threshold)
        {
            degMap[0].push_back(_vecSeed[i]);
            _vecSeed.pop_back();
        }
        else
        {
            break;
        }
    }

    // std::cout << "seedset size reaching 0.9 coverage: " << _vecSeed.size() << std::endl;
    // the following code is to select the nodes with large out-degree. 
    int seedSetSize = (degMap[0].size() > targetSize) ? targetSize : degMap[0].size();
    std::vector<std::pair<uint32_t, uint32_t>> vecHeap;

    for (int i = 0; i < seedSetSize; i++)
    {
        auto &node = degMap[0][i];
        vecHeap.push_back(std::make_pair(_vecOutDegree[node], node));
    }

    std::make_heap(vecHeap.begin(), vecHeap.end(), GreaterPair);
    const auto nodeNum = degMap[0].size();


    for (int i = seedSetSize; i < nodeNum; i++)
    {
        uint32_t node = degMap[0][i];
        uint32_t currDeg = _vecOutDegree[node];

        if (currDeg > vecHeap[0].first)
        {
            std::pop_heap(vecHeap.begin(), vecHeap.end());
            vecHeap.pop_back();
            vecHeap.push_back(std::make_pair(_vecOutDegree[node], node));
            std::push_heap(vecHeap.begin(), vecHeap.end());
        }
    }

    std::sort_heap(vecHeap.begin(), vecHeap.end(), GreaterPair);
    std::unordered_set<uint32_t> seedHashSet(_vecSeed.begin(), _vecSeed.end());

    for (auto &node : vecHeap)
    {
        if (_vecSeed.size() >= targetSize)
        {
            break;
        }

        if (seedHashSet.find(node.second) != seedHashSet.end())
        {
            std::cout << "node exist" << std::endl;
            continue;
        }

        _vecSeed.push_back(node.second);
        seedHashSet.insert(node.second);

        if (_vecVldtInf.size() < _vecSeed.size())
        {
            _vecVldtInf.push_back(sumInf);
        }
    }

    return 1.0 * _numV; // All RR sets are covered.
}

double Alg::MaxCoverTopK(const int targetSize)
{
    FRset coverage(_numV, 0);
    size_t maxDeg = 0;

    for (auto i = _numV; i--;)
    {
        const auto deg = _hyperGraph._FRsets[i].size();
        coverage[i] = deg;

        if (deg > maxDeg) maxDeg = deg;
    }

    RRsets degMap(maxDeg + 1); // degMap: map degree to the nodes with this degree

    for (auto i = _numV; i--;)
    {
        //if (coverage[i] == 0) continue;
        degMap[coverage[i]].push_back(i);
    }

    Nodelist sortedNode(_numV); // sortedNode: record the sorted nodes in ascending order of degree
    Nodelist nodePosition(_numV); // nodePosition: record the position of each node in the sortedNode
    Nodelist degreePosition(maxDeg + 2); // degreePosition: the start position of each degree in sortedNode
    uint32_t idxSort = 0;
    size_t idxDegree = 0;

    for (auto& nodes : degMap)
    {
        degreePosition[idxDegree + 1] = degreePosition[idxDegree] + (uint32_t)nodes.size();
        idxDegree++;

        for (auto& node : nodes)
        {
            nodePosition[node] = idxSort;
            sortedNode[idxSort++] = node;
        }
    }

    // check if an edge is removed
    std::vector<bool> edgeMark(_numRRsets, false);
    // record the total of top-k marginal gains
    size_t sumTopk = 0;

    for (auto deg = maxDeg + 1; deg--;)
    {
        if (degreePosition[deg] <= _numV - targetSize)
        {
            sumTopk += deg * (degreePosition[deg + 1] - (_numV - targetSize));
            break;
        }

        sumTopk += deg * (degreePosition[deg + 1] - degreePosition[deg]);
    }

    _boundMin = 1.0 * sumTopk;
    _vecSeed.clear();
    size_t sumInf = 0;

    /*
    * sortedNode: position -> node
    * nodePosition: node -> position
    * degreePosition: degree -> position (start position of this degree)
    * coverage: node -> degree
    * e.g., swap the position of a node with the start position of its degree
    * swap(sortedNode[nodePosition[node]], sortedNode[degreePosition[coverage[node]]])
    */
    for (auto k = targetSize; k--;)
    {
        const auto seed = sortedNode.back();
        sortedNode.pop_back();
        const auto newNumV = sortedNode.size();
        sumTopk += coverage[sortedNode[newNumV - targetSize]] - coverage[seed];
        sumInf += coverage[seed];
        _vecSeed.push_back(seed);
        coverage[seed] = 0;

        for (auto edgeIdx : _hyperGraph._FRsets[seed])
        {
            if (edgeMark[edgeIdx]) continue;

            edgeMark[edgeIdx] = true;

            for (auto nodeIdx : _hyperGraph._RRsets[edgeIdx])
            {
                if (coverage[nodeIdx] == 0) continue; // This node is seed, skip

                const auto currPos = nodePosition[nodeIdx]; // The current position
                const auto currDeg = coverage[nodeIdx]; // The current degree
                const auto startPos = degreePosition[currDeg]; // The start position of this degree
                const auto startNode = sortedNode[startPos]; // The node with the start position
                // Swap this node to the start position with the same degree, and update their positions in nodePosition
                std::swap(sortedNode[currPos], sortedNode[startPos]);
                nodePosition[nodeIdx] = startPos;
                nodePosition[startNode] = currPos;
                // Increase the start position of this degree by 1, and decrease the degree of this node by 1
                degreePosition[currDeg]++;
                coverage[nodeIdx]--;

                // If the start position of this degree is in top-k, reduce topk by 1
                if (startPos >= newNumV - targetSize) sumTopk--;
            }
        }

        _boundLast = 1.0 * (sumInf + sumTopk);

        if (_boundMin > _boundLast) _boundMin = _boundLast;
    }

    _boundMin *= 1.0 * _numV / _numRRsets;
    _boundLast *= 1.0 * _numV / _numRRsets;
    const auto finalInf = 1.0 * sumInf * _numV / _numRRsets;
    std::cout << "  >>>[greedy-topk] influence: " << finalInf << ", min-bound: " << _boundMin <<
              ", last-bound: " << _boundLast << '\n';
    return finalInf;
}

double Alg::MaxCover(const int targetSize)
{
    if (targetSize >= 1000) return MaxCoverTopK(targetSize);
    return MaxCoverVanilla(targetSize);
}

double Alg::ComputeFairMarginalGain(uint32_t u, bool usingVldtGraph)
{
    double fInf;
    if (usingVldtGraph)
    {
        fInf = _hyperGraphVldt.ComputeMarginalFairInf(u);
    }
    else
    {
        fInf = _hyperGraph.ComputeMarginalFairInf(u);
    }
    return fInf;
}

void Alg::UpdateCoverageStatus(uint32_t u, std::vector<bool>& needRecompute)
{
    _hyperGraph.UpdateCoverageStatus(u, needRecompute);
    return ;
}

double Alg::MaxCoverGreedyMarginalFairSentinelSet(const int targetSize, const int totalTargetSize)
{
    // Greedy algorithm with lazy evaluation
    _boundLast = DBL_MAX, _boundMin = DBL_MAX, _fairInf = 0.0;
    std::vector<double> marginalFair(_numV, 0.0); // the marginal fair gain of each node
    // std::vector<bool> edgeMark(_numRRsets, false); // check if an edge is removed
    std::vector<bool> selected(_numV, false); // check if a node is selected
    std::vector<bool> needRecompute(_numV, false); // the nodes need to be recomputed
    double updatedGain, generalMCFInf;

    _vecSeed.clear();
    _hyperGraph.InitGcoverage();

    for (int u = 0; u < _numV; ++u) 
    {
        marginalFair[u] = ComputeFairMarginalGain(u);
    }
    using GainNodePair = std::pair<double, uint32_t>;
    auto cmp = [](const GainNodePair& a, const GainNodePair& b) { return a.first < b.first; };
    std::priority_queue<GainNodePair, std::vector<GainNodePair>, decltype(cmp)> maxHeap(cmp);
    std::vector<GainNodePair> vecLargestMCFInf;

    for (int u = 0; u < _numV; ++u)
    {
        maxHeap.emplace(marginalFair[u], u);
    }

    for (int iter = 0; iter < targetSize && !maxHeap.empty(); ++iter)
    {
        while (!maxHeap.empty())
        {
            auto [estGain, u] = maxHeap.top(); maxHeap.pop();

            if (selected[u]) continue;

            // Recalculate the true marginal fairness
            if (needRecompute[u])
            {
                updatedGain = ComputeFairMarginalGain(u);
                needRecompute[u] = false;
            }
            else
            {
                updatedGain = estGain;
            }

            // If the valuation is lower than the current heap, put it back to the heap (lazy update)
            if (!maxHeap.empty() && updatedGain < maxHeap.top().first)
            {
                maxHeap.emplace(updatedGain, u);
                continue;
            }

            // Select the seed u
            selected[u] = true;
            _vecSeed.push_back(u);
            _fairInf += updatedGain;
            UpdateCoverageStatus(u, needRecompute);
            _vecSelfVldtGFCoverage.push_back(_fairInf * _numRRsets); //_vecSelfVldtGFCoverage is the general coverage
            // _vecSelfVldtGFCoverage.push_back(_fairInf);
            for (auto q = 0; q < _numMetric; q++) 
            {
                _vecFairInfDiffSize[q][iter] = _vecFairInf[q];
                _vecFairInfVldt[q] = _vecFairInf[q];
            }
            // Calculate the tighter upper bound
            generalMCFInf = 0.0;
            if (_hyperGraph.get_func_type() == FIM) // tighter upper bound for FIM
            {
                for (int i = 0; i < totalTargetSize; ++i)
                {
                    while (!maxHeap.empty())
                    {
                        auto [estGain, u] = maxHeap.top(); maxHeap.pop();
                        if (selected[u]) continue;
                        // Recalculate the true marginal fairness
                        if (needRecompute[u])
                        {
                            updatedGain = ComputeFairMarginalGain(u);
                            needRecompute[u] = false;
                        }
                        else
                        {
                            updatedGain = estGain;
                        }
                        // If the valuation is lower than the current heap, put it back to the heap
                        if (!maxHeap.empty() && updatedGain < maxHeap.top().first)
                        {
                            maxHeap.emplace(updatedGain, u);
                            continue;
                        }
                        generalMCFInf += updatedGain;
                        vecLargestMCFInf.push_back(std::make_pair(updatedGain, u));
                        break;
                    }
                }
                _boundLast = _fairInf + generalMCFInf;
                if (_boundMin > _boundLast) _boundMin = _boundLast;
                // Push back the temporarily removed top-k node candidates
                for (const auto& gainPair : vecLargestMCFInf) maxHeap.emplace(gainPair);
                vecLargestMCFInf.clear();
            }
            else // estimate the upper bound for RFIM
            {
                for (int i = 0; i < _k; ++i)
                {
                    while (!maxHeap.empty())
                    {
                        auto [estGain, u] = maxHeap.top(); maxHeap.pop();
                        if (selected[u]) continue;
                        // Recalculate the true marginal fairness
                        if (needRecompute[u])
                        {
                            updatedGain = ComputeFairMarginalGain(u);
                            needRecompute[u] = false;
                        }
                        else
                        {
                            updatedGain = estGain;
                        }
                        // If the valuation is lower than the current heap, put it back to the heap
                        if (!maxHeap.empty() && updatedGain < maxHeap.top().first)
                        {
                            maxHeap.emplace(updatedGain, u);
                            continue;
                        }
                        generalMCFInf += updatedGain;
                        vecLargestMCFInf.push_back(std::make_pair(updatedGain, u));
                        break;
                    }
                }
                _boundLast = _fairInf + generalMCFInf;
                if (_boundMin > _boundLast)
                {
                    _boundMin = _boundLast;
                    for (auto q = 0; q < _numMetric; q++)
                    {
                        _vecBoundMinFairInf[q] = _vecFairInf[q];
                        _hyperGraph.set_func_type(FIM, q);
                        for (const auto& gainPair : vecLargestMCFInf)
                        {
                            uint32_t u = gainPair.second;
                            double marginalFairGain = ComputeFairMarginalGain(u);
                            _vecBoundMinFairInf[q] += marginalFairGain;
                        }
                    }
                    _hyperGraph.set_func_type(RFIM);
                }
                // Push back the temporarily removed top-k node candidates
                for (const auto& gainPair : vecLargestMCFInf) maxHeap.emplace(gainPair);
                vecLargestMCFInf.clear();
            }
            break;
        }
    }
    return _fairInf;
}

double Alg::MaxCoverGreedyMarginalFairIMSentinel(std::vector<uint32_t> &seedSet, const int targetSize)
{
    std::unordered_set<uint32_t> subSeedSet(seedSet.begin(), seedSet.end());
    _boundLast = DBL_MAX, _boundMin = DBL_MAX, _fairInf = 0.0;
    std::vector<double> marginalFair(_numV, 0.0);
    std::vector<bool> selected(_numV, false);
    std::vector<bool> needRecompute(_numV, false);
    double updatedGain;

    _vecSeed.clear();
    _hyperGraph.InitGcoverage();

    for (auto u : seedSet)
    {
        _vecSeed.push_back(u);
        selected[u] = true;
        _fairInf += ComputeFairMarginalGain(u);
        UpdateCoverageStatus(u, needRecompute);
    }

    for (uint32_t u = 0; u < _numV; ++u)
    {
        if (!selected[u]) marginalFair[u] = ComputeFairMarginalGain(u);
    }

    using GainNodePair = std::pair<double, uint32_t>;
    auto cmp = [](const GainNodePair &a, const GainNodePair &b) { return a.first < b.first; };
    std::priority_queue<GainNodePair, std::vector<GainNodePair>, decltype(cmp)> maxHeap(cmp);

    for (int u = 0; u < _numV; ++u)
    {
        if (!selected[u]) maxHeap.emplace(marginalFair[u], u);
    }

    for (int iter = seedSet.size(); iter < targetSize && !maxHeap.empty(); ++iter)
    {
        while (!maxHeap.empty())
        {
            auto [estGain, u] = maxHeap.top();
            maxHeap.pop();
            if (selected[u]) continue;

            // Lazy update
            if (needRecompute[u])
            {
                updatedGain = ComputeFairMarginalGain(u);
                needRecompute[u] = false;
            }
            else
            {
                updatedGain = estGain;
            }

            if (!maxHeap.empty() && updatedGain < maxHeap.top().first)
            {
                maxHeap.emplace(updatedGain, u);
                continue;
            }

            // select a new seed
            selected[u] = true;
            _vecSeed.push_back(u);
            _fairInf += updatedGain;
            UpdateCoverageStatus(u, needRecompute);

            // update the estimation of upper bound
            std::vector<GainNodePair> tempTopK;
            double generalMCFInf = 0.0;

            if (_hyperGraph.get_func_type() == FIM)
            {
                for (int i = 0; i < targetSize; ++i)
                {
                    while (!maxHeap.empty())
                    {
                        auto [gain, cand] = maxHeap.top();
                        maxHeap.pop();
                        if (selected[cand]) continue;

                        if (needRecompute[cand])
                        {
                            updatedGain = ComputeFairMarginalGain(cand);
                            needRecompute[cand] = false;
                        }
                        else
                        {
                            updatedGain = gain;
                        }

                        if (!maxHeap.empty() && updatedGain < maxHeap.top().first)
                        {
                            maxHeap.emplace(updatedGain, cand);
                            continue;
                        }

                        generalMCFInf += updatedGain;
                        tempTopK.emplace_back(updatedGain, cand);
                        break;
                    }
                }
                _boundLast = _fairInf + generalMCFInf;
                if (_boundMin > _boundLast) _boundMin = _boundLast;

                for (const auto &[gain, node] : tempTopK) maxHeap.emplace(gain, node);
            }
            else //RFIM
            {
                for (int i = 0; i < _k; ++i)
                {
                    while (!maxHeap.empty())
                    {
                        auto [estGain, u] = maxHeap.top(); maxHeap.pop();
                        if (selected[u]) continue;
                        // Recalculate the true marginal fairness
                        if (needRecompute[u])
                        {
                            updatedGain = ComputeFairMarginalGain(u);
                            needRecompute[u] = false;
                        }
                        else
                        {
                            updatedGain = estGain;
                        }
                        // If the valuation is lower than the current heap, put it back to the heap
                        if (!maxHeap.empty() && updatedGain < maxHeap.top().first)
                        {
                            maxHeap.emplace(updatedGain, u);
                            continue;
                        }
                        generalMCFInf += updatedGain;
                        tempTopK.push_back(std::make_pair(updatedGain, u));
                        break;
                    }
                }
                _boundLast = _fairInf + generalMCFInf;
                if (_boundMin > _boundLast)
                {
                    _boundMin = _boundLast;
                    for (auto q = 0; q < _numMetric; q++)
                    {
                        _vecBoundMinFairInf[q] = _vecFairInf[q];
                        _hyperGraph.set_func_type(FIM, q);
                        for (const auto& gainPair : tempTopK)
                        {
                            uint32_t u = gainPair.second;
                            double marginalFairGain = ComputeFairMarginalGain(u);
                            _vecBoundMinFairInf[q] += marginalFairGain;
                        }
                    }
                    _hyperGraph.set_func_type(RFIM);
                }
                for (const auto& gainPair : tempTopK) maxHeap.emplace(gainPair);
                tempTopK.clear();
            }
            break;
        }
    }

    return _fairInf;
}

double Alg::AdjustGeneralBound(const double bound, const size_t numR)
{
    double adjustedFair = bound / numR;
    double consLB = _hyperGraph.ComputeConsLB(_vecSeed);
    double consUB = _hyperGraph.ComputeConsUB(numR);
    adjustedFair = (adjustedFair < consLB) ? consLB : adjustedFair;
    adjustedFair = (adjustedFair > consUB) ? consUB : adjustedFair;
    return adjustedFair * numR;
}

double Alg::AdjustGeneralInfluenceBound(const double gInfBound, const size_t numR)
{
    return AdjustGeneralBound(gInfBound * numR, numR) / numR;
}

void Alg::set_prob_dist(const ProbDist dist)
{
    _probDist = dist;
    _hyperGraph.set_prob_dist(dist);
    _hyperGraphVldt.set_prob_dist(dist);
}

void Alg::set_vanilla_sample(const bool isVanilla)
{
    if (isVanilla)
    {
        std::cout << "Vanilla sampling method is used" << std::endl;
    }

    _hyperGraph.set_vanilla_sample(isVanilla);
    _hyperGraphVldt.set_vanilla_sample(isVanilla);
}

double Alg::EfficInfVldtAlg()
{
    return EfficInfVldtAlg(_vecSeed);
}

double Alg::EfficInfVldtAlg(const Nodelist vecSeed)
{
    Timer EvalTimer("Inf. Eval.");
    std::cout << "  >>>Evaluating influence in [0.99,1.01]*EPT with prob. 99.9% ...\n";
    const auto inf = _hyperGraph.EfficInfVldtAlg(vecSeed);
    //const auto inf = _hyperGraphVldt.EfficInfVldtAlg(vecSeed);
    std::cout << "  >>>Down! influence: " << inf << ", time used (sec): " << EvalTimer.get_total_time() << '\n';
    return inf;
}

/// Evaluate the general fair influence of the seed set constructed
double Alg::EfficFairInfVldtAlg()
{
    return EfficFairInfVldtAlg(_vecSeed);
}

/// Evaluate the general fair influence of a given seed set
double Alg::EfficFairInfVldtAlg(const Nodelist vecSeed)
{
    return 0.0;
}

double Alg::estimateRRSize()
{
    const int sampleNum = 100;
    _hyperGraph.BuildRRsets(sampleNum);
    double avg= _hyperGraph.HyperedgeAvg();
    _hyperGraph.RefreshHypergraph();
    return avg;
}

double Alg::subsimOnly(const int targetSize, const double epsilon, const double delta)
{
    Timer timerSubsim("SUBSIM");
    const double e = exp(1);
    const double approx = 1 - 1.0 / e;
    const double alpha = sqrt(log(6.0 / delta));
    const double beta = sqrt((1 - 1 / e) * (logcnk(_numV, targetSize) + log(6.0 / delta)));
    const auto numRbase = size_t(2.0 * pow2((1 - 1 / e) * alpha + beta));
    const auto maxNumR = size_t(2.0 * _numV * pow2((1 - 1 / e) * alpha + beta) / targetSize / pow2(epsilon)) + 1;
    const auto numIter = (size_t)log2(maxNumR / numRbase) + 1;
    const double a1 = log(numIter * 3.0 / delta);
    const double a2 = log(numIter * 3.0 / delta);
    double time1 = 0.0, time2 = 0.0, time3 = 0.0;

    std::cout << std::endl;
    for (auto idx = 1; idx <= numIter; idx++)
    {
        const auto numR = numRbase << (idx-1);
        std::cout << "Iteration: " << idx << " RR set: " << numR << std::endl;
        timerSubsim.get_operation_time();
        _hyperGraph.BuildRRsets(numR); // R1
        _hyperGraphVldt.BuildRRsets(numR); // R2
        _numRRsets = _hyperGraph.get_RR_sets_size();
        time1 += timerSubsim.get_operation_time();
        const auto infSelf = MaxCover(targetSize);
        time2 += timerSubsim.get_operation_time();
        const auto infVldt = _hyperGraphVldt.CalculateInf(_vecSeed);

        const auto degVldt = infVldt * _numRRsets / _numV;
        auto upperBound = _boundMin;

        const auto upperDegOPT = upperBound * _numRRsets / _numV;
        const auto lowerSelect = pow2(sqrt(degVldt + a1 * 2.0 / 9.0) - sqrt(a1 / 2.0)) - a1 / 18.0;
        const auto upperOPT = pow2(sqrt(upperDegOPT + a2 / 2.0) + sqrt(a2 / 2.0));
        const auto currApprox = lowerSelect / upperOPT;
        std::cout << "lower bound: " << (lowerSelect * _numV / _numRRsets) << ", upperBound: " << (upperOPT * _numV / _numRRsets) << std::endl;
        std::cout << "-->SUBSIM (" << idx << "/" << numIter << ") approx. (max-cover): " << currApprox <<
                  " (" << infSelf / upperBound << "), #RR sets: " << _numRRsets << '\n';
        double avgSize = _hyperGraph.HyperedgeAvg();

        if (currApprox >= approx - epsilon)
        {
            _res.set_approximation(currApprox);
            _res.set_running_time(timerSubsim.get_total_time());
            _res.set_influence(infVldt);
            _res.set_influence_original(infSelf);
            _res.set_seed_vec(_vecSeed);
            _res.set_RR_sets_size(_numRRsets * 2);
            std::cout << "==>Influence via R2: " << infVldt << ", time: " << _res.get_running_time() << '\n';
            std::cout << "==>Time for RR sets and greedy: " << time1 << ", " << time2 << '\n';
            return 0;
        }
    }

    return 0.0;
}

int decideMultiple(int ratio, int numRRsets)
{
    int multiple = 1;

    if (numRRsets < 100)
    {
        return multiple;
    }

    if (ratio >= 32)
    {
        multiple = 8;
    }
    else if (ratio >= 16)
    {
        multiple = 4;
    }
    else if (ratio >= 4)
    {
        multiple = 2;
    }
    else
    {
        multiple = 1;
    }

    return multiple;
}

double Alg::IncreaseR2(std::unordered_set<uint32_t> &connSet, double a, double upperOPT, double targetAppr)
{
    if (_hyperGraph.get_func_type() == IM) // Code from SUBSIM
    {
        size_t vldtRRsets = _hyperGraphVldt.get_RR_sets_size();
        size_t R1RRsets = _hyperGraph.get_RR_sets_size();
        int multiple = 4;
        double estimateAppr;
        double lowerSelect;
        double infVldt = _hyperGraphVldt.CalculateInfEarlyStop();
        double degVldt = infVldt * vldtRRsets / _numV;
        lowerSelect = pow2(sqrt(degVldt * multiple + a * 2.0 / 9.0) - sqrt(a / 2.0)) - a / 18.0 ;
        estimateAppr = (lowerSelect / (multiple * vldtRRsets)) / (upperOPT / R1RRsets);

        if (estimateAppr < targetAppr)
        {
            return 0.0;
        }

        _hyperGraphVldt.BuildRRsetsEarlyStop(connSet, vldtRRsets * multiple);
        infVldt = _hyperGraphVldt.CalculateInfEarlyStop();
        vldtRRsets = _hyperGraphVldt.get_RR_sets_size();
        degVldt = infVldt * vldtRRsets / _numV;
        lowerSelect = pow2(sqrt(degVldt + a * 2.0 / 9.0) - sqrt(a / 2.0)) - a / 18.0;
        double newAppr = (lowerSelect / vldtRRsets) / (upperOPT / R1RRsets);
        return (newAppr > targetAppr) ? newAppr : 0.0;
    }
    else if (_hyperGraph.get_func_type() == FIM || _hyperGraph.get_func_type() == RFIM)
    {
        size_t vldtRRsets = _hyperGraphVldt.get_RR_sets_size();
        size_t R1RRsets = _hyperGraph.get_RR_sets_size();
        int multiple = 4;
        double estimateAppr;
        double lowerSelect;
        double infVldt = _hyperGraphVldt.CalculateInfEarlyStop();
        double degVldt = infVldt * vldtRRsets;
        _hyperGraphVldt.BuildRRsetsEarlyStop(connSet, vldtRRsets * multiple);
        infVldt = _hyperGraphVldt.CalculateInfEarlyStop();
        vldtRRsets = _hyperGraphVldt.get_RR_sets_size();
        degVldt = infVldt * vldtRRsets;
        _lowerDeg = degVldt;
        double lowerH = calculateLowerBound(a);
        double newAppr = lowerSelect / upperOPT;
        return (newAppr > targetAppr) ? newAppr : 0.0;
    }
    else
    {
        std::cout << "Error: Unknown function type" << std::endl;
        exit(1);
    }
}

double Alg::FindRemSet(const int targetSize, const double epsilon, const double targetEpsilon, const double delta)
{
    Timer timerSubsim("SUBSIM");
    size_t subSeedSetSize = _vecSeed.size();
    const double e = exp(1);
    const double approx = 1 - 1.0 / e;
    // delta for upper bound on the number of RR sets
    const double delta_upper = delta / 3.0;
    const double alpha = sqrt(log(3.0 / delta_upper));
    const double beta = sqrt((1 - 1 / e) * (logcnk(_numV, targetSize - subSeedSetSize) + log(3.0 / delta_upper)));
    //const auto numRbase = size_t(2.0 * pow2((1 - 1 / e) * alpha + beta));
    const auto numRbase = size_t(_baseNumRRsets);
    const auto maxNumR = size_t(2.0 * _numV * pow2(alpha + beta) / targetSize / pow2(epsilon)) + 1;
    const auto numIter = (size_t)log2(maxNumR / numRbase) + 1;
    const double a1 = log(numIter * 3.0 / delta);
    const double a2 = log(numIter * 3.0 / delta);
    double time1 = 0.0, time2 = 0.0, time3 = 0.0;
    double time4 = 0.0;
    double infVldt = 0.0;
    double currApprox = 0.0;
    double infSelf = 0.0;
    int multiple = 1;
    std::unordered_set<uint32_t> subSeedSet(_vecSeed.begin(), _vecSeed.end());
    std::vector<uint32_t> vecSubSeed(_vecSeed.begin(), _vecSeed.end());

    for (auto idx = 1; idx <= numIter; idx++)
    {
        const auto numR = numRbase << (idx-1);
        std::cout << "Iteration: " << idx << " RR set: " << numR << std::endl;
        timerSubsim.get_operation_time();
        _hyperGraph.BuildRRsetsEarlyStop(subSeedSet, numR); // R1
        _numRRsets = _hyperGraph.get_RR_sets_size();
        time1 += timerSubsim.get_operation_time();
        infSelf = MaxCoverIMSentinel(vecSubSeed, targetSize);
        time2 += timerSubsim.get_operation_time();
        std::unordered_set<uint32_t> connSet(_vecSeed.begin(), _vecSeed.end());
        infVldt = _hyperGraphVldt.EvalSeedSetInf(connSet, _numRRsets * multiple);
        time4 += timerSubsim.get_operation_time();
        const auto degVldt = infVldt * multiple * _numRRsets / _numV;
        auto upperBound = _boundMin;

        const auto upperDegOPT = upperBound * _numRRsets / _numV;
        const auto lowerSelect = (pow2(sqrt(degVldt + a2 * 2.0 / 9.0) - sqrt(a2 / 2.0)) - a2 / 18.0) / multiple;
        const auto upperOPT = pow2(sqrt(upperDegOPT + a2 / 2.0) + sqrt(a2 / 2.0));
        const auto currApprox = lowerSelect / upperOPT;

        std::cout << "lower bound: " << (lowerSelect * _numV / _numRRsets) << ", upperBound: " << (upperOPT * _numV / _numRRsets) << std::endl;
        std::cout << "-->SUBSIM (" << idx << "/" << numIter << ") approx. (max-cover): " << currApprox <<
                  " (" << infSelf / upperBound << "), #RR sets: " << _numRRsets << '\n';
        double fullRRSize = _hyperGraph.HyperedgeAvg();
        double truncRRSize = _hyperGraphVldt.EvalHyperedgeAvg();

        // if truncRRset is more efficient, increase the size of R2 in next iteration
        int ratio = fullRRSize / truncRRSize;
        multiple = decideMultiple(ratio, _numRRsets);

        if (currApprox >= approx - targetEpsilon)
        {
            _res.set_approximation(currApprox);
            _res.set_running_time(timerSubsim.get_total_time());
            _res.set_influence(infVldt);
            _res.set_influence_original(infSelf);
            _res.set_seed_vec(_vecSeed);
            _res.set_RR_sets_size(_numRRsets * 2);
            std::cout << "==>Time for full RR in IM-Sentinel phase: " << time1  << std::endl;
            std::cout << "==>Time for truncated RR in IM-Sentinel phase: " << time4 << std::endl;
            std::cout << "==>Time for greedy in IM-Sentinel phase: " << time2 << std::endl;
            std::cout << "==>Influence via R2 in IM-Sentinel phase: " << infVldt << ", time: " << _res.get_running_time() << '\n';
            return 0;
        }
    }

    return 0.0;
}

double Alg::FindDynamSub(const int totalTargetSize, const double epsilon, const double delta)
{
    Timer timerSubsim("SUBSIM");
    const double e = exp(1);
    const double x = (1.0 - 1.0 / totalTargetSize);
    const int minSubSize = ceil(log(1 - epsilon) / log(x));
    const double alpha = sqrt(log(6.0 / delta));
    const double beta = sqrt((1 - 1 / e) * (logcnk(_numV, totalTargetSize) + log(6.0 / delta)));
    const auto numRbasePrevious = size_t(2.0 * pow2((1 - 1 / e) * alpha + beta) / totalTargetSize);
    const auto numRbase = size_t(_baseNumRRsets);

    // the successful probability of at least 1-delta/3
    const auto maxNumR = size_t(2.0 * _numV * pow2((1 - 1 / e) * alpha + beta) / totalTargetSize / pow2(epsilon)) + 1;
    const auto numIter = (size_t)log2(maxNumR / numRbase) + 1;
    const double a1 = log(numIter * 3.0 / delta);
    const double a2 = log(numIter * 6.0 / delta);
    double time1 = 0.0, time2 = 0.0, time3 = 0.0;
    double time4 = 0.0;
    int multiple = 1;
    double infVldt = 0.0;
    bool firstRound = true;

    for (auto idx = 1; idx <= numIter; idx++)
    {
        const auto numR = numRbase << (idx-1);
        std::cout << "Iteration: " << idx << " RR set: " << numR << std::endl;
        timerSubsim.get_operation_time();
        //build R1
        _hyperGraph.BuildRRsets(numR); // R1
        _numRRsets = _hyperGraph.get_RR_sets_size();
        time1 += timerSubsim.get_operation_time();
        _vecVldtInf.clear();
        int targetSize = firstRound ? (totalTargetSize / 4) : (totalTargetSize / 8);
        firstRound = false;
        const auto infSelf = MaxCoverSentinelSet(targetSize, totalTargetSize);
        time2 += timerSubsim.get_operation_time();
        std::vector<double> vecAppro(_vecSeed.size());
        int lastPos = 0;
        bool found = false;

        if (_vecSeed.size() < targetSize)
        {
            //low influence
            continue;
        }

        double calcAppr = 0.0;

        for (int i = _vecSeed.size() - 1; i >= 0; i--)
        {
            infVldt = _vecVldtInf[i];
            double lowerDeg = infVldt;
            double upperDeg = _boundMin;
            double a = log(numIter * 6.0  / delta);
            double lower = pow2(sqrt(lowerDeg + a * 2.0 / 9.0) - sqrt(a / 2.0)) - a / 18.0;
            double upper = pow2(sqrt(upperDeg + a1 / 2.0) + sqrt(a1 / 2.0));
            upper = (upper > numR) ? numR : upper;
            vecAppro[i] = lower / upper;
            calcAppr = (1 - pow(x, i + 1) - epsilon) * 1.2;

            if (vecAppro[i] > calcAppr)
            {
                found = true;
                lastPos = i;
                break;
            }
        }


        if (!found)
        {
            lastPos = 0;
        }

        size_t setSize = (lastPos + 1);
        setSize = setSize > 10 ? setSize : 10;
        setSize = (setSize > targetSize) ? targetSize : setSize;
        setSize = (setSize < minSubSize) ? minSubSize : setSize;
        std::vector<uint32_t> dynSeedSet(_vecSeed.begin(), _vecSeed.begin() + setSize);
        _vecSeed.clear();
        _vecSeed.assign(dynSeedSet.begin(), dynSeedSet.end());
        std::unordered_set<uint32_t> connSet(_vecSeed.begin(), _vecSeed.end());
        _hyperGraphVldt.RefreshHypergraph();
        _hyperGraphVldt.BuildRRsetsEarlyStop(connSet, _numRRsets * multiple);
        infVldt = _hyperGraphVldt.CalculateInfEarlyStop();
        time4 += timerSubsim.get_operation_time();
        double degVldt = infVldt * multiple * _numRRsets / _numV;
        auto upperBound = _boundMin;

        double upperDegOPT = upperBound * _numRRsets / _numV;
        double lowerSelect = (pow2(sqrt(degVldt + a2 * 2.0 / 9.0) - sqrt(a2 / 2.0)) - a2 / 18.0) / multiple;

        if (lowerSelect < 0)
        {
            lowerSelect = 1.0 * _vecSeed.size() / _numV * _numRRsets * multiple;
        }

        double upperOPT = pow2(sqrt(upperDegOPT + a1 / 2.0) + sqrt(a1 / 2.0));
        upperOPT = (upperOPT > _numRRsets) ? _numRRsets : upperOPT;
        const auto currApprox = lowerSelect / upperOPT;
        std::cout << "lower bound: " << (lowerSelect * _numV / (_numRRsets)) << ", upperBound: " << (upperOPT * _numV / _numRRsets) << std::endl;
        std::cout << "-->SUBSIM (" << idx + 1 << "/" << numIter << ") approx. (max-cover): " << currApprox <<
                  " (" << infSelf / upperBound << "), #RR sets: " << _numRRsets << '\n';
        const double approx = 1 - pow(x, _vecSeed.size());
        double targetAppr = approx - epsilon;

        if (currApprox >= targetAppr)
        {
            goto succ;
        }

        if (_numRRsets < 100)
        {
            continue;
        }

        double fullRRSize = _hyperGraph.HyperedgeAvg();
        double truncRRSize = _hyperGraphVldt.HyperedgeAvg();

        if (fullRRSize / truncRRSize < 2)
        {
            continue;
        }

        double lowerThreshold = (upperOPT * _numV / _numRRsets) * targetAppr;

        if ((1.0 * infVldt / multiple) > lowerThreshold && lowerThreshold > 0)
        {
            double newAppr = IncreaseR2(connSet, a2, upperOPT, targetAppr);
            time4 += timerSubsim.get_operation_time();

            if (newAppr > targetAppr)
            {
                std::cout << "increase R2 successfully" << std::endl;
                infVldt = _hyperGraphVldt.CalculateInfEarlyStop();
                goto succ;
            }
        }
    }

succ:
    std::cout << "==>Time for full RR in SentinelSet phase: " << time1  << std::endl;
    std::cout << "==>Time for truncated RR in SentinelSet phase: " << time4 << std::endl;
    std::cout << "==>Time for greedy in SentinelSet phase: " << time2 << std::endl;
    std::cout << "==>size of sentinel set: " << _vecSeed.size() << ", inf: " << infVldt << std::endl;
    std::cout << "==>total time for SentinelSet phase: " << timerSubsim.get_total_time() << std::endl;
    return 0.0;
}

double Alg::subsimWithHIST(const int targetSize, const double epsilon, const double delta)
{
    Timer timerSubsim("SUBSIM");
    _baseNumRRsets = 3 * log(1 / delta);
    
    std::cout << std::endl;
    std::cout << "Sentinel Set Selection Phase" << std::endl;
    FindDynamSub(targetSize, epsilon / 2, delta / 2);
    _hyperGraph.RefreshHypergraph();
    _hyperGraphVldt.RefreshHypergraph();

    std::cout << std::endl;
    std::cout << "IM-Sentinel Phase" << std::endl;
    FindRemSet(targetSize, epsilon / 2, epsilon, delta / 2);
    _res.set_running_time(timerSubsim.get_total_time());
    return 0.0;
}

double Alg::calculateUpperBound(const double a1)
{
    size_t numR = _hyperGraph.get_RR_sets_size();
    double upperH = _hyperGraphVldt.get_isSaturateGreedy() ? 0.0 : _hyperGraphVldt.get_isWeightedAvgGreedy() ? 0.0 : DBL_MAX;
    if (_hyperGraph.get_func_type() == RFIM)
    {
        for (auto q = 0; q < _numMetric; q++)
        {
            double fairInf = _vecBoundMinFairInf[q];
            double upperDeg = fairInf * numR;
            double upper = (1.0 / numR) * pow2(sqrt(upperDeg + a1 / 2.0) + sqrt(a1 / 2.0));
            if (_hyperGraphVldt.get_isSaturateGreedy()) //S-HIST
            {
                upperH += std::min(upper / _vecOptFairInf[q], _c);
            }
            else
            {
                if (_hyperGraphVldt.get_isWeightedAvgGreedy()) // WeightedAvgGreedy
                    upperH += upper * _vecWeights[q];
                else // SingleGreedy
                    upperH = std::min(upperH, upper / _vecOptFairInf[q]);
            }
        }
    }
    else // FIM
    {
        upperH = (1.0 / numR) * pow2(sqrt(_boundMin * numR + a1 / 2.0) + sqrt(a1 / 2.0));
    }
    return upperH;
}

double Alg::calculateLowerBoundR1(const double a2, const int k)
{
    size_t numR = _hyperGraph.get_RR_sets_size();
    double lowerH = _hyperGraph.get_isSaturateGreedy() ? 0.0 : _hyperGraph.get_isWeightedAvgGreedy() ? 0.0 : DBL_MAX;
    if (_hyperGraph.get_func_type() == RFIM)
    {
        for (auto q = 0; q < _numMetric; q++)
        {
            double fairInf = _vecFairInfDiffSize[q][k];
            double lowerDeg = fairInf * numR;
            double lower = (1.0 / numR) * (pow2(sqrt(lowerDeg + a2 * 2.0 / 9.0) - sqrt(a2 / 2.0)) - a2 / 18.0);
            if (_hyperGraph.get_isSaturateGreedy()) //S-HIST
            {
                lowerH += std::min(lower / _vecOptFairInf[q], _c);
            }
            else
            {
                if (_hyperGraph.get_isWeightedAvgGreedy()) // WeightedAvgGreedy
                    lowerH += lower * _vecWeights[q];
                else // SingleGreedy
                    lowerH = std::min(lowerH, lower / _vecOptFairInf[q]);
            }
        }
    }
    else
    {
        double lowerDeg = _vecSelfVldtGFCoverage[k];
        lowerH = (1.0 / numR) * (pow2(sqrt(lowerDeg + a2 * 2.0 / 9.0) - sqrt(a2 / 2.0)) - a2 / 18.0);
    }
    return lowerH;
}

double Alg::calculateLowerBound(const double a2)
{
    size_t numR = _hyperGraphVldt.get_RR_sets_size();
    double lowerH = _hyperGraphVldt.get_isSaturateGreedy() ? 0.0 : _hyperGraphVldt.get_isWeightedAvgGreedy() ? 0.0 : DBL_MAX;

    if (_hyperGraphVldt.get_func_type() == RFIM)
    {
        for (auto q = 0; q < _numMetric; q++)
        {
            double fairInf = _vecFairInfVldt[q];
            double lowerDeg = fairInf * numR;
            double lower = (1.0 / numR) * (pow2(sqrt(lowerDeg + a2 * 2.0 / 9.0) - sqrt(a2 / 2.0)) - a2 / 18.0);
            if (_hyperGraphVldt.get_isSaturateGreedy()) //S-HIST
            {
                lowerH += std::min(lower / _vecOptFairInf[q], _c);
            }
            else
            {
                if (_hyperGraphVldt.get_isWeightedAvgGreedy()) // WeightedAvgGreedy
                    lowerH += lower * _vecWeights[q];
                else // SingleGreedy
                    lowerH = std::min(lowerH, lower / _vecOptFairInf[q]);
            }
        }
    }
    else //fim
    {
        lowerH = (1.0 / numR) * (pow2(sqrt(_lowerDeg + a2 * 2.0 / 9.0) - sqrt(a2 / 2.0)) - a2 / 18.0);
    }
    return lowerH;
}

double Alg::FindGeneralDynamSub(const int totalTargetSize, const double epsilon, const double delta)
{
    const double e = exp(1);
    double x = 1.0 - 1.0 / _k;
    double eps1 = epsilon, eps2 = epsilon, eps3 = epsilon, eps4 = epsilon;
    int minSubSize; size_t numRbase, maxNumR; size_t numIter;
    double alpha, beta, a1, a2;
    if (_hyperGraph.get_func_type() == RFIM && _hyperGraph.get_isSaturateGreedy())
    {
        eps1 = epsilon * _rFairInfMin / _numMetric;
        eps2 = epsilon * _rFairInfMin / _numMetric;
        minSubSize = ceil(log(1 - eps1) / log(x));
        alpha = sqrt(2 * log(6.0 * _numMetric / delta));
        beta = sqrt(((2 + 2.0 / 3 * eps1) * _fInfMax / _fInfMin) * (logcnk(_numV, _k) + log(6.0 * _numMetric / delta)));
        numRbase = size_t(_baseNumRRsets);
        numIter = (numIter > 30) ? 30 : numIter; // to acoid too many iterations, which is impractical
        maxNumR = size_t(2.0 * pow2(alpha + beta) / _fInfMin / pow2(eps1)) + 1;
        numIter = (size_t)log2(maxNumR / numRbase) + 1;
        a1 = log(numIter * 3.0 * _numMetric / delta);
        a2 = log(numIter * 6.0 * _numMetric / delta);
    }
    else //FIM
    {
        x = 1.0 - 1.0 / totalTargetSize;
        minSubSize = ceil(log(1 - epsilon) / log(x));
        alpha = sqrt(log(6.0 / delta));
        beta = sqrt((1 - 1 / e) * (logcnk(_numV, totalTargetSize) + log(6.0 / delta)));
        numRbase = size_t(_baseNumRRsets);
        maxNumR = size_t(2.0 * _numV * pow2((1 - 1 / e) * alpha + beta) / totalTargetSize / pow2(epsilon)) + 1;
        numIter = (size_t)log2(maxNumR / numRbase) + 1;
        a1 = log(numIter * 3.0 / delta);
        a2 = log(numIter * 6.0 / delta);
    }
    double infVldt = 0.0;
    for (auto idx = 1; idx <= numIter; idx++)
    {
        const auto numR = numRbase << (idx-1);
        std::cout << "Iteration: " << idx << " RR set: " << numR << std::endl;
        _hyperGraph.BuildRRsets(numR); // build R1
        _numRRsets = _hyperGraph.get_RR_sets_size();
        _vecSelfVldtGFCoverage.clear();
        int targetSize = (int)std::min(totalTargetSize / 4, _k); // aim for a small sentinel set
        for (auto q = 0; q < _numMetric; q++)
        {
            _vecFairInfDiffSize[q].clear();
            _vecFairInfDiffSize[q].resize(targetSize, 0.0);
        }
        double infSelf = MaxCoverGreedyMarginalFairSentinelSet(targetSize, totalTargetSize);
        std::vector<double> vecAppro(_vecSeed.size());
        int lastPos = 0;
        bool found = false;

        if (_vecSeed.size() < targetSize) continue;

        double calcAppr = 0.0;

        for (int i = _vecSeed.size() - 1; i >= 0; i--)
        {
            double upper = calculateUpperBound(a1);
            double lower = calculateLowerBoundR1(a2, i);
            vecAppro[i] = lower / upper;
            calcAppr = (1 - pow(x, i + 1) - epsilon) * 1.2; // slightly larger than target approximation
            if (vecAppro[i] > calcAppr)
            {
                found = true;
                lastPos = i;
                break;
            }
        }

        if (!found)
        {
            lastPos = 0;
        }

        size_t setSize = (lastPos + 1);
        setSize = setSize > 10 ? setSize : 10;
        setSize = (setSize > targetSize) ? targetSize : setSize;
        setSize = (setSize < minSubSize) ? minSubSize : setSize;
        std::vector<uint32_t> dynSeedSet(_vecSeed.begin(), _vecSeed.begin() + setSize);
        _vecSeed.clear();
        _vecSeed.assign(dynSeedSet.begin(), dynSeedSet.end());
        std::unordered_set<uint32_t> connSet(_vecSeed.begin(), _vecSeed.end());
        _hyperGraphVldt.RefreshHypergraph();
        _hyperGraphVldt.BuildRRsetsEarlyStop(connSet, _numRRsets);
        infVldt = _hyperGraphVldt.CalculateInfEarlyStop(); // get the fair influence by gHit
        double upperBound = _boundMin; // _boundMin is the general INFLUENCE

        _lowerDeg = infVldt * _numRRsets;
        double lowerSelect = calculateLowerBound(a2);
        double upperOPT = calculateUpperBound(a1);

        const auto currApprox = lowerSelect / upperOPT;
        std::cout << "lower bound: " << (lowerSelect) << ", upperBound: " << (upperOPT) << std::endl;
        std::cout << "-->SG-HIST (" << idx << "/" << numIter << ") approx. (max-cover): " << currApprox <<
                  " (" << infSelf / upperBound << "), #G-RR sets: " << _numRRsets << '\n';
        const double approx = 1 - pow(x, _vecSeed.size());
        double targetAppr = approx - epsilon;

        if (currApprox >= targetAppr) goto succ;
        if (_numRRsets < 100) continue;

        double fullRRSize = _hyperGraph.HyperedgeAvg();
        double truncRRSize = _hyperGraphVldt.HyperedgeAvg();
        if (fullRRSize / truncRRSize < 2) continue;

        double lowerThreshold = upperOPT * targetAppr;
        _res.set_RR_sets_size(_res.get_RRsets_size() + _hyperGraphVldt.get_RR_sets_size());
        if (infVldt > lowerThreshold && lowerThreshold > 0)
        {
            double newAppr = IncreaseR2(connSet, a2, upperOPT, targetAppr);

            if (newAppr > targetAppr)
            {
                std::cout << "increase R2 successfully" << std::endl;
                infVldt = _hyperGraphVldt.CalculateInfEarlyStop();
                infVldt = AdjustGeneralInfluenceBound(infVldt, _hyperGraphVldt.get_RR_sets_size());
                _res.set_RR_sets_size(_res.get_RRsets_size() + _hyperGraph.get_RR_sets_size());
                goto succ;
            }
        }
    }

succ:
    std::cout << "==>size of sentinel set: " << _vecSeed.size() << ", objective: " << infVldt << std::endl;
    return infVldt;
}

double Alg::FindGeneralRemSet(const int targetSize, const double epsilon, const double targetEpsilon, const double delta)
{
    size_t subSeedSetSize = _vecSeed.size();
    const double e = exp(1);
    double eps1 = epsilon, eps2 = epsilon, eps3 = epsilon, eps4 = epsilon;
    double approx, delta_upper, alpha, beta, a1, a2;
    size_t numRbase, maxNumR, numIter;
    if (_hyperGraph.get_func_type() == RFIM && _hyperGraph.get_isSaturateGreedy()) // s-hist
    {
        eps1 = epsilon * _rFairInfMin / _numMetric;
        eps2 = epsilon * _rFairInfMin / _numMetric;
        approx = 1 - exp(-_beta);
        alpha = sqrt(log(9.0 * _numMetric / delta));
        beta = sqrt((2 + 2.0 / 3 * eps2) * (logcnk(_numV - subSeedSetSize, targetSize - subSeedSetSize) + log(9.0 * _numMetric / delta)));
        numRbase = size_t(_baseNumRRsets);
        maxNumR = size_t(2.0 * _fInfMax * pow2(alpha + beta) / pow2(eps2 * _fInfMin)) + 1;
        numIter = (size_t)log2(maxNumR / numRbase) + 1;
        numIter = (numIter > 30) ? 30 : numIter; // to acoid too many iterations, which is impractical
        a1 = log(numIter * 3.0 * _numMetric / delta);
        a2 = log(numIter * 3.0 * _numMetric / delta);
    }
    else
    {
        approx = 1 - 1 / e; // FIM && other Baselines
        alpha = sqrt(log(3.0 / delta_upper));
        beta = sqrt((1 - 1 / e) * (logcnk(_numV - subSeedSetSize, targetSize - subSeedSetSize) + log(3.0 / delta_upper)));
        delta_upper = delta / 3.0;
        numRbase = size_t(_baseNumRRsets);
        maxNumR = size_t(2.0 * _numV * pow2(alpha + beta) / targetSize / pow2(epsilon)) + 1;
        numIter = (size_t)log2(maxNumR / numRbase) + 1;
        numIter = (numIter > 30) ? 30 : numIter; // to acoid too many iterations, which is impractical
        a1 = log(numIter * 3.0 / delta);
        a2 = log(numIter * 3.0 / delta);
    }

    double infVldt = 0.0;
    double currApprox = 0.0;
    double infSelf = 0.0;
    int multiple = 1;
    std::unordered_set<uint32_t> subSeedSet(_vecSeed.begin(), _vecSeed.end());
    std::vector<uint32_t> vecSubSeed(_vecSeed.begin(), _vecSeed.end());

    for (auto idx = 1; idx <= numIter; idx++)
    {
        const auto numR = numRbase << (idx-1);
        std::cout << "Iteration: " << idx << " RR set: " << numR << std::endl;
        _hyperGraph.BuildRRsetsEarlyStop(subSeedSet, numR); // R1
        _numRRsets = _hyperGraph.get_RR_sets_size();
        infSelf = MaxCoverGreedyMarginalFairIMSentinel(vecSubSeed, targetSize);
        std::unordered_set<uint32_t> connSet(_vecSeed.begin(), _vecSeed.end());
        infVldt = _hyperGraphVldt.EvalSeedSetInf(connSet, _numRRsets * multiple);
        double upperBound = _boundMin; // _boundMin is the minimum upper bound of general INFLUENCE
        _lowerDeg = infVldt * _numRRsets * multiple;
        double lowerSelect = calculateLowerBound(a2);
        double upperOPT = calculateUpperBound(a1);

        const auto currApprox = lowerSelect / upperOPT;

        std::cout << "lower bound: " << (lowerSelect) << ", upperBound: " << (upperOPT) << std::endl;
        std::cout << "-->SG-HIST (" << idx << "/" << numIter << ") approx. (max-cover): " << currApprox <<
                  " (" << infSelf / upperBound << "), #G-RR sets: " << _numRRsets * multiple << '\n';
        double fullRRSize = _hyperGraph.HyperedgeAvg();
        double truncRRSize = _hyperGraphVldt.EvalHyperedgeAvg();

        // if truncRRset is more efficient, increase the size of R2 in next iteration
        int ratio = fullRRSize / truncRRSize;
        multiple = decideMultiple(ratio, _numRRsets);

        if (upperOPT < _threshold && _hyperGraph.get_func_type() == RFIM && lowerSelect < _threshold && _hyperGraph.get_isSaturateGreedy() && idx >= 8)
        {
            _res.set_RR_sets_size(_res.get_RRsets_size() + _hyperGraph.get_RR_sets_size() + _hyperGraphVldt.get_RR_sets_size());
            std::cout << "==>Early stop since upper bound < threshold: " << upperOPT << " < " << _threshold << std::endl;
            return lowerSelect;
        }

        //EarlyStop, if currApprox is sufficiently large or we are using the weighted average baselines which do not have theoretical approxmation guarantees
        if (currApprox >= approx - targetEpsilon || ((currApprox >= 0.99 || _hyperGraph.get_isWeightedAvgGreedy()) && idx >= 8))
        {
            _res.set_approximation(currApprox);
            _res.set_influence(infVldt);
            _res.set_influence_original(infSelf);
            if (_hyperGraph.get_func_type() == RFIM)
            {
                if (_hyperGraph.get_isSaturateGreedy())
                {
                    _res.set_robustness_score(infVldt);
                    _res.set_robustness_score_ori(infSelf);
                    _res.set_fair_inf(lowerSelect / _numMetric);
                }
                else
                {
                    _res.set_fair_inf(infVldt);
                }
            }
            else if (_hyperGraph.get_func_type() == FIM)
            {
                _res.set_fair_inf(infVldt);
                _res.set_fair_inf_ori(infSelf);
            }
            _res.set_seed_vec(_vecSeed);
            _res.set_RR_sets_size(_res.get_RRsets_size() + _hyperGraph.get_RR_sets_size() + _hyperGraphVldt.get_RR_sets_size());
            if (_hyperGraph.get_func_type() == FIM) std::cout << "==>Generalized Fair Influence via R2 in Remaining Selection phase: " << infVldt << std::endl;
            else
            {
                if (_hyperGraph.get_isWeightedAvgGreedy()) // WeightedAvgGreedy
                    std::cout << "==>Weighted Average Fair Influence via R2 in Remaining Selection phase: " << infVldt << std::endl;
                else // S-HIST, SingleGreedy, AllGreedy
                    std::cout << "==>Robustness Score via R2 in Remaining Selection phase: " << infVldt << std::endl;
            }
            if (_hyperGraph.get_func_type() == FIM) return infVldt;
            else return lowerSelect;
        }
    }

    return infVldt / (1 + epsilon);
}

double Alg::SaturateGeneralizedHIST(const int targetSize, const double epsilon, const double delta)
{
    double res = 0.0;
    _baseNumRRsets = 3.0 * log(2.0 / delta);
    // if (_hyperGraph.get_func_type() == RFIM && !_hyperGraph.get_isWeightedAvgGreedy())
    if (_hyperGraph.get_func_type() == RFIM)
    {
        _baseNumRRsets = 3.0 * log(6.0 * _numMetric / delta);
        _baseNumRRsets = _baseNumRRsets < 1024 ? 1024 : _baseNumRRsets;  // To avoid useless iterations with few rrsets.
    }
    
    double epsilon0 = epsilon / 2;
    double delta0 = delta / 2;
    std::cout << std::endl;
    std::cout << "SG-HIST: Sentinel Set Selection Phase" << std::endl;
    FindGeneralDynamSub(targetSize, epsilon0, delta0);
    _hyperGraph.RefreshHypergraph();
    _hyperGraphVldt.RefreshHypergraph();

    std::cout << "SG-HIST: Remaining Set Selection Phase" << std::endl;
    res = FindGeneralRemSet(targetSize, epsilon0, epsilon, delta0);
    _hyperGraph.RefreshHypergraph();
    _hyperGraphVldt.RefreshHypergraph();
    return res;
}

double Alg::saturateGreedyWithHIST(const int targetSize, const double epsilon, const double delta, const double gamma, const double saturateRatio)
{
    Timer timerSGwH("saturateGreedyWithHIST");
    double c_min = 0.0, c_max = 2.0, accumulatedDelta = 0.0, e = exp(1);
    int seedSetSize;
    const double beta = (saturateRatio >= 1.0) ? saturateRatio : 1 + log(_numMetric) + log(6 / (e * epsilon));
    const double deltaHalf = 0.5 * delta; // if we treat the OPT of FIM of each metric is not given
    const double newDelta = deltaHalf / ceil(log(epsilon / (c_max - c_min)) / log(5.0 / 6));
    const double newEpsilon = epsilon / (6 * _numMetric);
    _k = targetSize;
    _beta = (saturateRatio >= 1.0) ? beta : beta;
    int saturateSeedSize = (int)ceil(beta * targetSize);
    double rFObjH = 0.0;
    // _fInfMin = targetSize * _alphaMin; // _alphaMin may be very small
    _fInfMin = targetSize; // roughly assume the minimum fair influence is at least k

    std::cout << "==> Precalculating the fair influence for each metric..." << std::endl;
    for(auto q = 0; q < _numMetric; q++)
    {
        double fairInf = 0.0;
        _hyperGraph.set_func_type(FIM, q);
        _hyperGraphVldt.set_func_type(FIM, q);
        _hyperGraphOri.set_func_type(FIM, q);
        fairInf = SaturateGeneralizedHIST(targetSize, epsilon / (6 * _numMetric), deltaHalf / _numMetric);
        _vecOptFairInf[q] = fairInf / (1.0 - 1.0 / e);
        std::cout << "---> The fair obj of metric " << q << " by G-HIST: " << fairInf << std::endl;
        _fInfMax = std::max(_fInfMax, _vecOptFairInf[q]);
    }

    std::cout << "\n\n==> Precalculating the fair influence for each metric done!" << std::endl;
    std::cout << "==> S-HIST: Start!" << std::endl;
    timerSGwH.refresh_time();
    _rFairInfMin = (_fInfMin / _fInfMax) * _numMetric; // we assume (_fInfMin / _fInfMax) < c for all c
    _hyperGraph.set_func_type(RFIM);
    _hyperGraphVldt.set_func_type(RFIM);
    _hyperGraphOri.set_func_type(RFIM);

    double c_min_last = 0.0; // for experiment c_max = 2.0
    bool if_update_c_min = false; // for experiment c_max = 2.0
    while (c_max - c_min >= epsilon)
    {
        _c = (c_min + c_max) / 2;
        double targetThreshold = _c * (_numMetric - epsilon / 3);
        _threshold = targetThreshold;
        rFObjH = SaturateGeneralizedHIST(saturateSeedSize, newEpsilon, newDelta);
        if (rFObjH <= targetThreshold)
        {
            c_max = _c;
            std::cout << "---> c_max <- " << _c << std::endl;
            if_update_c_min = false; // for experiment c_max = 2.0
        }
        else
        {
            c_min = _c * (1 - epsilon/3);
            _res.log_best_res();
            std::cout << "---> c_min <- " << _c * (1 - epsilon/3) << std::endl;
            if (if_update_c_min && c_min_last >= c_min) break; // for experiment c_max = 2.0
            if_update_c_min = true; // for experiment c_max = 2.0
            c_min_last = c_min; // for experiment c_max = 2.0
        }
    }
    // use the best _c to generate the final seeds
    _res.set_running_time(timerSGwH.get_total_time());
    _res.log_final_res();
    _vecSeed = _res.get_seed_vec();

    // MC simulation vertification
    std::cout << "==> Monte Carlo Simulation, Start!" << std::endl;
    _hyperGraphOri.set_isSaturateGreedy(false);
    std::unordered_set<uint32_t> seedSet(_vecSeed.begin(), _vecSeed.end());
    double fairInfMC = _hyperGraphOri.MonteCarloSimulate(seedSet, _numMC);
    _res.set_resMC(fairInfMC);

    return _res.get_fair_inf();
}

double Alg::singleGreedyWithHIST(const int targetSize, const double epsilon, const double delta, const double gamma, const double saturateRatio)
{
    Timer timerSG("singleGreedyWithHIST");
    const double e = exp(1);
    const double beta = (saturateRatio >= 1.0) ? saturateRatio : 1 + log(_numMetric) + log(6 / (e * epsilon));
    const double deltaHalf = 0.5 * delta; // if we treat the OPT of FIM of each metric is not given
    _k = targetSize;
    _beta = (saturateRatio >= 1.0) ? beta : beta;

    std::cout << "==> Precalculating the fair influence for each metric..." << std::endl;
    for(auto q = 0; q < _numMetric; q++)
    {
        double fairInf = 0.0;
        _hyperGraph.set_func_type(FIM, q);
        _hyperGraphVldt.set_func_type(FIM, q);
        _hyperGraphOri.set_func_type(FIM, q);
        fairInf = SaturateGeneralizedHIST(targetSize, epsilon, deltaHalf / _numMetric);
        _vecOptFairInf[q] = fairInf / (1.0 - 1.0 / e);
        std::cout << "---> The fair obj of metric " << q << " by G-HIST: " << fairInf << std::endl;
        // use MC simulations to estimate the greedy fair influence
        // std::unordered_set<uint32_t> seedSet(_vecSeed.begin(), _vecSeed.end());
        // double fairInfMC = _hyperGraphOri.MonteCarloSimulate(seedSet, _numMC);
        // std::cout << "---> The fair obj of metric " << q << " evaluated by MC simulation: " << fairInfMC << std::endl;
        // fairInfMC /= (1 - 1 / exp(1));
        // _vecOptFairInf[q] = fairInfMC;
    }
    std::cout << "\n\n==> Precalculating the fair influence for each metric done!" << std::endl;
    std::cout << "==> Single Greedy: Start!" << std::endl;
    timerSG.refresh_time();
    _hyperGraph.set_func_type(RFIM);
    _hyperGraphVldt.set_func_type(RFIM);
    _hyperGraph.set_isSaturateGreedy(false);
    _hyperGraphVldt.set_isSaturateGreedy(false);

    double fInf = SaturateGeneralizedHIST((int)ceil(beta * targetSize), epsilon, deltaHalf);
    _res.set_fair_inf(fInf);
    _res.set_running_time(timerSG.get_total_time());

    // MC simulation vertification
    std::cout << "==> Monte Carlo Simulation, Start!" << std::endl;
    _hyperGraphOri.set_func_type(RFIM);
    _hyperGraphOri.set_isSaturateGreedy(false);
    std::unordered_set<uint32_t> seedSet(_vecSeed.begin(), _vecSeed.end());
    double fairInfMC = _hyperGraphOri.MonteCarloSimulate(seedSet, _numMC);
    _res.set_resMC(fairInfMC);

    return _res.get_fair_inf();
}

double Alg::allGreedyWithHIST(const int targetSize, const double epsilon, const double delta, const double gamma, const double saturateRatio)
{
    Timer timerAGwH("allGreedyWithHIST");
    const double e = exp(1);
    // const double beta = (saturateRatio >= 1.0) ? saturateRatio : 1 + log(_numMetric) + log(3 / gamma);
    const double beta = (saturateRatio >= 1.0) ? saturateRatio : 1 + log(_numMetric) + log(6 / (e * epsilon));
    const double deltaHalf = 0.5 * delta; // if we treat the OPT of FIM of each metric is not given
    // const double deltaHalf = delta; // if we treat the OPT of FIM of each metric is given

    _k = targetSize;
    _beta = (saturateRatio >= 1.0) ? beta : beta;

    std::cout << "==> Precalculating the fair influence for each metric..." << std::endl;
    for(auto q = 0; q < _numMetric; q++)
    {
        double fairInf = 0.0;
        _hyperGraph.set_func_type(FIM, q);
        _hyperGraphVldt.set_func_type(FIM, q);
        _hyperGraphOri.set_func_type(FIM, q);
        fairInf = SaturateGeneralizedHIST(targetSize, epsilon, deltaHalf / _numMetric);
        _vecOptFairInf[q] = fairInf / (1.0 - 1.0 / e);
        std::cout << "---> The fair obj of metric " << q << " by G-HIST: " << fairInf << std::endl;
        // use MC simulations to estimate the greedy fair influence
        // std::unordered_set<uint32_t> seedSet(_vecSeed.begin(), _vecSeed.end());
        // double fairInfMC = _hyperGraphOri.MonteCarloSimulate(seedSet, _numMC);
        // std::cout << "---> The fair obj of metric " << q << " evaluated by MC simulation: " << fairInfMC << std::endl;
        // fairInfMC /= (1 - 1 / exp(1));
        // _vecOptFairInf[q] = fairInfMC;
    }
    std::cout << "\n\n==> Precalculating the fair influence for each metric done!" << std::endl;
    std::cout << "==> All Greedy: Start!" << std::endl;
    timerAGwH.refresh_time();
    const double newDelta = deltaHalf / (2 * _numMetric);
    double largestGFInf = 0.0;

    _hyperGraph.set_isSaturateGreedy(false);
    _hyperGraphVldt.set_isSaturateGreedy(false);

    Nodelist _vecCurrSeed;

    for(int q = 0; q < _numMetric; q++) // select optimal seeds under each metric
    {
        double fairInf = 0.0, gFairInf = 0.0;
        _hyperGraph.set_func_type(FIM, q);
        _hyperGraphVldt.set_func_type(FIM, q);
        fairInf = SaturateGeneralizedHIST((int)ceil(beta * targetSize), epsilon, newDelta);
        _vecCurrSeed = _res.get_seed_vec(); // find the opt seeds

        // Estimate the robust fair influence
        _hyperGraphOri.set_func_type(RFIM);
        _hyperGraphOri.set_isSaturateGreedy(false);
        std::unordered_set<uint32_t> connSet(_vecCurrSeed.begin(), _vecCurrSeed.end());
        gFairInf = _hyperGraphOri.MonteCarloSimulate(connSet, _numMC); // Vertified by MC simulation
        _hyperGraphOri.RefreshHypergraph();
        std::cout << " ==>All Greedy -> The objective rho of metric "<< q <<": " << gFairInf << std::endl;
        if (gFairInf > largestGFInf)
        {
            largestGFInf = gFairInf;
            _res.set_fair_inf(gFairInf);
            _res.log_best_res();
        }
    }

    _res.set_running_time(timerAGwH.get_total_time());
    _res.log_final_res();
    _vecSeed = _res.get_seed_vec();
    // MC simulation vertification
    std::cout << "==> Monte Carlo Simulation, Start!" << std::endl;
    _hyperGraphOri.set_isSaturateGreedy(false);
    std::unordered_set<uint32_t> seedSet(_vecSeed.begin(), _vecSeed.end());
    double fairInfMC = _hyperGraphOri.MonteCarloSimulate(seedSet, _numMC);
    _res.set_resMC(fairInfMC);
    return largestGFInf;
}

double Alg::wAvgGreedyWithHIST(const int targetSize, const double epsilon, const double delta, const double gamma, const double saturateRatio, const MethodType weightMode)
{
    Timer timerWAG("wAvgGreedyWithHIST");
    const double e = exp(1);
    const double beta = (saturateRatio >= 1.0) ? saturateRatio : 1 + log(_numMetric) + log(6 / (e * epsilon));
    const double deltaHalf = 0.5 * delta; // if we treat the OPT of FIM of each metric is not given
    _k = targetSize;
    _beta = (saturateRatio >= 1.0) ? beta : beta;

    std::cout << "==> Precalculating the fair influence for each metric..." << std::endl;
    for(auto q = 0; q < _numMetric; q++)
    {
        double fairInf = 0.0;
        _hyperGraph.set_func_type(FIM, q);
        _hyperGraphVldt.set_func_type(FIM, q);
        _hyperGraphOri.set_func_type(FIM, q);
        fairInf = SaturateGeneralizedHIST(targetSize, epsilon, deltaHalf / _numMetric);
        _vecOptFairInf[q] = fairInf / (1.0 - 1.0 / e);
        std::cout << "---> The fair obj of metric " << q << " by G-HIST: " << fairInf << std::endl;
        // use MC simulations to estimate the greedy fair influence
        // std::unordered_set<uint32_t> seedSet(_vecSeed.begin(), _vecSeed.end());
        // double fairInfMC = _hyperGraphOri.MonteCarloSimulate(seedSet, _numMC);
        // std::cout << "---> The fair obj of metric " << q << " evaluated by MC simulation: " << fairInfMC << std::endl;
        // fairInfMC /= (1 - 1 / exp(1));
        // _vecOptFairInf[q] = fairInfMC;
    }
    std::cout << "\n\n==> Precalculating the fair influence for each metric done!" << std::endl;
    if (weightMode == WEIGHTED_AVERAGE_GREEDY1)
    {
        std::cout << "==> Equal Weighted Average Greedy: Start!" << std::endl;
    }
    else if (weightMode == WEIGHTED_AVERAGE_GREEDY2)
    {
        std::cout << "==> Weighted Average Greedy with pre-calculated weights: Start!" << std::endl;
    }
    timerWAG.refresh_time();
    _hyperGraph.set_func_type(RFIM);
    _hyperGraphVldt.set_func_type(RFIM);
    _hyperGraph.set_isSaturateGreedy(false);
    _hyperGraphVldt.set_isSaturateGreedy(false);
    _hyperGraph.set_isWeightedAvgGreedy(true);
    _hyperGraphVldt.set_isWeightedAvgGreedy(true);

    double sumOptFairInf = 0.0;
    if (weightMode == WEIGHTED_AVERAGE_GREEDY1)
    {
        for(auto q = 0; q < _numMetric; q++) sumOptFairInf += _vecOptFairInf[q];
        for(auto q = 0; q < _numMetric; q++) _vecWeights[q] = 1.0 / sumOptFairInf;
    }
    else if (weightMode == WEIGHTED_AVERAGE_GREEDY2)
    {
        for(auto q = 0; q < _numMetric; q++) sumOptFairInf += _vecOptFairInf[q];
        for(auto q = 0; q < _numMetric; q++) _vecWeights[q] = _vecOptFairInf[q] / sumOptFairInf;
    }

    double fInf = SaturateGeneralizedHIST((int)ceil(beta * targetSize), epsilon, deltaHalf);
    _res.set_fair_inf(fInf);
    _res.set_running_time(timerWAG.get_total_time());

    // MC simulation vertification
    std::cout << "==> Monte Carlo Simulation, Start!" << std::endl;
    _hyperGraphOri.set_func_type(RFIM);
    _hyperGraphOri.set_isSaturateGreedy(false);
    _hyperGraphOri.set_isWeightedAvgGreedy(false);
    std::unordered_set<uint32_t> seedSet(_vecSeed.begin(), _vecSeed.end());
    double fairInfMC = _hyperGraphOri.MonteCarloSimulate(seedSet, _numMC);
    _res.set_resMC(fairInfMC);

    return _res.get_fair_inf();
}

