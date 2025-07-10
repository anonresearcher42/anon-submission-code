#pragma once

bool greater_first(const std::pair<uint32_t, float> x, const std::pair<uint32_t, float> y)
{
    if (x.second > y.second) return true;
    else if (x.second < y.second) return false;
    else if (x.first < y.first) return true;
    else return false;
}

class GraphBase
{
public:
    /// Format the input for future computing, which is much faster for loading. Vector serialization is used.
    /// Add commFilenames.
    static void FormatGraph(const std::string filename, ProbDist probDist,
        const float sum, const float prob, const std::string skewType,
        const std::string commFilenames = "")
    {
        size_t numV, numE;
        uint32_t srcId, dstId;
        float weight = 0.0;
        std::ifstream infile(filename);

        if (!infile.is_open())
        {
            std::cout << "The file \"" + filename + "\" can NOT be opened\n";
            return;
        }

        infile >> numV >> numE;
        Graph vecGRev(numV);
        Graph vecG(numV); //forward graph
        std::vector<size_t> vecInDeg(numV);

        for (auto i = numE; i--;)
        {
            if (probDist == WEIGHTS)
            {
                infile >> srcId >> dstId >> weight;
            }
            else
            {
                infile >> srcId >> dstId;
            }

            vecGRev[dstId].push_back(Edge(srcId, weight));
            vecG[srcId].push_back(Edge(dstId, weight));
        }

        infile.close();

        for (auto idx = 0; idx < numV; idx++)
        {
            vecInDeg[idx] = vecGRev[idx].size();
        }

        if (probDist == WC)
        {
            for (size_t i = 0; i < vecGRev.size(); i++)
            {
                if (vecGRev[i].size() == 0) continue;

                weight = sum / vecInDeg[i];

                for (size_t j = 0; j < vecGRev[i].size(); j++)
                {
                    vecGRev[i][j].second = weight;
                }
            }

            for (size_t i = 0; i < vecG.size(); i++)
            {
                if (vecG[i].size() == 0) continue;
                for (size_t j = 0; j < vecG[i].size(); j++)
                {
                    weight = vecGRev[vecG[i][j].first][0].second;
                    vecG[i][j].second = weight;
                }
            }
        }
        else if (probDist == UNIFORM)
        {
            // Uniform probability
            for (auto& nbrs : vecGRev)
            {
                for (auto& nbr : nbrs)
                {
                    nbr.second = prob;
                }
            }

            for (auto& nbrs : vecG)
            {
                for (auto& nbr : nbrs)
                {
                    nbr.second = prob;
                }
            }
        }
        // exponential distribution with lambada equal to its in-degree
        else if (probDist == SKEWED)
        {
            if (skewType == "weibull")
            {
                std::default_random_engine generator(time(NULL));
                double min_value = (1e-8 < 1.0/numV)? 1e-8: 1.0/numV;

                for (size_t i = 0; i< vecGRev.size(); i++)
                {

                    if (vecGRev[i].size() == 0) continue;
                    double sum = 0.0;
                    for (size_t j = 0; j < vecGRev[i].size(); j++)
                    {
                        // random number from (0, 10)
                        double a = dsfmt_gv_genrand_open_open() * 10;
                        double b = dsfmt_gv_genrand_open_open() * 10;
                        std::weibull_distribution<double> distribution(a,b);
                        auto weight =  distribution(generator);
                        vecGRev[i][j].second = weight;

                        sum += weight;
                    }

                    for (size_t j = 0; j < vecGRev[i].size(); j++)
                    {
                        auto weight = vecGRev[i][j].second/sum;
                        vecGRev[i][j].second = (weight > min_value)?weight:min_value; 
                    }
                    sort(vecGRev[i].begin(), vecGRev[i].end(), greater_first);
                }
            }
            else
            {
                double min_value = (1e-8 < 1.0/numV)? 1e-8: 1.0/numV;
                for (size_t i = 0; i<vecGRev.size(); i++)
                {
                    if (vecGRev[i].size() == 0) continue;
                    double sum = 0.0;
                    for (size_t j = 0; j < vecGRev[i].size(); j++)
                    {
                        // lambda = 1
                        auto weight = -log ( 1.0 - dsfmt_gv_genrand_open_open());
                        vecGRev[i][j].second = weight;
                        sum += weight;
                    }

                    for (size_t j = 0; j < vecGRev[i].size(); j++)
                    {
                        double weight = vecGRev[i][j].second/sum;
                        vecGRev[i][j].second = (weight > min_value)?weight:min_value;
                    }
                    sort(vecGRev[i].begin(), vecGRev[i].end(), greater_first);
                }
            }
        }
        else if (probDist == WEIGHTS)
        {
            for (size_t i = 0; i < vecGRev.size(); i++)
            {
                sort(vecGRev[i].begin(), vecGRev[i].end(), greater_first);
            }
        }

        std::cout << "probability distribution: " << probDist << std::endl;
        TIO::SaveGraphStruct(filename, vecGRev, true);
        TIO::SaveGraphProbDist(filename, (int)probDist);
        TIO::SaveGraphStruct(filename, vecG, false);

        if(commFilenames != "")
        {
            NodeComms vecComms;
            CommLists vecCommLists;
            CommFairList commFairList;
            double alpha = 0.5;
            int32_t commID = -1;
            std::string fairDistFilename;
            std::string commFilename;
            std::ifstream infile(commFilenames);
            int32_t commNum = -1;
            int32_t metricNum = -1;
            if (!infile.is_open())
            {
                std::cout << "The file \"" + commFilenames + "\" can NOT be opened\n";
                return;
            }
            infile >> metricNum;
            vecComms.resize(metricNum, NodeComm(numV, -1));
            vecCommLists.resize(metricNum);
            commFairList.resize(metricNum);
            for(auto j = 0; j < metricNum; j++) // the j-th community partition metric
            {
                infile >> commFilename >> fairDistFilename >> commNum;
                std::ifstream infileComm(commFilename);
                if (!infileComm.is_open())
                {
                    std::cout << "The sub-file \"" + commFilename + "\" can NOT be opened\n";
                    return;
                }
                std::ifstream infileFair(fairDistFilename);
                if (!infileFair.is_open())
                {
                    std::cout << "The sub-file \"" + fairDistFilename + "\" can NOT be opened\n";
                    return;
                }
                vecCommLists[j].resize(commNum);
                for(auto i = 0; i < numV; i++) // the i-th node
                {
                    infileComm >> commID; // the commID-th community
                    if (commID >= 0 && commID < commNum)
                    {
                        vecComms[j][i] = commID;
                        vecCommLists[j][commID].push_back(i);
                    }
                    else
                    {
                        std::cerr << "Error: commID " << commID << " is out of bounds for commNum " << commNum << std::endl;
                        return;
                    }
                }
                infileComm.close();
                for(auto i = 0; i < commNum; i++)
                {
                    infileFair >> alpha;
                    commFairList[j].push_back(alpha); // the j-th metric, i-th community's fairness parameter is alpha
                }
            }
            infile.close();
            TIO::SaveCommStruct(filename, vecComms, vecCommLists, commFairList);
            std::cout << "The community is formatted!" << std::endl;
        }
        std::cout << "The graph is formatted!" << std::endl;
    }

    /// Load graph via vector deserialization.
    static void LoadGraph(Graph &graph, const std::string graphName, bool isReverse = true)
    {
        TIO::LoadGraphStruct(graphName, graph, isReverse);
        return ;
    }

    static int LoadGraphProbDist(const std::string graphName)
    {
        return TIO::LoadGraphProbDist(graphName);
    }

    /// Load community structure via vector deserialization.
    //vecComms: metricID->nodeID->communityID
    //vecCommLists: metricID->communityID->nodeIDs
    static void LoadComm(NodeComms &vecComms, CommLists &vecCommLists, CommFairList &commFairList, const std::string graphName)
    {
        TIO::LoadCommStruct(graphName, vecComms, vecCommLists, commFairList);
        return ;
    }
};
