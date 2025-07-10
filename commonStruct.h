#pragma once
/// Node list
typedef std::vector<uint32_t> Nodelist;
/// Edge structure, neighbor id and the edge weight
typedef std::pair<uint32_t, float> Edge;
/// Edgelist structure from one source/target node
typedef std::vector<Edge> Edgelist;
/// Graph structure
typedef std::vector<Edgelist> Graph;
/// One forward reachable set
typedef std::vector<size_t> FRset;
/// A set of forward reachable sets
typedef std::vector<FRset> FRsets;
double decay_factor = 1.0;
/// One reverse reachable set
typedef std::vector<uint32_t> RRset;
/// A set of reverse reachable sets
typedef std::vector<RRset> RRsets;
bool optflag;
enum ProbDist {WEIGHTS, UNIFORM, WC, SKEWED, PROB_DIST_ERROR};
enum FuncType {FORMAT, IM, FIM, RFIM, FUNC_ERROR};
/// Different fairness calculation methods
enum FairType {LINEAR, SCALE_POWER, FAIR_ERROR};
/// Different baseline methods
enum MethodType {SATURATE_GREEDY_HIST, SINGLE_GREEDY, ALL_GREEDY, SATURATE_GREEDY_CELF, METHOD_ERROR};

enum EdgeStatus { UNCHECK, PASS, FAIL };
typedef std::vector<EdgeStatus> EdgeStatusList;
typedef std::vector<EdgeStatusList> GraphEdgeStatus;

/// New structure for RFIM
/// Community structure (a vector of the community id for each node)
typedef std::vector<int32_t> NodeComm;
/// A set of community partition metrics
typedef std::vector<NodeComm> NodeComms;

/// Community list (a vector of the vector of nodes' id for each community)
typedef std::vector<Nodelist> CommList;
/// A set of "CommList" under different community partition metrics
typedef std::vector<CommList> CommLists;

/// A set of fair parameter alpha for each community under one metric
typedef std::vector<double> CommFair;
/// A set of fair parameters under different community partition metrics
typedef std::vector<CommFair> CommFairList;

/// Node element with id and a property value
typedef struct NodeElement
{
    int id;
    double value;
} NodeEleType;

/// Smaller operation for node element
struct smaller
{
    bool operator()(const NodeEleType& Ele1, const NodeEleType& Ele2) const
    {
        return (Ele1.value < Ele2.value);
    }
};

/// Greater operation for node element
struct greater
{
    bool operator()(const NodeEleType& Ele1, const NodeEleType& Ele2) const
    {
        return (Ele1.value > Ele2.value);
    }
};
