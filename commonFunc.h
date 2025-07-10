#pragma once

#define ENUM_STR(x) #x

/// Log information
template <typename _Ty>
static inline void LogInfo(_Ty val)
{
    std::cout << val << std::endl;
}

/// Log information
template <typename _Ty>
static inline void LogInfo(const std::string title, _Ty val)
{
    std::cout << title << ": " << val << std::endl;
}

/// Math, pow2
static inline double pow2(const double t)
{
    return t * t;
}

/// Math, log2
static inline double log2(const size_t n)
{
    return log(n) / log(2);
}

/// Math, logcnk
static inline double logcnk(const size_t n, size_t k)
{
    k = k < n - k ? k : n - k;
    double res = 0;

    for (auto i = 1; i <= k; i++) res += log(double(n - k + i) / i);

    return res;
}

/// Make the vector to a min-heap.
inline void MakeMinHeap(FRset& vec)
{
    // Min heap
    const auto size = vec.size();

    if (2 <= size)
    {
        for (auto hole = (size + 1) / 2; hole--;)
        {
            const auto val = vec[hole];
            size_t i, child;

            for (i = hole; i * 2 + 1 < size; i = child)
            {
                // Find smaller child
                child = i * 2 + 2;

                if (child == size || vec[child - 1] < vec[child])
                {
                    // One child only or the left child is smaller than the right one
                    --child;
                }

                // Percolate one level
                if (vec[child] < val)
                {
                    vec[i] = vec[child];
                }
                else
                {
                    break;
                }
            }

            vec[i] = val;
        }
    }
}

/// Replace the value for the first element and down-heap this element.
inline void MinHeapReplaceMinValue(FRset& vec, const size_t& val)
{
    // Increase the value of the first element
    const auto size = vec.size();
    size_t i, child;

    for (i = 0; i * 2 + 1 < size; i = child)
    {
        // Find smaller child
        child = i * 2 + 2;

        if (child == size || vec[child - 1] < vec[child])
        {
            // One child only or the left child is smaller than the right one
            --child;
        }

        // Percolate one level
        if (vec[child] < val)
        {
            vec[i] = vec[child];
        }
        else
        {
            break;
        }
    }

    vec[i] = val;
}

/// Make the vector to a max-heap.
static inline void MakeMaxHeap(std::vector<std::pair<float, uint32_t>>& vec)
{
    // Max heap
    const auto size = vec.size();

    if (2 <= size)
    {
        for (auto hole = (size + 1) / 2; hole--;)
        {
            const auto val = vec[hole];
            size_t i, child;

            for (i = hole; i * 2 + 1 < size; i = child)
            {
                // Find smaller child
                child = i * 2 + 2;

                if (child == size || vec[child - 1] > vec[child])
                {
                    // One child only or the left child is greater than the right one
                    --child;
                }

                // Percolate one level
                if (vec[child] > val)
                {
                    vec[i] = vec[child];
                }
                else
                {
                    break;
                }
            }

            vec[i] = val;
        }
    }
}

/// Replace the value for the first element and down-heap this element.
static inline void MaxHeapReplaceMaxValue(std::vector<std::pair<float, uint32_t>>& vec, const float& val)
{
    // Increase the value of the first element
    const auto size = vec.size();
    size_t i, child;
    auto hole = vec[0];

    for (i = 0; i * 2 + 1 < size; i = child)
    {
        // Find smaller child
        child = i * 2 + 2;

        if (child == size || vec[child - 1] > vec[child])
        {
            // One child only or the left child is greater than the right one
            --child;
        }

        // Percolate one level
        if (vec[child].first > val)
        {
            vec[i] = vec[child];
        }
        else
        {
            break;
        }
    }

    hole.first = val;
    vec[i] = hole;
}

/// Generate one node with probabilities according to their weights for the LT cascade model
static inline size_t GenRandomNodeByWeightLT(const Edgelist& edges)
{
    const double weight = dsfmt_gv_genrand_open_close();
    size_t minIdx = 0, maxIdx = edges.size() - 1;

    if (weight < edges.front().second) return 0; // First element

    if (weight > edges.back().second) return edges.size() + 1; // No element

    while (maxIdx > minIdx)
    {
        const size_t meanIdx = (minIdx + maxIdx) / 2;
        const auto meanWeight = edges[meanIdx].second;

        if (weight <= meanWeight) maxIdx = meanIdx;
        else minIdx = meanIdx + 1;
    }

    return maxIdx;
}

/// Normalize the probabilities to a accumulative format, e.g., [0.2, 0.5, 0.3]->[0.2, 0.7, 1.0]
static inline void NormalizeAccumProb(Graph& vecGraph)
{
    for (auto& nbrs : vecGraph)
    {
        float accumVal = float(0.0);

        for (auto& nbr : nbrs)
        {
            accumVal += nbr.second;
            nbr.second = accumVal;
        }

        // Normalization
        for (auto& nbr : nbrs)
        {
            nbr.second /= accumVal;
        }
    }
}


bool SmallerPair(const std::pair<uint32_t, uint32_t> &x, const std::pair<uint32_t, uint32_t> &y)
{
    return x.first < y.first;
}

bool GreaterPair(const std::pair<uint32_t, uint32_t> &x, const std::pair<uint32_t, uint32_t> &y)
{
    return x.first > y.first;
}

inline double Logarithm(const double x)
{
    // log2f is 10% faster
    return log2f(x);
    //return log2(x);
    //return log(x);
}

// /* increase the sum of weights of all incoming edges to sumProb */
// void increaseWeights(Graph &graph, double sumProb)
// {
//     for (auto& nbrs : graph)
//     {
//         for (auto& nbr : nbrs)
//         {
//             nbr.second *= sumProb;
//         }
//     }
// }

/// Compute eta(n, alpha) in FIM
inline double eta(const size_t n, const double alpha)
{
    if (n == 1)
        return 1.0;

    double numerator = 1.0;
    for (size_t i = 1; i < n; ++i)
    {
        numerator *= static_cast<double>(i - alpha);
    }

    double denominator = 1.0;
    for (size_t i = 2; i <= n; ++i)
    {
        denominator *= static_cast<double>(i);
    }

    return numerator / denominator;
}

/// Compute the marginal fairness influence for SCALE_POWER model in one community and one metric
inline double computeMarginalScalePower(const size_t m1, const size_t m2,
    const double alpha, const size_t c)
{
    // m1: current number of covered RR sets (_numGCovered[q][j])
    // m2: new number of covered RR sets if v is added (_GCoverage[q][j][v])
    // alpha: fairness parameter
    // c: the total number of RR sets (_GRRsets[q][j].size() or _numRRsets)
    double res = 0.0, x1, x2 = 1.0, x3 = 1.0;
    for (auto n = 1; n < c - m1; ++n)
    {
        x1 = eta(n, alpha);
        for (auto i = 0; i < n; ++i)
        {
            x2 *= static_cast<double>(c - m1 - i) / static_cast<double>(c - i);
            x3 *= static_cast<double>(c - m1 - m2 - i) / static_cast<double>(c - i);
        }
        res += x1 * (x2 - x3);
    }
    return res;
}

/// Compute the fairness influence for SCALE_POWER model in one community and one metric
inline double computeScalePower(const size_t m1, const double alpha, const size_t c)
{
    // m1: current number of covered RR sets (_numGCovered[q][j] or _gHit[q][j])
    // alpha: fairness parameter
    // c: the total number of RR sets (_GRRsets[q][j].size() or _numRRsets)
    double res = 0.0, x1, x2 = 1.0;
    for (auto n = 1; n < c - m1; ++n)
    {
        x1 = eta(n, alpha);
        for (auto i = 0; i < n; ++i)
        {
            x2 *= static_cast<double>(c - m1 - i) / static_cast<double>(c - i);
        }
        res += x1 * x2;
    }
    return res;
}

// calcute the possible largest remaining number of iterations of Saturated Greedy in RIM.
inline int estRemIterCap(const double gamma, const double c_min, const double c_max)
{
    if (c_max - c_min <= gamma) return 0;
    else return 1 + estRemIterCap(gamma, (c_min + c_max) / 2 * (1 - gamma / 3), c_max);
}