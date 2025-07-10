#pragma once

class ResultInfo
{
private:
    double __RunningTime = -1.0;
    double __Influence = -1.0;
    double __InfluenceOriginal = -1.0;
    double __Approx = -1.0;
    int __SeedSize = 0;
    size_t __RRsetsSize = 0;
    Nodelist __VecSeed;

    double __ObjectiveMC = -1.0;

    // member for FIM and RFIM
    double __FairInfluence = -1.0;
    double __FairInfluenceOriginal = -1.0;
    double __RobustnessScore = -1.0;
    double __RobustnessScoreOriginal = -1.0;

    double __BestInfluence = -1.0;
    double __BestInfluenceOriginal = -1.0;
    double __BestApprox = -1.0;
    int __BestSeedSize = 0;
    size_t __BestRRsetsSize = 0;
    Nodelist __BestVecSeed;
    double __BestFairInfluence = -1.0;
    double __BestFairInfluenceOriginal = -1.0;
    double __BestRobustnessScore = -1.0;
    double __BestRobustnessScoreOriginal = -1.0;

public:
    ResultInfo()
    {
    }

    ~ResultInfo()
    {
    }

    /// Get running time.
    double get_running_time() const
    {
        return __RunningTime;
    }

    /// Get influence spread.
    double get_influence() const
    {
        return __Influence;
    }

    /// Get self-estimated influence spread.
    double get_influence_original() const
    {
        return __InfluenceOriginal;
    }

    /// Get approximation guarantee.
    double get_approximation() const
    {
        return __Approx;
    }

    /// Get seed sets.
    Nodelist get_seed_vec() const
    {
        return __VecSeed;
    }

    /// Get seed size.
    int get_seed_size() const
    {
        return __SeedSize;
    }

    /// Get the number of RR sets.
    size_t get_RRsets_size() const
    {
        return __RRsetsSize;
    }

    /// Set running time.
    void set_running_time(const double value)
    {
        __RunningTime = value;
    }

    /// Set influence spread.
    void set_influence(const double value)
    {
        __Influence = value;
    }

    /// Set self-estimated influence spread
    void set_influence_original(const double value)
    {
        __InfluenceOriginal = value;
    }

    /// Set approximation guarantee.
    void set_approximation(const double value)
    {
        __Approx = value;
    }

    /// Set seed sets.
    void set_seed_vec(Nodelist& vecSeed)
    {
        __VecSeed = vecSeed;
        set_seed_size((int)vecSeed.size());
    }

    /// Set seed size.
    void set_seed_size(const int value)
    {
        __SeedSize = value;
    }

    /// Set the number of RR sets.
    void set_RR_sets_size(const size_t value)
    {
        __RRsetsSize = value;
    }

    void set_fair_inf(const double value)
    {
        __FairInfluence = value;
    }

    double get_fair_inf() const
    {
        return __FairInfluence;
    }

    void set_fair_inf_ori(const double value)
    {
        __FairInfluenceOriginal = value;
    }

    double get_fair_inf_ori() const
    {
        return __FairInfluenceOriginal;
    }

    void set_robustness_score(const double value)
    {
        __RobustnessScore = value;
    }

    double get_robustness_score() const
    {
        return __RobustnessScore;
    }

    void set_robustness_score_ori(const double value)
    {
        __RobustnessScoreOriginal = value;
    }

    double get_robustness_score_ori() const
    {
        return __RobustnessScoreOriginal;
    }

    void set_resMC(const double value)
    {
        __ObjectiveMC = value;
    } 

    double get_resMC() const
    {
        return __ObjectiveMC;
    }

    void log_best_res() // not include ObjectiveMC
    {
        __BestInfluence = __Influence;
        __BestInfluenceOriginal = __InfluenceOriginal;
        __BestApprox = __Approx;
        __BestSeedSize = __SeedSize;
        __BestRRsetsSize = __RRsetsSize;
        __BestVecSeed = __VecSeed;
        __BestFairInfluence = __FairInfluence;
        __BestFairInfluenceOriginal = __FairInfluenceOriginal;
        __BestRobustnessScore = __RobustnessScore;
        __BestRobustnessScoreOriginal = __RobustnessScoreOriginal;
    }

    void log_final_res()
    {
        __Influence = __BestInfluence;
        __InfluenceOriginal = __BestInfluenceOriginal;
        __Approx = __BestApprox;
        __SeedSize = __BestSeedSize;
        __RRsetsSize = __BestRRsetsSize;
        __VecSeed = __BestVecSeed;
        __FairInfluence = __BestFairInfluence;
        __FairInfluenceOriginal = __BestFairInfluenceOriginal;
        __RobustnessScore = __BestRobustnessScore;
        __RobustnessScoreOriginal = __BestRobustnessScoreOriginal;
    }
};

typedef ResultInfo TResult;
typedef std::shared_ptr<ResultInfo> PResult;
