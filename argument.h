#pragma once


class Argument
{
public:
    // Function parameter.
    // format: format graph
    // im: influence maximization
    // fim: fair influence maximization
    // rfim: robust fair influence maximization
    std::string _funcStr = "rfim";
    FuncType _func = RFIM;

    // Method type
    // sgh: saturate greedy with hist
    // sg: single greedy with hist
    // ag: all greedy with hist
    std::string _methodStr = "sgh";
    MethodType _methodType = METHOD_ERROR;

    // Fair type
    // 1: Linear fair function
    // 0: Welfair fair function (currently not support for Welfair fair function)
    bool _isLinearFair = true;

    // The number of nodes to be selected. Default is 100.
    int _seedsize = 100;

    // For the uniform setting, every edge has the same diffusion probability.
    float _probEdge = float(0.1);

    // Error threshold 1-1/e-epsilon.
    double _eps = 0.1;

    // Precision gamma.
    double _gamma = 0.1;

    // Failure probability delta. Default is 1/#nodes.
    double _delta = -1.0;

    // Saturate ratio
    double _saturateRatio = -1.0;

    // Graph name. Default is "facebook".
    std::string _graphname = "facebook";

    // Probability distribution
    // weights: graph data with weights
    // wc: wc setting
    // uniform: uniform setting
    // skewed: skewed distribution
    std::string _probDistStr = "wc";
    ProbDist _probDist = WC;
    std::string _skewType = "exp";

    // Directory
    std::string _dir = "graphInfo";

    // Result folder
    std::string _resultFolder = "result";

    // File name of the result
    std::string _outFileName;

    // wc variant
    double _wcVar = 1.0;

    // sample RR set with the vanilla method
    bool _vanilla = true;

    Argument(int argc, char* argv[])
    {
        std::string param, value;

        for (int ind = 1; ind < argc; ind++)
        {
            if (argv[ind][0] != '-') break;

            std::stringstream sstr(argv[ind]);
            getline(sstr, param, '=');
            getline(sstr, value, '=');

            if (!param.compare("-func")) _funcStr = value;
            else if (!param.compare("-seedsize")) _seedsize = stoi(value);
            else if (!param.compare("-eps")) _eps = stod(value);
            else if (!param.compare("-delta")) _delta = stod(value);
            else if (!param.compare("-gamma")) _gamma = stod(value);
            else if (!param.compare("-graphname")) _graphname = value;
            else if (!param.compare("-dir")) _dir = value;
            else if (!param.compare("-outpath")) _resultFolder = value;
            else if (!param.compare("-pdist")) _probDistStr = value;
            else if (!param.compare("-pedge")) _probEdge = stof(value);
            else if (!param.compare("-wcvariant")) _wcVar = stod(value);
            else if (!param.compare("-skew")) _skewType = value;
            else if (!param.compare("-vanilla")) _vanilla = (value == "1");
            else if (!param.compare("-method")) _methodStr = value;
            else if (!param.compare("-sratio")) _saturateRatio = stod(value);
            else if (!param.compare("-linearfair")) _isLinearFair = (value == "1");
        }

        if (_wcVar <= 0)
        {
            //wrong input
            _wcVar = 1.0;
        }

        decode_func_type();
        decode_prob_dist();
        decode_method_type();
        if (_resultFolder == "result")
        {
            _resultFolder += "/" + _graphname;
        }
    }

    void build_outfilename(int seedSize, ProbDist dist, Graph& graph)
    {
        std::string distStr; 

        if (dist == WEIGHTS)
        {
            _probDistStr = "weights";
        }
        else if (dist == WC)
        {
            _probDistStr = "wc";
        }
        else if (dist == UNIFORM)
        {
            _probDistStr = "uniform";

            for (int i = 0; i < graph.size(); i++)
            {
                if (graph[i].size() > 0 )
                {
                    _probEdge = graph[i][0].second;
                    break;
                }
            }
        }
        else
        {
            _probDistStr = "skewed";
        }

        switch (_func)
        {
            case IM:
                _outFileName = TIO::BuildOutFileName(_graphname, "im", seedSize, _probDistStr, _probEdge);
                break;
            case FIM:
                _outFileName = TIO::BuildOutFileName(_graphname, "fim", seedSize, _probDistStr, _probEdge);
                break;
            case RFIM:
                switch (_methodType)
                {
                    case SATURATE_GREEDY_HIST:
                        _outFileName = TIO::BuildOutFileName(_graphname, "SHIST", seedSize, _probDistStr, _probEdge, _saturateRatio, _eps, _delta);
                        break;
                    case SINGLE_GREEDY:
                        _outFileName = TIO::BuildOutFileName(_graphname, "SingleGreedy", seedSize, _probDistStr, _probEdge, _saturateRatio, _eps, _delta);
                        break;
                    case ALL_GREEDY:
                        _outFileName = TIO::BuildOutFileName(_graphname, "AllGreedy", seedSize, _probDistStr, _probEdge, _saturateRatio, _eps, _delta);
                        break;
                    case SATURATE_GREEDY_CELF:
                        _outFileName = TIO::BuildOutFileName(_graphname, "SGwMC", seedSize, _probDistStr, _probEdge, _saturateRatio, _eps, _delta);
                        break;
                    default:
                        std::cerr<< "Error: Unknown method type" << std::endl;
                        break;
                }
                break;
            default:
                _outFileName = TIO::BuildOutFileName(_graphname, "im", seedSize, _probDistStr, _probEdge);
                break;
        }

        if (!_isLinearFair) _outFileName += "_welfair";

        return ;
    }

    void decode_prob_dist()
    {
        if (_probDistStr == "wc")
        {
            _probDist = WC;
        }
        else if (_probDistStr == "uniform")
        {
            _probDist = UNIFORM;
        }
        else if (_probDistStr == "skewed")
        {
            _probDist = SKEWED;
        }
        else if (_probDistStr == "weights")
        {
            _probDist = WEIGHTS;
        }
        else 
        {
            _probDist = PROB_DIST_ERROR;
        }
    }

    void decode_func_type()
    {
        if (_funcStr == "format")
        {
            _func = FORMAT;
        }
        else if (_funcStr == "im")
        {
            _func = IM;
        }
        else if (_funcStr == "fim")
        {
            _func = FIM;
        }
        else if (_funcStr == "rfim")
        {
            _func = RFIM;
        }
        else
        {
            _func = FUNC_ERROR;
        }
    }

    void decode_method_type()
    {
        if (_methodStr == "sgh")
        {
            _methodType = SATURATE_GREEDY_HIST;
        }
        else if (_methodStr == "sg")
        {
            _methodType = SINGLE_GREEDY;
        }
        else if (_methodStr == "ag")
        {
            _methodType = ALL_GREEDY;
        }
        else if (_methodStr == "sgc") // Currently, this code not support this baseline
        {
            _methodType = SATURATE_GREEDY_CELF;
        }
        else
        {
            _methodType = METHOD_ERROR;
        }
    }
};

using TArgument = Argument;
using PArgument = std::shared_ptr<TArgument>;
