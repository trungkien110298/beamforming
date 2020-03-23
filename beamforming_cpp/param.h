// Define class Param
#ifndef PARAM_H
#define PARAM_H

#include <Eigen/Dense>
using namespace Eigen;

class Param
{
public:
    long long rate;
    long long *freqRang;
    double tNoiseBeg;
    double tNoiseEnd;
    long long *fftsize;
    VectorXd awin;
    long long siglen;

    Param(long long _rate, long long *_freqRang, double _tNoiseBeg, double _tNoiseEnd, long long *_fftsize)
    {
        rate = _rate;
        freqRang = _freqRang;
        tNoiseBeg = _tNoiseBeg;
        tNoiseEnd = _tNoiseEnd;
        fftsize = _fftsize;
    };
};

#endif