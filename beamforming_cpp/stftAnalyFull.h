#include <iostream>
#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>
#include <tuple>
#include <math.h>
#include "Complex4DMatrix.h"
#include "Complex3DMatrix.h"
#include "param.h"
#include <complex>
#include <chrono>

using namespace std;
using namespace Eigen;
using namespace std::chrono;

VectorXd hanning(long long n)
{
    VectorXd result(n);
    for (long long i = 0; i < n; i++)
    {
        result(i) = 0.5 * (1 - cos(2 * M_PI * i / n));
    }
    return result;
}

MatrixXcd local_stft(VectorXd s, VectorXd awin, long long shift)
{
    //   short-time Fourier transform
    //   s: column vector
    //   S: matrix of num_freq * num_frame, num_freq = length(awin)/2 +1
    //   Sfull: matrix of num_freq * num_freq * num_frame
    FFT<double> fft;
    long long frameSize = awin.size();
    long long num_freq = frameSize / 2 + 1;

    long long sLength = s.size();
    long long num_frame = ceil(1.0 * (sLength - frameSize) / shift) + 1;
    num_frame = (num_frame > 1) ? num_frame : 1;

    long long newlen = frameSize + (num_frame - 1) * shift;
    VectorXd new_s(newlen);
    new_s << s, VectorXd::Zero(newlen - sLength);
    MatrixXcd S = MatrixXd::Zero(num_freq, num_frame);

    long long begin = 0;

    for (long long i = 0; i < num_frame; i++)
    {
        VectorXd winsig = new_s.segment(begin, frameSize).cwiseProduct(awin);
        VectorXcd fftS(frameSize);
        fft.fwd(fftS, winsig);
        S.col(i) = fftS.segment(0, num_freq);
        begin += shift;
    }

    return S;
}

tuple<Complex3DMatrix, Complex4DMatrix, Param> stftAnalyFull(MatrixXd &x, Param &param)
{

    //   stftAnaly: STFT analysis for x
    //
    //   [X, param] = stftAnaly(x, param)
    //
    //   x: cell, x{i}: array, sigLen x num_channels
    //   param: parameters
    //       param.awinsel: analysis window, 'hann' or 'sqrthann' (default)
    //   X: matrix of num_channels * num_frame * num_freq, num_freq = nfft/2 +1
    //   X_full: matrix of num_frame * num_freq * num_channels * num_channels

    long long nfft = param.fftsize[0];
    long long shift = param.fftsize[1];
    long long siglen = x.rows();
    int num_channels = x.cols();

    /* MATLAB CODE
        if isfield(param, 'awinsel')
            awinsel = param.awinsel;
        else
            awinsel = 'sqrthann';
        end

        switch awinsel
            case 'hann',
                awin = hanning(nfft, 'periodic');
            case 'sqrthann',
                awin = sqrt(hanning(nfft, 'periodic'));
        end
    */

    VectorXd awin = hanning(nfft).cwiseSqrt();
    param.awin = awin;
    param.siglen = siglen;
    long long num_frame = ceil(1.0 * (param.siglen - nfft) / shift) + 1;
    long long num_freq = nfft / 2 + 1;

    long long X_size[3] = {num_channels, num_frame, num_freq};
    Complex3DMatrix X = init_Complex3DMatrix(X_size);

    long long Xfull_size[4] = {num_frame, num_freq, num_channels, num_channels};
    Complex4DMatrix Xfull = init_Complex4DMatrix(Xfull_size);

    for (int mic = 0; mic < num_channels; mic++)
    {

        *X[mic] = local_stft(x.col(mic), param.awin, shift).transpose();
    }

    complex<double> sum_X = sum_Complex3DMatrix(X, X_size);

    for (long long k = 0; k < num_frame; k++)
    {
        for (long long f = 0; f < num_freq; f++)
        {
            MatrixXcd temp(num_channels, 1);
            for (int i = 0; i < num_channels; i++)
            {
                temp(i, 0) = (*X[i])(k, f);
            }
            // xff = temp * temp.transpose();
            *Xfull[k][f] = temp * temp.adjoint();
        }
    }

    return make_tuple(X, Xfull, param);
}