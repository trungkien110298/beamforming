#include <Eigen/Dense>
#include <tuple>
#include <math.h>
#include "param.h"
#include "Complex3DMatrix.h"
#include "Complex4DMatrix.h"
#include "Complex3DMatrix.h"
#include "Complex4DMatrix.h"

using namespace std;
using namespace Eigen;

tuple<long long, long long> freqBinRange(Param param)
{
    long long frame_size = param.fftsize[0];

    double freq_step = 1.0*param.rate / frame_size;
    long long initF = ceil(1.0 * param.freqRang[0] / freq_step);
    long long lastF = floor(1.0 * param.freqRang[1] / freq_step);
    return make_tuple(initF, lastF);
}

tuple<Complex3DMatrix, Complex4DMatrix> scaling(Complex3DMatrix X, Complex3DMatrix Y, Param param, int num_channels, int num_channels_out)
{
    // X: sensor observations, nMic * num_frame * num_freq
    // Y: separated signals, nOut * num_frame * num_freq
    // A: output, num_freq, num_channels, num_channels_out
    // YA: output, num_freq, num_channels_out, num_channels, num_frame (if needed)

    long long num_frame = (*X[0]).rows();
    long long num_freq = (*X[0]).cols();

    long long A_size[3] = {num_freq, num_channels, num_channels_out};
    Complex3DMatrix A = init_Complex3DMatrix(A_size);

    long long YA_size[4] = {num_freq, num_channels_out, num_channels, num_frame};
    Complex4DMatrix YA = init_Complex4DMatrix(YA_size);

    auto freqBinRange_output = freqBinRange(param);
    long long initF = get<0>(freqBinRange_output);
    long long lastF = get<1>(freqBinRange_output);

    for (long long f = initF; f <= lastF; f++)
    {
        MatrixXcd x(num_channels, num_frame);
        MatrixXcd y(num_channels_out, num_frame);

        //Copy data to x, y
        for (int c = 0; c < num_channels; c++)
        {
            for (long long frame = 0; frame < num_frame; frame++)
            {
                x(c, frame) = (*X[c])(frame, f);
            }
        }
        for (int c = 0; c < num_channels_out; c++)
        {
            for (long long frame = 0; frame < num_frame; frame++)
            {
                y(c, frame) = (*Y[c])(frame, f);
            }
        }

        *A[f] = (x * y.adjoint()) * (y * y.adjoint()).inverse();
        for (int c = 0; c < num_channels_out; c++)
        {
            VectorXcd y_row = y.row(c);
            VectorXcd Af_col = (*A[f]).col(c);
            *YA[f][c] = Af_col * y_row.transpose();
        }
    }
    return make_tuple(A, YA);
}
