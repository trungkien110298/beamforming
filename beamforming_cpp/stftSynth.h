#include <math.h>
#include <tuple>
#include <string>
#include <limits>
#include <Eigen/Core>
#include <unsupported/Eigen/FFT>

#include "Complex4DMatrix.h"
#include "Double3DMatrix.h"
#include "param.h"

using namespace std;
using namespace Eigen;

FFT<double> fft;

VectorXd local_istft(MatrixXcd Y, VectorXd swin, long long shift, long long siglen)
{
    long long frame_size = swin.size();
    long long num_freq = Y.rows();
    long long num_frame = Y.cols();

    VectorXd y = VectorXd::Zero(shift * (num_frame - 1) + frame_size);

    long long begin = 0;
    for (long long i = 0; i < num_frame; i++)
    {
        VectorXcd rv = Y.col(i).segment(1, num_freq - 2).reverse().conjugate();
        VectorXcd fftS(frame_size);

        fftS << Y.col(i), rv;
        VectorXd winsig(frame_size);
        fft.inv(winsig, fftS);
        y.segment(begin, frame_size) += winsig;
        begin = begin + shift;
    }
    return y.segment(0, siglen);
}

VectorXd local_synthwin(Param param, int center, int Nt)
{

    // DUET_SDW - Standard dual (biorthogonal) window for a frame.
    //
    // [swin,awin,A,B]=duet_sdw(kmax,timestep,numfreq,win,center,Nt)
    //

    int kmax = 1;
    VectorXd win = param.awin;
    long long timestep = param.fftsize[1];
    long long num_freq = param.fftsize[0];

    // KMAX * NUMFREQ is the size of returned dual window.
    // TIMESTEP is the # of samples between adjacent windows in time.
    // NUMFREQ is the # of frequency components per timestep.
    // WIN is the frame window.
    // CENTER [0], if 1, centers the input win for use as awin.
    // NT [256] is the number of t points in the (0,1) Zak domain.
    //
    // [*] denotes optional argument default value.
    //
    // SWIN is the (truncated) standard dual.
    // AWIN is win plus zero-padding to have length(swin).
    // A is the lower frame bound.
    // B is the upper frame bound.
    //
    // 24 July 2000 R. Balan & S. Rickard
    // (c) Siemens Corporate Research

    long long q = num_freq / timestep; // redundancy factor
    long long M = num_freq;
    long long N = timestep;

    MatrixXd g(M * Nt, 1);
    double numz = (1.0 * (kmax * num_freq) - win.size()) / 2;

    if (center && numz != 0)
    {
        VectorXd sub_1 = VectorXd::Zero(floor(numz));
        VectorXd sub_2 = win;
        VectorXd sub_3 = VectorXd::Zero(ceil(numz));
        VectorXd sub_4 = VectorXd::Zero(M * Nt - kmax * num_freq);
        VectorXd temp(M * Nt);
        temp << sub_1, sub_2, sub_3, sub_4;
        g.col(0) = temp;
    }
    else
    {
        VectorXd sub_1 = win;
        VectorXd sub_2 = VectorXd::Zero(M * Nt - win.size());
        VectorXd temp(M * Nt);

        temp << sub_1, sub_2;
        g.col(0) = temp;
    }

    //  zak
    MatrixXd g_reshape = (MatrixXd::Map(g.data(), M, Nt)).transpose();

    MatrixXd G(Nt, M);
    for (long long i = 0; i < M; i++)
    {
        VectorXd ifft_g(Nt);
        VectorXcd g_col = g_reshape.col(i);
        fft.inv(ifft_g, g_col);
        G.col(i) = ifft_g;
    }
    //fft.inv(G, g_reshape);
    G *= sqrt(M) * Nt;

    // // rescale the zak to create the standard dual zak
    VectorXi cind = VectorXi::LinSpaced((M - 1) / N + 1, 1, 1 + N * ((M - 1) / N));
    MatrixXd Gd = MatrixXd::Zero(Nt, M);
    double A = numeric_limits<double>::max();
    double B = 0;

    for (long long i = 0; i < N; i++)
    {
        VectorXd Gfactor = VectorXd::Zero(Nt);
        for (long long j = 0; j < Nt; j++)
        {
            for (long long k = 0; k < cind.size(); k++)
            {
                Gfactor(j) += pow(G(j, cind(k) + i - 1), 2);
            }
        }
        A = min(A, Gfactor.minCoeff());
        B = max(B, Gfactor.maxCoeff());

        // If we get divide by zero warnings, the window and
        // parameters do not form a frame. Nevertheless, we
        // might want to put the below line back in to avoid
        // the warnings. For now, leave it out.
        // Gfactor(Gfactor<eps) = 1;
        for (long long j = 0; j < q; j++)
        {
            Gd.col(j * N + i) = G.col(j * N + i).cwiseQuotient(Gfactor);
        }
    }

    MatrixXcd complex_swin(Nt, M);
    for (long long i = 0; i < M; i++)
    {
        VectorXcd fft_Gd(Nt);
        VectorXd Gd_col = Gd.col(i);
        fft.fwd(fft_Gd, Gd_col);
        complex_swin.col(i) = fft_Gd;
    }
    //fft.fwd(complex_swin, Gd);
    MatrixXd swin = complex_swin.real();
    VectorXd result = swin.row(0);
    result *= num_freq / (sqrt(M) * Nt);
    //VectorXd result = swin_reshape.col(0);
    return result;
}

tuple<Double3DMatrix, MatrixXd> stftSynth(Complex4DMatrix Y, Param param, int num_channels, int num_channels_out, long long num_frame, long long num_freq)
{
    // Input
    // Y: num_freq, num_channels_out, num_channels, num_frame

    long long nfft = param.fftsize[0];
    long long shift = param.fftsize[1];

    string swinsel = "dual";
    long long siglen = param.siglen;
    VectorXd swin(nfft);

    if (swinsel == "rect")
    {
        long long mult = nfft / (2 * shift);
        swin = VectorXd::Ones(nfft, 1) / mult;
    }
    else
    {
        swin = local_synthwin(param, 1, 2);
    }

    long long size[3] = {2, siglen, 2};
    Double3DMatrix y = init_Double3DMatrix(size);

    for (int i = 0; i < num_channels_out; i++)
    {
        for (int j = 0; j < num_channels; j++)
        {
            MatrixXcd Yji(num_freq, num_frame);
            for (long long f = 0; f < num_freq; f++)
            {
                Yji.row(f) = (*Y[f][i]).row(j);
            }

            (*y[i]).col(j) = local_istft(Yji, swin, shift, siglen);
        }
    }
    return make_tuple(y, swin);
}
