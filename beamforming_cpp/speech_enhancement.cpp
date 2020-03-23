#include <stdio.h>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/FFT>
#include <math.h>
#include <complex>
#include <iostream>
#include <vector>
#include <tuple>
#include <chrono>
#include <sndfile.hh>

#include "param.h"
#include "VAD.h"
#include "stftAnalyFull.h"
#include "stftSynth.h"
#include "Double3DMatrix.h"
#include "MVDRbeamf.h"
#include "scaling.h"
#include "wavfile.h"

using namespace std;
using namespace Eigen;
using namespace std::chrono;

#define NUM_CHANNEL 8

tuple<MatrixXd, MatrixXd> SE(MatrixXd &s, long long fs, double _noise_beg, double _noise_end, bool subnoisecov)
{
    long long num_samples = s.rows();
    int num_channels = s.cols();

    MatrixXd mvdr_out = MatrixXd::Zero(num_samples, 1);
    MatrixXd mvdr_derev_out = MatrixXd::Zero(num_samples, 1);

    long long freq_range[2] = {0, fs / 2};
    // long long fftsize[2] = {floor(pow(2.0,10)), floor(pow(2.0, 8))}

    long long fftsize[2] = {1024, 256};
    Param param(fs, freq_range, _noise_beg, _noise_end, fftsize);

    long long noise_beg = floor(param.rate * param.tNoiseBeg / param.fftsize[1]);
    long long noise_end = floor(param.rate * param.tNoiseEnd / param.fftsize[1]);

    //------------------------stftAnalyFull-----------------------------//
    auto start = high_resolution_clock::now();

    auto XnParam = stftAnalyFull(s, param);

    auto stop_sAF = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop_sAF - start);
    cout << "stftAnalyFull duration: " << duration.count() / 1000000.0 << endl;

    //-----------------------------------------------------//

    param = get<2>(XnParam);
    Complex3DMatrix X = get<0>(XnParam);
    Complex4DMatrix Xfull = get<1>(XnParam);
    
    //X's size (nMic, num_frame, num_freq)
    //Xfull size (num_frame, num_freq, nMic, nMic)

    long long num_frame = (*X[0]).rows();
    long long num_freq = (*X[0]).cols();

    long long size[3] = {num_freq, num_channels, num_channels};
    Complex3DMatrix PhiNN = init_Complex3DMatrix(size);
    Complex3DMatrix PhiXX = init_Complex3DMatrix(size);

    Complex3DMatrix VN = init_Complex3DMatrix(size);
    Complex3DMatrix VX = init_Complex3DMatrix(size);

    MatrixXcd DN(num_freq, num_channels);
    MatrixXcd DX(num_freq, num_channels);

    MatrixXcd VN_vec(num_freq, num_channels);
    MatrixXcd VX_vec(num_freq, num_channels);

    //Extract noise covariance matrix
    SelfAdjointEigenSolver<MatrixXcd> ces;

    for (int i = 0; i < num_freq; i++)
    {
        //----------- PhiNN ---------------//
        long long count = 0;
        for (int j = 0; j < noise_beg; j++)
        {
            *PhiNN[i] = *PhiNN[i] + (*Xfull[j][i]);
            count++;
        }

        //For noise period
        for (int j = num_frame - noise_end; j < num_frame; j++)
        {
            *PhiNN[i] = *PhiNN[i] + (*Xfull[j][i]);
            count++;
        }

        //Normalize
        *PhiNN[i] /= count;
        ces.compute(*PhiNN[i]);

        DN.row(i) = ces.eigenvalues();
        *VN[i] = ces.eigenvectors();

        //pick th comlumn which has largest eigen value
        int max_index = 0;
        for (int j = 1; j < num_channels; j++)
        {
            max_index = abs(DN(i, j)) > abs(DN(i, max_index)) ? j : max_index;
        }
        VN_vec.row(i) = (*VN[i]).col(max_index);

        //----------- PhiNN ---------------//

        count = 0;
        //For speech period
        for (int j = noise_beg; j < num_frame - noise_end; j++)
        {
            *PhiXX[i] = *PhiXX[i] + (*Xfull[j][i]);
            count++;
        }

        //Normalize
        *PhiXX[i] /= count;
        if (subnoisecov)
        {
            *PhiXX[i] = *PhiXX[i] - *PhiNN[i];
        }

        ces.compute(*PhiXX[i]);
        DX.row(i) = ces.eigenvalues();
        *VX[i] = ces.eigenvectors();

        //pick th comlumn which has largest eigen value
        max_index = 0;
        for (int j = 1; j < num_channels; j++)
        {
            max_index = abs(DX(i, j)) > abs(DX(i, max_index)) ? j : max_index;
        }
        VX_vec.row(i) = (*VX[i]).col(max_index);
    }

    auto stop_cal_vec = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop_cal_vec - stop_sAF);
    cout << "Calculate VN_vec, VX_vec duration: " << duration.count() / 1000000.0 << endl;

    long long speech_noise_MVDR_size[3] = {2, num_frame, num_freq};
    Complex3DMatrix speech_noise_MVDR = init_Complex3DMatrix(speech_noise_MVDR_size);



    //-------------------extract speech noise in the one utterance----------------------//
    
    *speech_noise_MVDR[0] = MVDRbeamf(VX_vec, PhiNN, X);
    *speech_noise_MVDR[1] = MVDRbeamf(VN_vec, PhiXX, X);

    auto stop_MVDR = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop_MVDR - stop_cal_vec);
    cout << "MVRDR duration: " << duration.count() / 1000000.0 << endl;

     //-----------------------------------------------------//

    auto scaling_output = scaling(X, speech_noise_MVDR, param, num_channels, 2);

    // A's size:      (num_freq, num_channels, num_channels_out)
    // Y_MVDR's size: (num_freq, num_channels_out, num_channels, num_frame)
    Complex3DMatrix A = get<0>(scaling_output);
    Complex4DMatrix Y_MVDR = get<1>(scaling_output);


    auto stop_scaling = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop_scaling - stop_MVDR);
    cout << "Scaling duration: " << duration.count() / 1000000.0 << endl;


    //-----------------------------------------------------//

    auto stftSynth_out = stftSynth(Y_MVDR, param, 2, 2, num_frame, num_freq);
    
    
    auto stop_sS = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop_sS - stop_scaling);
    cout << "stftSynth duration: " << duration.count() / 1000000.0 << endl;


    //-----------------------------------------------------//
    
    
    
    Double3DMatrix yt_mvdr = get<0>(stftSynth_out);
    mvdr_out.col(0) = (*yt_mvdr[0]).col(0);    
    return make_tuple(mvdr_out, mvdr_derev_out);
}

int main()
{
    // Input 8 channel audio data

    auto start = high_resolution_clock::now();

    int safe_region = 10;
    bool subnoisecov = true;
    bool derev = true;

    char file_name[] = "/home/kienpt/Documents/Beam/data/speaker0250-0020.wav";

    // -------- Read file ----------//

    SF_INFO sfinfo;
    SNDFILE *f = sf_open(file_name, SFM_READ, &sfinfo);
    MatrixXd audio = readfile(f, sfinfo);
    long long num_samples = sfinfo.frames;
    long long sample_rate = sfinfo.samplerate;
    int num_channels = sfinfo.channels;
    sf_close(f);

    // -------- Finished read file  ----------//

    double noise_begin = 100;
    double noise_end = 0;
    double threshold = 0.1;
    double win_dur = 0.05;
    double hop_dur = 0.025;
    int num_noise = 20;
    int argin = 1;

    MatrixXd sample(num_samples, 1);

    for (int i = 0; i < num_channels; i++)
    {
        sample.col(0) = audio.col(i);

        MatrixXi result = VAD(sample, sample_rate, threshold, win_dur, hop_dur, num_noise, argin);
        int left = 0;

        while (result(left, 0) != 1)
        {
            left++;
        }
        double left_side = max((left - safe_region + 1) * 0.025, 0.0);

        int right = result.rows() - 1;
        while (result(right, 0) != 1)
        {
            right--;
        }
        double right_side = min((right + safe_region + 1) * 0.025, num_samples * (1.0 / sample_rate));

        noise_begin = min(noise_begin, left_side);
        noise_end = max(noise_end, right_side);
    }

    auto stop_VAD = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop_VAD - start);
    cout << "VAD duration: " << duration.count() / 1000000.0 << endl;

    noise_end = num_samples * (1.0 / sample_rate) - noise_end;
    auto mvdr = SE(audio, sample_rate, noise_begin, noise_end, subnoisecov);
    MatrixXd mvdr_out = get<0>(mvdr);

    // -------- Write file ----------//

    char file_name_out[] = "/home/kienpt/Documents/Beam/data/speaker0250-0020_cpp_out.wav";
    SF_INFO sfinfo_out = sfinfo;
    sfinfo_out.channels = mvdr_out.cols();

    SNDFILE *f_out = sf_open(file_name_out, SFM_WRITE, &sfinfo_out);
    writefile(f_out, sfinfo_out, mvdr_out);
    sf_close(f_out);

    // -------- Finished write file  ----------//

    auto stop_SE = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop_SE - stop_VAD);
    cout << "SE duration: " << duration.count() / 1000000.0 << endl;

    return 0;
}
