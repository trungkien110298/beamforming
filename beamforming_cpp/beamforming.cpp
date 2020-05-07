#include <iostream>
#include <stdio.h>
#include <Eigen/Dense>
#include <eigen3-hdf5.hpp>
#include "beamforming.h"

using namespace std;
using namespace Eigen;

int main(int argc, char *argv[])
{
    // Input
    long long num_samples;
    long long num_channels;
    long long sample_rate;
    long long num_frames;
    long long window_order;

    cin >> num_samples >> num_channels >> sample_rate >> num_frames >> window_order;
    MatrixXd audio(num_samples, num_channels);
    MatrixXi vad(num_frames, num_channels);

    for (int i = 0; i < num_samples; i++)
    {
        for (int j = 0; j < num_channels; j++)
        {
            cin >> audio(i, j);
        }
    }
    
    for (int i = 0; i < num_frames; i++)
    {
        for (int j = 0; j < num_channels; j++)
        {
            cin >> vad(i, j);
        }
    }

    MatrixXd audio_enhenced = beamforming(audio, vad, sample_rate, window_order);

    // Output

    for (int i = 0; i < num_samples; i++)
    {
        if (isnan(audio_enhenced(i, 0))){
            audio_enhenced(i, 0) = 0;
        }
        cout << audio_enhenced(i, 0) << " ";
    }

    return 0;
}