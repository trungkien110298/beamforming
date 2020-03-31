#include <iostream>
#include <stdio.h>
#include <eigen3/Eigen>
#include "beamforming.h"


using namespace std;
using namespace Eigen;


int main(int argc, char *argv[]){
    // Input
    long long num_samples = argv[1];
    long long num_channels = argv[2];
    long long sample_rate = argv[3];
    long long num_frames = argv[4];

    MatrixXd audio(num_samples, num_channels);
    MatrixXi vad(num_frames, num_channels);

    for (int i = 0; i < num_channels; i++){
        for (int j = 0; j < num_samples; j++){
            cin >> audio(j, i);
        }
    }
    for (int i = 0; i < num_channels; i++){
        for (int j = 0; j < num_frames; j++){
            cin >> vad(j, i);
        }
    }

    
    MatrixXd audio_enhenced = beamforming(audio, vad, sample_rate);

    // Output

    for (int i = 0; i < num_channels; i++){
        for (int j = 0; j < num_samples; j++){
            cout << audio_enhenced(j, i) << " ";
        }
        cout << endl;
    }
    
    
    return 0;
}    