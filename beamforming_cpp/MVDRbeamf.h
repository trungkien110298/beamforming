#include <Eigen/Dense>
#include <iostream>
#include "Complex3DMatrix.h"

using namespace std;
using namespace Eigen;

MatrixXcd MVDRbeamf(MatrixXcd &FX, Complex3DMatrix PhiNN, Complex3DMatrix X)
{
    // input
    // FX = num_freq x num_mic
    // PhiXX = num_mic x num_mic x num_frame x num_freq
    // PhiNN = num_freq, num_channels, num_channels
    // X = num_mic X num_frame x num_freq

    int num_channels = FX.cols();
    long long num_frame = (*X[0]).rows();
    long long num_freq = (*X[0]).cols();

    
    MatrixXd channel_selection = MatrixXd::Zero(num_channels, 1);
    MatrixXd smooth_post_filter = MatrixXd::Zero(num_freq, num_frame);
    MatrixXcd y(num_frame, num_freq);

    MatrixXcd diag_matrix = (1e-4) * (VectorXcd::Ones(num_channels)).asDiagonal();

    for (long long i = 0; i < num_frame; i++)
    {
        for (long long j = 0; j < num_freq; j++)
        {
            MatrixXcd r_inv = (*PhiNN[j] + diag_matrix).inverse(); 
            VectorXcd fx = FX.row(j); 

            VectorXcd w = (r_inv * fx)/(fx.adjoint() * r_inv * fx);
            VectorXcd temp(num_channels);
            for (int c = 0; c < num_channels; c++)
            {
                temp(c) = (*X[c])(i, j);
            }

            y(i, j) = w.adjoint() * temp;
        }
    }
    complex sum_y = y.array().sum();
    return y;
}