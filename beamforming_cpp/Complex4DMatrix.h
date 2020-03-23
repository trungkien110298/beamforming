#ifndef COMPLEX4DMATRIX_H
#define COMPLEX4DMATRIX_H

#include <Eigen/Dense>
#include <complex>

using namespace std;
using namespace Eigen;

typedef MatrixXcd ***Complex4DMatrix;

Complex4DMatrix init_Complex4DMatrix(long long sizes[])
{
    Complex4DMatrix matrix = new MatrixXcd **[sizes[0]];
    for (long long d0 = 0; d0 < sizes[0]; d0++)
    {
        matrix[d0] = new MatrixXcd *[sizes[1]];
        for (long long d1 = 0; d1 < sizes[1]; d1++)
            matrix[d0][d1] = new (MatrixXcd)(MatrixXcd(sizes[2], sizes[3]));
    }
    return matrix;
}


complex<double> sum_Complex4DMatrix(Complex4DMatrix matrix, long long sizes[]){
    complex<double> sum = 0.;
    for (long long d0 = 0; d0 < sizes[0]; d0++){
        for (long long d1 = 0; d1 < sizes[1]; d1++)
            sum += (*matrix[d0][d1]).array().sum();
    }
    return sum;
}

#endif