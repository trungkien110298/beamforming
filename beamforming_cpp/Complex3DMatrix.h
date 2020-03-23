#ifndef COMPLEX3DMATRIX_H
#define COMPLEX3DMATRIX_H

#include <Eigen/Dense>
#include <complex>

using namespace std;
using namespace Eigen;

typedef MatrixXcd **Complex3DMatrix;

Complex3DMatrix init_Complex3DMatrix(long long sizes[])
{
    Complex3DMatrix matrix = new MatrixXcd *[sizes[0]];
    for (long long d0 = 0; d0 < sizes[0]; d0++)
    {
        matrix[d0] = new (MatrixXcd)(MatrixXcd::Zero(sizes[1], sizes[2]));
    }
    return matrix;
}

complex<double> sum_Complex3DMatrix(Complex3DMatrix matrix, long long sizes[]){
    complex<double> sum = 0.;
    for (long long d0 = 0; d0 < sizes[0]; d0++){
        sum += (*matrix[d0]).array().sum(); 
    }
    return sum;
}

#endif