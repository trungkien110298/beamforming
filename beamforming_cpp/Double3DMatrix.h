#ifndef DOUBLE3DMATRIX_H
#define DOUBLE3DMATRIX_H

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

typedef MatrixXd **Double3DMatrix;

Double3DMatrix init_Double3DMatrix(long long sizes[])
{
    Double3DMatrix matrix = new MatrixXd *[sizes[0]];
    for (long long d0 = 0; d0 < sizes[0]; d0++)
    {
        matrix[d0] = new (MatrixXd)(MatrixXd::Zero(sizes[1], sizes[2]));
    }
    return matrix;
}

#endif