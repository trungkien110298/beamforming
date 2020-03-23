#ifndef DOUBLE4DMATRIX_H
#define DOUBLE4DMATRIX_H

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

typedef MatrixXd ***Double4DMatrix;

Double4DMatrix init_Double4DMatrix(long long sizes[])
{
    Double4DMatrix matrix = new MatrixXd **[sizes[0]];
    for (long long d0 = 0; d0 < sizes[0]; d0++)
    {
        matrix[d0] = new MatrixXd *[sizes[1]];
        for (long long d1 = 0; d1 < sizes[1]; d1++)
            matrix[d0][d1] = new (MatrixXd)(MatrixXd::Zero(sizes[2], sizes[3]));
    }
    return matrix;
}

#endif