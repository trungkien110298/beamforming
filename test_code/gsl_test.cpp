#include <iostream>
#include <stdio.h>
#include <cmath>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_eigen.h>
using namespace std;
const int N = 4;

int main()
{

    gsl_eigen_hermv_workspace *workN = gsl_eigen_hermv_alloc(N);
    gsl_matrix_complex *A = gsl_matrix_complex_alloc(N, N);
    gsl_complex i = gsl_complex_rect(0.0, 1.0);
    gsl_complex ii = gsl_complex_rect(0.0, -1.0);
    gsl_vector *eval = gsl_vector_alloc(N);
    gsl_matrix_complex *evec = gsl_matrix_complex_alloc(N, N);

    double mTab[] = {
        0,
        0,
        1,
        0,
        0,
        0,
        5,
        0,
        1,
        0,
        0,
        0,
        5,
        0,
        0,
        0,
        0,
        0,
        5,
        0,
        0,
        0,
        1,
        0,
        5,
        0,
        0,
        0,
        1,
        0,
        0,
        0,
    };

    gsl_matrix_complex_view tmpM = gsl_matrix_complex_view_array(mTab, N, N);

    gsl_matrix_complex_memcpy(A, &tmpM.matrix);
    gsl_matrix_complex_set(A, 0, 3, ii);
    gsl_matrix_complex_set(A, 1, 2, ii);
    gsl_matrix_complex_set(A, 2, 1, i);
    gsl_matrix_complex_set(A, 3, 0, i);

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            gsl_complex z = gsl_matrix_complex_get(A, i, j);
            cout << GSL_REAL(z) << "+ i" << GSL_IMAG(z) << "    ";
        }
        cout << "\n";
    }

    gsl_eigen_hermv(A, eval, evec, workN);

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            gsl_complex z = gsl_matrix_complex_get(A, i, j);
            cout << GSL_REAL(z) << "+ i" << GSL_IMAG(z) << "    ";
        }
        cout << "\n";
    }

    cout << "\n";
    for (int i = 0; i < N; i++)
    {
        cout << gsl_vector_get(eval, i) << " ";
    }

    return 0;
}