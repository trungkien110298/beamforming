#include <blaze/Blaze.h>
#include <complex>
using namespace std;
int main()
{

    using blaze::DynamicMatrix;
    using blaze::DynamicVector;
    using blaze::rowMajor;
    using blaze::columnVector;
    DynamicMatrix<double, rowMajor> A(5, 5); // The general matrix A
    // ... Initialization

    DynamicVector<complex<double>, columnVector> w(5); // The vector for the complex eigenvalues
    DynamicMatrix<complex<double>, rowMajor> V(5, 5);  // The matrix for the left eigenvectors

    w = eigen(A);   // Computing only the eigenvalues of A (one argument)
    eigen(A, w);    // Computing only the eigenvalues of A (two arguments)
    eigen(A, w, V); // Computing both the eigenvalues and eigenvectors of A (three arguments)
}
