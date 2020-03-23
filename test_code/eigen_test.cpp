#include <iostream>
#include <complex>
#include <Eigen/Eigenvalues> 

int main()
{
  const int n = 2;
  Eigen::MatrixXcd a(n, n);
  typedef std::complex<double> C;
  a <<
    C(1.0, 2.0), C(2.0, 1.0), C(3.0, -1.0), C(4.0, -2.0);
  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ces;
  ces.compute(a);
  std::cout << "The eigenvalues of a are:" << std::endl << ces.eigenvalues() << std::endl;

  std::cout << "The eigenvectors of a are:" << std::endl << ces.eigenvectors() << std::endl;
  return 0;
}