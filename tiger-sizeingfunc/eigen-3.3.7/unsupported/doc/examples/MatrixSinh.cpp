#include <spdlog/spdlog.h> 
 #include <unsupported/Eigen/MatrixFunctions>
#include <iostream>

using namespace Eigen;

int main()
{
  MatrixXf A = MatrixXf::Random(3,3);
  spdlog::info("A = \n" << A << "\n\n");

  MatrixXf sinhA = A.sinh();
  spdlog::info("sinh(A) = \n" << sinhA << "\n\n");

  MatrixXf coshA = A.cosh();
  spdlog::info("cosh(A) = \n" << coshA << "\n\n");
  
  // The matrix functions satisfy cosh^2(A) - sinh^2(A) = I, 
  // like the scalar functions.
  spdlog::info("cosh^2(A) - sinh^2(A) = \n" << coshA*coshA - sinhA*sinhA << "\n\n");
}
