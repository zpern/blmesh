#include <spdlog/spdlog.h> 
 #include <unsupported/Eigen/MatrixFunctions>
#include <iostream>

using namespace Eigen;

int main()
{
  const double pi = std::acos(-1.0);

  MatrixXd A(3,3);
  A << 0,    -pi/4, 0,
       pi/4, 0,     0,
       0,    0,     0;
  spdlog::info("The matrix A is:\n" << A << "\n\n");
  spdlog::info("The matrix exponential of A is:\n" << A.exp() << "\n\n");
}
