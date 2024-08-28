#include <spdlog/spdlog.h> 
 #include <unsupported/Eigen/MatrixFunctions>
#include <iostream>

using namespace Eigen;

int main()
{
  const double pi = std::acos(-1.0);

  MatrixXd A(2,2);
  A << cos(pi/3), -sin(pi/3), 
       sin(pi/3),  cos(pi/3);
  spdlog::info("The matrix A is:\n" << A << "\n\n");
  spdlog::info("The matrix square root of A is:\n" << A.sqrt() << "\n\n");
  spdlog::info("The square of the last matrix is:\n" << A.sqrt() * A.sqrt() );
}
