#include <spdlog/spdlog.h> 
 #include <iostream>
#include <Eigen/Dense>

using namespace Eigen;

int main()
{
  MatrixXd m(2,5);
  m.resize(4,3);
  std::cout << "The matrix m is of size "
            << m.rows() << "x" << m.cols() << std::endl;
  spdlog::info("It has " << m.size() << " coefficients");
  VectorXd v(2);
  v.resize(5);
  spdlog::info("The vector v is of size " << v.size());
  std::cout << "As a matrix, v is of size "
            << v.rows() << "x" << v.cols() << std::endl;
}
