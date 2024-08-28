#include <spdlog/spdlog.h> 
 #include <iostream>
#include <Eigen/Dense>

using namespace Eigen;

int main()
{
  Matrix2d a;
  a << 1, 2,
       3, 4;
  MatrixXd b(2,2);
  b << 2, 3,
       1, 4;
  spdlog::info("a + b =\n" << a + b);
  spdlog::info("a - b =\n" << a - b);
  spdlog::info("Doing a += b;");
  a += b;
  spdlog::info("Now a =\n" << a);
  Vector3d v(1,2,3);
  Vector3d w(1,0,0);
  spdlog::info("-v + w - v =\n" << -v + w - v);
}
