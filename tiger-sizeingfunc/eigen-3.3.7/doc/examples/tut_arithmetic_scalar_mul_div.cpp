#include <spdlog/spdlog.h> 
 #include <iostream>
#include <Eigen/Dense>

using namespace Eigen;

int main()
{
  Matrix2d a;
  a << 1, 2,
       3, 4;
  Vector3d v(1,2,3);
  spdlog::info("a * 2.5 =\n" << a * 2.5);
  spdlog::info("0.1 * v =\n" << 0.1 * v);
  spdlog::info("Doing v *= 2;");
  v *= 2;
  spdlog::info("Now v =\n" << v);
}
