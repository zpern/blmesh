#include <spdlog/spdlog.h> 
 #include <iostream>
#include <Eigen/Dense>

using namespace Eigen;
int main()
{
  Matrix2d mat;
  mat << 1, 2,
         3, 4;
  Vector2d u(-1,1), v(2,0);
  spdlog::info("Here is mat*mat:\n" << mat*mat);
  spdlog::info("Here is mat*u:\n" << mat*u);
  spdlog::info("Here is u^T*mat:\n" << u.transpose()*mat);
  spdlog::info("Here is u^T*v:\n" << u.transpose()*v);
  spdlog::info("Here is u*v^T:\n" << u*v.transpose());
  spdlog::info("Let's multiply mat by itself");
  mat = mat*mat;
  spdlog::info("Now mat is mat:\n" << mat);
}
