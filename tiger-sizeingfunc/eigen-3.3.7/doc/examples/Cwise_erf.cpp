#include <spdlog/spdlog.h> 
 #include <Eigen/Core>
#include <unsupported/Eigen/SpecialFunctions>
#include <iostream>
using namespace Eigen;
int main()
{
  Array4d v(-0.5,2,0,-7);
  spdlog::info(v.erf());
}
