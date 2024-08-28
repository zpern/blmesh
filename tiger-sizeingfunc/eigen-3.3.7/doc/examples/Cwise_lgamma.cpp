#include <spdlog/spdlog.h> 
 #include <Eigen/Core>
#include <unsupported/Eigen/SpecialFunctions>
#include <iostream>
using namespace Eigen;
int main()
{
  Array4d v(0.5,10,0,-1);
  spdlog::info(v.lgamma());
}
