#include <spdlog/spdlog.h> 
 #include <iostream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;
int main()
{
  MatrixXf mat(2,4);
  mat << 1, 2, 6, 9,
         3, 1, 7, 2;
  
  MatrixXf::Index   maxIndex;
  float maxNorm = mat.colwise().sum().maxCoeff(&maxIndex);
  
  spdlog::info("Maximum sum at position " << maxIndex);

  spdlog::info("The corresponding vector is: ");
  spdlog::info(mat.col( maxIndex ));
  spdlog::info("And its sum is is: " << maxNorm);
}
