#include <spdlog/spdlog.h> 
 #include <Eigen/Dense>
#include <iostream>

using namespace Eigen;

void copyUpperTriangularPart(MatrixXf& dst, const MatrixXf& src)
{
  dst.triangularView<Upper>() = src.triangularView<Upper>();
}

int main()
{
  MatrixXf m1 = MatrixXf::Ones(4,4);
  MatrixXf m2 = MatrixXf::Random(4,4);
  spdlog::info("m2 before copy:");
  spdlog::info(m2 << std::endl);
  copyUpperTriangularPart(m2, m1);
  spdlog::info("m2 after copy:");
  spdlog::info(m2 << std::endl);
}
