#include <spdlog/spdlog.h> 
 #include <Eigen/Dense>
#include <iostream>

using namespace Eigen;

template <typename Derived1, typename Derived2>
void copyUpperTriangularPart(MatrixBase<Derived1>& dst, const MatrixBase<Derived2>& src)
{
  /* Note the 'template' keywords in the following line! */
  dst.template triangularView<Upper>() = src.template triangularView<Upper>();
}

int main()
{
  MatrixXi m1 = MatrixXi::Ones(5,5);
  MatrixXi m2 = MatrixXi::Random(4,4);
  spdlog::info("m2 before copy:");
  spdlog::info(m2 << std::endl);
  copyUpperTriangularPart(m2, m1.topLeftCorner(4,4));
  spdlog::info("m2 after copy:");
  spdlog::info(m2 << std::endl);
}
