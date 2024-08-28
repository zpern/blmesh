#include <spdlog/spdlog.h> 
 #include <Eigen/Core>
#include <iostream>

using namespace Eigen;

// [functor]
template<class ArgType, class RowIndexType, class ColIndexType>
class indexing_functor {
  const ArgType &m_arg;
  const RowIndexType &m_rowIndices;
  const ColIndexType &m_colIndices;
public:
  typedef Matrix<typename ArgType::Scalar,
                 RowIndexType::SizeAtCompileTime,
                 ColIndexType::SizeAtCompileTime,
                 ArgType::Flags&RowMajorBit?RowMajor:ColMajor,
                 RowIndexType::MaxSizeAtCompileTime,
                 ColIndexType::MaxSizeAtCompileTime> MatrixType;

  indexing_functor(const ArgType& arg, const RowIndexType& row_indices, const ColIndexType& col_indices)
    : m_arg(arg), m_rowIndices(row_indices), m_colIndices(col_indices)
  {}

  const typename ArgType::Scalar& operator() (Index row, Index col) const {
    return m_arg(m_rowIndices[row], m_colIndices[col]);
  }
};
// [functor]

// [function]
template <class ArgType, class RowIndexType, class ColIndexType>
CwiseNullaryOp<indexing_functor<ArgType,RowIndexType,ColIndexType>, typename indexing_functor<ArgType,RowIndexType,ColIndexType>::MatrixType>
indexing(const Eigen::MatrixBase<ArgType>& arg, const RowIndexType& row_indices, const ColIndexType& col_indices)
{
  typedef indexing_functor<ArgType,RowIndexType,ColIndexType> Func;
  typedef typename Func::MatrixType MatrixType;
  return MatrixType::NullaryExpr(row_indices.size(), col_indices.size(), Func(arg.derived(), row_indices, col_indices));
}
// [function]


int main()
{
  spdlog::info("[main1]\n");
  Eigen::MatrixXi A = Eigen::MatrixXi::Random(4,4);
  Array3i ri(1,2,1);
  ArrayXi ci(6); ci << 3,2,1,0,0,2;
  Eigen::MatrixXi B = indexing(A, ri, ci);
  spdlog::info("A =");
  spdlog::info(A << std::endl);
  spdlog::info("A([" << ri.transpose() << "], [" << ci.transpose() << "]) =");
  spdlog::info(B);
  spdlog::info("[main1]\n");

  spdlog::info("[main2]\n");
  B =  indexing(A, ri+1, ci);
  spdlog::info("A(ri+1,ci) =");
  spdlog::info(B << std::endl);
#if __cplusplus >= 201103L
  B =  indexing(A, ArrayXi::LinSpaced(13,0,12).unaryExpr([](int x){return x%4;}), ArrayXi::LinSpaced(4,0,3));
  spdlog::info("A(ArrayXi::LinSpaced(13,0,12).unaryExpr([](int x){return x%4;}), ArrayXi::LinSpaced(4,0,3)) =");
  spdlog::info(B << std::endl);
#endif
  spdlog::info("[main2]\n");
}

