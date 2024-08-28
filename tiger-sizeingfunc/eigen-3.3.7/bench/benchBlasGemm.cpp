// g++ -O3 -DNDEBUG -I.. -L /usr/lib64/atlas/ benchBlasGemm.cpp -o benchBlasGemm -lrt -lcblas
// possible options:
//    -DEIGEN_DONT_VECTORIZE
//    -msse2

// #define EIGEN_DEFAULT_TO_ROW_MAJOR
#define _FLOAT

#include <iostream>

#include <Eigen/Core>
#include <spdlog/spdlog.h> 
 #include "BenchTimer.h"

// include the BLAS headers
extern "C" {
#include <cblas.h>
}
#include <string>

#ifdef _FLOAT
typedef float Scalar;
#define CBLAS_GEMM cblas_sgemm
#else
typedef double Scalar;
#define CBLAS_GEMM cblas_dgemm
#endif


typedef Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> MyMatrix;
void bench_eigengemm(MyMatrix& mc, const MyMatrix& ma, const MyMatrix& mb, int nbloops);
void check_product(int M, int N, int K);
void check_product(void);

int main(int argc, char *argv[])
{
  // disable SSE exceptions
  #ifdef __GNUC__
  {
    int aux;
    asm(
    "stmxcsr   %[aux]           \n\t"
    "orl       $32832, %[aux]   \n\t"
    "ldmxcsr   %[aux]           \n\t"
    : : [aux] "m" (aux));
  }
  #endif

  int nbtries=1, nbloops=1, M, N, K;

  if (argc==2)
  {
    if (std::string(argv[1])=="check")
      check_product();
    else
      M = N = K = atoi(argv[1]);
  }
  else if ((argc==3) && (std::string(argv[1])=="auto"))
  {
    M = N = K = atoi(argv[2]);
    nbloops = 1000000000/(M*M*M);
    if (nbloops<1)
      nbloops = 1;
    nbtries = 6;
  }
  else if (argc==4)
  {
    M = N = K = atoi(argv[1]);
    nbloops = atoi(argv[2]);
    nbtries = atoi(argv[3]);
  }
  else if (argc==6)
  {
    M = atoi(argv[1]);
    N = atoi(argv[2]);
    K = atoi(argv[3]);
    nbloops = atoi(argv[4]);
    nbtries = atoi(argv[5]);
  }
  else
  {
    spdlog::info("Usage: " << argv[0] << " size  \n");
    spdlog::info("Usage: " << argv[0] << " auto size\n");
    spdlog::info("Usage: " << argv[0] << " size nbloops nbtries\n");
    spdlog::info("Usage: " << argv[0] << " M N K nbloops nbtries\n");
    spdlog::info("Usage: " << argv[0] << " check\n");
    spdlog::info("Options:\n");
    spdlog::info("    size       unique size of the 2 matrices (integer)\n");
    spdlog::info("    auto       automatically set the number of repetitions and tries\n");
    spdlog::info("    nbloops    number of times the GEMM routines is executed\n");
    spdlog::info("    nbtries    number of times the loop is benched (return the best try)\n");
    spdlog::info("    M N K      sizes of the matrices: MxN  =  MxK * KxN (integers)\n");
    spdlog::info("    check      check eigen product using cblas as a reference\n");
    exit(1);
  }

  double nbmad = double(M) * double(N) * double(K) * double(nbloops);

  if (!(std::string(argv[1])=="auto"))
    spdlog::info(M << " x " << N << " x " << K );

  Scalar alpha, beta;
  MyMatrix ma(M,K), mb(K,N), mc(M,N);
  ma = MyMatrix::Random(M,K);
  mb = MyMatrix::Random(K,N);
  mc = MyMatrix::Random(M,N);

  Eigen::BenchTimer timer;

  // we simply compute c += a*b, so:
  alpha = 1;
  beta = 1;

  // bench cblas
  // ROWS_A, COLS_B, COLS_A, 1.0,  A, COLS_A, B, COLS_B, 0.0, C, COLS_B);
  if (!(std::string(argv[1])=="auto"))
  {
    timer.reset();
    for (uint k=0 ; k<nbtries ; ++k)
    {
        timer.start();
        for (uint j=0 ; j<nbloops ; ++j)
              #ifdef EIGEN_DEFAULT_TO_ROW_MAJOR
              CBLAS_GEMM(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, alpha, ma.data(), K, mb.data(), N, beta, mc.data(), N);
              #else
              CBLAS_GEMM(CblasColMajor, CblasNoTrans, CblasNoTrans, M, N, K, alpha, ma.data(), M, mb.data(), K, beta, mc.data(), M);
              #endif
        timer.stop();
    }
    if (!(std::string(argv[1])=="auto"))
      spdlog::info("cblas: " << timer.value() << " (" << 1e-3*floor(1e-6*nbmad/timer.value()) << " GFlops/s)\n");
    else
        spdlog::info(M << " : " << timer.value() << " ; " << 1e-3*floor(1e-6*nbmad/timer.value()) );
  }

  // clear
  ma = MyMatrix::Random(M,K);
  mb = MyMatrix::Random(K,N);
  mc = MyMatrix::Random(M,N);

  // eigen
//   if (!(std::string(argv[1])=="auto"))
  {
      timer.reset();
      for (uint k=0 ; k<nbtries ; ++k)
      {
          timer.start();
          bench_eigengemm(mc, ma, mb, nbloops);
          timer.stop();
      }
      if (!(std::string(argv[1])=="auto"))
        spdlog::info("eigen : " << timer.value() << " (" << 1e-3*floor(1e-6*nbmad/timer.value()) << " GFlops/s)\n");
      else
        spdlog::info(M << " : " << timer.value() << " ; " << 1e-3*floor(1e-6*nbmad/timer.value()) );
  }

  spdlog::info("l1: " << Eigen::l1CacheSize());
  spdlog::info("l2: " << Eigen::l2CacheSize());
  

  return 0;
}

using namespace Eigen;

void bench_eigengemm(MyMatrix& mc, const MyMatrix& ma, const MyMatrix& mb, int nbloops)
{
  for (uint j=0 ; j<nbloops ; ++j)
      mc.noalias() += ma * mb;
}

#define MYVERIFY(A,M) if (!(A)) { \
    spdlog::info("FAIL: " << M ); \
  }
void check_product(int M, int N, int K)
{
  MyMatrix ma(M,K), mb(K,N), mc(M,N), maT(K,M), mbT(N,K), meigen(M,N), mref(M,N);
  ma = MyMatrix::Random(M,K);
  mb = MyMatrix::Random(K,N);
  maT = ma.transpose();
  mbT = mb.transpose();
  mc = MyMatrix::Random(M,N);

  MyMatrix::Scalar eps = 1e-4;

  meigen = mref = mc;
  CBLAS_GEMM(CblasColMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1, ma.data(), M, mb.data(), K, 1, mref.data(), M);
  meigen += ma * mb;
  MYVERIFY(meigen.isApprox(mref, eps),". * .");

  meigen = mref = mc;
  CBLAS_GEMM(CblasColMajor, CblasTrans, CblasNoTrans, M, N, K, 1, maT.data(), K, mb.data(), K, 1, mref.data(), M);
  meigen += maT.transpose() * mb;
  MYVERIFY(meigen.isApprox(mref, eps),"T * .");

  meigen = mref = mc;
  CBLAS_GEMM(CblasColMajor, CblasTrans, CblasTrans, M, N, K, 1, maT.data(), K, mbT.data(), N, 1, mref.data(), M);
  meigen += (maT.transpose()) * (mbT.transpose());
  MYVERIFY(meigen.isApprox(mref, eps),"T * T");

  meigen = mref = mc;
  CBLAS_GEMM(CblasColMajor, CblasNoTrans, CblasTrans, M, N, K, 1, ma.data(), M, mbT.data(), N, 1, mref.data(), M);
  meigen += ma * mbT.transpose();
  MYVERIFY(meigen.isApprox(mref, eps),". * T");
}

void check_product(void)
{
  int M, N, K;
  for (uint i=0; i<1000; ++i)
  {
    M = internal::random<int>(1,64);
    N = internal::random<int>(1,768);
    K = internal::random<int>(1,768);
    M = (0 + M) * 1;
    spdlog::info(M << " x " << N << " x " << K );
    check_product(M, N, K);
  }
}

