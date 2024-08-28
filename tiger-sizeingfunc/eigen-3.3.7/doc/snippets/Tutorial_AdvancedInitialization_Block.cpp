#include <spdlog/spdlog.h> 
 MatrixXf matA(2, 2);
matA << 1, 2, 3, 4;
MatrixXf matB(4, 4);
matB << matA, matA/10, matA/10, matA;
spdlog::info(matB);
