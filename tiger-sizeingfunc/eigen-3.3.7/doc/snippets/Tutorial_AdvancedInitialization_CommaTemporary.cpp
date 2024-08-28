#include <spdlog/spdlog.h> 
 MatrixXf mat = MatrixXf::Random(2, 3);
spdlog::info(mat << std::endl);
mat = (MatrixXf(2,2) << 0, 1, 1, 0).finished() * mat;
spdlog::info(mat);
