#include <spdlog/spdlog.h> 
 spdlog::info("A fixed-size array:\n");
Array33f a1 = Array33f::Zero();
spdlog::info(a1 << "\n\n");


spdlog::info("A one-dimensional dynamic-size array:\n");
ArrayXf a2 = ArrayXf::Zero(3);
spdlog::info(a2 << "\n\n");


spdlog::info("A two-dimensional dynamic-size array:\n");
ArrayXXf a3 = ArrayXXf::Zero(3, 4);
spdlog::info(a3 );
