#include <spdlog/spdlog.h> 
 MatrixXf a(2,2);
spdlog::info("a is of size " << a.rows() << "x" << a.cols());
MatrixXf b(3,3);
a = b;
spdlog::info("a is now of size " << a.rows() << "x" << a.cols());
