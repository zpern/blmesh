#include <spdlog/spdlog.h> 
 RowVectorXd vec1(3);
vec1 << 1, 2, 3;
spdlog::info("vec1 = " << vec1);

RowVectorXd vec2(4);
vec2 << 1, 4, 9, 16;
spdlog::info("vec2 = " << vec2);

RowVectorXd joined(7);
joined << vec1, vec2;
spdlog::info("joined = " << joined);
