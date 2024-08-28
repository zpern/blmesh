#include <spdlog/spdlog.h> 
 ArrayXXf table(10, 4);
table.col(0) = ArrayXf::LinSpaced(10, 0, 90);
table.col(1) = M_PI / 180 * table.col(0);
table.col(2) = table.col(1).sin();
table.col(3) = table.col(1).cos();
spdlog::info("  Degrees   Radians      Sine    Cosine\n");
spdlog::info(table);
