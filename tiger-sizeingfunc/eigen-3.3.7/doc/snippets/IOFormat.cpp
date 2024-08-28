#include <spdlog/spdlog.h> 
 std::string sep = "\n----------------------------------------\n";
Matrix3d m1;
m1 << 1.111111, 2, 3.33333, 4, 5, 6, 7, 8.888888, 9;

IOFormat CommaInitFmt(StreamPrecision, DontAlignCols, ", ", ", ", "", "", " << ", ";");
IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
IOFormat OctaveFmt(StreamPrecision, 0, ", ", ";\n", "", "", "[", "]");
IOFormat HeavyFmt(FullPrecision, 0, ", ", ";\n", "[", "]", "[", "]");

spdlog::info(m1 << sep);
spdlog::info(m1.format(CommaInitFmt) << sep);
spdlog::info(m1.format(CleanFmt) << sep);
spdlog::info(m1.format(OctaveFmt) << sep);
spdlog::info(m1.format(HeavyFmt) << sep);
