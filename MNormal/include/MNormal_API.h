
#include <string>
#include <vector>

#include "blpre.h"


namespace MNormal {
void generateFirstLayer(blpreConfig& blcf,
                        std::string &f,
                        std::vector<std::array<double, 3>> &points,
                        ControlVolume &cv1,bool isALM = false);

void generateFullLayer(blpreConfig& blcf,
              std::string &f,
              std::vector<std::array<double, 3>> &points,
              ControlVolume &cv1,
              ControlVolume &cv2);

} // namespace MNormal
