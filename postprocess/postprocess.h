#ifndef POSTPROCESS_H_
#define POSTPROCESS_H_
#include <iostream>
#include <exception>
#include <vector>
#include <array>
#include <cmath>
#include <set>
#include "../BLMesh/include/BLVector.h"
#include <limits>

//return true if a extra optimization is necessary
bool MeshOptimize(std::vector<BLVector> &coordinate, std::vector < std::vector < int>> &connector, std::set<int>& boundary_point);
#endif // DATA_H
