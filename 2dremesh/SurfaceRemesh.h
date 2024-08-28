#pragma once

#include <functional>
#include <array>
#include <set>
#include "Types.h"
#include "PolyMesh.h"

namespace BLMESH
{
namespace zju
{

class SurfaceRemesh
{

public:

  SurfaceRemesh(Mesh &xyz_mesh);

  void iso_remesh(const int &iternum = 10);

  void set_sizing_function(std::function<double(double, double, double)  >& function) { sizingfunction_ = function; }
  void set_pro_function(std::function<std::array<double, 3>(double, double, double) >& function) { projectfunction_ = function; }
  void set_fix_id(std::set<int> s) { fix_id_ = s; }

private:
  

  void split_long_edges();

  void collapse_short_edges();

  void equalize_valences();

  void tangential_relaxation();

  double calculateTargetEdgeLength();
  
  double calculate_target_lenth_shortest();

private:


	std::function<double(double, double,double) > sizingfunction_;
	std::function < std::array<double, 3>(double, double, double)> projectfunction_;
	std::set<int> fix_id_;

  PolyMesh* mesh;
  Mesh &xyz_mesh_;
  

};

} // namespace zju

} // namespace VGModel

