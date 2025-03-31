#pragma once
#ifndef _BLPRE_H
#define _BLPRE_H
#include<sstream>
#include <string>
#include <vector>
#include <array>
#include <map>
using namespace std;
struct blpreConfig {
	int n;
	double len;
	double Ro;
	std::vector<double> len_vec;
	std::vector<int> layer_vec;
	vector<int> symm;
	vector<int> box;
	vector<int> match;
	vector<int> per;
	vector<int> adjacent;
	bool use_multiple_normals;
};
struct ControlVolume{
	std::vector<std::array<double, 3>> v;
	std::vector<std::vector<int>> f;
	int lower_point_num;
	int add_point_num;
};
namespace PRE {
	static ControlVolume empty;
	std::tuple<std::string, double*, int*, int*, std::vector<double>> blpre(std::string& f, blpreConfig cf, std::vector < std::array< double, 3> > points = std::vector < std::array< double, 3> >(), ControlVolume& cv = empty);
}
#endif
