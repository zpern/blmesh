#pragma once
#ifndef _BLPRE_H
#define _BLPRE_H
#include <array>
#include <map>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

struct blpreConfig {
    int n;
    double len;
    double Ro;
    std::vector<double> len_vec;
    std::vector<int> layer_vec;
    bool use_multiple_normals;
    int max_layer_diff;
    double clearance;
    vector<double> max_skewness;
    vector<double> max_orth;

    vector<int> box;
    vector<int> wall;
    vector<int> symm;
    vector<int> match;
    vector<int> per;
    vector<int> adjacent;
};

struct ControlVolume {
    std::vector<std::array<double, 3>> v;
    std::vector<std::vector<int>> f;
    std::vector<std::array<int, 3>> s;
    int lower_point_num;
    int add_point_num;
};


namespace PRE {
std::tuple<
    std::string,
    std::vector<std::array<double, 3>>,
    std::vector<std::array<int, 4>>,
    std::vector<int>,
    std::vector<double>> 
    blpre(blpreConfig cf,std::string &f,std::vector<std::array<double, 3>> points);

std::tuple<std::string, double *, int *, int *, std::vector<double>>temptransform(std::tuple<std::string,
                                                                                  std::vector<std::array<double, 3>>,
                                                                                  std::vector<std::array<int, 4>>,
                                                                                  std::vector<int>,
                                                                                  std::vector<double>> &in);
} // namespace PRE
#endif
