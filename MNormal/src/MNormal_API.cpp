#include <stdio.h>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <spdlog/spdlog.h> 
#include "MNormal_API.h"
#include "MNormalMesh.h"
using namespace std;

namespace MNormal {
void generateFirstLayer(blpreConfig& blcf,
                        std::string &f,
                        std::vector<std::array<double, 3>> &points,
                        ControlVolume &cv1, bool isALM)
{
    // read step_length
    double len = blcf.len;
    std::map<std::array<double, 3>, double> point_to_length;

    if (blcf.point_length_vec.empty()) {
        for (int i = 0; i < points.size(); i++) {
            point_to_length[points[i]] = len;
        }
    }
    else{
        for (int i = 0; i < points.size(); i++) {
            point_to_length[points[i]] = blcf.point_length_vec[i];
        }
    }
    vector<int> box = blcf.box;
    vector<int> wall = blcf.wall;
    vector<int> symm = blcf.symm;
    vector<int> match = blcf.match;
    vector<int> per_face = blcf.per;
    vector<int> adjacent_face = blcf.adjacent;
    int faceCount = symm.size() + wall.size() + box.size() + match.size() + adjacent_face.size();

    std::vector<std::array<double, 3>> points_multiply;
    std::vector<std::array<double, 3>> points_nonwall;
    std::string f_nonwall;
    std::string f_multiply;
    splite_by_faceID(points, points_multiply, points_nonwall, f, f_multiply, f_nonwall, wall);

    spdlog::info("begin splite.");

    // for fallback
    bool multiplySuccess = true;

    MNormalMesh chamfer;
    chamfer.number_of_layer = 1;
    chamfer.point_to_length = point_to_length;

    chamfer.ReadPlsBuf(f_multiply, points_multiply);
    chamfer.FixedLength();
    chamfer.CalculateMultiNormal();
    chamfer.BuildTopo(faceCount);
    chamfer.Generate_preVol(cv1.v, cv1.f,isALM);

    if (chamfer.multiplySuccess) {
        chamfer.WriteSurMesh(f_multiply, points_multiply);
        spdlog::info("PreJob Finished.");
    } else {
        multiplySuccess = chamfer.multiplySuccess;
        spdlog::info("mulltiply failed.");
    }

    // 建一个 coordinate → length 的查找表
    std::map<std::array<double, 3>, double> coord2length;
    if (blcf.len_vec.size()) {
        for (int i = 0; i < points_multiply.size(); ++i) {
            std::array<double, 3> key = {
                points_multiply[i][0],
                points_multiply[i][1],
                points_multiply[i][2]};
            coord2length[key] = chamfer.length[i];
        }
    }


    combine_by_faceID(points, points_multiply, points_nonwall, f, f_multiply, f_nonwall);

    
    std::vector<double> points_length(points.size(), 0.0);
    if (blcf.len_vec.size()) {
        // 给最终 points 赋值
        for (int i = 0; i < points.size(); ++i) {
            auto it1 = coord2length.find(points[i]);
            auto it2 = point_to_length.find(points[i]);
            if (it2 != point_to_length.end()) {
                points_length[i] = it2->second;
            } else if (it1 != coord2length.end()) {
                points_length[i] = it1->second; // 原点 → 继承 length
            } else {
                points_length[i] = 0.0;        // 新点 → 0
            }
        }

        blcf.point_length_vec = std::move(points_length);
    }
    spdlog::info("Job Finished.");
}

void generateFullLayer(blpreConfig blcf,
              std::string &f,
              std::vector<std::array<double, 3>> &points,
              ControlVolume &cv1,
              ControlVolume &cv2)
{
    generateFirstLayer(blcf, f, points, cv1);

    // read step_length
    double len = blcf.len;
    std::map<std::array<double, 3>, double> point_to_length;

    if (blcf.point_length_vec.empty()) {
        for (int i = 0; i < points.size(); i++) {
            point_to_length[points[i]] = len;
        }
    }
    else{
        for (int i = 0; i < points.size(); i++) {
            point_to_length[points[i]] = blcf.point_length_vec[i];
        }
    }
    vector<int> box = blcf.box;
    vector<int> wall = blcf.wall;
    vector<int> symm = blcf.symm;
    vector<int> match = blcf.match;
    vector<int> per_face = blcf.per;
    vector<int> adjacent_face = blcf.adjacent;
    int faceCount = symm.size() + wall.size() + box.size() + match.size() + adjacent_face.size();
    wall.push_back(faceCount+1);

    std::vector<std::array<double, 3>> points_multiply;
    std::vector<std::array<double, 3>> points_nonwall;
    std::string f_nonwall;
    std::string f_multiply;
    splite_by_faceID(points, points_multiply, points_nonwall, f, f_multiply, f_nonwall, wall);

    spdlog::info("begin Full_layer.");

    // for fallback
    bool multiplySuccess = true;

    MNormalMesh chamfer;
    chamfer.number_of_layer = blcf.n;
    chamfer.point_to_length = point_to_length;
    chamfer.firstLayer = true;
    chamfer.ReadPlsBuf(f_multiply, points_multiply);
    chamfer.FixedLength();
    chamfer.CalculateMultiNormal();
    chamfer.BuildTopo(faceCount);
    chamfer.Generate_Vol(cv2.v, cv2.f);

    if (chamfer.multiplySuccess) {
        chamfer.WriteSurMesh(f_multiply, points_multiply);
        spdlog::info("PreJob Finished.");
    } else {
        multiplySuccess = chamfer.multiplySuccess;
        spdlog::info("mulltiply failed.");
    }

    combine_by_faceID(points, points_multiply, points_nonwall, f, f_multiply, f_nonwall);
    spdlog::info("Job Finished.");
}
} // namespace MNormal
