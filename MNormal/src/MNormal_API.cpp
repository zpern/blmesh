#include "MNormal_API.h"
#include "MNormalMesh.h"
#include <fstream>
#include <iostream>
#include <spdlog/spdlog.h>
#include <sstream>
#include <stdio.h>
#include <string>
#include <vector>
using namespace std;

namespace MNormal {
void generateFirstLayer(blpreConfig &blcf, std::string &f, std::vector<std::array<double, 3>> &points, ControlVolume &cv1, bool isALM)
{
    // read step_length
    double len = blcf.len;
    std::map<std::array<double, 3>, double> point_to_length;
    std::vector<double> original_point_length(points.size(), len);

    if (blcf.point_length_vec.empty()) {
        for (int i = 0; i < points.size(); i++) {
            point_to_length[points[i]] = len;
        }
    } else {
        for (int i = 0; i < points.size(); i++) {
            point_to_length[points[i]] = blcf.point_length_vec[i];
            original_point_length[i] = blcf.point_length_vec[i];
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
    std::vector<int> point_multiply_to_original;
    std::vector<int> point_nonwall_to_original;
    splite_by_faceID(points, points_multiply, points_nonwall, f, f_multiply, f_nonwall, wall, &point_multiply_to_original,
                     &point_nonwall_to_original);

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

    chamfer.node_original_id_.assign(chamfer.real_node_id_.size(), -1);
    for (int i = 0; i < chamfer.real_node_id_.size(); ++i) {
        int local_id = chamfer.real_node_id_[i];
        if (0 <= local_id && local_id < point_multiply_to_original.size()) {
            chamfer.node_original_id_[i] = point_multiply_to_original[local_id];
        }
    }

    chamfer.Generate_preVol(cv1.v, cv1.f, isALM);

    std::vector<double> current_multiply_length(chamfer.real_node_id_.size(), len);
    for (int i = 0; i < chamfer.node_original_id_.size(); ++i) {
        int original_id = chamfer.node_original_id_[i];
        if (0 <= original_id && original_id < original_point_length.size()) {
            current_multiply_length[i] = original_point_length[original_id];
        }
        if (i < chamfer.length.size()) {
            current_multiply_length[i] = chamfer.length[i];
        }
    }

    if (chamfer.multiplySuccess) {
        chamfer.WriteSurMesh(f_multiply, points_multiply);
        spdlog::info("PreJob Finished.");
    } else {
        multiplySuccess = chamfer.multiplySuccess;
        spdlog::info("mulltiply failed.");
    }

    std::vector<double> nonwall_length(point_nonwall_to_original.size(), len);
    for (int i = 0; i < point_nonwall_to_original.size(); ++i) {
        int original_id = point_nonwall_to_original[i];
        if (0 <= original_id && original_id < original_point_length.size()) {
            nonwall_length[i] = original_point_length[original_id];
        }
    }

    std::vector<int> point_to_original;
    std::vector<double> points_length;
    combine_by_faceID(points, points_multiply, points_nonwall, f, f_multiply, f_nonwall, &chamfer.node_original_id_,
                      &point_nonwall_to_original, &point_to_original, &current_multiply_length, &nonwall_length, &points_length);

    blcf.point_original_id_vec = std::move(point_to_original);

    if (blcf.len_vec.size()) {
        blcf.point_length_vec = std::move(points_length);
    }
    spdlog::info("Job Finished.");
}

void generateFullLayer(blpreConfig &blcf,
                       std::string &f,
                       std::vector<std::array<double, 3>> &points,
                       ControlVolume &cv1,
                       ControlVolume &cv2)
{
    generateFirstLayer(blcf, f, points, cv1);
    std::vector<int> previous_point_to_original = blcf.point_original_id_vec;

    // read step_length
    double len = blcf.len;
    std::map<std::array<double, 3>, double> point_to_length;

    if (blcf.point_length_vec.empty()) {
        for (int i = 0; i < points.size(); i++) {
            point_to_length[points[i]] = len;
        }
    } else {
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
    wall.push_back(faceCount + 1);

    std::vector<std::array<double, 3>> points_multiply;
    std::vector<std::array<double, 3>> points_nonwall;
    std::string f_nonwall;
    std::string f_multiply;
    std::vector<int> point_multiply_to_original;
    std::vector<int> point_nonwall_to_original;
    splite_by_faceID(points, points_multiply, points_nonwall, f, f_multiply, f_nonwall, wall, &point_multiply_to_original,
                     &point_nonwall_to_original);

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

    chamfer.node_original_id_.assign(chamfer.real_node_id_.size(), -1);
    for (int i = 0; i < chamfer.real_node_id_.size(); ++i) {
        int local_id = chamfer.real_node_id_[i];
        if (0 <= local_id && local_id < point_multiply_to_original.size()) {
            int current_point_id = point_multiply_to_original[local_id];
            if (0 <= current_point_id && current_point_id < previous_point_to_original.size()) {
                chamfer.node_original_id_[i] = previous_point_to_original[current_point_id];
            } else {
                spdlog::warn("missing original id mapping for multiply point {}", current_point_id);
            }
        }
    }

    chamfer.Generate_Vol(cv2.v, cv2.f);

    if (chamfer.multiplySuccess) {
        chamfer.WriteSurMesh(f_multiply, points_multiply);
        spdlog::info("PreJob Finished.");
    } else {
        multiplySuccess = chamfer.multiplySuccess;
        spdlog::info("mulltiply failed.");
    }

    std::vector<int> current_nonwall_to_original(point_nonwall_to_original.size(), -1);
    for (int i = 0; i < point_nonwall_to_original.size(); ++i) {
        int current_point_id = point_nonwall_to_original[i];
        if (0 <= current_point_id && current_point_id < previous_point_to_original.size()) {
            current_nonwall_to_original[i] = previous_point_to_original[current_point_id];
        } else {
            spdlog::warn("missing original id mapping for nonwall point {}", current_point_id);
        }
    }

    std::vector<int> point_to_original;
    combine_by_faceID(points, points_multiply, points_nonwall, f, f_multiply, f_nonwall, &chamfer.node_original_id_,
                      &current_nonwall_to_original, &point_to_original);
    blcf.point_original_id_vec = std::move(point_to_original);
    spdlog::info("Job Finished.");
}
} // namespace MNormal
