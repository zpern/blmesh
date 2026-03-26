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
void generateFirstLayer(blpreConfig blcf,
                        std::string &f,
                        std::vector<std::array<double, 3>> &points,
                        ControlVolume &cv1)
{

    // read step_length
    double len = blcf.len;
    std::map<std::array<double, 3>, double> point_to_length;

    if (blcf.len_vec.empty()) {
        for (int i = 0; i < points.size(); i++) {
            point_to_length[points[i]] = len;
        }
    }
    else{
        for (int i = 0; i < points.size(); i++) {
            point_to_length[points[i]] = blcf.len_vec[i];
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


    chamfer.pre_WriteVol(cv1.v, cv1.f, cv1.lower_point_num, cv1.add_point_num);

    if (chamfer.multiplySuccess) {
        chamfer.WriteMesh(f_multiply, points_multiply, blcf.len);
        spdlog::info("PreJob Finished.");
    } else {
        multiplySuccess = chamfer.multiplySuccess;
    }
    combine_by_faceID(points, points_multiply, points_nonwall, f, f_multiply, f_nonwall);
    spdlog::info("Job Finished.");
}

//void generateFullLayer(blpreConfig blcf,
//              std::string &f,
//              std::vector<std::array<double, 3>> &points,
//              ControlVolume &cv1,
//              ControlVolume &cv2)
//{
//    int n = blcf.n;
//    double len = blcf.len;
//
//    vector<int> box = blcf.box;
//    vector<int> wall = blcf.wall;
//    vector<int> symm = blcf.symm;
//    vector<int> match = blcf.match;
//    vector<int> per_face = blcf.per;
//    vector<int> adjacent_face = blcf.adjacent;
//    int faceCount = symm.size() + wall.size() + box.size() + match.size() + adjacent_face.size();
//
//    spdlog::info("begin splite.");
//
//    // for fallback
//    bool multiplySuccess = true;
//
//    // --- Step 2: map point ¡ú length ---
//    std::map<std::array<double, 3>, double> point_to_length;
//    if (blcf.len_vec.empty()) {
//        for (int i = 0; i < points.size(); i++) {
//            point_to_length[points[i]] = len;
//        }
//    }
//    else{
//        for (int i = 0; i < points.size(); i++) {
//            point_to_length[points[i]] = blcf.len_vec[i];
//        }
//    }
//
//    // --- Step 3: prism pre-layer processing ---
//    if (blcf.preMultiple) {
//        MNormalMesh chamfer;
//
//        chamfer.number_of_layer = 1;
//        chamfer.step_of_length = blcf.len;
//        chamfer.point_to_length = point_to_length;
//
//        chamfer.ReadPlsBuf(f, points);
//
//        chamfer.CalculateMultiNormal();
//        chamfer.BuildTopo(faceCount);
//
//        spdlog::info("Handling output mesh!");
//
//        chamfer.pre_WriteVol(cv2.v, cv2.f, cv2.lower_point_num, cv2.add_point_num);
//
//        if (chamfer.multiplySuccess) {
//            chamfer.WriteMesh(f, points, blcf.len);
//            spdlog::info("PreJob Finished.");
//        } else {
//            multiplySuccess = chamfer.multiplySuccess;
//        }
//    }
//
//    // --- Step 4: actual multilayer extrusion ---
//    if (blcf.n > 0 && multiplySuccess) {
//        ChamferBehavior behavior;
//        MNormalMesh chamfer;
//
//        chamfer.number_of_layer = blcf.n;
//        chamfer.step_of_length = blcf.len;
//        chamfer.point_to_length = point_to_length;
//        chamfer.exist_prism = blcf.preMultiple;
//
//        chamfer.SetBehavior(behavior);
//        chamfer.ReadPlsBuf(f, points);
//
//        spdlog::info("Done!");
//
//        chamfer.CalculateMultiNormal();
//        chamfer.BuildTopo(faceCount);
//
//        spdlog::info("Handling output mesh!");
//        chamfer.WriteVol(cv1.v, cv1.f, cv1.lower_point_num, cv1.add_point_num);
//
//        if (chamfer.multiplySuccess) {
//            chamfer.WriteMesh(f, points, blcf.len);
//        }
//    } else {
//        spdlog::info("Skipping Chamfer processing because number_of_layer is 0");
//    }
//
//    spdlog::info("Job Finished.");
//}
} // namespace MNormal
