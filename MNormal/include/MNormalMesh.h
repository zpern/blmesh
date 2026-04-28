#pragma once
#ifndef _MNORMAL_MESH_H_
#define _MNORMAL_MESH_H_
#include "./geometry_check.hpp" // Assumes it contains the declaration for tri_tri_intersect.
#include "./hexa_octree.hpp"    // Uses the previously defined binary tree.
#include "./hexa_tag.h"         // Provides HexaTag and HexaTagHash.
#include "./intersection_check.hpp"
#include "complexnode.h"
#include <array>
#include <string>
#include <vector>

class MNormalMesh {
public:
    MNormalMesh() {}

public:
    unsigned int number_of_origin_node;
    unsigned int number_of_node;
    std::vector<BLVector> coordinate;
    std::vector<BLVector> point_normals;
    std::vector<int> real_node_id_; /// mark as virutal tags, id=the coorespoding real node id;

    unsigned int number_of_origin_triangles;
    unsigned int number_of_triangles;
    std::vector<std::array<int, 3>> connector;
    std::vector<int> attribute;
    
    vector<ComplexNode> node_array;
    std::vector<std::vector<int>> point_neighbors_;

    int number_of_layer;
    std::vector<double> length;
    std::map<std::array<double, 3>, double> point_to_length; /*建立点——步长对应关系 **/

    bool firstLayer = false;
    bool multiplySuccess = true;
    std::vector<int> avoid_spliteNode;
    /*生成边界层的准备数据 **/
    struct GenContext {
        std::vector<std::array<int, 3>> lower_ids; // 每个 connector 三角片在底层唯一点编号中的映射
        std::vector<std::array<double, 3>> bottom_points; // 底层唯一点坐标
        std::set<int> duplicate_lower_ids;                // connector 中多法向点对应的多法向点号
        int lower_num = 0;
    }; 
enum class CellType
{
    TET,
    PRISM
};

    struct TempCell {
        CellType type = CellType::TET;
        std::vector<int> verts;         // 单元顶点编号
        std::vector<int> source_points; // 这个单元关联的原始表面点编号
        int tri_id = -1;                // 来自哪个 connector 三角片
        int layer = -1;                 // 位于第几层之间
    };

    struct CandidateVolume {
        std::vector<std::array<double, 3>> vertices;
        std::vector<TempCell> cells;
    };
    struct CheckResult {
        bool ok = true;
        std::set<int> bad_points;   // 需要回缩步长的原始点
        std::vector<int> bad_cells; // 出问题的单元编号
        std::vector<int> bad_faces; // 出问题的面编号
    };

double SignedTetVolume(const std::array<double, 3>& a,
                       const std::array<double, 3>& b,
                       const std::array<double, 3>& c,
                       const std::array<double, 3>& d) const;

bool MNormalMesh::IsBadPrismByScaledJacobian(const BLVector &p0,
                                             const BLVector &p1,
                                             const BLVector &p2,
                                             const BLVector &p3,
                                             const BLVector &p4,
                                             const BLVector &p5) const;
    void ReadPlsBuf(std::string f, std::vector<std::array<double, 3>> &points);
    void FixedLength();
    void CalculateMultiNormal();
    void SmoothNormalsField(int itertion);
    void BuildTopo(int faceCount);
    void MNormalMesh::RebuildPointNeighbors();

    void Generate_preVol(std::vector<std::array<double, 3>> &v,
                         std::vector<std::vector<int>> &f,bool isALM = false);
    void Generate_Vol(std::vector<std::array<double, 3>> &v,
                     std::vector<std::vector<int>> &f);
    void WriteSurMesh(std::string &f, std::vector<std::array<double, 3>> &points);

    std::vector<double> MNormalMesh::InitLengthField(const GenContext &ctx) const;
    bool MNormalMesh::ResolveLengthField(const GenContext &ctx,
                                         std::vector<double> &length,
                                         CandidateVolume &final_candidate,bool IsPrism = false,bool isALM = false) const;
    GenContext BuildPreGenContext() const;
    void MNormalMesh::SmoothLengthField(std::vector<double> &length) const;
    std::vector<std::array<double, 3>> MNormalMesh::BuildLayerPoints(
        const GenContext &ctx,
        const std::vector<double> &length) const;
    int MNormalMesh::Idx(int base, int lower_num, int layer = 0) const;
    CandidateVolume MNormalMesh::BuildCandidateVolume(const GenContext &ctx,
                                                      const std::vector<double> &length) const;
    MNormalMesh::CandidateVolume MNormalMesh::BuildPrismVolume(
        const GenContext &ctx,
        const std::vector<double> &length) const;
    void MNormalMesh::CutPrisms(const GenContext &ctx,
                                             std::vector<std::vector<int>> &f) const;
    CheckResult MNormalMesh::CheckSurfaceIntersection(const GenContext &ctx,
                                                      const std::vector<double> &length) const;
    void MNormalMesh::CommitCandidate(const GenContext &ctx,
                                      const CandidateVolume &candidate,
                                      std::vector<std::array<double, 3>> &v,
                                      std::vector<std::vector<int>> &f) const;
    CheckResult MNormalMesh::CheckCellJacobian(const CandidateVolume &candidate) const;
    void MNormalMesh::ShrinkLengthField(std::vector<double> &length,
                                        const std::set<int> &bad_points,
                                        double ratio) const
    {
        for (int p : bad_points) {
            length[p] *= ratio;
        }
    }
    void MNormalMesh::ZeroLengthField(std::vector<double> &length,
                                      const std::set<int> &bad_points) const
    {
        for (int p : bad_points) {
            length[p] = 0.0;
        }
    }

protected:
    void CaculateFrontNormal();
    void CalculateNodeNormal();
};
void splite_by_faceID(std::vector<std::array<double, 3>> &point,
                      std::vector<std::array<double, 3>> &point_multiply,
                      std::vector<std::array<double, 3>> &point_nonwall,
                      std::string &f,
                      std::string &f_multiply,
                      std::string &f_nonwall,
                      std::vector<int> srufaceID);
void combine_by_faceID(std::vector<std::array<double, 3>> &points,
                       std::vector<std::array<double, 3>> points_multiply,
                       std::vector<std::array<double, 3>> points_nonwall,
                       std::string &f,
                       std::string f_multiply,
                       std::string f_nonwall);
#endif
