#ifndef INTERSECTION_CHECK_HPP
#define INTERSECTION_CHECK_HPP
#include "./hexa_octree.hpp"   // Uses the previously defined binary tree.
#include "./hexa_tag.h"        // Provides HexaTag and HexaTagHash.
#include "./geometry_check.hpp"  // Assumes it contains the declaration for tri_tri_intersect.
#include <vector>
#include <unordered_map>
#include <memory>
#include <iostream>
#include <array>
#include <functional>
#include <Eigen/Dense>
#include <limits>
#include <algorithm>
#include <fstream>
#include "../geom/geom_func.h"
// Type alias for a triangle simplex (each triangle is a 3x3 matrix with 3 vertices in 3D).
using Simplex = Eigen::Matrix<double, 3, 3>;

// Type alias for a bounding box.
// Format: [min_x, min_y, min_z, max_x, max_y, max_z]
using BoundingBox = std::array<double, 6>;

/**
 * @brief Structure representing an element (polygon) with its triangulated simplices,
 *        its computed bounding box, and the identifiers returned by the spatial tree.
 *
 * For a triangle input, the vector 'simplices' will contain one simplex.
 * For a quadrilateral input, it will contain two simplices (obtained by triangulation).
 */
struct Element {
    std::vector<Simplex> simplices;  // Triangle simplices that represent the polygon.
    BoundingBox bbox;                // Axis-aligned bounding box computed from the polygon.
    std::vector<int> tree_ids;       // Identifiers returned by the BinaryTree for each simplex.
    std::vector<int> ids;

    Element()
        : simplices(), bbox({ 0, 0, 0, 0, 0, 0 }), tree_ids(), ids(){
    }
    Element(const std::vector<Simplex>& s, const BoundingBox& b, const std::vector<int>& treeids, const std::vector<int>& oids)
        : simplices(s), bbox(b), tree_ids(treeids) , ids(oids){
    }
};

/**
 * @brief IntersecChecker uses a spatial BinaryTree (from hexa_octree.hpp) to store
 *        polygon elements (triangles or quadrilaterals) and quickly check for intersections.
 *
 * In this implementation each element (a polygon assumed to be either a triangle or a quadrilateral)
 * is inserted into the binary tree after being triangulated. In addition, the element is stored
 * in a vector while an unordered_map maps a HexaTag to its index in the vector.
 */
class IntersecChecker {
public:
    IntersecChecker();
    ~IntersecChecker();

    void init(const BoundingBox box);

    // Adds a group of points and returns the starting index (in the points_ vector).
    int addPoint(const std::vector<BLVector>& coords);

    // Adds an element (polygon) with the specified tag and point indices.
    // Expects exactly three (triangle) or four (quadrilateral) indices.
    bool addElement(const HexaTag& tag, const std::vector<int>& point_ids);
    bool addElement(const std::array<HexaTag, 2>& tags, const std::vector<int>& point_ids);

    // Adds an extra tag to an existing element specified by the origin_tag.
    // Returns true if successful; returns false if the origin_tag is not found.
    bool add_extra_tag(const HexaTag& origin_tag, const HexaTag& extra_tag);

    // Batch addition of elements.
    void addElements(const std::vector<std::pair<HexaTag, std::vector<int>>>& elems);

    // Removes an element identified by its tag.
    void removeElement(const HexaTag& tag);

    // Checks if the given element (polygon specified by point indices) intersects any existing element.
    bool checkIntersect(const std::vector<int>& point_ids) const;

    // Retrieves an element by its tag.
    const Element* getElement(const HexaTag& tag) const;

    void exportToVtk(const std::string& filename) const;
private:
    // Computes the bounding box for a given polygon (vector of points).
    BoundingBox computeBoundingBox(const std::vector<BLVector>& poly) const;

    // Determines if two bounding boxes intersect.
    bool bboxIntersect(const BoundingBox& a, const BoundingBox& b) const;

    // Helper: Triangulates a polygon.
    // For a triangle (3 vertices), returns one simplex.
    // For a quadrilateral (4 vertices), returns two simplices: (0,1,2) and (0,2,3).
    std::vector<Simplex> triangulate(const std::vector<BLVector>& poly) const;

    // Stores all points.
    std::vector<BLVector> points_;

    int tree_id;
    // The actual element storage, where each Element object is stored.
    std::vector<Element> element_storage_;

    std::vector<bool> isable_;

    // Hash map that maps a HexaTag to an index in element_storage_.
    std::unordered_map<HexaTag, size_t, HexaTagHash> tag_to_index_;

    // Pointer to the binary tree used for spatial partitioning.
    // The BinaryTree is templated with Simplex so that each inserted ¡°data¡±
    // contains a triangle and its computed bounding box.
    std::unique_ptr<HEXA::BinaryTree<Simplex>> tree_;
};

//
// Implementation
//

// Adds an extra tag by mapping extra_tag to the same element index as origin_tag.
inline bool IntersecChecker::add_extra_tag(const HexaTag& origin_tag, const HexaTag& extra_tag) {
    auto it = tag_to_index_.find(origin_tag);
    if (it == tag_to_index_.end())
        return false; // origin_tag not found

    if (tag_to_index_.find(extra_tag) != tag_to_index_.end())
        return false; // extra_tag already exists

    tag_to_index_[extra_tag] = it->second;
    return true;
}

// Constructor: Initializes the binary tree with a sufficiently large initial bounding box.
inline IntersecChecker::IntersecChecker() {
}

inline void IntersecChecker::init(const BoundingBox box) {
    GEOM_FUNC::exactinit();
    tree_ = std::make_unique<HEXA::BinaryTree<Simplex>>(box);
    tree_->kMaxAllowedTreeHeight = 12;
    tree_->kMaxAllowedNodeDataSize = 80;
    tree_id = 0;
}

// Destructor.
inline IntersecChecker::~IntersecChecker() {
    // unique_ptr automatically releases memory.
}

// Adds a group of points and returns the starting index.
inline int IntersecChecker::addPoint(const std::vector<BLVector>& coords) {
    int start_id = static_cast<int>(points_.size());
    points_.insert(points_.end(), coords.begin(), coords.end());
    return start_id;
}

// Computes the axis-aligned bounding box for a given polygon (vector of points).
inline BoundingBox IntersecChecker::computeBoundingBox(const std::vector<BLVector>& poly) const {
    if (poly.empty()) return { 0, 0, 0, 0, 0, 0 };
    BoundingBox bbox = { poly[0].x, poly[0].y, poly[0].z, poly[0].x, poly[0].y, poly[0].z };
    for (const auto& p : poly) {
        bbox[0] = std::min(bbox[0], p.x);
        bbox[1] = std::min(bbox[1], p.y);
        bbox[2] = std::min(bbox[2], p.z);
        bbox[3] = std::max(bbox[3], p.x);
        bbox[4] = std::max(bbox[4], p.y);
        bbox[5] = std::max(bbox[5], p.z);
    }
    return bbox;
}

// Determines if two bounding boxes intersect.
inline bool IntersecChecker::bboxIntersect(const BoundingBox& a, const BoundingBox& b) const {
    return (a[0] <= b[3] && a[3] >= b[0]) &&
        (a[1] <= b[4] && a[4] >= b[1]) &&
        (a[2] <= b[5] && a[5] >= b[2]);
}

// Triangulates a polygon. For 3 vertices, returns one triangle; for 4 vertices, returns two triangles.
inline std::vector<Simplex> IntersecChecker::triangulate(const std::vector<BLVector>& poly) const {
    std::vector<Simplex> triangles;
    if (poly.size() == 3) {
        Simplex s;
        for (int i = 0; i < 3; ++i) {
            s.row(i) = Eigen::RowVector3d(poly[i].x, poly[i].y, poly[i].z);
        }
        triangles.push_back(s);
    }
    else if (poly.size() == 4) {
        // First triangle: vertices 0, 1, 2.
        Simplex s1;
        s1.row(0) = Eigen::RowVector3d(poly[0].x, poly[0].y, poly[0].z);
        s1.row(1) = Eigen::RowVector3d(poly[1].x, poly[1].y, poly[1].z);
        s1.row(2) = Eigen::RowVector3d(poly[2].x, poly[2].y, poly[2].z);
        triangles.push_back(s1);
        // Second triangle: vertices 0, 2, 3.
        Simplex s2;
        s2.row(0) = Eigen::RowVector3d(poly[0].x, poly[0].y, poly[0].z);
        s2.row(1) = Eigen::RowVector3d(poly[2].x, poly[2].y, poly[2].z);
        s2.row(2) = Eigen::RowVector3d(poly[3].x, poly[3].y, poly[3].z);
        triangles.push_back(s2);
    }
    return triangles;
}

// Adds an element (polygon) given a tag and point indices.
// The polygon may have 3 vertices (triangle) or 4 vertices (quadrilateral).
inline bool IntersecChecker::addElement(const HexaTag& tag, const std::vector<int>& point_ids) {
    if (point_ids.size() != 3 && point_ids.size() != 4) {
        std::cerr << "Element must have exactly 3 (triangle) or 4 (quadrilateral) point indices.\n";
        return false;
    }
    if (tag_to_index_.find(tag) != tag_to_index_.end()) {
        std::cerr << "Element with the given tag already exists.\n";
        return false;
    }
    // Build the polygon (vector of points) from the given point indices.
    std::vector<BLVector> poly;
    for (int pid : point_ids) {
        if (pid < 0 || pid >= static_cast<int>(points_.size())) {
            std::cerr << "Invalid point index.\n";
            return false;
        }
        poly.push_back(points_[pid]);
    }
    // Triangulate the polygon.
    std::vector<Simplex> triangles = triangulate(poly);
    if (triangles.empty()) {
        std::cerr << "Failed to triangulate the polygon.\n";
        return false;
    }
    // Compute the bounding box for the polygon.
    BoundingBox bbox = computeBoundingBox(poly);
    // Insert each triangle into the spatial tree and record its tree ID.
    std::vector<int> tree_ids;
    
    for (const auto& s : triangles) {
        int id = tree_->insert(s, tree_id++);  // The second parameter is a dummy id; modify as needed.
        tree_ids.push_back(id);
    }
    // Create and store the element.
    Element elem(triangles, bbox, tree_ids, point_ids);
    size_t index = element_storage_.size();
    element_storage_.push_back(elem);
    tag_to_index_[tag] = index;
    isable_.push_back(true);
    
    //if (checkIntersect(point_ids))
    //    std::cerr << " intersection check error";


    return true;
}

inline bool IntersecChecker::addElement(const std::array<HexaTag, 2>& tags, const std::vector<int>& point_ids) {
    return addElement(tags[0], point_ids) && add_extra_tag(tags[0], tags[1]);
}

// Batch addition of elements.
inline void IntersecChecker::addElements(const std::vector<std::pair<HexaTag, std::vector<int>>>& elems) {
    for (const auto& pair : elems) {
        addElement(pair.first, pair.second);
    }
}

// Removes an element identified by its tag.
inline void IntersecChecker::removeElement(const HexaTag& tag) {
    auto it = tag_to_index_.find(tag);
    if (it != tag_to_index_.end()) {
        // Remove each triangle from the binary tree.
        for (int id : element_storage_[it->second].tree_ids) {
            tree_->remove(id);
        }
        isable_[it->second] = false;
        // Note: Here only the mapping for the provided tag is removed.
        // If one element corresponds to multiple tags, additional tags must be cleared.
        tag_to_index_.erase(it);
        // If removal from element_storage_ is required, one may use swap-and-pop strategies
        // and update other tags accordingly.
        
    }
    else {
        std::cerr << "Element with the given tag does not exist.\n";
    }
}
// Checks if the given element (polygon specified by point indices) intersects any existing element.
inline bool IntersecChecker::checkIntersect(const std::vector<int>& point_ids) const {
    if (point_ids.size() != 3 && point_ids.size() != 4) {
        std::cerr << "Element must have exactly 3 (triangle) or 4 (quadrilateral) point indices.\n";
        return false;
    }
    // Build the query polygon.
    std::vector<BLVector> poly;
    for (int pid : point_ids) {
        if (pid < 0 || pid >= static_cast<int>(points_.size())) {
            std::cerr << "Invalid point index.\n";
            return false;
        }
        poly.push_back(points_[pid]);
    }
    // Triangulate the query polygon.
    std::vector<Simplex> queryTriangles = triangulate(poly);
    if (queryTriangles.empty()) {
        std::cerr << "Failed to triangulate the query polygon.\n";
        return false;
    }
    // Compute the query polygon's bounding box.
    BoundingBox bbox = computeBoundingBox(poly);
    // Query the binary tree for candidate indices whose bounding boxes intersect.
    std::vector<size_t> candidate_indices;
    tree_->query(bbox, candidate_indices);
    // For each candidate, if the bounding boxes intersect, perform precise triangle-triangle intersection tests.
    for (size_t idx : candidate_indices) {
        const auto& candidate_box = tree_->tree_data_array_[idx].box;
        if (bboxIntersect(bbox, candidate_box)) {
            // Extract the candidate triangle simplex.
            const Simplex& candidateSimplex = tree_->tree_data_array_[idx].simplex;
            // Convert the candidate simplex vertices to BLVector type.
            BLVector candV[3], queryV[3];
            candV[0] = {candidateSimplex(0, 0), candidateSimplex(0, 1), candidateSimplex(0, 2)};
            candV[1] = { candidateSimplex(1, 0), candidateSimplex(1, 1), candidateSimplex(1, 2) };
            candV[2] = { candidateSimplex(2, 0), candidateSimplex(2, 1), candidateSimplex(2, 2) };
            // Check each query triangle against the candidate triangle.
            for (const auto& s : queryTriangles) {
                queryV[0] = {s(0, 0), s(0, 1), s(0, 2)};
                queryV[1] = { s(1, 0), s(1, 1), s(1, 2) };
                queryV[2] = { s(2, 0), s(2, 1), s(2, 2) };
                // Call tri_tri_intersect (expected signature:
                // bool tri_tri_intersect(const BLVector&, const BLVector&, const BLVector&, const BLVector&, const BLVector&, const BLVector&))
                /*if (intersection3d::tri_tri_intersect(queryV0, queryV1, queryV2, candV0, candV1, candV2))
                    return true;*/
                int same_count = 0;
                double p[6][3];
                int inter1, inter2;
                for (int k = 0; k < 3; k++) {
                    for (int l = 0; l < 3; l++) {
                        if (candV[k] == queryV[l]) {
                            same_count++;
                            inter1 = k; inter2 = l;
                        }
                        p[k][l] = candV[k][l];
                        p[k+3][l] = queryV[k][l];
                    }
                    
                }
                if (same_count == 3) {
                    continue;
                }
                if (same_count==0&&GEOM_FUNC::tri_tri_overlap_test_3d(p[0], p[1], p[2], p[3], p[4], p[5]))
                    return true;
               // TiGER_GEOM_FUNC::lin_tri_intersect3d();
                double linep[2][3]; double facep[3][3];
                int intTyp; int intCod;
                double intPnt[3]; bool bEpsilon=false;
                if (same_count == 1) {
                    for (int k = 0; k < 3; k++) {
                        linep[0][k] = p[(inter1 + 1) % 3][k];
                        linep[1][k] = p[(inter1 + 2) % 3][k];
                        facep[0][k] = p[3][k];
                        facep[1][k] = p[4][k];
                        facep[2][k] = p[5][k];
                    }
                    if(GEOM_FUNC::lin_tri_intersect3d(linep, facep, &intTyp, &intCod,intPnt, bEpsilon))
                        return true;
                    for (int k = 0; k < 3; k++) {
                        linep[0][k] = p[(inter2 + 1) % 3+3][k];
                        linep[1][k] = p[(inter2 + 2) % 3+3][k];
                        facep[0][k] = p[0][k];
                        facep[1][k] = p[1][k];
                        facep[2][k] = p[2][k];
                    }
                    if (GEOM_FUNC::lin_tri_intersect3d(linep, facep, &intTyp, &intCod, intPnt, bEpsilon))
                        return true;
                }
               
            }
        }
    }
    return false;
}

// Retrieves an element by its tag. Returns a pointer to the element or nullptr if not found.
inline const Element* IntersecChecker::getElement(const HexaTag& tag) const {
    auto it = tag_to_index_.find(tag);
    if (it != tag_to_index_.end() && it->second < element_storage_.size())
        return &element_storage_[it->second];
    return nullptr;
}
inline void IntersecChecker::exportToVtk(const std::string& filename) const {
    std::ofstream vtk_file(filename);

    if (!vtk_file.is_open()) {
        std::cerr << "Error opening file for writing: " << filename << '\n';
        return;
    }

    // Write VTK file header
    vtk_file << "# vtk DataFile Version 3.0\n";
    vtk_file << "Polygon data exported from IntersecChecker\n";
    vtk_file << "ASCII\n";
    vtk_file << "DATASET UNSTRUCTURED_GRID\n";

    // Write points
    vtk_file << "POINTS " << points_.size() << " float\n";
    for (const auto& point : points_) {
        vtk_file << point.x << " " << point.y << " " << point.z << "\n";
    }

    // Count the number of valid cells and total indices
    int cell_count = 0;
    int index_count = 0;
    for (size_t i = 0; i < isable_.size(); i++) {
        if (isable_[i]) {
            cell_count++;
            index_count += element_storage_[i].ids.size();
        }
    }

    // Write cells (each cell consists of its size followed by point indices)
    vtk_file << "CELLS " << cell_count << " " << (cell_count + index_count) << "\n";
    for (size_t i = 0; i < isable_.size(); i++) {
        if (isable_[i]) {
            vtk_file << element_storage_[i].ids.size() << " ";
            for (auto k : element_storage_[i].ids) {
                vtk_file << k << " ";
            }
            vtk_file << "\n";
        }
    }

    // Write cell types (assuming all are polygons or triangles)
    vtk_file << "CELL_TYPES " << cell_count << "\n";
    for (size_t i = 0; i < isable_.size(); i++) {
        if (isable_[i]) {
            size_t vertices = element_storage_[i].ids.size();
            if (vertices == 3) {
                vtk_file << "5\n";  // VTK_TRIANGLE
            }
            else if (vertices == 4) {
                vtk_file << "9\n";  // VTK_QUAD
            }
            else {
                vtk_file << "7\n";  // VTK_POLYGON (for >4 vertices)
            }
        }
    }

    vtk_file.close();
    std::cout << "Successfully exported to VTK: " << filename << '\n';
}

#endif // INTERSECTION_CHECK_HPP
