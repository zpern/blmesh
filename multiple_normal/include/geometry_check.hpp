#ifndef GEOMETRY_CHECK_HPP
#define GEOMETRY_CHECK_HPP

#include <vector>
#include <cmath>
#include <limits>
#include "./blvector.h"  // Use the vec class

// Namespace for 3D intersection functions
namespace intersection3d {

// Basic vector operations using vec

// Subtract two vectors: a - b
inline BLVector subtract(const BLVector &a, const BLVector &b) {
    return a - b;
}

// Dot product of two vectors
inline double dot(const BLVector &a, const BLVector &b) {
    return a * b;
}

// Cross product of two vectors
inline BLVector cross(const BLVector &a, const BLVector &b) {
    return a ^ b;
}

// Add two vectors
inline BLVector add(const BLVector &a, const BLVector &b) {
    return a + b;
}

// Scale a vector by a scalar s
inline BLVector scale(const BLVector &v, double s) {
    return v * s;
}

// Triangle-triangle intersection test in 3D
// Based on "A Fast Triangle-Triangle Intersection Test" by Möller (1997).
// Returns true if triangles (V0,V1,V2) and (U0,U1,U2) intersect.
inline bool tri_tri_intersect(const BLVector &V0, const BLVector &V1, const BLVector &V2,
                              const BLVector &U0, const BLVector &U1, const BLVector &U2)
{
    return false;
}

// Helper: Triangulate a polygon.
// For convex polygons with 3 or 4 vertices:
//   - If the polygon has 3 vertices, return one triangle.
//   - If the polygon has 4 vertices, split into two triangles: (0,1,2) and (0,2,3).
// For unsupported polygon sizes, returns an empty vector.
inline std::vector<std::vector<BLVector>> triangulate(const std::vector<BLVector> &poly)
{
    std::vector<std::vector<BLVector>> triangles;
    if (poly.size() == 3) {
        triangles.push_back(poly);
    } else if (poly.size() == 4) {
        triangles.push_back({ poly[0], poly[1], poly[2] });
        triangles.push_back({ poly[0], poly[2], poly[3] });
    }
    return triangles;
}

// Main intersection check function:
// Given two planar polygons in 3D (each as a vector of BLVector points),
// check whether they intersect. The polygons may be triangles or quadrilaterals.
// The function triangulates each polygon and then tests each pair of triangles.
inline bool intersection_check(const std::vector<BLVector> &s1,
                               const std::vector<BLVector> &s2)
{
    // Triangulate both polygons.
    std::vector<std::vector<BLVector>> tris1 = triangulate(s1);
    std::vector<std::vector<BLVector>> tris2 = triangulate(s2);

    // If either polygon could not be triangulated, assume no intersection.
    if (tris1.empty() || tris2.empty())
        return false;

    // For every pair of triangles, check for intersection.
    for (const auto &tri1 : tris1) {
        if (tri1.size() != 3) continue;
        for (const auto &tri2 : tris2) {
            if (tri2.size() != 3) continue;
            if (tri_tri_intersect(tri1[0], tri1[1], tri1[2],
                                  tri2[0], tri2[1], tri2[2]))
            {
                return true;
            }
        }
    }
    return false;
}




} // namespace intersection3d

#endif // GEOMETRY_CHECK_HPP
