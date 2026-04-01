#pragma once
#include <array>
#include <cstdint>
#include <map>
#include <set>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

// #define NULL 0
#define DIM2      2
#define DIM3      3
#define NEIG_NULL -1
#define MAX_CONN  3

enum MeshType {
    MESH_2D = 0,
    MESH_3D = 1
};

struct Elm {
    int nconn;
    int conn[MAX_CONN];
    int neig[MAX_CONN];
    int igeom;
};

class MeshInfo {
public:
    MeshInfo(MeshType mtType)
        : m_mtType(mtType)
        , m_nElm(0)
        , m_nPt(0)
    {}

    ~MeshInfo() = default;

    void Initialize(int nPt,
                    int nElm,
                    const std::vector<std::array<double, 3>> &point_array,
                    const std::vector<std::array<int, 4>> &Elm);
    void CheckDuplicateFacesOrThrow();
    void buildPntElm();
    void buildElmNeig();

    std::vector<int> MeshInfo::getBoundaryPoints() const;
    std::vector<std::array<int, 2>> MeshInfo::getBoundaryEdges() const;
    std::vector<std::array<int, 2>> MeshInfo::getBoundaryEdges(int fidx) const;

    static uint64_t edgeKey(int a, int b)
    {
        if (a > b) {
            std::swap(a, b);
        }
        return (static_cast<uint64_t>(static_cast<uint32_t>(a)) << 32) | static_cast<uint32_t>(b);
    }

private:
    MeshType m_mtType;

    int m_nPt;
    int m_nElm;
    std::vector<std::array<double, 3>> m_ppoint;
    std::vector<Elm> m_pElm;

    std::vector<std::vector<int>> m_pointToElms;

    std::set<int> m_setBdryPt;
    std::multimap<int, int> m_mpaBdryElm;
};
