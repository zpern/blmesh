#include "MeshInfo.h"
#include <spdlog/spdlog.h> 
 #include <vector>
#include <iostream>
#include <stdexcept>
#include <sstream>

void MeshInfo::Initialize(int nPt,
                          int nElm,
                          const std::vector<std::array<double, 3>> &point_array,
                          const std::vector<std::array<int, 4>> &Elem)
{
    m_nPt = nPt;
    m_nElm = nElm;
    m_ppoint = point_array;
    m_pElm.resize(m_nElm);

    const int nconn = (m_mtType == MeshType::MESH_2D) ? 2 : 3;

    for (int i = 0; i < m_nElm; ++i)
    {
        m_pElm[i].nconn = nconn;

        for (int j = 0; j < nconn; ++j)
        {
            m_pElm[i].conn[j] = Elem[i][j];
            m_pElm[i].neig[j] = -1;
        }

        m_pElm[i].igeom = Elem[i][nconn];
    }

    CheckDuplicateFacesOrThrow();
    buildPntElm();
    buildElmNeig();
}
void MeshInfo::CheckDuplicateFacesOrThrow()
{
    std::map<std::array<int,3>, int> seen; // key -> first elem index

    for (int i = 0; i < m_nElm; ++i) {
        // 只针对三角形面：conn[0..2]
        int a = m_pElm[i].conn[0];
        int b = m_pElm[i].conn[1];
        int c = m_pElm[i].conn[2];

        // 可选：退化三角形也报错
        if (a == b || b == c || a == c) {
            throw std::logic_error("degenerate triangle at elem " + std::to_string(i));
        }

		// 重叠三角形报错
        std::array<int,3> key{a,b,c};
        std::sort(key.begin(), key.end());  // 与顺序无关：{1,2,3} 和 {1,3,2} 会变成同一个 key

        auto it = seen.find(key);
        if (it != seen.end()) {
            int j = it->second;
            // 更详细也可以 spdlog 输出 i/j 的原始顺序、igeom 等
            throw std::logic_error(
                "overlapped faces: elem " + std::to_string(j) +
                " and elem " + std::to_string(i) +
                " share same vertices (" +
                std::to_string(key[0]) + "," +
                std::to_string(key[1]) + "," +
                std::to_string(key[2]) + ")"
            );
        }
        seen.emplace(key, i);
    }
}
void MeshInfo::buildPntElm()
{
    const int dim = (m_mtType == MeshType::MESH_3D) ? DIM3 : DIM2;

    m_pointToElms.resize(m_nPt);

    for (int i = 0; i < m_nElm; ++i)
    {
        for (int j = 0; j < dim; ++j)
        {
            int p = m_pElm[i].conn[j];
            m_pointToElms[p].push_back(i);
        }
    }

    for (int p = 0; p < m_nPt; ++p)
    {
        if (m_pointToElms[p].empty())
        {
            throw std::logic_error(
                "Node " + std::to_string(p) + " is not connected with any element!");
        }
    }
}

void MeshInfo::buildElmNeig()
{

    const int dim = (m_mtType == MeshType::MESH_3D) ? DIM3 : DIM2;


    // 边 -> 所有关联的(单元id, 局部边号)
    std::unordered_map<uint64_t, std::vector<std::pair<int, int>>> edge2inc;
    edge2inc.reserve(static_cast<size_t>(m_nElm) * dim);

    for (int i = 0; i < m_nElm; ++i)
    {
        for (int k = 0; k < dim; ++k)
        {
            int a, b;
            if (dim == 3)
            {
                a = m_pElm[i].conn[k % 3];
                b = m_pElm[i].conn[(k +1) % 3];
            }
            else
            {
                a = m_pElm[i].conn[k];
                b = m_pElm[i].conn[(k + 1) % 2];
            }

            edge2inc[edgeKey(a, b)].push_back({i, k});
        }
    }

    // 根据每条边的关联单元数，建立邻接
    for (const auto& kv : edge2inc)
    {
        const auto& inc = kv.second;

        if (inc.size() == 1)
        {
            // 只有一个单元占用这条边，说明开口，不闭合
            int ei = inc[0].first;
            int lk = inc[0].second;
            spdlog::info("open edge found at elm {}, local edge {}", ei, lk);
            throw std::logic_error("non-closed input!");
        }
        else if (inc.size() == 2)
        {
            int e0 = inc[0].first, k0 = inc[0].second;
            int e1 = inc[1].first, k1 = inc[1].second;

            m_pElm[e0].neig[k0] = e1;
            m_pElm[e1].neig[k1] = e0;
        }
        else
        {
            // 一条边被 3 个及以上单元共享，非流形
            throw std::logic_error("non-manifold edge!");
        }
    }

    // 如果不同igeom之间不算邻居，就在这里断开
    for (int ei = 0; ei < m_nElm; ++ei)
    {
        for (int k = 0; k < dim; ++k)
        {
            int nj = m_pElm[ei].neig[k];
            if (nj < 0) continue;

            if (m_pElm[ei].igeom != m_pElm[nj].igeom)
            {
                m_pElm[ei].neig[k] = -1;
            }
        }
    }
}


std::vector<int> MeshInfo::getBoundaryPoints() const
{
    const int dim = (m_mtType == MeshType::MESH_3D) ? DIM3 : DIM2;
    std::set<int> boundaryPoints;

    for (int i = 0; i < m_nElm; ++i)
    {
        for (int j = 0; j < dim; ++j)
        {
            if (m_pElm[i].neig[j] != NEIG_NULL)
                continue;

            if (dim == 2)
            {
                boundaryPoints.insert(m_pElm[i].conn[j  % dim]);
            }
            else
            {
                boundaryPoints.insert(m_pElm[i].conn[j  % dim]);
                boundaryPoints.insert(m_pElm[i].conn[(j + 1) % dim]);
            }
        }
    }

    return std::vector<int>(boundaryPoints.begin(), boundaryPoints.end());
}

//获取每个面的几何边界
std::vector<std::array<int, 2>> MeshInfo::getBoundaryEdges() const
{
    const int dim = (m_mtType == MeshType::MESH_3D) ? DIM3 : DIM2;
    std::set<std::pair<int,int>> edgeSet;

    for (int i = 0; i < m_nElm; ++i)
    {
        for (int j = 0; j < dim; ++j)
        {
            if (m_pElm[i].neig[j] != NEIG_NULL)
                continue;

            int a = m_pElm[i].conn[j  % dim];
            int b = m_pElm[i].conn[(j + 1) % dim];
            if (a > b) std::swap(a, b);

            edgeSet.insert({a, b});
        }
    }

    std::vector<std::array<int, 2>> edges;
    edges.reserve(edgeSet.size());

    for (const auto& e : edgeSet)
    {
        edges.push_back({e.first, e.second});
    }

    return edges;
}

//获取指定几何面的边界
std::vector<std::array<int, 2>> MeshInfo::getBoundaryEdges(int fidx) const
{
    const int dim = (m_mtType == MeshType::MESH_3D) ? DIM3 : DIM2;
    std::set<std::pair<int,int>> edgeSet;

    for (int i = 0; i < m_nElm; ++i)
    {
        if (m_pElm[i].igeom != fidx)
            continue;

        for (int j = 0; j < dim; ++j)
        {
            if (m_pElm[i].neig[j] != NEIG_NULL)
                continue;

            int a = m_pElm[i].conn[j  % dim];
            int b = m_pElm[i].conn[(j + 1) % dim];
            if (a > b) std::swap(a, b);

            edgeSet.insert({a, b});
        }
    }

    std::vector<std::array<int, 2>> edges;
    edges.reserve(edgeSet.size());

    for (const auto& e : edgeSet)
    {
        edges.push_back({e.first, e.second});
    }

    return edges;
}
