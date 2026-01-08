#include "MeshInfo.h"
#include <spdlog/spdlog.h> 
 #include <vector>
#include <iostream>
#include <stdexcept>
#include <sstream>

MeshInfo::~MeshInfo(void)
{
	delete[] m_pPntElm;
	delete[] m_pPntidx;
	if (m_pElm) {
		delete[]m_pElm;
		m_pElm = nullptr;
	}
}

void MeshInfo::RemoveBarePoint(int nPt, int nElm, int * Elm)
{
	std::vector<int> count(nPt,0);
	std::vector<int> map_point(nPt,-1);
	for (int i = 0; i < nElm; i++) {
		for (int k = 0; k < 3; k++)
		count[Elm[4 * i + k]]++;
	}
	int invalid_count = 0;
	for (int i = 0; i < nPt; i++) {
		if (!count[i]) {
			invalid_count++;
		}
		else {
			for (int k = 0; k < 3; k++)
				m_ppoint[3 * (i - invalid_count) + k] = m_ppoint[3 * (i) + k];
			map_point[i] = i - invalid_count;
		}
	}
	for (int i = 0; i < nElm; i++) {
		for (int k = 0; k < 3; k++)
			if (map_point[Elm[4 * i + k]]>=0)
				Elm[4 * i + k] = map_point[Elm[4 * i + k]];
	}
	for (int i = 0; i < nElm; i++) {
		for(int k=0;k< m_pElm[i].nconn;k++)
			if (map_point[m_pElm[i].conn[k]] >= 0)
				m_pElm[i].conn[k] = map_point[m_pElm[i].conn[k]];
			
	}
	m_nPt -= invalid_count;
	if (invalid_count) {
		spdlog::info("{} points removed",invalid_count);

	}

}
/*!
 * 
 * 
 * \param point_array
 * \param nPt
 * \param nElm
 * \param Elem
 */
void MeshInfo::Initialize(double* point_array,int nPt, int nElm, int *Elem)
{
	int i, j, nconn;
	m_nPt = nPt;
	m_nElm = nElm;
	m_ppoint = point_array;
	m_pElm = new Elm[m_nElm];
	for (i=0; i<nElm; i++)
	{
		if(m_mtType == MeshType::MESH_2D)
			nconn = 2;
		else
			nconn = 3;

		m_pElm[i].nconn = nconn;

		for (j=0; j<nconn; j++)
		{
			m_pElm[i].conn[j] = Elem[(nconn+1)*i + j];
			m_pElm[i].neig[j] = -1;
		}
		m_pElm[i].igeom = Elem[(nconn+1)*i + nconn];
	}
	RemoveBarePoint(nPt,nElm,Elem);
	//CheckDuplicateFacesOrThrow();
	CalPntElm();
	CalElmNeig();
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
void MeshInfo::CalPntElm()
{
	int i, j, start, *idx, *pelem;
	int npoints, nelm, dim;
	
	nelm = m_nElm;
	npoints = m_nPt;

	if(m_mtType == MeshType::MESH_3D)
		dim = DIM3;
	else
		dim = DIM2;

	idx = new int[npoints + 1];
	pelem = new int[nelm*dim];

	m_pPntidx = idx;
	m_pPntElm = pelem;
	
	for(i = 0; i <= npoints; i++)
		idx[i] = 0;
	
	/*----------------------------------------------------------------------------
	| Count the number of elements for each point:
	----------------------------------------------------------------------------*/
	for (i=0; i<nelm; i++)
	{
		for(j=0; j<dim; j++)
			idx[m_pElm[i].conn[j]]++;
	}
	
	/*----------------------------------------------------------------------------
	| From the numbers of elements compute the startindex for each point.
	| The elements of point 'i' are then stored in pelem[idx[i]] to
	| pelem[idx[i+1]]-1.
	----------------------------------------------------------------------------*/
	start = 0;
	for(i = 0; i <= npoints; i++)
	{
		int count = idx[i];
		if(count == 0 && i < npoints)
		{
			spdlog::info("Node {} is not connected with any element!\n", i);
			return ;
		}
		
		idx[i] = start;
		start += count;
	}
	
	/*----------------------------------------------------------------------------
	| Store the elements for each point. Thereby the 'idx' field is modified
	| and has to be restored afterwards.
	----------------------------------------------------------------------------*/
	for(i = 0; i < nelm; i++)
	{
		for(j = 0; j < dim; j++)
		{
			int pnt = m_pElm[i].conn[j];
			pelem[idx[pnt]] = i;
			idx[pnt]++;
		}
	}
	
	/*----------------------------------------------------------------------------
	| Restore 'idx' field:
	----------------------------------------------------------------------------*/
	start = 0;
	for(i = 0; i < npoints; i++)
	{
		int nstart = idx[i];
		idx[i] = start;
		start = nstart;
	}

}

void MeshInfo::CalElmNeig()
{
	int i, j, k, nsurf, beg, end;
	int pidx1, pidx2, eidx;
	int dim, nElm, *pntidx, *pntelm;

	nElm = m_nElm;
	pntidx = m_pPntidx;
	pntelm = m_pPntElm;

	if(m_mtType == MeshType::MESH_3D)
		dim = DIM3;
	else
		dim = DIM2;

	//遍历每个曲面建立网格单元邻接关系
	for (i=0; i<nElm; i++)
	{
		for (k=0; k<dim; k++)
		{
			if (m_pElm[i].neig[k] == -1)
			{
				//网格点(k+1)%dim的共享单元中找出k对应的邻接单元
				pidx1 = m_pElm[i].conn[(k+1)%dim];
				pidx2 = m_pElm[i].conn[(k+2)%dim];
				int pbeg = pntidx[pidx1];
				int pend = pntidx[pidx1+1];

				//遍历共享网格点pidx1的网格单元
				for (int j1=pbeg; j1<pend; j1++)
				{
					int peidx = pntelm[j1], lidx;

					if (peidx != i && IsElmEdge(peidx, pidx1, pidx2, &lidx))
					{
						m_pElm[i].neig[k] = peidx;
						m_pElm[peidx].neig[lidx] = i;
						break;
					}
				}
			}
		}
	}
	if (dim == DIM3)
{
    auto edge_key = [](uint32_t a, uint32_t b) -> uint64_t {
        if (a > b) std::swap(a, b);           // 无向边：排序
        return (uint64_t(a) << 32) | uint64_t(b);
    };

    // 每条边 -> 关联的(单元id, 该边在单元中的局部边号k)
    std::unordered_map<uint64_t, std::vector<std::pair<int,int>>> edge2inc;
    edge2inc.reserve(static_cast<size_t>(nElm) * 3);

    for (int ei = 0; ei < nElm; ++ei)
    {
        for (int k = 0; k < 3; ++k)
        {
            // 与你原逻辑保持一致：k 对应“对边”，边端点是 (k+1)%3 与 (k+2)%3
            uint32_t a = static_cast<uint32_t>(m_pElm[ei].conn[(k + 1) % 3]);
            uint32_t b = static_cast<uint32_t>(m_pElm[ei].conn[(k + 2) % 3]);
            edge2inc[edge_key(a, b)].push_back({ei, k});
        }
    }

    for (const auto& kv : edge2inc)
    {
        const uint64_t key = kv.first;
        const auto& inc = kv.second;

        if (inc.size() >= 3)
        {
            uint32_t p0 = static_cast<uint32_t>(key >> 32);
            uint32_t p1 = static_cast<uint32_t>(key & 0xffffffffu);

            std::ostringstream oss;
            oss << "[";
            for (size_t t = 0; t < inc.size(); ++t)
            {
                int eid = inc[t].first;
                oss << eid;
                if (t + 1 < inc.size()) oss << ",";
            }
            oss << "]";

            printf("error: non-manifold edge (%u, %u), incident elms = %s\n",
                   p0, p1, oss.str().c_str());
            spdlog::info("non-manifold edge ({}, {}), incident elms = {}\n",
                         p0, p1, oss.str());

            // 额外打印每个关联三角形，便于定位
            for (auto [eid, lk] : inc)
            {
                spdlog::info("  elm {} (igeom={}) conn=({}, {}, {}) localEdgeK={}\n",
                    eid, m_pElm[eid].igeom,
                    m_pElm[eid].conn[0], m_pElm[eid].conn[1], m_pElm[eid].conn[2], lk);
            }

            throw std::logic_error("non-manifold edge!");
        }
    }
}
	int ineig;
	for (i=0; i<nElm; i++)
	{
		for (j=0; j<dim; j++)
		{
			ineig = m_pElm[i].neig[j];
			if(ineig == NEIG_NULL)
			{
				printf("error neigh of input surface id = %d!\n",i);
				spdlog::info("the three point indexes of this triangle are {} {} {}\n", m_pElm[i].conn[0], m_pElm[i].conn[1], m_pElm[i].conn[2]);
				throw std::logic_error("non-closed input!");
				//exit(0);
			}
		}
	}

	ineig;
	for (i=0; i<nElm; i++)
	{
		for (j=0; j<dim; j++)
		{
			ineig = m_pElm[i].neig[j];
			if(ineig == NEIG_NULL)
				continue;

			for (k=0; k<dim; k++)
			{
				if(m_pElm[ineig].neig[k] == i)
					break;
			}

			if(m_pElm[i].igeom != m_pElm[ineig].igeom)
				m_pElm[i].neig[j] = m_pElm[ineig].neig[k] = -1;
		}
	}
}

bool MeshInfo::IsElmEdge(int eidx, int idx1, int idx2, int *pidx)
{
	int i, p1, p2;

	if (m_mtType == MeshType::MESH_2D)
	{
		for (i=0; i<DIM2; i++)
		{
			p1 = m_pElm[eidx].conn[i];
			*pidx = (i+1)%2;	//第三个点的局部编号

			if(p1 == idx1)
				return true;
		}
	} 
	else if(m_mtType == MeshType::MESH_3D)
	{
		for (i=0; i<DIM3; i++)
		{
			p1 = m_pElm[eidx].conn[i];
			p2 = m_pElm[eidx].conn[(i+1)%3];
			*pidx = (i+2)%3;	//第三个点的局部编号

			if((p1 == idx1 && p2 == idx2) || (p1 == idx2 && p2 == idx1))
				return true;
		}
	}
	
	return false;
}

void MeshInfo::GetBdryPt(int *npt, int **pt)
{
	int i, j, nElm, dim, ineig;
	int idx1, idx2;
	nElm = m_nElm;

	if(m_mtType == MeshType::MESH_3D)
		dim = DIM3;
	else
		dim = DIM2;

	for (i=0; i<nElm; i++)
	{
		for (j=0; j<dim; j++)
		{
			ineig = m_pElm[i].neig[j];
			if(ineig == NEIG_NULL)
			{
				if (m_mtType == MeshType::MESH_2D)
				{
					idx1 = m_pElm[i].conn[(j+1)%dim];

					m_setBdryPt.insert(idx1);
				} 
				else if(m_mtType == MeshType::MESH_3D)
				{
					idx1 = m_pElm[i].conn[(j+1)%dim];
					idx2 = m_pElm[i].conn[(j+2)%dim];

					m_setBdryPt.insert(idx1);
					m_setBdryPt.insert(idx2);
				}
			}
		}
	}

	*npt = m_setBdryPt.size();
	*pt = new int[*npt];

	i = 0;
	std::set<int>::iterator it = m_setBdryPt.begin();
	while(it != m_setBdryPt.end())
	{
		(*pt)[i++] = *it;
		++it;
	}
}

void MeshInfo::GetBdryElm(int *nbdry, int **bdry)
{
	int i, j, nElm, dim, ineig;
	int idx1, idx2, min, max;
	bool bexist;
	nElm = m_nElm;

	if(m_mtType == MeshType::MESH_3D)
		dim = DIM3;
	else
		dim = DIM2;

	m_mpaBdryElm.clear();

	for (i=0; i<nElm; i++)
	{
		for (j=0; j<dim; j++)
		{
			ineig = m_pElm[i].neig[j];
			if(ineig == NEIG_NULL)
			{
				idx1 = m_pElm[i].conn[(j+1)%dim];
				idx2 = m_pElm[i].conn[(j+2)%dim];
				min = idx1;
				max = idx2;
				if(idx1 > idx2)
				{ min = idx2; max = idx1; }
				bexist = false;
				std::multimap<int, int>::iterator beg, end;
				beg = m_mpaBdryElm.lower_bound(min);
				end = m_mpaBdryElm.upper_bound(min);
				while (beg != end)
				{
					if(beg->second == max)
					{
						bexist = true;
						break;
					}
					++beg;
				}

				if(!bexist)
					m_mpaBdryElm.insert(std::make_pair(min, max));
			}
		}
	}

	*nbdry = m_mpaBdryElm.size();
	*bdry = new int[(*nbdry)*2];

	i = 0;
	std::multimap<int, int>::iterator it = m_mpaBdryElm.begin();
	while(it != m_mpaBdryElm.end())
	{
		(*bdry)[i*2 + 0] = it->first;
		(*bdry)[i*2 + 1] = it->second;
		++i;
		++it;
	}
}

void MeshInfo::GetBdryElm(int *nbdry, int **bdry, int fidx)
{
	int i, j, nElm, dim, ineig;
	int idx1, idx2, min, max;
	bool bexist;
	nElm = m_nElm;

	if(m_mtType == MeshType::MESH_3D)
		dim = DIM3;
	else
		dim = DIM2;

	m_mpaBdryElm.clear();

	for (i=0; i<nElm; i++)
	{
		if(m_pElm[i].igeom != fidx)
			continue;

		for (j=0; j<dim; j++)
		{
			ineig = m_pElm[i].neig[j];
			if(ineig == NEIG_NULL)
			{
				idx1 = m_pElm[i].conn[(j+1)%dim];
				idx2 = m_pElm[i].conn[(j+2)%dim];
				min = idx1;
				max = idx2;
				if(idx1 > idx2)
				{ min = idx2; max = idx1; }
				bexist = false;
				std::multimap<int, int>::iterator beg, end;
				beg = m_mpaBdryElm.lower_bound(min);
				end = m_mpaBdryElm.upper_bound(min);
				while (beg != end)
				{
					if(beg->second == max)
					{
						bexist = true;
						break;
					}
					++beg;
				}

				if(!bexist)
					m_mpaBdryElm.insert(std::make_pair(min, max));
			}
		}
	}

	*nbdry = m_mpaBdryElm.size();
	*bdry = new int[(*nbdry)*2];

	i = 0;
	std::multimap<int, int>::iterator it = m_mpaBdryElm.begin();
	while(it != m_mpaBdryElm.end())
	{
		(*bdry)[i*2 + 0] = it->first;
		(*bdry)[i*2 + 1] = it->second;
		++i;
		++it;
	}
}
int MeshInfo::GetNumPoint()
{  return m_nPt; 
}
