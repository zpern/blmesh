#include "iso3d_multidomain.h"

#include <spdlog/spdlog.h> 
 #include <iostream>
#include <vector>
#include <queue>
#include <map>
#include <set>
#include <array>
#include <algorithm>
#include "../third/tetgen/tetgen.h"
using namespace std;
using TETGEN_BLMESH_::tetgenio;
using TETGEN_BLMESH_::tetrahedralize;
using std::vector;
int splitMultiDomain(
	int bndPtNum, double bndPts[],
	int bndFctNum, int bndFcts[],
	int *domainNum,
	int *domainBndPtNum[], double* domainBndPts[],			//DIM of *domainBndPts = sum(*domainBndPtNum[i]) * 3
	int *domainBndFctNum[], int* domainBndFcts[],			//DIM of *domainBndFcts = sum(*domainBndFctNum[i]) * 3
	int *ltowBndPtMap[]									//DIM of *ltowBndPtMap = domainNum,  DIM of *ltowBndPtMap[domainIndex] = domainBndPtNum[domainIndex]
)
{
	TETGEN_BLMESH_::tetgenbehavior behavior;
	tetgenio surface;
	tetgenio volume;
	behavior.plc = 1;
	behavior.supsteiner_level = 0;
	behavior.steinerleft = 0;
	behavior.nobisect = 1;
	//behavior.nobisect_nomerge = 1;
	behavior.neighout = 1;
	behavior.nomergefacet = 1;
	behavior.nomergevertex = 1;
	vector<int> BndFcts, PtMap;
	vector<double> BndPts;

	double* coord = new double[bndPtNum * 3];
	for (int i = 0; i < 3 * bndPtNum; i++) {
		coord[i] = bndPts[i];
	}
	surface.pointlist = coord;
	surface.numberofpoints = bndPtNum;
	surface.numberoffacets = bndFctNum;
	surface.facetlist = new tetgenio::facet[bndFctNum];

	for (int i = 0; i < bndFctNum; i++) {
		auto f = &surface.facetlist[i];
		tetgenio::init(f);
		f->numberofpolygons = 1;
		f->polygonlist = new tetgenio::polygon[1];
		auto p = &f->polygonlist[0];
		tetgenio::init(p);
		// Read the number of vertices, it should be greater than 0.
		p->numberofvertices = 3;
		p->vertexlist = new int[3];
		//cout << i + 1 << " ";
		for (int k = -0; k < 3; k++) {
			p->vertexlist[k] = bndFcts[3 * i + k];
		//	cout << bndFcts[3 * i + k] << " ";
		}
		//cout << endl;
	}
	tetrahedralize(&behavior, &surface, &volume);


	//build connective graph
	std::vector<std::map<std::pair<int, int>,int>> ban_graph(volume.numberofpoints);
	std::vector<std::map<std::pair<int, int>, int>> graph(volume.numberofpoints);
	for (int i = 0; i < bndFctNum; i++) {
		std::array<int, 3> p{ bndFcts[3 * i + 0] ,bndFcts[3 * i + 1] ,bndFcts[3 * i + 2] };
		std::sort(p.begin(), p.end());
		ban_graph[p[0]][std::pair<int, int>(p[1], p[2])]=i;
	}
	vector<vector<int>> neigbhour(volume.numberoftetrahedra);
	
	const int facelist[4][3] = { {0,1,2},{0,1,3},{0,2,3},{1,2,3} };
	
	for (int i = 0; i < volume.numberoftetrahedra; i++) {
		std::array<int, 4> p{volume.tetrahedronlist[4*i+0],volume.tetrahedronlist[4 * i + 1] ,volume.tetrahedronlist[4 * i + 2] ,volume.tetrahedronlist[4 * i + 3] };
		std::sort(p.begin(), p.end());
		for (int k = 0; k < 4; k++) {
			int key = p[facelist[k][0]];
			std::pair<int, int> spair = { p[facelist[k][1]], p[facelist[k][2]] };
			if (ban_graph[key].find(spair) != ban_graph[key].end()) {
				continue;
			}
			if (graph[key].find(spair) == graph[key].end() )
				graph[key][spair] = i;
			else {
				neigbhour[i].push_back(graph[key][spair]);
				neigbhour[graph[key][spair]].push_back(i);
			}
		}
	}

	std::vector<std::vector<int>> ans;
	//bfs
	int before = -1;
	std::vector<bool> isread(volume.numberoftetrahedra, false);
	while (true) {
		int i  ;
		for (i = before + 1; i < volume.numberoftetrahedra; i++) {
			if (!isread[i])
				break;
		}
		before = i;
		if (i == volume.numberoftetrahedra)
			break;
		std::queue<int> q;
		q.push(i);
		ans.push_back(std::vector<int>());
		isread[i] = true;
		while (!q.empty()) {
			int c = q.front();
			q.pop();
			ans.rbegin()->push_back(c);
			for (auto j : neigbhour[c]) {
				if (!isread[j]) {
					q.push(j);
					isread[j] = true;
				}
			}
		}
	}

	*domainNum=ans.size();
	*domainBndFctNum = new int[ans.size()];
	*domainBndPtNum = new int[ans.size()];

	for (int i = 0; i < ans.size(); i++) {
		map<int, int> loc_map;
		int count = 0;
		vector<std::array<int,3>> face_list;

		for (auto j : ans[i]) {
			std::array<int, 4> arr{ volume.tetrahedronlist[4 * j + 0],volume.tetrahedronlist[4 * j + 1],volume.tetrahedronlist[4 * j + 2],volume.tetrahedronlist[4 * j + 3] };
			sort(arr.begin(), arr.end());
			for (int k = 0; k < 4; k++) {
				std::pair<int, int> spair = { arr[facelist[k][1]], arr[facelist[k][2]] };
				if (ban_graph[arr[facelist[k][0]]].find(spair) != ban_graph[arr[facelist[k][0]]].end()) {
					std::array<int, 3> tri{ arr[facelist[k][0]], arr[facelist[k][1]], arr[facelist[k][2]] };
					face_list.push_back(tri);
				}
			}
		}


		for (auto j : face_list) {
			for (int k =0; k < 3; k++) {
				if (loc_map.find(j[k]) == loc_map.end()) {
					loc_map[j[k]] = count++;
				}
			}
		}
		(*domainBndPtNum)[i] = loc_map.size();
		int old_size = BndPts.size();
		BndPts.resize(BndPts.size() + 3*loc_map.size());
		int old_p_size = PtMap.size();
		PtMap.resize(PtMap.size() + loc_map.size());
		
		for (auto j : loc_map) {
			PtMap[old_p_size +j.second] = j.first;
			for(int k=0;k<3;k++)
				BndPts[old_size+3*j.second+k] = volume.pointlist[3*j.first+k];
		}
		int face_count = 0;
		for (auto j : face_list) {
			BndFcts.push_back(loc_map[j[0]]);
			BndFcts.push_back(loc_map[j[1]]);
			BndFcts.push_back(loc_map[j[2]]);
		}

		(*domainBndFctNum)[i] = face_list.size();
	}

	*domainBndFcts = new int [BndFcts.size()];
	for (int i = 0; i < BndFcts.size(); i++)
		(*domainBndFcts)[i] = BndFcts[i];

	*ltowBndPtMap = new int[PtMap.size()];
	for (int i = 0; i < PtMap.size(); i++)
		(*ltowBndPtMap)[i] = PtMap[i];

	*domainBndPts = new double[BndPts.size()];
	for (int i = 0; i < BndPts.size(); i++)
		(*domainBndPts)[i] = BndPts[i];
	return 0;
}