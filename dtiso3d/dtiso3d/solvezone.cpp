#include "solvezone.h"
#include <spdlog/spdlog.h> 
 #include <map>
#include <array>
#include <queue>
#include <algorithm>
#include <set>
#include "iso3d.h"
#include "fileio.h"

using namespace std;

int ZoneSolver::divideZone(const std::vector<int>& bdTris, const std::vector<std::vector<int>>& trisToBlock)
{
	int nElem = elems.size() / 4;
	int nPt = pts.size() / 3;
	int nTri = bdTris.size() / 3;

	std::vector<std::map<std::pair<int, int>, int>> bd_graph(nPt);
	for (int i = 0; i < nTri; i++) {
		std::array<int, 3> p{ bdTris[3 * i + 0] ,bdTris[3 * i + 1] ,bdTris[3 * i + 2] };
		std::sort(p.begin(), p.end());
		bd_graph[p[0]][std::pair<int, int>(p[1], p[2])] = i;
	}

	//找到边界角点，角边界三角单元
	auto cmpPt = [&](int p1, int p2) {
		if (pts[p1 * 3] != pts[p2 * 3])
			return pts[p1 * 3] > pts[p2 * 3];
		else if (pts[p1 * 3 + 1] != pts[p2 * 3 + 1])
			return pts[p1 * 3 + 1] > pts[p2 * 3 + 1];
		else
			return pts[p1 * 3 + 2] > pts[p2 * 3 + 2];
	};
	int cornerPt = 0;
	for (int i = 1; i < nPt; i++) {
		if (cmpPt(i, cornerPt))
			cornerPt = i;
	}
	int cornerTri = std::find(bdTris.begin(), bdTris.end(), cornerPt) - bdTris.begin();
	cornerTri /= 3;	
	struct ZoneSeed
	{
		int iElem;
		set<int> blocks;
	};
	std::queue<ZoneSeed> seeds;

	//建立邻接关系
	vector<vector<int>> neighbour(elems.size() / 4);
	vector<vector<int>> neighbourBnd(elems.size() / 4); //邻居对是否经过边界，是则值为边界索引，否则为-1
	std::vector<std::map<std::pair<int, int>, int>> graph(nPt);
	const int facelist[4][3] = { {0,1,2},{0,1,3},{0,2,3},{1,2,3} };
	for (int i = 0; i < nElem; i++) {
		std::array<int, 4> p{ elems[4 * i + 0],elems[4 * i + 1] ,elems[4 * i + 2] ,elems[4 * i + 3] };
		std::sort(p.begin(), p.end());
		for (int k = 0; k < 4; k++) {
			int key = p[facelist[k][0]];
			std::pair<int, int> spair = { p[facelist[k][1]], p[facelist[k][2]] };
			if (bd_graph[key].find(spair) != bd_graph[key].end() && bd_graph[key][spair] == cornerTri) {
				//找到角边界单元，设定该单元为初始种子
				ZoneSeed initsd;
				initsd.iElem = i;
				initsd.blocks.insert(trisToBlock[cornerTri].begin(), trisToBlock[cornerTri].end());
				seeds.push(initsd);
			}
			if (graph[key].find(spair) == graph[key].end())
				graph[key][spair] = i;
			else {
				neighbour[i].push_back(graph[key][spair]);
				neighbour[graph[key][spair]].push_back(i);
				if (bd_graph[key].find(spair) == bd_graph[key].end()) {
					neighbourBnd[i].push_back(-1);
					neighbourBnd[graph[key][spair]].push_back(-1);
				}
				else {
					neighbourBnd[i].push_back(bd_graph[key][spair]);
					neighbourBnd[graph[key][spair]].push_back(bd_graph[key][spair]);
				}
			}
		}
	}

	int count = 0;

	// bfs
	std::vector<bool> visited(nElem, false);
	std::vector<bool> seeded(nElem, false);
	seeded[seeds.front().iElem] = true;
	
	int zoneId = 0;
	elemToZone.resize(elems.size() / 4);
	int ct = 0;

	while (!seeds.empty())
	{
		ZoneSeed cursd = seeds.front();
		seeds.pop();
		if(visited[cursd.iElem]) continue;
		visited[cursd.iElem] = true;
		zoneToBlock.emplace_back(cursd.blocks.begin(), cursd.blocks.end());

		std::queue<int> q;
		q.push(cursd.iElem);
		while (!q.empty())
		{
			int curElem = q.front();
			q.pop();
			elemToZone[curElem] = zoneId;
			for (int j = 0; j < neighbour[curElem].size(); j++) {
				int jElem = neighbour[curElem][j];
				if (!visited[jElem]) {
					int coBnd = neighbourBnd[curElem][j];
					if (coBnd == -1) {
						//相邻单元不穿过边界面，为区域内部单元
						q.push(jElem);
						visited[jElem] = true;
					}
					else if(!seeded[jElem]){
						seeded[jElem] = true;
						ZoneSeed newsd;
						newsd.iElem = jElem;
						newsd.blocks = cursd.blocks;
						//如果已在block中，则边界为出边界，否则为入边界
						for (const auto& b : trisToBlock[coBnd]) {
							if (newsd.blocks.find(b) == newsd.blocks.end())
								newsd.blocks.insert(b);
							else
								newsd.blocks.erase(b);
						}
						seeds.push(newsd);
					}
				}
			}
		}
		zoneId++;
	}

	for (const auto v : visited) {
		if (!v) {
			throw std::logic_error("Error in zone dividing!");
		}
	}

	return 0;
}

int ZoneSolver::deleteBlock(int block)
{
	std::vector<int> newElems;
	std::vector<int> newElemToZone;
	std::vector<int> newSeq;
	std::vector<std::vector<int>> newZoneToBlock;

	std::vector<int> zoneToDelete;
	std::vector<int> zoneIndexTable;

	int deletedZoneNum = 0;
	for (int i = 0; i < this->zoneToBlock.size(); i++)
	{
		int flag = 0;
		for (int j = 0; j < zoneToBlock[i].size(); j++)
		{
			if (zoneToBlock[i][j] == block)
			{
				zoneToDelete.push_back(i);
				flag = 1;
				break;
			}
		}
		if (!flag)
		{
			newZoneToBlock.push_back(zoneToBlock[i]);
			zoneIndexTable.push_back(i - deletedZoneNum);
		}
		else
		{
			zoneIndexTable.push_back(-1);
			deletedZoneNum++;
		}
	}

	for (int i = 0; i < this->elemToZone.size(); i++)
	{
		int flag = 0;
		for (int j = 0; j < zoneToDelete.size(); j++)
		{
			if (zoneToDelete[j] == elemToZone[i])
			{
				flag = 1;
				break;
			}
		}
		if (!flag)
		{
			newElemToZone.push_back(elemToZone[i]);
			newElems.push_back(elems[4 * i]);
			newElems.push_back(elems[4 * i + 1]);
			newElems.push_back(elems[4 * i + 2]);
			newElems.push_back(elems[4 * i + 3]);
		}
	}
	for (int i = 0; i < newElemToZone.size(); i++)
	{
		newElemToZone[i] = zoneIndexTable[newElemToZone[i]];
	}


	elems = newElems;
	elemToZone = newElemToZone;
	zoneToBlock = newZoneToBlock;
	return 0;
}

ZoneSolver::ZoneSolver(BlockedArray_Node* pnodes, BlockedArray_Elem* pelems)
{
	pts.resize(pnodes->getArraySize() * 3);
	elems.resize(pelems->getArraySize() * 4);
	for (int i = 0; i < pnodes->getArraySize(); i++) {
		pts[3 * i] = (*pnodes)[i].pt[0];
		pts[3 * i + 1] = (*pnodes)[i].pt[1];
		pts[3 * i + 2] = (*pnodes)[i].pt[2];
	}
	for (int i = 0; i < pelems->getArraySize(); i++) {
		elems[4 * i] = (*pelems)[i].form[0];
		elems[4 * i + 1] = (*pelems)[i].form[1];
		elems[4 * i + 2] = (*pelems)[i].form[2];
		elems[4 * i + 3] = (*pelems)[i].form[3];
	}
}

/*
 * 在网格中去除没有用到的点
 */
int ZoneSolver::getShrinkedMesh(std::vector<double>& pts_out, std::vector<int>& elems_out, std::vector<int>& ptshash)
{
	std::vector<double> newPts;
	elems_out = elems;

	// record which pt needs to be deleted
	std::vector<int> ptsTag(pts.size() / 3, 0);

	// record the old and new pt index
	std::vector<int> ptsIndexTable;

	for (int i = 0; i < elems.size(); i++)
	{
		ptsTag[elems[i]] = -1;
	}
	// j is the number of the points which have been deleted
	int j = 0;
	for (int i = 0; i < ptsTag.size(); i++)
	{
		if (ptsTag[i] < 0) {
			newPts.push_back(pts[3*i]);
			newPts.push_back(pts[3*i + 1]);
			newPts.push_back(pts[3*i + 2]);
			ptsIndexTable.push_back(i - j);
		}
		else
		{
			ptsIndexTable.push_back(-1);
			j++;
		}
	}

	// convert the old to the new pts index in elems
	for (int j = 0; j < elems_out.size(); j++)
		elems_out[j] = ptsIndexTable[elems_out[j]];

	pts_out = newPts;
	ptshash = ptsIndexTable;
	return 0;
}

std::map<int, double> ZoneAttachInfo::getFaceSizes()
{
	std::map<int, double> ret;
	for (int i = 0; i < blockSize.size(); i++) {
		if (blockSize[i] > 0) {
			double sz = blockSize[i];
			for (int fc = block2FaceOffsets[i]; fc < block2FaceOffsets[i + 1]; fc++) {
				int face = block2Faces[fc];
				if (ret.find(face) == ret.end()) {
					ret[face] = blockSize[i];
				}
				else {
					ret[face] = min(blockSize[i], ret[face]);
				}
			}
		}
	}
	return ret;
}
