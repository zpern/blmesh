#include "tiger_volmtet.h"
#include <spdlog/spdlog.h> 
 #include <cstdio>
#include <array>
#include "tools/EdgeGrapher.h"
//#include "MeshMetric_dll.h"

TigerVolmTet::TigerVolmTet()
{
	nPnt = 0;
	nTet = 0;
	coord = nullptr;
	tets = nullptr;
}


TigerVolmTet::~TigerVolmTet()
{
}

int TigerVolmTet::setCoords(int nP, double * pntCoord)
{
	nPnt = nP;

	if (coord != nullptr)
	{
		delete[] coord;
		coord = nullptr;
	}

	if (nPnt > 0)
	{
		coord = new double[nPnt * 3];
	}
	else
	{
		spdlog::info("Wrong point number\n");
		return 1;
	}

	for (int i = 0; i < nPnt; i++)
	{
		coord[3 * i] = pntCoord[3 * i];
		coord[3 * i + 1] = pntCoord[3 * i + 1];
		coord[3 * i + 2] = pntCoord[3 * i + 2];		
	}

	return 0;
}

int TigerVolmTet::setTets(int nT, int * tetConn)
{
	nTet = nT;
	if (tets != nullptr)
	{
		delete[] tets;
		tets = nullptr;
	}

	if (nTet > 0)
	{
		tets = new int[nTet * 4];
	}
	else
	{
		//spdlog::info("Wrong tet number\n");
		return 1;
	}

	for (int i = 0; i < nTet; i++)
	{
		tets[4 * i] = tetConn[4 * i];
		tets[4 * i + 1] = tetConn[4 * i + 1];
		tets[4 * i + 2] = tetConn[4 * i + 2];
		tets[4 * i + 3] = tetConn[4 * i + 3];
	}

	return 0;
}




int TigerVolmTet::allocCoordMem()
{
	if (coord != nullptr)
	{
		delete[] coord;
		coord = nullptr;
	}

	if (nPnt > 0)
	{
		coord = new double[nPnt*3];
	}		
	else
	{
		spdlog::info("Wrong point number\n");
		return 1;
	}


	return 0;
}

int TigerVolmTet::allocTetMem()
{
	if (tets != nullptr)
	{
		delete[] tets;
		tets = nullptr;
	}

	if (nTet > 0)
	{
		tets = new int[nTet * 4];
	}
	else
	{
		//spdlog::info("Wrong tet number\n");
		return 1;
	}

	return 0;
}

void TigerVolmTet::freeMemory()
{
	if (tets != nullptr)
	{
		delete[] tets;
		tets = nullptr;
	}

	if (coord != nullptr)
	{
		delete[] coord;
		coord = nullptr;
	}
}

int TigerVolmTet::initPnt(int nP)
{	
	nPnt = nP;
	return allocCoordMem();
}

int TigerVolmTet::initTet(int nT)
{
	nTet = nT;	
	return allocTetMem();
}

int TigerVolmTet::getPntNum()
{
	return nPnt;
}

int TigerVolmTet::getTetNum()
{
	return nTet;
}

int TigerVolmTet::getCoord(int idx, double &x, double &y, double &z)
{
	if (idx <= 0 || idx > nPnt)
	{
		spdlog::info("Point index over size\n");
		return 1;
	}

	x = coord[3 * (idx - 1)];
	y = coord[3 * (idx - 1) + 1];
	z = coord[3 * (idx - 1) + 2];

	return 0;
}

int TigerVolmTet::getTet(int idx, int &a, int &b, int &c, int &d)
{
	if (idx <= 0 || idx > nTet)
	{
		spdlog::info("Tet index over size\n");
		return 1;
	}

	a = tets[4 * (idx - 1)];
	b = tets[4 * (idx - 1) + 1];
	c = tets[4 * (idx - 1) + 2];
	d = tets[4 * (idx - 1) + 3];
	
	return 0;
}

int TigerVolmTet::setCoord(int idx, double x, double y, double z)
{
	if (idx <= 0 || idx > nPnt)
	{
		spdlog::info("Point index over size\n");
		return 1;
	}
	coord[3 * (idx - 1)] = x;
	coord[3 * (idx - 1) + 1] = y;
	coord[3 * (idx - 1) + 2] = z;
	
	return 0;
}

int TigerVolmTet::setTet(int idx, int a, int b, int c, int d)
{
	if (idx <= 0 || idx > nTet)
	{
		spdlog::info("Tet index over size\n");
		return 1;
	}

	tets[4 * (idx - 1)] = a;
	tets[4 * (idx - 1) + 1] = b;
	tets[4 * (idx - 1) + 2] = c;
	tets[4 * (idx - 1) + 3] = d;
	
	return 0;
}

TigerMesh::TigerMesh()
	:TigerVolmTet()
{
	nSurTri = 0;
	surTris = nullptr;	

	nMidPoint = 0;
	midPointCoord = nullptr;
	tet2MidPoint = nullptr;

	zoneinfo = std::make_shared<ZoneAttachInfo>();
	zoneinfo->valid = false;
}

TigerMesh::~TigerMesh()
{
	if (surTris) delete[] surTris;
	if (midPointCoord) delete[] midPointCoord;
	if (tet2MidPoint) delete[] tet2MidPoint;
}

void TigerMesh::setZoneInfo(std::shared_ptr<ZoneAttachInfo> azoneinfo)
{
	mutex.lock();
	zoneinfo = azoneinfo;
	mutex.unlock();
}

void TigerMesh::setZoneOut(std::shared_ptr<ZoneOutInfo> azoneout)
{
	zoneout = azoneout;
}

int TigerMesh::getSurTriNum()
{
	return nSurTri;
}

int* TigerMesh::getTris()
{
	return surTris;
}

int TigerMesh::getTri(int idx, int &x, int &y, int &z)
{
	x = surTris[idx * 3 + 0];
	y = surTris[idx * 3 + 1];
	z = surTris[idx * 3 + 2];
	return 0;
}

int TigerMesh::initTri(int nT)
{
	nSurTri = nT;
	surTris = new int[nT*3];
	return 0;
}

int TigerMesh::setTri(int nT, int *tetConn)
{
	surTris[nT * 3 + 0] = tetConn[0];
	surTris[nT * 3 + 1] = tetConn[1];
	surTris[nT * 3 + 2] = tetConn[2];
	return 0;
}

int TigerMesh::setTri(int idx, int n0, int n1, int n2){
	surTris[idx * 3 + 0] = n0;
	surTris[idx * 3 + 1] = n1;
	surTris[idx * 3 + 2] = n2;
	return 0;
}

int TigerMesh::initMidPoint()
{
	// 计算高阶中点信息
	tet2MidPoint = new int[nTet * 6];
	std::vector<std::array<double, 3>> cps;
	const int iedge[6][2] = {
		{0,1},
		{1,2},
		{2,0},
		{0,3},
		{1,3},
		{2,3} };
	EdgeGrapher<int> eg;
	for (int i = 0; i < nTet; i++) {
		for (int e = 0; e < 6; e++) {
			int n1 = tets[4 * i + iedge[e][0]];
			int n2 = tets[4 * i + iedge[e][1]];
			if (eg.hasEdge(n1, n2)) continue;
			std::array<double, 3> cp;
			for (int d = 0; d < 3; d++) {
				cp[d] = (coord[3 * n1 + d] + coord[3 * n2 + d]) / 2;
			}
			eg.addEdgeOverride(n1, n2, cps.size());
			cps.push_back(cp);
		}
	}
	nMidPoint = cps.size();
	midPointCoord = new double[nMidPoint * 3];
	for (int i = 0; i < nMidPoint; i++) {
		midPointCoord[i * 3 + 0] = cps[i][0];
		midPointCoord[i * 3 + 1] = cps[i][1];
		midPointCoord[i * 3 + 2] = cps[i][2];
	}

	for (int i = 0; i < nTet; i++) {
		for (int e = 0; e < 6; e++) {
			int n1 = tets[4 * i + iedge[e][0]];
			int n2 = tets[4 * i + iedge[e][1]];
			int icp = eg.getEdgeVal(n1, n2);
			tet2MidPoint[i * 6 + e] = icp;
		}
	}
	return 0;
}

int TigerMesh::getPntNum(bool withMidPoint)
{
	return withMidPoint?(nMidPoint + nPnt):(nPnt);
}

double* TigerVolmTet::getCoords()
{
	return coord;
}

int *TigerVolmTet::getTets(){
	return tets;
}

int TigerVolmTet::testQuality()
{
	return 0;// meshmetric_tetra(nPnt, coord, nTet, tets);
}

int TigerMesh::getCoord(int idx, double &x, double &y, double &z)
{
	if (idx < nPnt) {
		x = coord[idx * 3 + 0];
		y = coord[idx * 3 + 1];
		z = coord[idx * 3 + 2];
		return 0;
	}
	if (idx < nMidPoint + nPnt) {
		idx -= nPnt;
		x = midPointCoord[idx * 3 + 0];
		y = midPointCoord[idx * 3 + 1];
		z = midPointCoord[idx * 3 + 2];
		return 0;
	}
	return -1;
}

int TigerMesh::getTetMP(int idx, int pt[6])
{
	for (int i = 0; i < 6; i++) {
		pt[i] = tet2MidPoint[idx * 6 + i] + nPnt;
	}
	return 0;
}
