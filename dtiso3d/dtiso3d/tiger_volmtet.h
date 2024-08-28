#pragma once
#include <mutex>
#include <memory>
#include "solvezone.h"

class TigerVolmTet
{
public:
	TigerVolmTet();
	~TigerVolmTet();

	void lock() { mutex.lock(); }
	void unlock() { mutex.unlock(); }

	void freeMemory();

	int initPnt(int nP);
	int initTet(int nT);
	int setCoords(int nP, double *pntCoord);
	int setTets(int nT, int *tetConn);
	virtual int getPntNum();
	int getTetNum();
	virtual double* getCoords();
	virtual int getCoord(int idx, double &x, double &y, double &z);
	int getTet(int idx, int &a, int &b, int &c, int &d);
	int setCoord(int idx, double x, double y, double z);
	int setTet(int idx, int a, int b, int c, int d);
	int *getTets();
	int testQuality();
protected:
	int allocCoordMem();
	int allocTetMem();
	

protected:
	std::mutex mutex;
	int nPnt;
	int nTet;
	double *coord;
	int* tets;
};

class TigerMesh: public TigerVolmTet
{
public:
	TigerMesh();
	~TigerMesh();
	void setZoneInfo(std::shared_ptr<ZoneAttachInfo> zoneinfo);
	std::shared_ptr<ZoneAttachInfo> getZoneInfo() { return zoneinfo; }
	void setZoneOut(std::shared_ptr<ZoneOutInfo> zoneout);
	std::shared_ptr<ZoneOutInfo> getZoneOut() { return zoneout; }

	int getSurTriNum();
	int* getTris();
	int getTri(int idx, int &x, int &y, int &z);
	int initTri(int nT);
	int setTri(int nT, int *tetConn);
	int setTri(int idx, int n0, int n1, int n2);

	int initMidPoint();
	int getPntNum(bool withMidPoint);
	virtual int getCoord(int idx, double &x, double &y, double &z);
	int getTetMP(int idx, int pt[6]);

protected:
	int nSurTri;
	int* surTris;

	int nMidPoint;
	double* midPointCoord;
	int* tet2MidPoint;

	std::shared_ptr<ZoneAttachInfo> zoneinfo;
	std::shared_ptr<ZoneOutInfo> zoneout;
};

