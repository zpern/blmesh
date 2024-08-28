#pragma once


class TigerVolmTet
{
public:
	TigerVolmTet();
	~TigerVolmTet();

	int initPnt(int nP);
	int initTet(int nT);
	int setCoord(int nP, double *pntCoord);
	int setTet(int nT, int *tetConn);
	int setPntNum(int nP);
	int setTetNum(int nT);
	int getPntNum();
	int getTetNum();
	int getCoord(int idx, double &x, double &y, double &z);
	int getTet(int idx, int &a, int &b, int &c, int &d);
	int setCoord(int idx, double x, double y, double z);
	int setTet(int idx, int a, int b, int c, int d);
private:
	int allocCoordMem();
	int allocTetMem();


private:
	int nPnt;
	int nTet;
	double *coord;
	int* tets;
};

