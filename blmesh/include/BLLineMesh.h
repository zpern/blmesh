#pragma once
#include <cstdio>
#include <cstdlib>
class BLLineMesh
{
public:
	BLLineMesh(void);
	~BLLineMesh(void);

	void GenMeshLine(int nbpt, double *bpt, double *spa, int *npt, double **pt, int *nElm, int **pElm);

private:
	double m_pEndPt[2];
	double m_pSpacing[2];
};

