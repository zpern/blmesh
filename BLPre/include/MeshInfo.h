#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <set>
#include <map>

//#define NULL 0
#define DIM3 3
#define DIM2 2
#define NEIG_NULL -1
#define MAX_CONN 4

enum MeshType
{
	MESH_2D = 0,
	MESH_3D = 1
};

typedef struct
{
	int nconn;
	int conn[MAX_CONN];
	int neig[MAX_CONN];
	int igeom;
}Elm;

class MeshInfo
{
public:
	MeshInfo(MeshType mtType)
	{
		m_mtType = mtType;
		m_pPntidx = NULL;
		m_pPntElm = NULL;
		m_pElm = NULL;
	}

	~MeshInfo(void);
	void RemoveBarePoint(int nPt, int nElm, int *Elm);
	void Initialize(double* point_array,int nPt, int nElm, int *Elm);
	void CalPntElm();

	void CalElmNeig();
	bool IsElmEdge(int eidx, int idx1, int idx2, int *pidx);
	void GetBdryPt(int *npt, int **pt);
	void GetBdryElm(int *nbdry, int **bdry);
	void GetBdryElm(int *nbdry, int **bdry, int fidx);
	int GetNumPoint();
private:
	MeshType m_mtType;
	int *m_pPntidx;
	int *m_pPntElm;
	Elm *m_pElm;
	double *m_ppoint;
	int m_nElm;
	int m_nPt;
	std::set<int> m_setBdryPt;
	std::multimap<int, int> m_mpaBdryElm;
};

