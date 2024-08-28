#pragma once
#include "BLMesh_define.h"
#include "BLMultipole_define.h"
#include "FMM.h"
#include "PotentialBEM3D.h"

class CMultipole;
class BLDomain
{
public:
	BLDomain(void);
	~BLDomain(void);

	int DomainInit(PotentialBEM3D* pPotenialBEM, int nNods, Node* pNods, int nElms, int *pElms, double *pUv, double *pQv, double* pNorm);
	void SurfValBEM(double *pUv, double *pQv, double **pBc, double **pU);
	int ComputeMoment();
	double GetResultPt(double *pt, double *val, int ilyer = -1);
	void PrintMmn();
	void PrintLmn();

private:
	int SetPotentialFlux(int nElm, double *u, double *q);
	void allocateMatrices(complex<double>**** INT_RmnU, complex<double>**** INT_RmnQ);
	void freeMatrices(complex<double>**** INT_RmnU, complex<double>**** INT_RmnQ);
	int getRegularPatches(int *nTriPtch, CTrglPatch *TriPtchAry);
	int getRegularINT(CTrglPatch *pPatch,complex<double>*** INT_RmnU, complex<double>*** INT_RmnQ, D3POINT lfcenter, int INT_sequenceNo, int iIntElemt);
	void getPatchLocalCoef_Constant(double w,  D3POINT Dpt,double *vShp, D3NORMAL normal,complex<double>*** INT_RmnU, complex<double>*** INT_RmnQ, D3POINT LeafCenter, int INT_sequenceNo);
	int getMatriceRow(complex<double>*** INT_RmnU, complex<double>*** INT_RmnQ, D3POINT lfcenter, int INT_sequenceNo, int iIntElemt);
	int getIntGeoDataAndShape(int num, CEPOINT *vCpt, D3POINT *vDpt, D3NORMAL *vNml, double *vJcb, double *vShp, int iIntElemt);
	int getIntGeoData(int num, CEPOINT *vCpt, D3POINT *vDpt, D3NORMAL *vNml, double *vJcb, int iIntElemt);

public:
	Node m_dmCen;
	double m_dLength;
	int tlNode;
	int tlUnknown;

	int m_nElms;
	int m_nNods;
	int* m_pFEle;
	Node *m_pFNod;
	Node *m_pElmCen;

	double *m_pBC;
	double *m_pUv;
	double *m_pNorm;

	//
	complex<double>* m_pU;
	complex<double>* m_pQ;

	//
	CMultipole *pMultipole;
	int p;
	/*
	D3POINT LeafCenter;
	int INT_nodeNumRow;
	int INT_nodeNumCol;
	long* NodeSequence;
	complex<double>*** INT_RmnU;
	complex<double>*** INT_RmnQ;
	int INT_sequenceNo;
	int iIntElemt;*/
	double** localU;
	double** localQ;
	int nRglrGauss;

	//
	//PotentialBEM3D* m_pPotentialBEM;
};

