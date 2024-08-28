#ifndef _POTENTIAL3D_
#define _POTENTIAL3D_
#pragma once
#include "PotentialBEM.h"
#include <cmath>
#include "BLVector.h"
#ifndef __CUDACC__

#define PI_FLOAT     3.1415926f
#define PIBY2_FLOAT  1.5707963f


#define sqrt14 sqrt

class PotentialBEM3D :
	public PotentialBEM
{
public:
	PotentialBEM3D(void);
	~PotentialBEM3D(void);

	void initatan2();

public:

	virtual void BEMCodeInfo(){};

	/*
	 * @nElm:  number of boundary elements (in)
	 * @nNod:  number of boundary nodes (in)
	 * @pElm:  connectivity of boundary elements (in)
	 * @pPt:   coordinates of boundary nodes (in)
	 * @pBc:   boundary condition types and potential/flux values of boundary elements (in)
	 * @pNorm: normal of boundary elements (in)
	 * @pU:    potential/flux values of boundary elements (out)
	 */
	virtual void SurfEvalBEM(int nElm, int nNod, int* pElm, double* pPt, 
		double* pBc, double* pNorm, double* pU);

	void SurfEvalBEM(int nElm, double* pBc, double* pU, char* filename);

	/*
	 * @nField: number of field points (in)
	 * @pFpt:   coordinates of boundary nodes (in)
	 * @f:      potential values of field points (out)
	 * @g:      potential gradient values of field points (out)
	 */
	virtual void FieldEvalBEM(int nField, double* pFpt, double* f, double* g, int nElm,
		int nBpt, double* pBpt, double* pBc, int* pElm, double* pNorm, double* pU, bool bVector = false, int iNod = -1);

	void FieldEvalBEM(int nField, double* pFpt, double* fldvl, double* g, int nElm,
		int *elms, double* pBc, double* pU, bool bVector, int iNod=-1);

	//useful routine
	inline double gammai(double p1, double p2, double qet, double cr1, double cr2, double zn, double zd);
	double pt_tria_tst(double w, double u, double v, const double *theta, int *theta_id);
	inline double asym_chi(double p1, double p2, double d);
	void ortho_comp_basis(const double *t1, const double *t2, double *e1, double *e2, double *e3);

	double atan22(double y, double x)
	{
		if (x >= 0)
		{
			if (y >= 0)
			{
				if (x >= y)
					return ATAN2_TABLE_PPY[(int)(SIZE * y / x + 0.5)];
				else
					return ATAN2_TABLE_PPX[(int)(SIZE * x / y + 0.5)];
			}
			else
			{
				if (x >= -y)
					return ATAN2_TABLE_PNY[(int)(EZIS * y / x + 0.5)];
				else
					return ATAN2_TABLE_PNX[(int)(EZIS * x / y + 0.5)];
			}
		}
		else
		{
			if (y >= 0)
			{
				if (-x >= y)
					return ATAN2_TABLE_NPY[(int)(EZIS * y / x + 0.5)];
				else
					return ATAN2_TABLE_NPX[(int)(EZIS * x / y + 0.5)];
			}
			else
			{
				if (x <= y) // (-x >= -y)
					return ATAN2_TABLE_NNY[(int)(SIZE * y / x + 0.5)];
				else
					return ATAN2_TABLE_NNX[(int)(SIZE * x / y + 0.5)];
			}
		}
	}

	//fast_atan2f
	double atan23( double/*double*/ y, double/*double*/ x )
	{
		if ( x == 0.0f )
		{
			if ( y > 0.0f ) return PIBY2_FLOAT;
			if ( y == 0.0f ) return 0.0f;
			return -PIBY2_FLOAT;
		}
		double atan;
		double z = y/x;
		if ( fabsf( z ) < 1.0f )
		{
			atan = z/(1.0f + 0.28f*z*z);
			if ( x < 0.0f )
			{
				if ( y < 0.0f ) return atan - PI_FLOAT;
				return atan + PI_FLOAT;
			}
		}
		else
		{
			atan = PIBY2_FLOAT - z/(z*z + 0.28f);
			if ( y < 0.0f ) return atan - PI_FLOAT;
		}
		return atan;
	}

	double fabs( double f ) {
		int tmp = * ( int * ) &f;
		tmp &= 0x7FFFFFFF;
		return * ( double * ) &tmp;
	}

	virtual void SurfElmInfo(int nElm, double* pBpt, int* pElm);
	void AllocSurfElmInfo(int nElm);
	void InitSurfElmInfo(int nElm, double* pBpt, int* pElm);
	void FreeSurfElmInfo(int nElm);

private:
	//
	int m_nElm;
	double *m_pE1, *m_pE2, *m_pE3;
	double *m_pXeQ;
	double *m_pAQ, *m_pBQ, *m_pCQ;
	double *m_pTheta;
	double *m_pCscf, *m_pSncf;

	//for atan2
	int           SIZE;
	double        STRETCH;
	// Output will swing from -STRETCH to STRETCH (default: Math.PI)
	// Useful to change to 1 if you would normally do "atan2(y, x) / Math.PI"

	// Inverse of SIZE
	int        EZIS;
	double*    ATAN2_TABLE_PPY;
	double*    ATAN2_TABLE_PPX;
	double*    ATAN2_TABLE_PNY;
	double*    ATAN2_TABLE_PNX;
	double*    ATAN2_TABLE_NPY;
	double*    ATAN2_TABLE_NPX;
	double*    ATAN2_TABLE_NNY;
	double*    ATAN2_TABLE_NNX;
};
#else
class PotentialBEM3D :
	public PotentialBEM
{
public:
	PotentialBEM3D(void);
	~PotentialBEM3D(void);

	void initatan2();

public:

	virtual void BEMCodeInfo(){};

	/*
	 * @nElm:  number of boundary elements (in)
	 * @nNod:  number of boundary nodes (in)
	 * @pElm:  connectivity of boundary elements (in)
	 * @pPt:   coordinates of boundary nodes (in)
	 * @pBc:   boundary condition types and potential/flux values of boundary elements (in)
	 * @pNorm: normal of boundary elements (in)
	 * @pU:    potential/flux values of boundary elements (out)
	 */
	virtual void SurfEvalBEM(int nElm, int nNod, int* pElm, double* pPt, 
		double* pBc, double* pNorm, double* pU);

	void SurfEvalBEM(int nElm, double* pBc, double* pU, char* filename);

	/*
	 * @nField: number of field points (in)
	 * @pFpt:   coordinates of boundary nodes (in)
	 * @f:      potential values of field points (out)
	 * @g:      potential gradient values of field points (out)
	 */
	virtual void FieldEvalBEM(int nField, double* pFpt, double* f, double* g, int nElm,
		int nBpt, double* pBpt, double* pBc, int* pElm, double* pNorm, double* pU);


	virtual void SurfElmInfo(int nElm, double* pBpt, int* pElm);

	void ReleaseMemory();

private:
	//
	int m_nElm;
	double *m_pE1, *m_pE2, *m_pE3;
	double *m_pXeQ;
	double *m_pAQ, *m_pBQ, *m_pCQ;
	double *m_pTheta;
	double *m_pCscf, *m_pSncf;

	//for atan2
	int           SIZE;
	double        STRETCH;
	// Output will swing from -STRETCH to STRETCH (default: Math.PI)
	// Useful to change to 1 if you would normally do "atan2(y, x) / Math.PI"

	// Inverse of SIZE
	int        EZIS;
	double*    ATAN2_TABLE_PPY;
	double*    ATAN2_TABLE_PPX;
	double*    ATAN2_TABLE_PNY;
	double*    ATAN2_TABLE_PNX;
	double*    ATAN2_TABLE_NPY;
	double*    ATAN2_TABLE_NPX;
	double*    ATAN2_TABLE_NNY;
	double*    ATAN2_TABLE_NNX;
};

//useful routine
// inline double gammai(double p1, double p2, double qet, double cr1, double cr2, double zn, double zd);
// double pt_tria_tst(double w, double u, double v, const double *theta, int *theta_id);
// inline double asym_chi(double p1, double p2, double d);
// void ortho_comp_basis(const double *t1, const double *t2, double *e1, double *e2, double *e3);
#endif

#endif
