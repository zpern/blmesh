#pragma once
#include "PotentialBEM.h"

extern "C" {
	extern void  POTENTIAL(int*, int*, int*, double*, double*, double*, double*);
	extern void  DOMAIN_GRADIENT(int*, double*, double*, double*, int *,			//nfield,xfield,f,g,n
		double*, double*, int*, double*, double*, int*);		//y,bc,node,dnorm,u,bvec
}

class PotentialBEM2D :
	public PotentialBEM
{
public:
	PotentialBEM2D(void);
	~PotentialBEM2D(void);

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

	/*
	 * @nField: number of field points (in)
	 * @pFpt:   coordinates of boundary nodes (in)
	 * @f:      potential values of field points (out)
	 * @g:      potential gradient values of field points (out)
	 */
	virtual void FieldEvalBEM(int nField, double* pFpt, double* f, double* g, int nElm,
		int nBpt, double* pBpt, double* pBc, int* pElm, double* pNorm, double* pU, bool bVector = false, int iNod = -1);

	virtual void SurfElmInfo(int nElm, double* pBpt, int* pElm){}

	virtual void ReleaseMemory(){}
};

