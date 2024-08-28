#pragma once
class PotentialBEM
{
public:
	PotentialBEM(void);
	virtual ~PotentialBEM(void);

public:
	virtual void BEMCodeInfo() = 0;

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
							double* pBc, double* pNorm, double* pU) = 0;

	/*
	 * @nField: number of field points (in)
	 * @pFpt:   coordinates of boundary nodes (in)
	 * @f:      potential values of field points (out)
	 * @g:      potential gradient values of field points (out)
	 */
	virtual void FieldEvalBEM(int nField, double* pFpt, double* f, double* g, int nElm,
		int nBpt, double* pBpt, double* pBc, int* pElm, double* pNorm, double* pU, bool bVector=false, int iNod = -1) = 0;

	virtual void SurfElmInfo(int nElm, double* pBpt, int* pElm) = 0;

#ifdef __CUDACC__
	virtual void ReleaseMemory() = 0;
#endif
};

