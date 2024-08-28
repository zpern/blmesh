#include <spdlog/spdlog.h> 
 #include "PotentialBEM2D.h"


PotentialBEM2D::PotentialBEM2D(void)
{
}


PotentialBEM2D::~PotentialBEM2D(void)
{
}

void PotentialBEM2D::SurfEvalBEM(int nElm, int nNod, int* pElm, double* pPt, 
							double* pBc, double* pNorm, double* pU)
{
	//POTENTIAL(&nElm, &nNod, pElm, pPt, pBc, pNorm, pU);
}

void PotentialBEM2D::FieldEvalBEM(int nField, double* pFpt, double* f, double* g, int nElm,
	int nBpt, double* pBpt, double* pBc, int* pElm, double* pNorm, double* pU, bool bVector, int iNod)
{
	int nvec = bVector?1:0;

	//DOMAIN_GRADIENT(&nField, pFpt, f, g, &nElm, pBpt, pBc, pElm, pNorm, pU, &nvec);
}