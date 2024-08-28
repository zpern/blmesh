#include "BLDomain.h"
#include "algebra.h"
#include <spdlog/spdlog.h> 
 #include <time.h>

BLDomain::BLDomain(void)
{
	m_pElmCen = NULL;
	m_pFEle = NULL;
	m_pFNod = NULL;

	m_pQ = NULL;
	m_pU = NULL;

	nRglrGauss = 4;

	//m_pPotentialBEM = NULL;
}


BLDomain::~BLDomain(void)
{
	if (m_pElmCen)
	{
		delete []m_pElmCen;
		m_pElmCen = NULL;
	}

	if (m_pFEle)
	{
		delete []m_pFEle;
		m_pFEle = NULL;
	}

	if (m_pFNod)
	{
		delete []m_pFNod;
		m_pFNod = NULL;
	}
}

int BLDomain::getRegularPatches(int *nTriPtch, CTrglPatch *TriPtchAry)
{
	*nTriPtch=1;
	TriPtchAry[0].x1=1.0;
	TriPtchAry[0].y1=0.0;
	TriPtchAry[0].x2=0.0;
	TriPtchAry[0].y2=1.0;
	TriPtchAry[0].x3=0.0;
	TriPtchAry[0].y3=0.0;
	return 0;
}

////////////////////////////////////////////////////////////////////
void BLDomain::allocateMatrices(complex<double>**** INT_RmnU, complex<double>**** INT_RmnQ)
{
	int j, k;
	//allocate memory
/*	localU = new double*[pMultipole->max_nodes];
	for(j=0; j<pMultipole->max_nodes; j++) 
		localU[j] = new double[pMultipole->max_nearNodes];

	localQ = new double*[pMultipole->max_nodes];
	for(j=0; j<pMultipole->max_nodes; j++) 
		localQ[j] = new double[pMultipole->max_nearNodes];*/

	*INT_RmnU = new  complex<double>**[pMultipole->max_nodes];
	*INT_RmnQ = new  complex<double>**[pMultipole->max_nodes];
	for(k=0; k<pMultipole->max_nodes; k++)
	{
		(*INT_RmnU)[k] = new  complex<double>*[p+1];
		(*INT_RmnQ)[k] = new  complex<double>*[p+1];
		for(j=0;j<=p;j++)
		{
			(*INT_RmnU)[k][j] = new  complex<double>[j+1];
			(*INT_RmnQ)[k][j] = new  complex<double>[j+1];
		}
	}
}
////////////////////////////////////////////////////////////////////
void BLDomain::freeMatrices(complex<double>**** INT_RmnU, complex<double>**** INT_RmnQ)
{
	int i, j;
	//freeing memory
/*	for(i=0; i<pMultipole->max_nodes; i++) delete[] localU[i];
	delete[] localU;
	for(i=0; i<pMultipole->max_nodes; i++) delete[] localQ[i];
	delete[] localQ;*/

	for(i=0; i<pMultipole->max_nodes; i++)
	{
		for(j=0;j<=p;j++)
		{
			delete[] (*INT_RmnU)[i][j];
			delete[] (*INT_RmnQ)[i][j];
		}
		delete[] (*INT_RmnU)[i];
		delete[] (*INT_RmnQ)[i];
	}
	delete[] (*INT_RmnU);
	delete[] (*INT_RmnQ);
}

int BLDomain::DomainInit(PotentialBEM3D* pPotenialBEM, int nNods, Node* pNods, int nElms, int *pElms, double *pUv, double *pQv, double* pNorm)
{
	int i, j, idx[DIM3];
	double minx[DIM3], maxx[DIM3], mina, maxa, lenx = 0.0;
	double tmbeg, tmend;

	tmbeg = clock();

	m_nNods = nNods;
	m_nElms = nElms;
	tlNode = nElms;
	tlUnknown = 0;

	m_pFNod = new Node[nNods];
	if(!m_pFNod)
	{
		spdlog::info("Error: Not enough memory!\n");
		return OUT_OF_MEMORY;
	}
	m_pFEle = new int[nElms*DIM3];
	if(!m_pFEle)
	{
		spdlog::info("Error: Not enough memory!\n");
		return OUT_OF_MEMORY;
	}

	m_pElmCen = new Node[nElms];
	if(!m_pElmCen)
	{
		spdlog::info("Error: Not enough memory!\n");
		return OUT_OF_MEMORY;
	}

	m_pNorm = new double[nElms*DIM3];
	if(!m_pNorm)
	{
		spdlog::info("Error: Not enough memory!\n");
		return OUT_OF_MEMORY;
	}
	m_pBC = new double[(DIM3+1)*nElms];
	m_pUv = new double[DIM3*nElms];

	minx[0] = minx[1] = minx[2] = mina = 1.0e20;
	maxx[0] = maxx[1] = maxx[2] = maxa = -1.0e20;
	for (i=0; i<nNods; i++)
	{
		m_pFNod[i] = pNods[i];

		for (j=0; j<DIM3; j++)
		{
			if(m_pFNod[i].coord[j] < minx[j])
				minx[j] = m_pFNod[i].coord[j];

			if(m_pFNod[i].coord[j] > maxx[j])
				maxx[j] = m_pFNod[i].coord[j];
		}
	}

#if 0
	for (j=0; j<DIM3; j++)
	{
		if(minx[j] < mina)
			mina = minx[j];

		if(maxx[j] > maxa)
			maxa = maxx[j];
	}
	for (j=0; j<DIM3; j++)
	{
		minx[j] = mina;

		maxx[j] = maxa;
	}
#endif

	for (i=0; i<DIM3; i++)
	{
		minx[i] -= 1/*0.0001*/;
		maxx[i] += 1/*0.0001*/;

		m_dmCen.coord[i] = (minx[i]+maxx[i])/2;

		if(lenx < (maxx[i]-minx[i]))
			lenx = maxx[i]-minx[i];
	}
	m_dLength = lenx;

	for (i=0; i<nElms; i++)
	{
		idx[0] = m_pFEle[i*DIM3 + 0] = pElms[i*DIM3 + 0]-1;
		idx[1] = m_pFEle[i*DIM3 + 1] = pElms[i*DIM3 + 1]-1;
		idx[2] = m_pFEle[i*DIM3 + 2] = pElms[i*DIM3 + 2]-1;

		for(j=0; j<DIM3; j++)
			m_pElmCen[i].coord[j] = (m_pFNod[idx[0]].coord[j] + m_pFNod[idx[1]].coord[j] + m_pFNod[idx[2]].coord[j])/3;

		m_pNorm[i*DIM3 + 0] = pNorm[i*DIM3 + 0];
		m_pNorm[i*DIM3 + 1] = pNorm[i*DIM3 + 1];
		m_pNorm[i*DIM3 + 2] = pNorm[i*DIM3 + 2];
	}
	
	//***create octatree***//
	pMultipole=new CMultipole;
	if(pMultipole->createMultipole(this))
	{
		spdlog::info("Error: Create multipole failed for domain!\n");
		return BLM_FAILED;
	}
	p=pMultipole->p;//number of truncated terms
	pMultipole->m_pPotentialBEM = pPotenialBEM;
	//allocateMatrices();

	SetPotentialFlux(nElms, pUv, pQv);
	SurfValBEM(pUv, pQv, &m_pBC, &m_pUv);

	//3D only right now
// 	m_pPotentialBEM = new PotentialBEM3D;
// 	m_pPotentialBEM->AllocSurfElmInfo(pMultipole->max_nodes*pMultipole->max_nearNodes);

	tmend = clock();
	printf("finished FMM initialization, wall-clock time used: %fs\n", (tmend - tmbeg)/CLOCKS_PER_SEC);

	return BLM_SUCCESS;
}

void BLDomain::SurfValBEM(double *pUv, double *pQv, double **pBc, double **pU)
{
	int i;
	for (i=0; i<m_nElms; i++)
	{
		(*pBc)[i*(DIM3+1)+0] = 1;
		(*pBc)[i*(DIM3+1)+1] = pUv[i*DIM3+0];
		(*pBc)[i*(DIM3+1)+2] = pUv[i*DIM3+1];
		(*pBc)[i*(DIM3+1)+3] = pUv[i*DIM3+2];

		(*pU)[i*DIM3 + 0] = pQv[i*DIM3 + 0];
		(*pU)[i*DIM3 + 1] = pQv[i*DIM3 + 1];
		(*pU)[i*DIM3 + 2] = pQv[i*DIM3 + 2];
	}
}

int BLDomain::SetPotentialFlux(int nElm, double *u, double *q)
{
	int i, j;
	if (nElm != m_nElms)
	{
		spdlog::info("Error: inconsistent counts of elements and potentials");
		return BLM_FAILED;
	}
#if 1
	m_pU = new complex<double>[nElm*DIM3];
	m_pQ = new complex<double>[nElm*DIM3];

	for (i=0; i<m_nElms; i++)
	{
		for (j=0; j<DIM3; j++)
		{
			m_pU[i*DIM3+j].real(0.0);
			m_pU[i*DIM3+j].imag(0.0);
			m_pQ[i*DIM3+j].real(0.0);
			m_pQ[i*DIM3+j].imag(0.0);

			m_pU[i*DIM3+j] = u[i*DIM3+j];
			m_pQ[i*DIM3+j] = q[i*DIM3+j];
		}
	}
#else
	m_pU = new complex<double>[nElm];
	m_pQ = new complex<double>[nElm];

	for (i=0; i<m_nElms; i++)
	{
		m_pU[i].real(0.0);
		m_pU[i].imag(0.0);
		m_pQ[i].real(0.0);
		m_pQ[i].imag(0.0);

		m_pU[i] = u[i*DIM3+0];
		m_pQ[i] = q[i*DIM3+0];
	}
#endif
}

int BLDomain::ComputeMoment()
{
	int i, j, k, l;
	int INT_nodeNumRow, INT_nodeNumCol, *NodeSequence = NULL, INT_sequenceNo, iIntElemt;
	D3POINT LeafCenter;
	complex<double> ***INT_RmnU = NULL, ***INT_RmnQ = NULL;
	double tmbeg, tmmnt, tmtrfm;

	allocateMatrices(&INT_RmnU, &INT_RmnQ);
	tmbeg = clock();
	#ifndef DEBUG
#endif

	//#pragma omp parallel for /*private(INT_nodeNumRow, INT_nodeNumCol, NodeSequence, LeafCenter, INT_RmnU, INT_RmnQ, INT_sequenceNo)*/
	for(k=0; k<pMultipole->tlLeaf; k++)
	{
		//allocateMatrices(&INT_RmnU, &INT_RmnQ);
		//prepare initial data for integration
		INT_nodeNumRow=pMultipole->arryLeafPoint[k]->tlNode;
		INT_nodeNumCol=pMultipole->arryLeafPoint[k]->tlLocalNode;
		NodeSequence=pMultipole->arryLeafPoint[k]->sequence;
		LeafCenter.x=pMultipole->arryLeafPoint[k]->cube.center.x;
		LeafCenter.y=pMultipole->arryLeafPoint[k]->cube.center.y;
		LeafCenter.z=pMultipole->arryLeafPoint[k]->cube.center.z;

		//initializing
		for(i=0; i<INT_nodeNumRow; i++)
		{
			for(j=0; j<=p; j++)
			{
				for(l=0; l<=j; l++)
				{
					INT_RmnU[i][j][l]=0.0;
					INT_RmnQ[i][j][l]=0.0;
				}
			}
		}

		//begin integration
		for(i=0; i<INT_nodeNumRow; i++)
		{
			//iIntElemt=m_pFEle[NodeSequence[i]];
			iIntElemt=NodeSequence[i];
			INT_sequenceNo=i;

			if(getMatriceRow(INT_RmnU, INT_RmnQ, LeafCenter,INT_sequenceNo, iIntElemt) != BLM_SUCCESS)
			{
				spdlog::info("Error: getMatrices fail!\n");
				//return BLM_FAILED;
			}
		}

		pMultipole->initALeaf_PotentialAndFlux(k, m_pU, m_pQ, INT_RmnU, INT_RmnQ);
		//freeMatrices(&INT_RmnU, &INT_RmnQ);
	}
	tmmnt = clock();
	printf("finished moment initialization, wall-clock time used: %fs\n", (tmmnt - tmbeg)/CLOCKS_PER_SEC);

	pMultipole->ComputeTransform();
	tmtrfm = clock();

	printf("finished moment transformation, wall-clock time used: %fs\n", (tmtrfm - tmmnt)/CLOCKS_PER_SEC);
	//PrintMmn();
	//PrintLmn();
	freeMatrices(&INT_RmnU, &INT_RmnQ);
	return BLM_SUCCESS;
}

int BLDomain::getMatriceRow(complex<double>*** INT_RmnU, complex<double>*** INT_RmnQ, D3POINT lfcenter, int INT_sequenceNo, int iIntElemt)
{
	int i,j, nTrglPtch;
	CTrglPatch vTrglPtch[1];

	getRegularPatches(&nTrglPtch, vTrglPtch);
	
	for(i=0; i<nTrglPtch; i++)
	{
		if(getRegularINT(&vTrglPtch[i], INT_RmnU, INT_RmnQ, lfcenter, INT_sequenceNo,iIntElemt) != BLM_SUCCESS) 
		{
			spdlog::info("Error: regular integration failed!\n");
			return BLM_FAILED;
		}
	}

	return BLM_SUCCESS;
}

int BLDomain::getRegularINT(CTrglPatch *pPatch, complex<double>*** INT_RmnU, complex<double>*** INT_RmnQ, D3POINT lfcenter, int INT_sequenceNo, int iIntElemt)
{
	int i,j, tlPt=1;
	CEPOINT gsElmtCEPt[7];
	double gsElemtJcb[7], gsIntJcbQ[7], gsIntShpQ[7];
	D3POINT gsIntD3PtQ[7];
	D3NORMAL gsIntNmlQ[7];

	double *vShp;
	vShp=new double [tlPt];

	//获取高斯点及权重(除以2)
	pPatch->getElemtCoordAndJb_Full(gsElmtCEPt, gsElemtJcb);

	if(getIntGeoDataAndShape(nRglrGauss, gsElmtCEPt, gsIntD3PtQ, gsIntNmlQ, gsIntJcbQ, gsIntShpQ, iIntElemt) != BLM_SUCCESS)
	{
		spdlog::info("Error: getRegularINT for Full failed!\n");
		return BLM_FAILED;
	}
	VectorXVector(nRglrGauss, gsIntJcbQ, gsElemtJcb);

	for(i=0;i<nRglrGauss;i++)
	{
		//for(j=0;j<ElemtArray[NodeSequence[SRC_ElemtNo]]->tlSrcPt;j++)
		for(j=0;j<1;j++)
			vShp[j]=gsIntShpQ[j*nRglrGauss+i];
		
		getPatchLocalCoef_Constant(gsIntJcbQ[i], gsIntD3PtQ[i],vShp, gsIntNmlQ[i], INT_RmnU, INT_RmnQ, lfcenter, INT_sequenceNo);
	}

	if(vShp) delete [] vShp; vShp=NULL; 
	return 0;
}

void BLDomain::getPatchLocalCoef_Constant(double w,  D3POINT Dpt,double *vShp, D3NORMAL normal,
	complex<double>*** INT_RmnU, complex<double>*** INT_RmnQ, D3POINT LeafCenter, int INT_sequenceNo)//计算原始的Mmn,对应的项为3.53
{
	int m, n,is;
	complex<double> RSmn0[21][21];	//suppose p<=10

	//pMultipole->getExpansionR(LeafCenter, Dpt);
	pMultipole->getExpansionR_Local(LeafCenter, Dpt, RSmn0);

	for(n=0; n<=p; n++)
	{
		for(m=0; m<=n; m++)
		{
			INT_RmnU[INT_sequenceNo][n][m]+=RSmn0[n][m]*w/**vShp[0]*/;
		}
	}

	//Partial derivates of R0,0 are all zero
	//w*=(-1.0);
	complex<double> dR1, dR2;
	double dR3;

	//for(is=0;is<ElemtArray[NodeSequence[INT_ElemtNo]]->tlSrcPt;is++)
	for(is=0;is<1;is++)
	{	
		dR1.real(normal.x*w*vShp[is]/2.0);
		dR1.imag(normal.y*w*vShp[is]/2.0);
// 		dR1.real(normal.x*w/2.0);
// 		dR1.imag(normal.y*w/2.0);
		dR2.real(-dR1.real());
		dR2.imag(dR1.imag());
		dR3=normal.z*w*vShp[is];
//		dR3=normal.z*w;

		//Partial derivates of R0,0 are all zero
		for(n=1; n<=pMultipole->p; n++)
		{
			if(n > 1)
				INT_RmnQ[INT_sequenceNo][n][0]+=(-1.0)*conj(RSmn0[n-1][1])*dR1
				+RSmn0[n-1][1]*dR2;
			INT_RmnQ[INT_sequenceNo][n][0]+=RSmn0[n-1][0]*dR3;

			for(m=1; m<=n; m++)
			{
				INT_RmnQ[INT_sequenceNo][n][m]+=RSmn0[n-1][m-1]*dR1;
				if(n-1 >= m+1)
					INT_RmnQ[INT_sequenceNo][n][m]+=RSmn0[n-1][m+1]*dR2;
				if(n-1 >= m)
					INT_RmnQ[INT_sequenceNo][n][m]+=RSmn0[n-1][m]*dR3;
			}
		}
	}
}

int BLDomain::getIntGeoDataAndShape(int num, CEPOINT *vCpt, D3POINT *vDpt, D3NORMAL *vNml, double *vJcb, double *vShp, int iIntElemt)
{
	if(getIntGeoData(num, vCpt, vDpt, vNml, vJcb, iIntElemt) != BLM_SUCCESS)
		return BLM_FAILED;

	for(int i=0; i<num; i++)
		vShp[i]=1.0;

	return BLM_SUCCESS;
}

int BLDomain::getIntGeoData(int num, CEPOINT *vCpt, D3POINT *vDpt, D3NORMAL *vNml, double *vJcb, int iIntElemt)
{
	int i;
	D3POINT vpt[DIM3];
	for (i=0; i<DIM3; i++)
	{
		vpt[i].x = m_pFNod[m_pFEle[iIntElemt*DIM3 + i]].coord[0];
		vpt[i].y = m_pFNod[m_pFEle[iIntElemt*DIM3 + i]].coord[1];
		vpt[i].z = m_pFNod[m_pFEle[iIntElemt*DIM3 + i]].coord[2];
	}
	vNml[0]=CrossProduct(vpt[0] - vpt[2], vpt[1] - vpt[2]);
	vJcb[0]=sqrt(vNml[0].x*vNml[0].x+vNml[0].y*vNml[0].y+vNml[0].z*vNml[0].z);
	//ASSERT(vJcb[0] > 1.0e-30);
	vNml[0]=vNml[0]/vJcb[0];
	
	for(i=0; i<num; i++)
	{
		vDpt[i].x=vpt[0].x*vCpt[i].s+vpt[1].x*vCpt[i].t+vpt[2].x*(1.0-vCpt[i].s-vCpt[i].t);
		vDpt[i].y=vpt[0].y*vCpt[i].s+vpt[1].y*vCpt[i].t+vpt[2].y*(1.0-vCpt[i].s-vCpt[i].t);
		vDpt[i].z=vpt[0].z*vCpt[i].s+vpt[1].z*vCpt[i].t+vpt[2].z*(1.0-vCpt[i].s-vCpt[i].t);
		if(i == 0)continue;

		vNml[i]=vNml[0];
		vJcb[i]=vJcb[0];
	}

	return BLM_SUCCESS;
}

void BLDomain::PrintMmn()
{
	int i, j, k, l, m, n;
	FILE *fout = NULL;
	if ((fout =fopen("Mmn.out","w")) == NULL )
	{
		spdlog::info("Error: can not open file %s!\n");
		return;
	}

	for(k=0; k<pMultipole->tlLeaf; k++)
	{
		fprintf(fout, "%d \n", k+1);
		for(n=0; n<=p; n++)
		{
			for(m=0; m<=n; m++)
				fprintf(fout, "%20.16e ", std::real(pMultipole->arryLeafPoint[k]->Mmn[0][n][m]));
		}
		fprintf(fout, "\n");
	}

	fclose(fout);
	fout = NULL;
}

void BLDomain::PrintLmn()
{
	int i, j, k, l, m, n;
	FILE *fout = NULL;
	if ((fout =fopen("Lmn.out","w")) == NULL )
	{
		spdlog::info("Error: can not open file %s!\n");
		return;
	}

	for(k=0; k<pMultipole->tlLeaf; k++)
	{
		fprintf(fout, "%d \n", k+1);
		for(n=0; n<=p; n++)
		{
			for(m=0; m<=n; m++)
				fprintf(fout, "%lf ", std::real(pMultipole->arryLeafPoint[k]->Lmn[0][n][m]));
		}
		fprintf(fout, "\n");
	}

	fclose(fout);
	fout = NULL;
}

double BLDomain::GetResultPt(double *pt, double *val, int ilayer)
{
	//double val;
	pMultipole->GetResultPt(pt, val, ilayer);
	return val[0];
}
