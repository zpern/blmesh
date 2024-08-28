#ifndef __CUDACC__
#include <spdlog/spdlog.h> 
 #include "PotentialBEM3D.h"
#include <cmath>
#include <stdio.h>
#include <cstring>
#include <cstdlib>
#include "Vector.h"

double pi = 4.0*atan(1.0); /* pi = 3.1415926535897932 */
double pix2 = 2.0*pi;
double fac = 1.0/(4.0*pi); /* fac == Factor of the Green's function */

PotentialBEM3D::PotentialBEM3D(void)
{
	//initatan2();
	m_pAQ = NULL; m_pBQ = NULL; m_pCQ = NULL;
	m_pXeQ = NULL;
	m_pE1 = NULL;
	m_pE2 = NULL;
	m_pE3 = NULL;
	m_pTheta = NULL;
	m_pCscf = NULL;
	m_pSncf = NULL;
}


PotentialBEM3D::~PotentialBEM3D(void)
{
}

void PotentialBEM3D::initatan2()
{
	SIZE                 = /*1024*/4096;
	STRETCH            = pi;
	// Output will swing from -STRETCH to STRETCH (default: Math.PI)
	// Useful to change to 1 if you would normally do "atan2(y, x) / Math.PI"

	// Inverse of SIZE
	EZIS            = -SIZE;
	ATAN2_TABLE_PPY    = new double[SIZE + 1];
	ATAN2_TABLE_PPX    = new double[SIZE + 1];
	ATAN2_TABLE_PNY    = new double[SIZE + 1];
	ATAN2_TABLE_PNX    = new double[SIZE + 1];
	ATAN2_TABLE_NPY    = new double[SIZE + 1];
	ATAN2_TABLE_NPX    = new double[SIZE + 1];
	ATAN2_TABLE_NNY    = new double[SIZE + 1];
	ATAN2_TABLE_NNX    = new double[SIZE + 1];

	for (int i = 0; i <= SIZE; i++)
	{
		double f = (double)i / SIZE;
		ATAN2_TABLE_PPY[i] = (double)(atan(f) * STRETCH / pi);
		ATAN2_TABLE_PPX[i] = STRETCH * 0.5f - ATAN2_TABLE_PPY[i];
		ATAN2_TABLE_PNY[i] = -ATAN2_TABLE_PPY[i];
		ATAN2_TABLE_PNX[i] = ATAN2_TABLE_PPY[i] - STRETCH * 0.5f;
		ATAN2_TABLE_NPY[i] = STRETCH - ATAN2_TABLE_PPY[i];
		ATAN2_TABLE_NPX[i] = ATAN2_TABLE_PPY[i] + STRETCH * 0.5f;
		ATAN2_TABLE_NNY[i] = ATAN2_TABLE_PPY[i] - STRETCH;
		ATAN2_TABLE_NNX[i] = -STRETCH * 0.5f - ATAN2_TABLE_PPY[i];
	}
}

void PotentialBEM3D::SurfEvalBEM(int nElm, int nNod, int* pElm, double* pPt, 
	double* pBc, double* pNorm, double* pU)
{
	int i, tmp;
	FILE *fout = NULL, *fin = NULL;
	char str[128];

	//for testing
	int nfield = 0/*20*/;
	double *fpnt = NULL, *fldv, *g;

	fldv = new double[nfield];
	fpnt = new double[nfield*3];
	g = new double[nfield*3];
#if 1
	fout = fopen("input.dat", "w");
	if (!fout)
	{
		spdlog::info("Error: cann't open file!\n");
		return;
	}

	fprintf(fout, "BLMesh Version 1.0\n");
	fprintf(fout, "\t%d \t%d \t%d\n", nElm, nNod, nfield);
	fprintf(fout, "\t1.0 \t-1.0\n");
	fprintf(fout, "#Nodes:\n");

	for (i=0; i<nNod; i++)
		fprintf(fout, "\t%d \t%e \t%e \t%e\n", i+1, pPt[i*3+0], pPt[i*3+1], pPt[i*3+2]);

	fprintf(fout, "#Elements:\n");

	for (i=0; i<nElm; i++)
	{
		int bct = pBc[2*i+0];
		fprintf(fout, "\t%10d %10d %10d %10d %10d \t%25.15e\n", i+1, pElm[i*3+0], pElm[i*3+1], pElm[i*3+2], bct, pBc[2*i+1]);
	}

	//for testing (box)
	if(nfield > 0)
	{
		fprintf(fout, "#Field Points Inside Domain\n");
		double itv = /*1.0/(nfield-1)*/0.01;
		for (i=0; i<nfield; i++)
		{
			fpnt[i*3+0] = 1.0+itv*(i);
			fpnt[i*3+1] = 0.52;
			fpnt[i*3+2] = 0.52;

			fprintf(fout, "\t%10d \t%e \t%e \t%e\n", i+1, fpnt[i*3+0], fpnt[i*3+1], fpnt[i*3+2]);
		}
	}

	fprintf(fout, "#End of the File\n");
	fclose(fout);
	fout = NULL;
#endif

	system("3D_Potential_CHBIE_FMM_64.exe");

	fin = fopen("output.dat", "r");
	if (!fin)
	{
		spdlog::info("Error: cann't open file!\n");
		return;
	}
	fscanf(fin, "%s %s %s %s %s\n", str, str, str, str, str);
	fscanf(fin, "%s %s %s %s\n", str, str, str, str);
	for (i=0; i<nElm; i++)
	{
		double fv, gv;
		fscanf(fin, "%d %lf %lf\n", &tmp, &fv, &gv);

		double av = pBc[i*2+1];
		if (pBc[i*2+0] == 1)
		{
			pU[i] = gv;
		}
		else
		{
			pU[i] = fv;
		}
	}

	fclose(fin);
	fin = NULL;

	//for testing (box)
	if(nfield > 0 && 0)
	{
		FieldEvalBEM(nfield, fpnt, fldv, g, nElm, nNod, pPt, pBc, pElm, pNorm, pU);
		fout = fopen("filed.out", "w");
		for (i=0; i<nfield; i++)
		{
			BLVector vec(0.0, 0.0, 0.0);
			vec.x = g[3*i+0];
			vec.y = g[3*i+1];
			vec.z = g[3*i+2];

			vec.normalize();

			fprintf(fout, "%d \t%e (%lf %lf %lf)\n", i+1, fldv[i], vec.x, vec.y, vec.z);
		}
		fclose(fout);
		fout = NULL;
	}
}

void PotentialBEM3D::SurfEvalBEM(int nElm, double* pBc, double* pU, char* filename)
{
	FILE *fin = NULL;
	char str[32];
	int i, tmp;

	fin = fopen(filename, "r");
	if (!fin)
	{
		spdlog::info("Error: cann't open file!\n");
		return;
	}
	fscanf(fin, "%s %s %s %s %s\n", str, str, str, str, str);
	fscanf(fin, "%s %s %s %s\n", str, str, str, str);
	for (i=0; i<nElm; i++)
	{
		double fv, gv;
		fscanf(fin, "%d %lf %lf\n", &tmp, &fv, &gv);

		double av = pBc[i*2+1];
		if (pBc[i*2+0] == 1)
		{
			pU[i] = gv;
		}
		else
		{
			pU[i] = fv;
		}
	}

	fclose(fin);
	fin = NULL;
}




void PotentialBEM3D::SurfElmInfo(int nElm, double* pBpt, int* pElm)
{
	int elQ, k, mx_nd_elt = 3, i;
	double aQ,bQ,cQ,bmc, xeQ[3][3], vij[3][3], aQs, theta[3], e1[3], e2[3], e3[3];
	double alpha1, alpha2, alpha3, kQ3, cscf1, sncf1, cscf2, sncf2, cscf3, sncf3;

	//allocate memeory
	m_nElm = nElm;
	m_pAQ = new double[nElm]; m_pBQ = new double[nElm]; m_pCQ = new double[nElm];
	m_pXeQ = new double[nElm*3];
	m_pE1 = new double[nElm*3];
	m_pE2 = new double[nElm*3];
	m_pE3 = new double[nElm*3];
	m_pTheta = new double[nElm*3];
	m_pCscf = new double[nElm*3];

	m_pSncf = new double[nElm*3];

	for (elQ=0; elQ<nElm; elQ++)
	{
		for (k=0; k<mx_nd_elt; k++) /* Loop over all vertices of element elQ */
			for (i=0; i<3; i++) /* Loop over all components */
				xeQ[k][i] = pBpt[(pElm[elQ*3+k] - 1)*3 + i];
		//xeQ[k][i] = bnd[belt[elQ][k]][i]; /* xeQ[k][0:2] == Cartesian coordinates of the k-th vertex on the triangle elQ */
		/* belt[elQ][k] == Global index of the k-th node on the boundary element elQ   */

		for (k=0; k<3; k++) /* Loop over all components */
		{
			vij[0][k] = xeQ[1][k] - xeQ[0][k]; /* y2-y1 <== Side 1 (Edge1) of element elQ */
			vij[1][k] = xeQ[2][k] - xeQ[0][k]; /* y3-y1 <== Side 3 (Edge3) of element elQ */
		}

		ortho_comp_basis(vij[0], vij[1],e1,e2,e3); /* Deternime the orthonormal companion basis associated with element elQ */

		bQ = vij[0][0]*e1[0] + vij[0][1]*e1[1] + vij[0][2]*e1[2]; /* bQ>0.0 always */
		aQ = vij[1][0]*e2[0] + vij[1][1]*e2[1] + vij[1][2]*e2[2]; /* aQ>0.0 always */
		cQ = vij[1][0]*e1[0] + vij[1][1]*e1[1] + vij[1][2]*e1[2];

		bmc = bQ-cQ; aQs = aQ*aQ;

		theta[0] = acos(cQ/sqrt14(cQ*cQ + aQs)); theta[1] = acos(bmc/sqrt14(bmc*bmc + aQs)); theta[2] = pi - (theta[0]+theta[1]);

		alpha1 = 0.0;	//added by xzf
		alpha2 = pi-theta[1]; alpha3 = pi+theta[0]; /* alpha1 = 0.0 */

		/*       kQ2 = -bmc/aQ;  kQ2 == Slope of Side 2 of element elQ */
		kQ3 = cQ/aQ; /* kQ3 == Slope of Side 3 of element elQ */

		cscf1 = cos(alpha1); sncf1 = sin(alpha1); cscf2 = cos(alpha2); sncf2 = sin(alpha2); cscf3 = cos(alpha3); sncf3 = sin(alpha3);

		m_pAQ[elQ] = aQ; m_pBQ[elQ] = bQ; m_pCQ[elQ] = cQ;
		for (k=0; k<3; k++)
		{
			m_pXeQ[elQ*3+k] = xeQ[0][k];
			m_pE1[elQ*3+k] = e1[k];
			m_pE2[elQ*3+k] = e2[k];
			m_pE3[elQ*3+k] = e3[k];

			m_pTheta[elQ*3+k] = theta[k];
		}
		m_pCscf[elQ*3+0] = cscf1; m_pSncf[elQ*3+0] = sncf1;
		m_pCscf[elQ*3+1] = cscf2; m_pSncf[elQ*3+1] = sncf2;
		m_pCscf[elQ*3+2] = cscf3; m_pSncf[elQ*3+2] = sncf3;
	}
}

void PotentialBEM3D::AllocSurfElmInfo(int nElm)
{
	//allocate memeory
	m_pAQ = new double[nElm]; m_pBQ = new double[nElm]; m_pCQ = new double[nElm];
	m_pXeQ = new double[nElm*3];
	m_pE1 = new double[nElm*3];
	m_pE2 = new double[nElm*3];
	m_pE3 = new double[nElm*3];
	m_pTheta = new double[nElm*3];
	m_pCscf = new double[nElm*3];
	m_pSncf = new double[nElm*3];
}

void PotentialBEM3D::FreeSurfElmInfo(int nElm)
{
	//free memeory
	if(m_pAQ) { delete []m_pAQ; m_pAQ = NULL;}
	if(m_pBQ) { delete []m_pBQ; m_pBQ = NULL;}
	if(m_pCQ) { delete []m_pCQ; m_pCQ = NULL;}

	if(m_pXeQ) { delete []m_pXeQ; m_pXeQ = NULL;}
	if(m_pE1) { delete []m_pE1; m_pE1 = NULL;}
	if(m_pE2) { delete []m_pE2; m_pE2 = NULL;}
	if(m_pE3) { delete []m_pE3; m_pE3 = NULL;}
	if(m_pTheta) { delete []m_pTheta; m_pTheta = NULL;}
	if(m_pCscf) { delete []m_pCscf; m_pCscf = NULL;}
	if(m_pSncf) { delete []m_pSncf; m_pSncf = NULL;}
}

void PotentialBEM3D::InitSurfElmInfo(int nElm, double* pBpt, int* pElm)
{
	int elQ, k, mx_nd_elt = 3, i;
	double aQ,bQ,cQ,bmc, xeQ[3][3], vij[3][3], aQs, theta[3], e1[3], e2[3], e3[3];
	double alpha1, alpha2, alpha3, kQ3, cscf1, sncf1, cscf2, sncf2, cscf3, sncf3;

	m_nElm = nElm;

	for (elQ=0; elQ<nElm; elQ++)
	{
		for (k=0; k<mx_nd_elt; k++) /* Loop over all vertices of element elQ */
			for (i=0; i<3; i++) /* Loop over all components */
				xeQ[k][i] = pBpt[(pElm[elQ*3+k] - 1)*3 + i];
		//xeQ[k][i] = bnd[belt[elQ][k]][i]; /* xeQ[k][0:2] == Cartesian coordinates of the k-th vertex on the triangle elQ */
		/* belt[elQ][k] == Global index of the k-th node on the boundary element elQ   */

		for (k=0; k<3; k++) /* Loop over all components */
		{
			vij[0][k] = xeQ[1][k] - xeQ[0][k]; /* y2-y1 <== Side 1 (Edge1) of element elQ */
			vij[1][k] = xeQ[2][k] - xeQ[0][k]; /* y3-y1 <== Side 3 (Edge3) of element elQ */
		}

		ortho_comp_basis(vij[0], vij[1],e1,e2,e3); /* Deternime the orthonormal companion basis associated with element elQ */

		bQ = vij[0][0]*e1[0] + vij[0][1]*e1[1] + vij[0][2]*e1[2]; /* bQ>0.0 always */
		aQ = vij[1][0]*e2[0] + vij[1][1]*e2[1] + vij[1][2]*e2[2]; /* aQ>0.0 always */
		cQ = vij[1][0]*e1[0] + vij[1][1]*e1[1] + vij[1][2]*e1[2];

		bmc = bQ-cQ; aQs = aQ*aQ;

		theta[0] = acos(cQ/sqrt14(cQ*cQ + aQs)); theta[1] = acos(bmc/sqrt14(bmc*bmc + aQs)); theta[2] = pi - (theta[0]+theta[1]);

		alpha1 = 0.0;	//added by xzf
		alpha2 = pi-theta[1]; alpha3 = pi+theta[0]; /* alpha1 = 0.0 */

		/*       kQ2 = -bmc/aQ;  kQ2 == Slope of Side 2 of element elQ */
		kQ3 = cQ/aQ; /* kQ3 == Slope of Side 3 of element elQ */

		cscf1 = cos(alpha1); sncf1 = sin(alpha1); cscf2 = cos(alpha2); sncf2 = sin(alpha2); cscf3 = cos(alpha3); sncf3 = sin(alpha3);

		m_pAQ[elQ] = aQ; m_pBQ[elQ] = bQ; m_pCQ[elQ] = cQ;
		for (k=0; k<3; k++)
		{
			m_pXeQ[elQ*3+k] = xeQ[0][k];
			m_pE1[elQ*3+k] = e1[k];
			m_pE2[elQ*3+k] = e2[k];
			m_pE3[elQ*3+k] = e3[k];

			m_pTheta[elQ*3+k] = theta[k];
		}
		m_pCscf[elQ*3+0] = cscf1; m_pSncf[elQ*3+0] = sncf1;
		m_pCscf[elQ*3+1] = cscf2; m_pSncf[elQ*3+1] = sncf2;
		m_pCscf[elQ*3+2] = cscf3; m_pSncf[elQ*3+2] = sncf3;
	}
}


void printinfo(int inod, int elm, double flx, double phi, double flxpa, double thetgam, double omega, double gam, double et, double fldvl)
{
	FILE *fout = NULL;
	char filename[256];
	
	memset(filename, 0, sizeof(filename));
	sprintf(filename, "info_%d.data", inod+1);

	fout = fopen(filename, "a");
	if(fout == NULL)
	{
		spdlog::info("Error: can not open file %s\n", filename);
		return;
	}
	if(elm % 40 == 0)
		fprintf(fout, "%3s %10s %10s %10s %11s %12s %13s %14s %15s\n", "elQ", "omega", "gam", "thetgam", "et", "flx", "phi", "fldvl", "fldvl_add");
	
	fprintf(fout, "%d %13.10lf %13.10lf %13.10lf %13.10lf %13.10lf %13.10lf %13.10lf %13.10lf\n", elm+1, omega, gam, thetgam, et, flx, phi, fldvl, (omega - et*thetgam)*flx + thetgam*phi);

	fclose(fout);
	fout = NULL;
}

void PotentialBEM3D::FieldEvalBEM(int nField, double* pFpt, double* fldvl, double* g, int nElm,
	int nBpt, double* pBpt, double* pBc, int* pElm, double* pNorm, double* pU, bool bVector, int iNod)
{
#if 1
	/* ===purpose=== */

	/**********************************************************************************************************************************/
	/********************************  COLLOCATION BEM for 3D LAPLACE EQUATION: SINGULAR FORMULATION  *********************************/
	/***************************************  3-NODED TRIANGULAR ELEMENT: ANALYTIC INTEGRATION  ***************************************/
	/***********************************************  PIECEWISE CONSTANT APPROXIMATION  ***********************************************/
	/********************************  EXact INTegration with Uniform Distribution on a flat triangle  ********************************/
	/*********************************  Given an array of field points ipts[i][0:2], i=0,1,...,nip-1  *********************************/
	/**********  Determine the potential at all field points ipts[i][0:2] due to boundary distribution of potential and flux  *********/
	/**********************************************************************************************************************************/

	/* Developer: Sylvain Nintcheu Fata */
	/* Version:   V.1.3                 */
	/* Released:  Jun 14, 2012          */


	/* typical triangular element with vertices {y1,y2,y3}

	*           y3
	*	    *
	*	   / \
	*	  /   \
	*	 /     \
	*       /       \
	*      /	 \
	*     /		  \
	*    *-------------*
	*  y1               y2

	* Oriented edge L1 = [y1,y2] <== Edge1
	* Oriented edge L2 = [y2,y3] <== Edge2
	* Oriented edge L3 = [y3,y1] <== Edge3
	*/


	/* Note: fldvl[i] == Potential at the i-th interior point */
	double gm_eps = 1.0E-8/*1.0E-10*/, f_eps = 1.0E-8/*1.0E-10*/;	/* gm_eps == Tolerance on the geometry */
	/* f_eps == Tolerance on the argument of a function */

	const int mx_nd_elt = 3;						/* mx_nd_elt == Maximum number of nodes per element */
	//double *phi = NULL, *flx = NULL;				/* phi == BLVector containing boundary potentials */
	/* flx == BLVector containing boundary fluxes */


	int k,i,th_id,elQ;							/* elQ == Global index of the current boundary element */
	double e1[3],e2[3],e3[3];					/* Orthonormal companion basis associated with element elQ */
	double xeQ[mx_nd_elt][3];					/* xeQ[k][0:2] == Cartesian coordinates of the k-th vertex of the triangle elQ */
	double aQ,bQ,cQ;								/* Geometric parameters of element elQ: bQ = base, aQ = height, */
	/* cQ = relative position of the 3rd node of element elQ        */
	double r1[3],xi,zt,eth,et,ets,q[3],qs[3],rho[3],chi[3],gamma[3],gam;
	double d[3];		//added by xzf
	double vij[2][3],aQs,bmc,theta[3],alpha1, alpha2,alpha3,Rh[3],Rhs[3],ch11,ch12,ch22,ch23,ch33,ch31;
	double kQ3; /* kQ3 == Slope of Edge3 of element elQ */
	double shc[mx_nd_elt],sncf1,sncf2,sncf3,cscf1,cscf2,cscf3,qetx2[3];
	double p11,p12,p22,p23,p33,p31,x3,z3,rh1s,rh2s,rh3s,p11s,p12s,p33s,p11rh1,p31rh1,omega;
	double cr11,cr12,zn1,zd1,cr22,cr23,zn2,zd2,cr31,cr33,zn3,zd3,p12rh2,p22rh2,p23rh3,p33rh3,theta0,ThetGam;

	//added by xzf
	double prho1, prho2, prho3, qd[3], dsin[3], dcos[3], delta, deltas, deltac;
	double omes, omec, i3, i3et, i3e, i3d, i5ets, i5e, i5d;

	/*** Note: Below, the Projection of the collocation point ipts[i][0:2] means the Projection of the collocation point ipts[i][0:2] on
	*** Note: the plane of element elQ in the direction e3. e3 is the unit normal to element elQ */

	for (k=0; k<nField; k++) 
	{
		if(!bVector)
			fldvl[k] = 0.0; /* Initialize the potential fldvl[k] at all interior points */
		else
			fldvl[k*3+0] = fldvl[k*3+1] = fldvl[k*3+2] = 0.0;
		g[k*3+0] = g[k*3+1] = g[k*3+2] = 0.0; /* Initialize the gradients g[k] at all interior points */
	}

	for (elQ=0; elQ<nElm; elQ++) /* Loop over all boundary elements */
	{
#if 0
		for (k=0; k<mx_nd_elt; k++) /* Loop over all vertices of element elQ */
			for (i=0; i<3; i++) /* Loop over all components */
				xeQ[k][i] = pBpt[(pElm[elQ*3+k] - 1)*3 + i];
		//xeQ[k][i] = bnd[belt[elQ][k]][i]; /* xeQ[k][0:2] == Cartesian coordinates of the k-th vertex on the triangle elQ */
		/* belt[elQ][k] == Global index of the k-th node on the boundary element elQ   */

		for (k=0; k<3; k++) /* Loop over all components */
		{
			vij[0][k] = xeQ[1][k] - xeQ[0][k]; /* y2-y1 <== Side 1 (Edge1) of element elQ */
			vij[1][k] = xeQ[2][k] - xeQ[0][k]; /* y3-y1 <== Side 3 (Edge3) of element elQ */
		}

		ortho_comp_basis(vij[0], vij[1],e1,e2,e3); /* Deternime the orthonormal companion basis associated with element elQ */

		bQ = vij[0][0]*e1[0] + vij[0][1]*e1[1] + vij[0][2]*e1[2]; /* bQ>0.0 always */
		aQ = vij[1][0]*e2[0] + vij[1][1]*e2[1] + vij[1][2]*e2[2]; /* aQ>0.0 always */
		cQ = vij[1][0]*e1[0] + vij[1][1]*e1[1] + vij[1][2]*e1[2];

		bmc = bQ-cQ; aQs = aQ*aQ;

		theta[0] = acos(cQ/sqrt(cQ*cQ + aQs)); theta[1] = acos(bmc/sqrt(bmc*bmc + aQs)); theta[2] = pi - (theta[0]+theta[1]);

		alpha1 = 0.0;	//added by xzf
		alpha2 = pi-theta[1]; alpha3 = pi+theta[0]; /* alpha1 = 0.0 */

		/*       kQ2 = -bmc/aQ;  kQ2 == Slope of Side 2 of element elQ */
		kQ3 = cQ/aQ; /* kQ3 == Slope of Side 3 of element elQ */

		cscf1 = cos(alpha1); sncf1 = sin(alpha1); cscf2 = cos(alpha2); sncf2 = sin(alpha2); cscf3 = cos(alpha3); sncf3 = sin(alpha3);
#else 
		aQ = m_pAQ[elQ], bQ = m_pBQ[elQ], cQ = m_pCQ[elQ];
		kQ3 = cQ/aQ;
		for (k=0; k<3; k++)
		{
			xeQ[0][k] = m_pXeQ[elQ*3+k];
			e1[k] = m_pE1[elQ*3+k];
			e2[k] = m_pE2[elQ*3+k];
			e3[k] = m_pE3[elQ*3+k];
		}
		theta[0] = m_pTheta[elQ*3+0], theta[1] = m_pTheta[elQ*3+1], theta[2] = m_pTheta[elQ*3+2];
		cscf1 = m_pCscf[elQ*3+0], sncf1 = m_pSncf[elQ*3+0];
		cscf2 = m_pCscf[elQ*3+1], sncf2 = m_pSncf[elQ*3+1];
		cscf3 = m_pCscf[elQ*3+2], sncf3 = m_pSncf[elQ*3+2];
#endif

		for (i=0; i<nField; i++) /* Loop over all interior points */
		{
			for (k=0; k<3; k++) /* Loop over all components */
				r1[k] = pFpt[i*3+k] - xeQ[0][k]; /* r1 = ipts[i][0:2]-y1 */

			//r1(x,y1),r1*el, r1*e2, r1*e3 (xzf)
			xi = r1[0]*e1[0] + r1[1]*e1[1] + r1[2]*e1[2];
			zt = r1[0]*e2[0] + r1[1]*e2[1] + r1[2]*e2[2];
			eth = r1[0]*e3[0] + r1[1]*e3[1] + r1[2]*e3[2];

			//the coordinates of three vertices of elQ, y1:(p11,q11), y2:(p12,q12), y3(p13,q13) and q11 = q12 (xzf)
			p11 = -xi; p12 = bQ-xi; q[0] = -zt; et = -eth;

			if (fabs(et) < gm_eps) et = 0.0;	 //the point ipt[i] is on the plane defined by elQ

			//here, p13 = x3, q13 = z3 (xzf)
			x3 = cQ+p11; z3 = aQ+q[0];
			p22 = p12*cscf2 + q[0]*sncf2; p23 = x3*cscf2 + z3*sncf2; q[1] = q[0]*cscf2 - p12*sncf2;
			p31 = p11*cscf3 + q[0]*sncf3; p33 = x3*cscf3 + z3*sncf3; q[2] = q[0]*cscf3 - p11*sncf3;

			ets = et*et; p11s = p11*p11; p12s = p12*p12; p33s = p33*p33; qs[0] = q[0]*q[0]; qs[1] = q[1]*q[1]; qs[2] = q[2]*q[2];

			//Rhs[0~2], the square of distance between the original point and three vertices (xzf)
			//rh1s, rh2s, rh3s, the square of distance between the source point and three vertices (xzf)
			Rhs[0] = p11s + qs[0]; Rhs[1] = p12s + qs[0]; Rhs[2] = p33s + qs[2]; /* Rhs(3) = x3*x3 + z3*z3 */
			Rh[0] = sqrt14(Rhs[0]); Rh[1] = sqrt14(Rhs[1]); Rh[2] = sqrt14(Rhs[2]);
			rh1s = Rhs[0] + ets; rh2s = Rhs[1] + ets; rh3s = Rhs[2] + ets;
			rho[0] = sqrt14(rh1s); rho[1] = sqrt14(rh2s); rho[2] = sqrt14(rh3s);

			shc[2] = -q[0]/aQ; shc[1] = (kQ3*q[0] - p11)/bQ; 
			shc[0] = 1.0 - (shc[2]+shc[1]); /* shc[0] = 1.0 + (p11 - kQ2*q[0])/bQ */

			if ( fabs(shc[0]) < gm_eps ) { shc[0] = 0.0; q[1] = 0.0; }

			if ( fabs(shc[1]) < gm_eps ) { shc[1] = 0.0; q[2] = 0.0; }
			else if ( fabs(shc[1]-1.0) < gm_eps ) shc[1] = 1.0;

			if ( fabs(shc[2]) < gm_eps ) { shc[2] = 0.0; q[0] = 0.0; }
			else if ( fabs(shc[2]-1.0) < gm_eps ) shc[2] = 1.0;

			theta0 = pt_tria_tst(shc[0],shc[1],shc[2],theta,&th_id); /* Point in triangle test */

			if ( th_id == 4 ) /* if the Projection of the collocation point ipts[i][0:2] is at vertex y1 */
			{
				Rh[0] = 0.0;
				if (et == 0.0) rho[0] = 0.0; /* ipts[i][0:2] = y1 */
			}
			else if ( th_id == 5 ) /* if the Projection of the collocation point ipts[i][0:2] is at vertex y2 */
			{
				Rh[1] = 0.0;
				if (et == 0.0) rho[1] = 0.0; /* ipts[i][0:2] = y2 */
			}
			else if ( th_id == 6 ) /* if the Projection of the collocation point ipts[i][0:2] is at vertex y3 */
			{
				Rh[2] = 0.0;
				if (et == 0.0) rho[2] = 0.0; /* ipts[i][0:2] = y3 */
			}


			qetx2[0] = 2.0*q[0]*et; qetx2[1] = 2.0*q[1]*et; qetx2[2] = 2.0*q[2]*et;

			p11rh1 = p11*rho[0]; p12rh2 = p12*rho[1];
			p22rh2 = p22*rho[1]; p23rh3 = p23*rho[2];
			p33rh3 = p33*rho[2]; p31rh1 = p31*rho[0];

			//added by xzf
			d[0] = qs[0]+ets; d[1] = qs[1]+ets; d[2] = qs[2]+ets;

			if ( rho[0] == 0.0 ) /* if The collocation point ipts[i][0:2] is at vertex y1, i.e., ipts[i][0:2] = y1 */
			{
				gam = 0.0;
				omega = q[1]*log( (p22+rho[1])/(p23+rho[2]) ); /* chi[1] = log( (p22+rho[1])/(p23+rho[2]) ) */
			}
			else if ( rho[1] == 0.0 ) /* if The collocation point ipts[i][0:2] is at vertex y2, i.e., ipts[i][0:2] = y2 */
			{
				gam = 0.0;
				omega = q[2]*log( (p33+rho[2])/(p31+rho[0]) ); /* chi[2] = log( (p33+rho[2])/(p31+rho[0]) ) */
			}
			else if ( rho[2] == 0.0 ) /* if The collocation point ipts[i][0:2] is at vertex y3, i.e., ipts[i][0:2] = y3 */
			{
				gam = 0.0;
				omega = q[0]*log( (p11+rho[0])/(p12+rho[1]) ); /* chi[0] = log( (p11+rho[0])/(p12+rho[1]) ) */
			}
			else
			{
				if ( Rh[0] == 0.0 ) /* if the Projection of the collocation point ipts[i][0:2] is at vertex y1 */
				{
					cr22 = qs[1]*rh2s-p22*p22*ets; cr23 = qs[1]*rh3s-p23*p23*ets;
					zn2 = qetx2[1]*(p23rh3*cr22 - p22rh2*cr23); zd2 = cr22*cr23 + qetx2[1]*qetx2[1]*p22rh2*p23rh3;
					gam = gammai(p22,p23,qetx2[1],cr22,cr23,zn2,zd2); /* gam = gamma[1]; gamma[0] = -gamma[2] = +pi or -pi */
				}
				else if ( Rh[1] == 0.0 ) /* if the Projection of the collocation point ipts[i][0:2] is at vertex y2 */
				{
					cr33 = qs[2]*rh3s - p33s*ets; cr31 = qs[2]*rh1s - p31*p31*ets;
					zn3 = qetx2[2]*(p31rh1*cr33 - p33rh3*cr31); zd3 = cr31*cr33 + qetx2[2]*qetx2[2]*p31rh1*p33rh3;
					gam = gammai(p33,p31,qetx2[2],cr33,cr31,zn3,zd3); /* gam = gamma[2]; gamma[0] = -gamma[1] = +pi or -pi */
				}
				else if ( Rh[2] == 0.0 ) /* if the Projection of the collocation point ipts[i][0:2] is at vertex y3 */
				{
					cr11 = qs[0]*rh1s - p11s*ets; cr12 = qs[0]*rh2s - p12s*ets;
					zn1 = qetx2[0]*(p12rh2*cr11 - p11rh1*cr12); zd1 = cr11*cr12 + qetx2[0]*qetx2[0]*p11rh1*p12rh2;
					gam = gammai(p11,p12,qetx2[0],cr11,cr12,zn1,zd1); /* gam = gamma[0]; gamma[1] = -gamma[2] = +pi or -pi */
				}
				else if ( q[0] == 0.0 ) /* if the Projection of the collocation point ipts[i][0:2] is on the line parallel to edge [y1,y2] */
				{
					cr22 = qs[1]*rh2s-p22*p22*ets; cr23 = qs[1]*rh3s-p23*p23*ets;
					zn2 = qetx2[1]*(p23rh3*cr22 - p22rh2*cr23); zd2 = cr22*cr23 + qetx2[1]*qetx2[1]*p22rh2*p23rh3;

					cr33 = qs[2]*rh3s - p33s*ets; cr31 = qs[2]*rh1s - p31*p31*ets;
					zn3 = qetx2[2]*(p31rh1*cr33 - p33rh3*cr31); zd3 = cr31*cr33 + qetx2[2]*qetx2[2]*p31rh1*p33rh3;

					gamma[1] = gammai(p22,p23,qetx2[1],cr22,cr23,zn2,zd2);
					gamma[2] = gammai(p33,p31,qetx2[2],cr33,cr31,zn3,zd3);
					gam = gamma[1]+gamma[2];
				}
				else if ( q[1] == 0.0 ) /* if the Projection of the collocation point ipts[i][0:2] is on the line parallel to edge [y2,y3] */
				{
					cr11 = qs[0]*rh1s - p11s*ets; cr12 = qs[0]*rh2s - p12s*ets;
					zn1 = qetx2[0]*(p12rh2*cr11 - p11rh1*cr12); zd1 = cr11*cr12 + qetx2[0]*qetx2[0]*p11rh1*p12rh2;

					cr33 = qs[2]*rh3s - p33s*ets; cr31 = qs[2]*rh1s - p31*p31*ets;
					zn3 = qetx2[2]*(p31rh1*cr33 - p33rh3*cr31); zd3 = cr31*cr33 + qetx2[2]*qetx2[2]*p31rh1*p33rh3;

					gamma[0] = gammai(p11,p12,qetx2[0],cr11,cr12,zn1,zd1);
					gamma[2] = gammai(p33,p31,qetx2[2],cr33,cr31,zn3,zd3);
					gam = gamma[0]+gamma[2];
				}
				else if ( q[2] == 0.0 ) /* if the Projection of the collocation point ipts[i][0:2] is on the line parallel to edge [y3,y1] */
				{
					cr11 = qs[0]*rh1s - p11s*ets; cr12 = qs[0]*rh2s - p12s*ets;
					zn1 = qetx2[0]*(p12rh2*cr11 - p11rh1*cr12); zd1 = cr11*cr12 + qetx2[0]*qetx2[0]*p11rh1*p12rh2;

					cr22 = qs[1]*rh2s-p22*p22*ets; cr23 = qs[1]*rh3s-p23*p23*ets;
					zn2 = qetx2[1]*(p23rh3*cr22 - p22rh2*cr23); zd2 = cr22*cr23 + qetx2[1]*qetx2[1]*p22rh2*p23rh3;

					gamma[0] = gammai(p11,p12,qetx2[0],cr11,cr12,zn1,zd1);
					gamma[1] = gammai(p22,p23,qetx2[1],cr22,cr23,zn2,zd2);
					gam = gamma[0]+gamma[1];
				}
				else
				{
					// alpha - beta = arctan[(x-y)/(1+xy)]   (xzf)
					cr11 = qs[0]*rh1s - p11s*ets; cr12 = qs[0]*rh2s - p12s*ets;
					zn1 = qetx2[0]*(p12rh2*cr11 - p11rh1*cr12); zd1 = cr11*cr12 + qetx2[0]*qetx2[0]*p11rh1*p12rh2;

					cr22 = qs[1]*rh2s-p22*p22*ets; cr23 = qs[1]*rh3s-p23*p23*ets;
					zn2 = qetx2[1]*(p23rh3*cr22 - p22rh2*cr23); zd2 = cr22*cr23 + qetx2[1]*qetx2[1]*p22rh2*p23rh3;

					cr33 = qs[2]*rh3s - p33s*ets; cr31 = qs[2]*rh1s - p31*p31*ets;
					zn3 = qetx2[2]*(p31rh1*cr33 - p33rh3*cr31); zd3 = cr31*cr33 + qetx2[2]*qetx2[2]*p31rh1*p33rh3;

					gamma[0] = gammai(p11,p12,qetx2[0],cr11,cr12,zn1,zd1);
					gamma[1] = gammai(p22,p23,qetx2[1],cr22,cr23,zn2,zd2);
					gamma[2] = gammai(p33,p31,qetx2[2],cr33,cr31,zn3,zd3);
					gam = gamma[0]+gamma[1]+gamma[2];
				}


				ch11 = p11+rho[0]; ch12 = p12+rho[1];
				ch22 = p22+rho[1]; ch23 = p23+rho[2];
				ch33 = p33+rho[2]; ch31 = p31+rho[0];

				if ( (ch11 >= f_eps) && (ch12 >= f_eps) )
					chi[0] = log(ch11/ch12);
				else
					chi[0] = asym_chi(p11,p12,qs[0] + ets);

				if ( (ch22 >= f_eps) && (ch23 >= f_eps) )
					chi[1] = log(ch22/ch23);
				else
					chi[1] = asym_chi(p22,p23,qs[1] + ets);

				if ( (ch33 >= f_eps) && (ch31 >= f_eps) )
					chi[2] = log(ch33/ch31);
				else
					chi[2] = asym_chi(p33,p31,qs[2] + ets);

				omega = q[0]*chi[0] + q[1]*chi[1] + q[2]*chi[2];


				//added by xzf
				prho1 = p11/rho[0]-p12/rho[1];
				prho2 = p22/rho[1]-p23/rho[2];
				prho3 = p33/rho[2]-p31/rho[0];
				
				qd[0] = q[0]/d[0];
				qd[1] = q[1]/d[1];
				qd[2] = q[2]/d[2];
			 
				dsin[0] = sncf1/d[0];
				dsin[1] = sncf2/d[1];
				dsin[2] = sncf3/d[2];
				 
				dcos[0] = cscf1/d[0];
				dcos[1] = cscf2/d[1];
				dcos[2] = cscf3/d[2];
				 
				delta = qd[0]*prho1 + qd[1]*prho2 + qd[2]*prho3;
				deltas = dsin[0]*prho1 + dsin[1]*prho2 + dsin[2]*prho3;
				deltac = dcos[0]*prho1 + dcos[1]*prho2 + dcos[2]*prho3;
				 
				omes = sncf1*chi[0] + sncf2*chi[1] + sncf3*chi[2];
				omec = cscf1*chi[0] + cscf2*chi[1] + cscf3*chi[2];

			} /* if The collocation point ipts[i][0:2] is at vertex y1, i.e., ipts[i][0:2] = y1 */


			if ( et >= 0.0 )
				ThetGam = 0.5*gam + theta0;
			else
				ThetGam = 0.5*gam - theta0;

			double flx, phi;
			if(!bVector)
			{
				if (pBc[elQ*2+0] == 1)
				{
					phi = pBc[elQ*2+1];
					flx = pU[elQ];
				}
				else
				{
					phi = pU[elQ];
					flx = pBc[elQ*2+1];
				}

				fldvl[i] = fldvl[i] + (omega - et*ThetGam)*flx + ThetGam*phi;

				//added by xzf
				i3et = ThetGam;
				if(et == 0.0)
					i3 = 0/*i3et*/;
				else
					i3 = (1.0/et)*ThetGam;

				i3e = omes;
				i3d = (-1.0)*omec;

				i5e = deltas;
				if(et == 0.0)
					i5ets = delta;
				else
					i5ets = delta + (1.0/et)*ThetGam;

				i5d = (-1.0)*deltac;

				//1
				// 	  g[i*3+0] = g[i*3+0] + i3e*flx[elQ] + i5e*phi[elQ];
				// 	  g[i*3+1] = g[i*3+1] + i3d*flx[elQ] + i5d*phi[elQ];
				// 	  g[i*3+2] = g[i*3+2] - i3et*flx[elQ] - i5ets*phi[elQ];

				//2
				// 	  g[i*3+0] = g[i*3+0] + flx*(e1[0]*i3e + e2[0]*i3d - e3[0]*i3et) + phi*(et*(e1[0]*i5e + e2[0]*i5d) - e3[0]*i5ets);
				// 	  g[i*3+1] = g[i*3+1] + flx*(e1[1]*i3e + e2[1]*i3d - e3[1]*i3et) + phi*(et*(e1[1]*i5e + e2[1]*i5d) - e3[1]*i5ets);
				// 	  g[i*3+2] = g[i*3+2] + flx*(e1[2]*i3e + e2[2]*i3d - e3[2]*i3et) + phi*(et*(e1[2]*i5e + e2[2]*i5d) - e3[2]*i5ets);

				//3
				if(et == 0.0)
				{
					g[i*3+0] = g[i*3+0] + flx*(e1[0]*i3e + e2[0]*i3d + e3[0]*i3et);
					g[i*3+1] = g[i*3+1] + flx*(e1[1]*i3e + e2[1]*i3d + e3[1]*i3et);
					g[i*3+2] = g[i*3+2] + flx*(e1[2]*i3e + e2[2]*i3d + e3[2]*i3et);
				}
				else
				{
					// 		  g[i*3+0] = g[i*3+0] + flx*(e1[0]*i3e + e2[0]*i3d + e3[0]*i3et) - phi*(-i3*e3[0] - (et*(e1[0]*i5e + e2[0]*i5d) + e3[0]*i5ets));
					// 		  g[i*3+1] = g[i*3+1] + flx*(e1[1]*i3e + e2[1]*i3d + e3[1]*i3et) - phi*(-i3*e3[1] - (et*(e1[1]*i5e + e2[1]*i5d) + e3[1]*i5ets));
					// 		  g[i*3+2] = g[i*3+2] + flx*(e1[2]*i3e + e2[2]*i3d + e3[2]*i3et) - phi*(-i3*e3[2] - (et*(e1[2]*i5e + e2[2]*i5d) + e3[2]*i5ets));

					g[i*3+0] = g[i*3+0] + flx*(e1[0]*i3e + e2[0]*i3d + e3[0]*i3et) - phi*(i3*e3[0] - (et*(e1[0]*i5e + e2[0]*i5d) + e3[0]*i5ets));
					g[i*3+1] = g[i*3+1] + flx*(e1[1]*i3e + e2[1]*i3d + e3[1]*i3et) - phi*(i3*e3[1] - (et*(e1[1]*i5e + e2[1]*i5d) + e3[1]*i5ets));
					g[i*3+2] = g[i*3+2] + flx*(e1[2]*i3e + e2[2]*i3d + e3[2]*i3et) - phi*(i3*e3[2] - (et*(e1[2]*i5e + e2[2]*i5d) + e3[2]*i5ets));
				}
			}
			else
			{
				for (int ndim=0; ndim<3; ndim++)
				{
					if (pBc[elQ*4+0] == 1)
					{
						phi = pBc[elQ*4+1+ndim];
						flx = pU[elQ*3+ndim];
					}
					else
					{
						phi = pU[elQ*3+ndim];
						flx = pBc[elQ*4+1+ndim];
					}
					fldvl[i*3+ndim] = fldvl[i*3+ndim] + (omega - et*ThetGam)*flx + ThetGam*phi;

#if 0
					//if(ndim == 0 && (iNod == 5286 || iNod == 5287 || iNod == 5512 || iNod == 5514))
					if(ndim == 2 && (iNod == 390 || iNod == 391 /*|| iNod == 5151*/))
						printinfo(iNod, elQ, flx, phi, omega - et*ThetGam, ThetGam, omega, gam, et, fldvl[i*3+ndim]);
					int ttttt = 0;
#endif
				}
			}

		} /* Loop over all interior points */

	} /* Loop over all boundary elements */

	for (k=0; k<nField; k++)
	{
		if(!bVector)
		{
			fldvl[k] *= fac; /* Correct the potential at all interior points by the factor fac */

			g[k*3+0] *= fac;
			g[k*3+1] *= fac;
			g[k*3+2] *= fac;
		}
		else
		{
			fldvl[k*3+0] *= fac; /* Correct the potential at all interior points by the factor fac */
			fldvl[k*3+1] *= fac;
			fldvl[k*3+2] *= fac;
		}
	}
#endif
}

void PotentialBEM3D::FieldEvalBEM(int nField, double* pFpt, double* fldvl, double* g, int nElm,
	int *elms, double* pBc, double* pU, bool bVector, int iNod)
{
#if 1
	/* ===purpose=== */

	/**********************************************************************************************************************************/
	/********************************  COLLOCATION BEM for 3D LAPLACE EQUATION: SINGULAR FORMULATION  *********************************/
	/***************************************  3-NODED TRIANGULAR ELEMENT: ANALYTIC INTEGRATION  ***************************************/
	/***********************************************  PIECEWISE CONSTANT APPROXIMATION  ***********************************************/
	/********************************  EXact INTegration with Uniform Distribution on a flat triangle  ********************************/
	/*********************************  Given an array of field points ipts[i][0:2], i=0,1,...,nip-1  *********************************/
	/**********  Determine the potential at all field points ipts[i][0:2] due to boundary distribution of potential and flux  *********/
	/**********************************************************************************************************************************/

	/* Developer: Sylvain Nintcheu Fata */
	/* Version:   V.1.3                 */
	/* Released:  Jun 14, 2012          */


	/* typical triangular element with vertices {y1,y2,y3}

	*           y3
	*	    *
	*	   / \
	*	  /   \
	*	 /     \
	*       /       \
	*      /	 \
	*     /		  \
	*    *-------------*
	*  y1               y2

	* Oriented edge L1 = [y1,y2] <== Edge1
	* Oriented edge L2 = [y2,y3] <== Edge2
	* Oriented edge L3 = [y3,y1] <== Edge3
	*/


	/* Note: fldvl[i] == Potential at the i-th interior point */
	double gm_eps = 1.0E-8/*1.0E-10*/, f_eps = 1.0E-8/*1.0E-10*/;	/* gm_eps == Tolerance on the geometry */
	/* f_eps == Tolerance on the argument of a function */

	const int mx_nd_elt = 3;						/* mx_nd_elt == Maximum number of nodes per element */
	//double *phi = NULL, *flx = NULL;				/* phi == BLVector containing boundary potentials */
	/* flx == BLVector containing boundary fluxes */


	int k,i,th_id,elQ;							/* elQ == Global index of the current boundary element */
	double e1[3],e2[3],e3[3];					/* Orthonormal companion basis associated with element elQ */
	double xeQ[mx_nd_elt][3];					/* xeQ[k][0:2] == Cartesian coordinates of the k-th vertex of the triangle elQ */
	double aQ,bQ,cQ;								/* Geometric parameters of element elQ: bQ = base, aQ = height, */
	/* cQ = relative position of the 3rd node of element elQ        */
	double r1[3],xi,zt,eth,et,ets,q[3],qs[3],rho[3],chi[3],gamma[3],gam;
	double d[3];		//added by xzf
	double vij[2][3],aQs,bmc,theta[3],alpha1, alpha2,alpha3,Rh[3],Rhs[3],ch11,ch12,ch22,ch23,ch33,ch31;
	double kQ3; /* kQ3 == Slope of Edge3 of element elQ */
	double shc[mx_nd_elt],sncf1,sncf2,sncf3,cscf1,cscf2,cscf3,qetx2[3];
	double p11,p12,p22,p23,p33,p31,x3,z3,rh1s,rh2s,rh3s,p11s,p12s,p33s,p11rh1,p31rh1,omega;
	double cr11,cr12,zn1,zd1,cr22,cr23,zn2,zd2,cr31,cr33,zn3,zd3,p12rh2,p22rh2,p23rh3,p33rh3,theta0,ThetGam;

	//added by xzf
	double prho1, prho2, prho3, qd[3], dsin[3], dcos[3], delta, deltas, deltac;
	double omes, omec, i3, i3et, i3e, i3d, i5ets, i5e, i5d;

	/*** Note: Below, the Projection of the collocation point ipts[i][0:2] means the Projection of the collocation point ipts[i][0:2] on
	*** Note: the plane of element elQ in the direction e3. e3 is the unit normal to element elQ */

	for (k=0; k<nField; k++) 
	{
		if(!bVector)
			fldvl[k] = 0.0; /* Initialize the potential fldvl[k] at all interior points */
		else
			fldvl[k*3+0] = fldvl[k*3+1] = fldvl[k*3+2] = 0.0;
		g[k*3+0] = g[k*3+1] = g[k*3+2] = 0.0; /* Initialize the gradients g[k] at all interior points */
	}

	for (int eidx=0; eidx<nElm; eidx++) /* Loop over all boundary elements */
	{
		elQ = elms[eidx];
#if 0
		for (k=0; k<mx_nd_elt; k++) /* Loop over all vertices of element elQ */
			for (i=0; i<3; i++) /* Loop over all components */
				xeQ[k][i] = pBpt[(pElm[elQ*3+k] - 1)*3 + i];
		//xeQ[k][i] = bnd[belt[elQ][k]][i]; /* xeQ[k][0:2] == Cartesian coordinates of the k-th vertex on the triangle elQ */
		/* belt[elQ][k] == Global index of the k-th node on the boundary element elQ   */

		for (k=0; k<3; k++) /* Loop over all components */
		{
			vij[0][k] = xeQ[1][k] - xeQ[0][k]; /* y2-y1 <== Side 1 (Edge1) of element elQ */
			vij[1][k] = xeQ[2][k] - xeQ[0][k]; /* y3-y1 <== Side 3 (Edge3) of element elQ */
		}

		ortho_comp_basis(vij[0], vij[1],e1,e2,e3); /* Deternime the orthonormal companion basis associated with element elQ */

		bQ = vij[0][0]*e1[0] + vij[0][1]*e1[1] + vij[0][2]*e1[2]; /* bQ>0.0 always */
		aQ = vij[1][0]*e2[0] + vij[1][1]*e2[1] + vij[1][2]*e2[2]; /* aQ>0.0 always */
		cQ = vij[1][0]*e1[0] + vij[1][1]*e1[1] + vij[1][2]*e1[2];

		bmc = bQ-cQ; aQs = aQ*aQ;

		theta[0] = acos(cQ/sqrt(cQ*cQ + aQs)); theta[1] = acos(bmc/sqrt(bmc*bmc + aQs)); theta[2] = pi - (theta[0]+theta[1]);

		alpha1 = 0.0;	//added by xzf
		alpha2 = pi-theta[1]; alpha3 = pi+theta[0]; /* alpha1 = 0.0 */

		/*       kQ2 = -bmc/aQ;  kQ2 == Slope of Side 2 of element elQ */
		kQ3 = cQ/aQ; /* kQ3 == Slope of Side 3 of element elQ */

		cscf1 = cos(alpha1); sncf1 = sin(alpha1); cscf2 = cos(alpha2); sncf2 = sin(alpha2); cscf3 = cos(alpha3); sncf3 = sin(alpha3);
#else 
		aQ = m_pAQ[elQ], bQ = m_pBQ[elQ], cQ = m_pCQ[elQ];
		kQ3 = cQ/aQ;
		for (k=0; k<3; k++)
		{
			xeQ[0][k] = m_pXeQ[elQ*3+k];
			e1[k] = m_pE1[elQ*3+k];
			e2[k] = m_pE2[elQ*3+k];
			e3[k] = m_pE3[elQ*3+k];
		}
		theta[0] = m_pTheta[elQ*3+0], theta[1] = m_pTheta[elQ*3+1], theta[2] = m_pTheta[elQ*3+2];
		cscf1 = m_pCscf[elQ*3+0], sncf1 = m_pSncf[elQ*3+0];
		cscf2 = m_pCscf[elQ*3+1], sncf2 = m_pSncf[elQ*3+1];
		cscf3 = m_pCscf[elQ*3+2], sncf3 = m_pSncf[elQ*3+2];
#endif

		for (i=0; i<nField; i++) /* Loop over all interior points */
		{
			for (k=0; k<3; k++) /* Loop over all components */
				r1[k] = pFpt[i*3+k] - xeQ[0][k]; /* r1 = ipts[i][0:2]-y1 */

			//r1(x,y1),r1*el, r1*e2, r1*e3 (xzf)
			xi = r1[0]*e1[0] + r1[1]*e1[1] + r1[2]*e1[2];
			zt = r1[0]*e2[0] + r1[1]*e2[1] + r1[2]*e2[2];
			eth = r1[0]*e3[0] + r1[1]*e3[1] + r1[2]*e3[2];

			//the coordinates of three vertices of elQ, y1:(p11,q11), y2:(p12,q12), y3(p13,q13) and q11 = q12 (xzf)
			p11 = -xi; p12 = bQ-xi; q[0] = -zt; et = -eth;

			if (fabs(et) < gm_eps) et = 0.0;	 //the point ipt[i] is on the plane defined by elQ

			//here, p13 = x3, q13 = z3 (xzf)
			x3 = cQ+p11; z3 = aQ+q[0];
			p22 = p12*cscf2 + q[0]*sncf2; p23 = x3*cscf2 + z3*sncf2; q[1] = q[0]*cscf2 - p12*sncf2;
			p31 = p11*cscf3 + q[0]*sncf3; p33 = x3*cscf3 + z3*sncf3; q[2] = q[0]*cscf3 - p11*sncf3;

			ets = et*et; p11s = p11*p11; p12s = p12*p12; p33s = p33*p33; qs[0] = q[0]*q[0]; qs[1] = q[1]*q[1]; qs[2] = q[2]*q[2];

			//Rhs[0~2], the square of distance between the original point and three vertices (xzf)
			//rh1s, rh2s, rh3s, the square of distance between the source point and three vertices (xzf)
			Rhs[0] = p11s + qs[0]; Rhs[1] = p12s + qs[0]; Rhs[2] = p33s + qs[2]; /* Rhs(3) = x3*x3 + z3*z3 */
			Rh[0] = sqrt14(Rhs[0]); Rh[1] = sqrt14(Rhs[1]); Rh[2] = sqrt14(Rhs[2]);
			rh1s = Rhs[0] + ets; rh2s = Rhs[1] + ets; rh3s = Rhs[2] + ets;
			rho[0] = sqrt14(rh1s); rho[1] = sqrt14(rh2s); rho[2] = sqrt14(rh3s);

			shc[2] = -q[0]/aQ; shc[1] = (kQ3*q[0] - p11)/bQ; 
			shc[0] = 1.0 - (shc[2]+shc[1]); /* shc[0] = 1.0 + (p11 - kQ2*q[0])/bQ */

			if ( fabs(shc[0]) < gm_eps ) { shc[0] = 0.0; q[1] = 0.0; }

			if ( fabs(shc[1]) < gm_eps ) { shc[1] = 0.0; q[2] = 0.0; }
			else if ( fabs(shc[1]-1.0) < gm_eps ) shc[1] = 1.0;

			if ( fabs(shc[2]) < gm_eps ) { shc[2] = 0.0; q[0] = 0.0; }
			else if ( fabs(shc[2]-1.0) < gm_eps ) shc[2] = 1.0;

			theta0 = pt_tria_tst(shc[0],shc[1],shc[2],theta,&th_id); /* Point in triangle test */

			if ( th_id == 4 ) /* if the Projection of the collocation point ipts[i][0:2] is at vertex y1 */
			{
				Rh[0] = 0.0;
				if (et == 0.0) rho[0] = 0.0; /* ipts[i][0:2] = y1 */
			}
			else if ( th_id == 5 ) /* if the Projection of the collocation point ipts[i][0:2] is at vertex y2 */
			{
				Rh[1] = 0.0;
				if (et == 0.0) rho[1] = 0.0; /* ipts[i][0:2] = y2 */
			}
			else if ( th_id == 6 ) /* if the Projection of the collocation point ipts[i][0:2] is at vertex y3 */
			{
				Rh[2] = 0.0;
				if (et == 0.0) rho[2] = 0.0; /* ipts[i][0:2] = y3 */
			}


			qetx2[0] = 2.0*q[0]*et; qetx2[1] = 2.0*q[1]*et; qetx2[2] = 2.0*q[2]*et;

			p11rh1 = p11*rho[0]; p12rh2 = p12*rho[1];
			p22rh2 = p22*rho[1]; p23rh3 = p23*rho[2];
			p33rh3 = p33*rho[2]; p31rh1 = p31*rho[0];

			//added by xzf
			d[0] = qs[0]+ets; d[1] = qs[1]+ets; d[2] = qs[2]+ets;

			if ( rho[0] == 0.0 ) /* if The collocation point ipts[i][0:2] is at vertex y1, i.e., ipts[i][0:2] = y1 */
			{
				gam = 0.0;
				omega = q[1]*log( (p22+rho[1])/(p23+rho[2]) ); /* chi[1] = log( (p22+rho[1])/(p23+rho[2]) ) */
			}
			else if ( rho[1] == 0.0 ) /* if The collocation point ipts[i][0:2] is at vertex y2, i.e., ipts[i][0:2] = y2 */
			{
				gam = 0.0;
				omega = q[2]*log( (p33+rho[2])/(p31+rho[0]) ); /* chi[2] = log( (p33+rho[2])/(p31+rho[0]) ) */
			}
			else if ( rho[2] == 0.0 ) /* if The collocation point ipts[i][0:2] is at vertex y3, i.e., ipts[i][0:2] = y3 */
			{
				gam = 0.0;
				omega = q[0]*log( (p11+rho[0])/(p12+rho[1]) ); /* chi[0] = log( (p11+rho[0])/(p12+rho[1]) ) */
			}
			else
			{
				if ( Rh[0] == 0.0 ) /* if the Projection of the collocation point ipts[i][0:2] is at vertex y1 */
				{
					cr22 = qs[1]*rh2s-p22*p22*ets; cr23 = qs[1]*rh3s-p23*p23*ets;
					zn2 = qetx2[1]*(p23rh3*cr22 - p22rh2*cr23); zd2 = cr22*cr23 + qetx2[1]*qetx2[1]*p22rh2*p23rh3;
					gam = gammai(p22,p23,qetx2[1],cr22,cr23,zn2,zd2); /* gam = gamma[1]; gamma[0] = -gamma[2] = +pi or -pi */
				}
				else if ( Rh[1] == 0.0 ) /* if the Projection of the collocation point ipts[i][0:2] is at vertex y2 */
				{
					cr33 = qs[2]*rh3s - p33s*ets; cr31 = qs[2]*rh1s - p31*p31*ets;
					zn3 = qetx2[2]*(p31rh1*cr33 - p33rh3*cr31); zd3 = cr31*cr33 + qetx2[2]*qetx2[2]*p31rh1*p33rh3;
					gam = gammai(p33,p31,qetx2[2],cr33,cr31,zn3,zd3); /* gam = gamma[2]; gamma[0] = -gamma[1] = +pi or -pi */
				}
				else if ( Rh[2] == 0.0 ) /* if the Projection of the collocation point ipts[i][0:2] is at vertex y3 */
				{
					cr11 = qs[0]*rh1s - p11s*ets; cr12 = qs[0]*rh2s - p12s*ets;
					zn1 = qetx2[0]*(p12rh2*cr11 - p11rh1*cr12); zd1 = cr11*cr12 + qetx2[0]*qetx2[0]*p11rh1*p12rh2;
					gam = gammai(p11,p12,qetx2[0],cr11,cr12,zn1,zd1); /* gam = gamma[0]; gamma[1] = -gamma[2] = +pi or -pi */
				}
				else if ( q[0] == 0.0 ) /* if the Projection of the collocation point ipts[i][0:2] is on the line parallel to edge [y1,y2] */
				{
					cr22 = qs[1]*rh2s-p22*p22*ets; cr23 = qs[1]*rh3s-p23*p23*ets;
					zn2 = qetx2[1]*(p23rh3*cr22 - p22rh2*cr23); zd2 = cr22*cr23 + qetx2[1]*qetx2[1]*p22rh2*p23rh3;

					cr33 = qs[2]*rh3s - p33s*ets; cr31 = qs[2]*rh1s - p31*p31*ets;
					zn3 = qetx2[2]*(p31rh1*cr33 - p33rh3*cr31); zd3 = cr31*cr33 + qetx2[2]*qetx2[2]*p31rh1*p33rh3;

					gamma[1] = gammai(p22,p23,qetx2[1],cr22,cr23,zn2,zd2);
					gamma[2] = gammai(p33,p31,qetx2[2],cr33,cr31,zn3,zd3);
					gam = gamma[1]+gamma[2];
				}
				else if ( q[1] == 0.0 ) /* if the Projection of the collocation point ipts[i][0:2] is on the line parallel to edge [y2,y3] */
				{
					cr11 = qs[0]*rh1s - p11s*ets; cr12 = qs[0]*rh2s - p12s*ets;
					zn1 = qetx2[0]*(p12rh2*cr11 - p11rh1*cr12); zd1 = cr11*cr12 + qetx2[0]*qetx2[0]*p11rh1*p12rh2;

					cr33 = qs[2]*rh3s - p33s*ets; cr31 = qs[2]*rh1s - p31*p31*ets;
					zn3 = qetx2[2]*(p31rh1*cr33 - p33rh3*cr31); zd3 = cr31*cr33 + qetx2[2]*qetx2[2]*p31rh1*p33rh3;

					gamma[0] = gammai(p11,p12,qetx2[0],cr11,cr12,zn1,zd1);
					gamma[2] = gammai(p33,p31,qetx2[2],cr33,cr31,zn3,zd3);
					gam = gamma[0]+gamma[2];
				}
				else if ( q[2] == 0.0 ) /* if the Projection of the collocation point ipts[i][0:2] is on the line parallel to edge [y3,y1] */
				{
					cr11 = qs[0]*rh1s - p11s*ets; cr12 = qs[0]*rh2s - p12s*ets;
					zn1 = qetx2[0]*(p12rh2*cr11 - p11rh1*cr12); zd1 = cr11*cr12 + qetx2[0]*qetx2[0]*p11rh1*p12rh2;

					cr22 = qs[1]*rh2s-p22*p22*ets; cr23 = qs[1]*rh3s-p23*p23*ets;
					zn2 = qetx2[1]*(p23rh3*cr22 - p22rh2*cr23); zd2 = cr22*cr23 + qetx2[1]*qetx2[1]*p22rh2*p23rh3;

					gamma[0] = gammai(p11,p12,qetx2[0],cr11,cr12,zn1,zd1);
					gamma[1] = gammai(p22,p23,qetx2[1],cr22,cr23,zn2,zd2);
					gam = gamma[0]+gamma[1];
				}
				else
				{
					// alpha - beta = arctan[(x-y)/(1+xy)]   (xzf)
					cr11 = qs[0]*rh1s - p11s*ets; cr12 = qs[0]*rh2s - p12s*ets;
					zn1 = qetx2[0]*(p12rh2*cr11 - p11rh1*cr12); zd1 = cr11*cr12 + qetx2[0]*qetx2[0]*p11rh1*p12rh2;

					cr22 = qs[1]*rh2s-p22*p22*ets; cr23 = qs[1]*rh3s-p23*p23*ets;
					zn2 = qetx2[1]*(p23rh3*cr22 - p22rh2*cr23); zd2 = cr22*cr23 + qetx2[1]*qetx2[1]*p22rh2*p23rh3;

					cr33 = qs[2]*rh3s - p33s*ets; cr31 = qs[2]*rh1s - p31*p31*ets;
					zn3 = qetx2[2]*(p31rh1*cr33 - p33rh3*cr31); zd3 = cr31*cr33 + qetx2[2]*qetx2[2]*p31rh1*p33rh3;

					gamma[0] = gammai(p11,p12,qetx2[0],cr11,cr12,zn1,zd1);
					gamma[1] = gammai(p22,p23,qetx2[1],cr22,cr23,zn2,zd2);
					gamma[2] = gammai(p33,p31,qetx2[2],cr33,cr31,zn3,zd3);
					gam = gamma[0]+gamma[1]+gamma[2];
				}


				ch11 = p11+rho[0]; ch12 = p12+rho[1];
				ch22 = p22+rho[1]; ch23 = p23+rho[2];
				ch33 = p33+rho[2]; ch31 = p31+rho[0];

				if ( (ch11 >= f_eps) && (ch12 >= f_eps) )
					chi[0] = log(ch11/ch12);
				else
					chi[0] = asym_chi(p11,p12,qs[0] + ets);

				if ( (ch22 >= f_eps) && (ch23 >= f_eps) )
					chi[1] = log(ch22/ch23);
				else
					chi[1] = asym_chi(p22,p23,qs[1] + ets);

				if ( (ch33 >= f_eps) && (ch31 >= f_eps) )
					chi[2] = log(ch33/ch31);
				else
					chi[2] = asym_chi(p33,p31,qs[2] + ets);

				omega = q[0]*chi[0] + q[1]*chi[1] + q[2]*chi[2];


				//added by xzf
				prho1 = p11/rho[0]-p12/rho[1];
				prho2 = p22/rho[1]-p23/rho[2];
				prho3 = p33/rho[2]-p31/rho[0];
				
				qd[0] = q[0]/d[0];
				qd[1] = q[1]/d[1];
				qd[2] = q[2]/d[2];
			 
				dsin[0] = sncf1/d[0];
				dsin[1] = sncf2/d[1];
				dsin[2] = sncf3/d[2];
				 
				dcos[0] = cscf1/d[0];
				dcos[1] = cscf2/d[1];
				dcos[2] = cscf3/d[2];
				 
				delta = qd[0]*prho1 + qd[1]*prho2 + qd[2]*prho3;
				deltas = dsin[0]*prho1 + dsin[1]*prho2 + dsin[2]*prho3;
				deltac = dcos[0]*prho1 + dcos[1]*prho2 + dcos[2]*prho3;
				 
				omes = sncf1*chi[0] + sncf2*chi[1] + sncf3*chi[2];
				omec = cscf1*chi[0] + cscf2*chi[1] + cscf3*chi[2];

			} /* if The collocation point ipts[i][0:2] is at vertex y1, i.e., ipts[i][0:2] = y1 */


			if ( et >= 0.0 )
				ThetGam = 0.5*gam + theta0;
			else
				ThetGam = 0.5*gam - theta0;

			double flx, phi;
			if(!bVector)
			{
				if (pBc[elQ*2+0] == 1)
				{
					phi = pBc[elQ*2+1];
					flx = pU[elQ];
				}
				else
				{
					phi = pU[elQ];
					flx = pBc[elQ*2+1];
				}

				fldvl[i] = fldvl[i] + (omega - et*ThetGam)*flx + ThetGam*phi;

				//added by xzf
				i3et = ThetGam;
				if(et == 0.0)
					i3 = 0/*i3et*/;
				else
					i3 = (1.0/et)*ThetGam;

				i3e = omes;
				i3d = (-1.0)*omec;

				i5e = deltas;
				if(et == 0.0)
					i5ets = delta;
				else
					i5ets = delta + (1.0/et)*ThetGam;

				i5d = (-1.0)*deltac;

				//1
				// 	  g[i*3+0] = g[i*3+0] + i3e*flx[elQ] + i5e*phi[elQ];
				// 	  g[i*3+1] = g[i*3+1] + i3d*flx[elQ] + i5d*phi[elQ];
				// 	  g[i*3+2] = g[i*3+2] - i3et*flx[elQ] - i5ets*phi[elQ];

				//2
				// 	  g[i*3+0] = g[i*3+0] + flx*(e1[0]*i3e + e2[0]*i3d - e3[0]*i3et) + phi*(et*(e1[0]*i5e + e2[0]*i5d) - e3[0]*i5ets);
				// 	  g[i*3+1] = g[i*3+1] + flx*(e1[1]*i3e + e2[1]*i3d - e3[1]*i3et) + phi*(et*(e1[1]*i5e + e2[1]*i5d) - e3[1]*i5ets);
				// 	  g[i*3+2] = g[i*3+2] + flx*(e1[2]*i3e + e2[2]*i3d - e3[2]*i3et) + phi*(et*(e1[2]*i5e + e2[2]*i5d) - e3[2]*i5ets);

				//3
				if(et == 0.0)
				{
					g[i*3+0] = g[i*3+0] + flx*(e1[0]*i3e + e2[0]*i3d + e3[0]*i3et);
					g[i*3+1] = g[i*3+1] + flx*(e1[1]*i3e + e2[1]*i3d + e3[1]*i3et);
					g[i*3+2] = g[i*3+2] + flx*(e1[2]*i3e + e2[2]*i3d + e3[2]*i3et);
				}
				else
				{
					// 		  g[i*3+0] = g[i*3+0] + flx*(e1[0]*i3e + e2[0]*i3d + e3[0]*i3et) - phi*(-i3*e3[0] - (et*(e1[0]*i5e + e2[0]*i5d) + e3[0]*i5ets));
					// 		  g[i*3+1] = g[i*3+1] + flx*(e1[1]*i3e + e2[1]*i3d + e3[1]*i3et) - phi*(-i3*e3[1] - (et*(e1[1]*i5e + e2[1]*i5d) + e3[1]*i5ets));
					// 		  g[i*3+2] = g[i*3+2] + flx*(e1[2]*i3e + e2[2]*i3d + e3[2]*i3et) - phi*(-i3*e3[2] - (et*(e1[2]*i5e + e2[2]*i5d) + e3[2]*i5ets));

					g[i*3+0] = g[i*3+0] + flx*(e1[0]*i3e + e2[0]*i3d + e3[0]*i3et) - phi*(i3*e3[0] - (et*(e1[0]*i5e + e2[0]*i5d) + e3[0]*i5ets));
					g[i*3+1] = g[i*3+1] + flx*(e1[1]*i3e + e2[1]*i3d + e3[1]*i3et) - phi*(i3*e3[1] - (et*(e1[1]*i5e + e2[1]*i5d) + e3[1]*i5ets));
					g[i*3+2] = g[i*3+2] + flx*(e1[2]*i3e + e2[2]*i3d + e3[2]*i3et) - phi*(i3*e3[2] - (et*(e1[2]*i5e + e2[2]*i5d) + e3[2]*i5ets));
				}
			}
			else
			{
				for (int ndim=0; ndim<3; ndim++)
				{
					if (pBc[elQ*4+0] == 1)
					{
						phi = pBc[elQ*4+1+ndim];
						flx = pU[elQ*3+ndim];
					}
					else
					{
						phi = pU[elQ*3+ndim];
						flx = pBc[elQ*4+1+ndim];
					}
					fldvl[i*3+ndim] = fldvl[i*3+ndim] + (omega - et*ThetGam)*flx + ThetGam*phi;

#if 0
					//if(ndim == 0 && (iNod == 5286 || iNod == 5287 || iNod == 5512 || iNod == 5514))
					if(ndim == 2 && (iNod == 390 || iNod == 391 /*|| iNod == 5151*/))
						printinfo(iNod, elQ, flx, phi, omega - et*ThetGam, ThetGam, omega, gam, et, fldvl[i*3+ndim]);
					int ttttt = 0;
#endif
				}
			}

		} /* Loop over all interior points */

	} /* Loop over all boundary elements */

	for (k=0; k<nField; k++)
	{
		if(!bVector)
		{
			fldvl[k] *= fac; /* Correct the potential at all interior points by the factor fac */

			g[k*3+0] *= fac;
			g[k*3+1] *= fac;
			g[k*3+2] *= fac;
		}
		else
		{
			fldvl[k*3+0] *= fac; /* Correct the potential at all interior points by the factor fac */
			fldvl[k*3+1] *= fac;
			fldvl[k*3+2] *= fac;
		}
	}
#endif
}


double PotentialBEM3D::gammai(double p1, double p2, double qet, double cr1, double cr2, double zn, double zd)
{

	/* ===purpose=== */

	/**********************************************************************************************************************************/
	/***************************************************  Determine gamma_i  **********************************************************/
	/**********************************************************************************************************************************/

	/* Developer: Sylvain Nintcheu Fata, Oak Ridge National Laboratory */
	/* Version:   V.1.2                                                */
	/* Released:  Apr 14, 2011                                         */


	extern double pix2; /* pix2 = 2*pi, pi = 3.1415926535897932 */



	if ( (p1<0.0) && (p2>0.0) ) /* 1-----x----------2 <==> if Projection of collocation point falls in side k */
	{
		if ( qet<0.0 )
		{
			if ( (cr1>0.0) && (cr2<0.0) && (zn>0.0) && (zd<0.0) )
				return atan2(zn,zd) - pix2;
			else if ( (cr1<0.0) && (cr2>0.0) && (zn>0.0) && (zd<0.0) )
				return atan2(zn,zd) - pix2;
			else if ( (cr1<0.0) && (cr2<0.0) && (zn>0.0) )
				return atan2(zn,zd) - pix2;
			else
				return atan2(zn,zd);
		}
		else if ( qet>0.0 )
		{
			if ( (cr1>0.0) && (cr2<0.0) && (zn<0.0) && (zd<0.0) )
				return atan2(zn,zd) + pix2;
			else if ( (cr1<0.0) && (cr2>0.0) && (zn<0.0) && (zd<0.0) )
				return atan2(zn,zd) + pix2;
			else if ( (cr1<0.0) && (cr2<0.0) && (zn<0.0) )
				return atan2(zn,zd) + pix2;
			else
				return atan2(zn,zd);
		}
		else
			return 0.0;    /* gammai = atan2(zn,zd) = 0.0 */
	}
	else /* x 1----------------2 or 1----------------2 x <==> if Projection of collocation point does not fall in side k */
		return atan2(zn,zd);
}


double PotentialBEM3D::pt_tria_tst(double w, double u, double v, const double *theta, int *theta_id)
{

	/* ===purpose=== */

	/**********************************************************************************************************************************/
	/*************************************************  Point in Triangle Test  *******************************************************/
	/***********************  Given a point with barycentric coordinates (u,v) on a plane containing a triangle  **********************/
	/************  Determine whether the point (u,v) is outside, or inside, or on an edge, or at a vertex of the triangle  ************/
	/**************************  Set the flag theta_id characterizing the position of the query point (u,v)  **************************/
	/**********************************************************************************************************************************/

	/* Developer: Sylvain Nintcheu Fata, Oak Ridge National Laboratory */
	/* Version:   V.1.2                                                */
	/* Released:  Apr 14, 2011                                         */


	/* Triangle with vertices {y1,y2,y3}

	*           y3
	*	    *
	*	   / \
	*	  /   \
	*	 /     \
	*       /       \
	*      /	 \
	*     /		  \
	*    *-------------*
	*  y1               y2

	* Oriented edge L1 = [y1,y2] <== Edge1
	* Oriented edge L2 = [y2,y3] <== Edge2
	* Oriented edge L3 = [y3,y1] <== Edge3
	*/


	/* Note: (u,v) == Barycentric coordinates of the query point on the plane containing the triangle
	* Note: w = 1-(u+v)
	* Note: The length of vector theta is 3, i.e. theta[3] */
	extern double pi,pix2; /* pi = 3.1415926535897932, pix2 = 2*pi */



	*theta_id = -1;

	if ( u<0.0 || v<0.0 || w<0.0 )
	{
		*theta_id = 7; return 0.0; /* Query point is Outside the triangle */
	}
	else if ( u==0.0 )
	{
		if ( v==0.0 )
		{
			*theta_id = 4; return theta[0]; /* Query point is at Vertex1 */
		}
		else if ( v==1.0 )
		{
			*theta_id = 6; return theta[2]; /* Query point is at Vertex3 */
		}
		else if ( v<0.0 || v>1.0)
		{
			*theta_id = 7; return 0.0; /* Query point is Outside the triangle */
		}
		else /* v>0.0 && v<1.0 */
		{
			*theta_id = 3; return  pi; /* Query point is on Edge3 of the triangle */
		}   
	}
	else if ( v==0.0 )
	{
		if ( u==1.0 )
		{
			*theta_id = 5; return theta[1]; /* Query point is at Vertex2 */
		}
		else if ( u<0.0 || u>1.0 )
		{
			*theta_id = 7; return 0.0; /* Query point is Outside the triangle */
		}
		else /* u>0.0 && u<1.0 */
		{
			*theta_id = 1; return pi; /* Query point is on Edge1 of the triangle */
		}
	}
	else if ( w==0.0 )
	{
		if ( u<0.0 || v<0.0 )
		{
			*theta_id = 7; return 0.0; /* Query point is Outside the triangle */
		}
		else /* u>0.0 && v>0.0 */
		{
			*theta_id = 2; return pi; /* Query point is on Edge2 of the triangle */
		}
	}
	else /* u>0.0 && v>0.0 && w>0.0 */
	{
		*theta_id = 0; return pix2; /* Query point is Inside the triangle */
	}

}


double PotentialBEM3D::asym_chi(double p1, double p2, double d)
{

	/* ===purpose=== */

	/**********************************************************************************************************************************/
	/************************************  return ASYMTOTIC VALUE of chi when q^2 + et^2 = d << 1  ************************************/
	/**********************************************************************************************************************************/

	/* Developer: Sylvain Nintcheu Fata, Oak Ridge National Laboratory */
	/* Version:   V.1.1                                                */
	/* Released:  Mar 21, 2011                                         */


	double z1,z2;



	if ( fabs(p1)>0.0 && fabs(p2)>0.0 )
	{
		z1 = (d/p1)/p1; z2 = (d/p2)/p2;
		return log( fabs(p2/p1)*( (1.0-0.25*z1*(1.0-0.5*z1))/(1.0-0.25*z2*(1.0-0.5*z2)) ) );
	}
	else return 0.0;

}


void PotentialBEM3D::ortho_comp_basis(const double *t1, const double *t2, double *e1, double *e2, double *e3)
{

	/* ===purpose=== */

	/**********************************************************************************************************************************/
	/************************************  Given two linearly independent vectors t1 and t2 in 3D  ************************************/
	/** Compute the orthonormal triad {e1,e2,e3} such that e1=t1/norm(t1), norm(e2) = 1 and e2 is orthogonal to e1, and e3 = e1 x e2 **/
	/**********************************************************************************************************************************/

	/* Developer: Sylvain Nintcheu Fata, Oak Ridge National Laboratory */
	/* Version:   V.1.2                                                */
	/* Released:  Apr 14, 2011                                         */


	/* Note: The length of vectors t1, t2, e1, e2, e3 is 3 */
	int i;
	double al,snrm_te,nrm_te,nrm_tx;



	/*********  Orthonormal local basis: e1[3],e2[3],e3[3]  ***************************************************************************/
	snrm_te = t1[0]*t1[0] + t1[1]*t1[1] + t1[2]*t1[2];
	nrm_te = sqrt14(snrm_te);

	al = (t1[0]*t2[0] + t1[1]*t2[1] + t1[2]*t2[2])/snrm_te;

	for (i=0; i<3; i++) /* Loop over all components */
		e2[i] = t2[i] - al*t1[i];

	nrm_tx = sqrt14(e2[0]*e2[0] + e2[1]*e2[1] + e2[2]*e2[2]);

	for (i=0; i<3; i++) /* Loop over all components */
	{   
		e1[i] = t1[i]/nrm_te; e2[i] = e2[i]/nrm_tx;
	}   

	/*** Unit Normal vector to the plane spanned by {e1,e2}: e3 = e1 x e2 */
	e3[0] = e1[1]*e2[2] - e2[1]*e1[2];
	e3[1] = e1[2]*e2[0] - e2[2]*e1[0];
	e3[2] = e1[0]*e2[1] - e2[0]*e1[1];

}

#if 0
double PotentialBEM3D::atan2(double y, double x)
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
#endif




#if 0
/**
* ATAN2
*/

double atan2(double y, double x)
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
#endif
#else
#include "PotentialBEM3D.h"
#include <math.h>
#include <stdio.h>
#include <windows.h>
#include "BLVector.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

// double pi = 4.0*atan(1.0); /* pi = 3.1415926535897932 */
// double pix2 = 2.0*pi;
// double fac = 1.0/(4.0*pi); /* fac == Factor of the Green's function */
#define pi 4.0*atan(1.0)
#define pix2 2.0*pi
#define fac 1.0/(4.0*pi)

__device__ double gammai(double/*double*/ p1, double/*double*/ p2, double/*double*/ qet, double/*double*/ cr1, double/*double*/ cr2, double/*double*/ zn, double/*double*/ zd)
{

	/* ===purpose=== */

	/**********************************************************************************************************************************/
	/***************************************************  Determine gamma_i  **********************************************************/
	/**********************************************************************************************************************************/

	/* Developer: Sylvain Nintcheu Fata, Oak Ridge National Laboratory */
	/* Version:   V.1.2                                                */
	/* Released:  Apr 14, 2011                                         */


	//   extern double pix2; /* pix2 = 2*pi, pi = 3.1415926535897932 */



	if ( (p1<0.0) && (p2>0.0) ) /* 1-----x----------2 <==> if Projection of collocation point falls in side k */
	{
		if ( qet<0.0 )
		{
			if ( (cr1>0.0) && (cr2<0.0) && (zn>0.0) && (zd<0.0) )
				return atan2f(zn,zd) - pix2;
			else if ( (cr1<0.0) && (cr2>0.0) && (zn>0.0) && (zd<0.0) )
				return atan2f(zn,zd) - pix2;
			else if ( (cr1<0.0) && (cr2<0.0) && (zn>0.0) )
				return atan2f(zn,zd) - pix2;
			else
				return atan2f(zn,zd);
		}
		else if ( qet>0.0 )
		{
			if ( (cr1>0.0) && (cr2<0.0) && (zn<0.0) && (zd<0.0) )
				return atan2f(zn,zd) + pix2;
			else if ( (cr1<0.0) && (cr2>0.0) && (zn<0.0) && (zd<0.0) )
				return atan2f(zn,zd) + pix2;
			else if ( (cr1<0.0) && (cr2<0.0) && (zn<0.0) )
				return atan2f(zn,zd) + pix2;
			else
				return atan2f(zn,zd);
		}
		else
			return 0.0;    /* gammai = atan2(zn,zd) = 0.0 */
	}
	else /* x 1----------------2 or 1----------------2 x <==> if Projection of collocation point does not fall in side k */
		return atan2f(zn,zd);
}


__device__ double pt_tria_tst(double/*double*/ w, double/*double*/ u, double/*double*/ v, const double *theta, int *theta_id)
{

	/* ===purpose=== */

	/**********************************************************************************************************************************/
	/*************************************************  Point in Triangle Test  *******************************************************/
	/***********************  Given a point with barycentric coordinates (u,v) on a plane containing a triangle  **********************/
	/************  Determine whether the point (u,v) is outside, or inside, or on an edge, or at a vertex of the triangle  ************/
	/**************************  Set the flag theta_id characterizing the position of the query point (u,v)  **************************/
	/**********************************************************************************************************************************/

	/* Developer: Sylvain Nintcheu Fata, Oak Ridge National Laboratory */
	/* Version:   V.1.2                                                */
	/* Released:  Apr 14, 2011                                         */


	/* Triangle with vertices {y1,y2,y3}

	*           y3
	*	    *
	*	   / \
	*	  /   \
	*	 /     \
	*       /       \
	*      /	 \
	*     /		  \
	*    *-------------*
	*  y1               y2

	* Oriented edge L1 = [y1,y2] <== Edge1
	* Oriented edge L2 = [y2,y3] <== Edge2
	* Oriented edge L3 = [y3,y1] <== Edge3
	*/


	/* Note: (u,v) == Barycentric coordinates of the query point on the plane containing the triangle
	* Note: w = 1-(u+v)
	* Note: The length of vector theta is 3, i.e. theta[3] */
	//   extern double pi,pix2; 
	// double pi = 3.1415926535897932, pix2 = 2*pi;



	*theta_id = -1;

	if ( u<0.0 || v<0.0 || w<0.0 )
	{
		*theta_id = 7; return 0.0; /* Query point is Outside the triangle */
	}
	else if ( u==0.0 )
	{
		if ( v==0.0 )
		{
			*theta_id = 4; return theta[0]; /* Query point is at Vertex1 */
		}
		else if ( v==1.0 )
		{
			*theta_id = 6; return theta[2]; /* Query point is at Vertex3 */
		}
		else if ( v<0.0 || v>1.0)
		{
			*theta_id = 7; return 0.0; /* Query point is Outside the triangle */
		}
		else /* v>0.0 && v<1.0 */
		{
			*theta_id = 3; return  pi; /* Query point is on Edge3 of the triangle */
		}   
	}
	else if ( v==0.0 )
	{
		if ( u==1.0 )
		{
			*theta_id = 5; return theta[1]; /* Query point is at Vertex2 */
		}
		else if ( u<0.0 || u>1.0 )
		{
			*theta_id = 7; return 0.0; /* Query point is Outside the triangle */
		}
		else /* u>0.0 && u<1.0 */
		{
			*theta_id = 1; return pi; /* Query point is on Edge1 of the triangle */
		}
	}
	else if ( w==0.0 )
	{
		if ( u<0.0 || v<0.0 )
		{
			*theta_id = 7; return 0.0; /* Query point is Outside the triangle */
		}
		else /* u>0.0 && v>0.0 */
		{
			*theta_id = 2; return pi; /* Query point is on Edge2 of the triangle */
		}
	}
	else /* u>0.0 && v>0.0 && w>0.0 */
	{
		*theta_id = 0; return pix2; /* Query point is Inside the triangle */
	}

}


__device__ double asym_chi(double/*double*/ p1, double/*double*/ p2, double/*double*/ d)
{

	/* ===purpose=== */

	/**********************************************************************************************************************************/
	/************************************  return ASYMTOTIC VALUE of chi when q^2 + et^2 = d << 1  ************************************/
	/**********************************************************************************************************************************/

	/* Developer: Sylvain Nintcheu Fata, Oak Ridge National Laboratory */
	/* Version:   V.1.1                                                */
	/* Released:  Mar 21, 2011                                         */


	double z1,z2;



	if ( fabs(p1)>0.0 && fabs(p2)>0.0 )
	{
		z1 = (d/p1)/p1; z2 = (d/p2)/p2;
		return log( fabs(p2/p1)*( (1.0-0.25*z1*(1.0-0.5*z1))/(1.0-0.25*z2*(1.0-0.5*z2)) ) );
	}
	else return 0.0;

}


__host__ __device__ void ortho_comp_basis(const double *t1, const double *t2, double *e1, double *e2, double *e3)
{

	/* ===purpose=== */

	/**********************************************************************************************************************************/
	/************************************  Given two linearly independent vectors t1 and t2 in 3D  ************************************/
	/** Compute the orthonormal triad {e1,e2,e3} such that e1=t1/norm(t1), norm(e2) = 1 and e2 is orthogonal to e1, and e3 = e1 x e2 **/
	/**********************************************************************************************************************************/

	/* Developer: Sylvain Nintcheu Fata, Oak Ridge National Laboratory */
	/* Version:   V.1.2                                                */
	/* Released:  Apr 14, 2011                                         */


	/* Note: The length of vectors t1, t2, e1, e2, e3 is 3 */
	int i;
	double al,snrm_te,nrm_te,nrm_tx;



	/*********  Orthonormal local basis: e1[3],e2[3],e3[3]  ***************************************************************************/
	snrm_te = t1[0]*t1[0] + t1[1]*t1[1] + t1[2]*t1[2];
	nrm_te = sqrt(snrm_te);

	al = (t1[0]*t2[0] + t1[1]*t2[1] + t1[2]*t2[2])/snrm_te;

	for (i=0; i<3; i++) /* Loop over all components */
		e2[i] = t2[i] - al*t1[i];

	nrm_tx = sqrt(e2[0]*e2[0] + e2[1]*e2[1] + e2[2]*e2[2]);

	for (i=0; i<3; i++) /* Loop over all components */
	{   
		e1[i] = t1[i]/nrm_te; e2[i] = e2[i]/nrm_tx;
	}   

	/*** Unit Normal vector to the plane spanned by {e1,e2}: e3 = e1 x e2 */
	e3[0] = e1[1]*e2[2] - e2[1]*e1[2];
	e3[1] = e1[2]*e2[0] - e2[2]*e1[0];
	e3[2] = e1[0]*e2[1] - e2[0]*e1[1];

}

__global__ void PotentialBEMGPU(int nField, double* pFpt, double* fldvl, double* g, int nElm,
	double* pBpt, double* pBc, int* pElm, double* pNorm, double* pU/*, double *pAQ, double* pBQ,
																   double *pCQ, double *pXeQ, double *pE1, double* pE2, double *pE3, double *pTheta,
																   double *pCscf, double *pSncf*/)
{
#if 1
	/* ===purpose=== */

	/**********************************************************************************************************************************/
	/********************************  COLLOCATION BEM for 3D LAPLACE EQUATION: SINGULAR FORMULATION  *********************************/
	/***************************************  3-NODED TRIANGULAR ELEMENT: ANALYTIC INTEGRATION  ***************************************/
	/***********************************************  PIECEWISE CONSTANT APPROXIMATION  ***********************************************/
	/********************************  EXact INTegration with Uniform Distribution on a flat triangle  ********************************/
	/*********************************  Given an array of field points ipts[i][0:2], i=0,1,...,nip-1  *********************************/
	/**********  Determine the potential at all field points ipts[i][0:2] due to boundary distribution of potential and flux  *********/
	/**********************************************************************************************************************************/

	/* Developer: Sylvain Nintcheu Fata */
	/* Version:   V.1.3                 */
	/* Released:  Jun 14, 2012          */


	/* typical triangular element with vertices {y1,y2,y3}

	*           y3
	*	    *
	*	   / \
	*	  /   \
	*	 /     \
	*       /       \
	*      /	 \
	*     /		  \
	*    *-------------*
	*  y1               y2

	* Oriented edge L1 = [y1,y2] <== Edge1
	* Oriented edge L2 = [y2,y3] <== Edge2
	* Oriented edge L3 = [y3,y1] <== Edge3
	*/


	/* Note: fldvl[i] == Potential at the i-th interior point */
	double gm_eps = 1.0E-8/*1.0E-10*/, f_eps = 1.0E-8/*1.0E-10*/;	/* gm_eps == Tolerance on the geometry */
	/* f_eps == Tolerance on the argument of a function */

	const int mx_nd_elt = 3;						/* mx_nd_elt == Maximum number of nodes per element */
	//double *phi = NULL, *flx = NULL;				/* phi == BLVector containing boundary potentials */
	/* flx == BLVector containing boundary fluxes */


	int k,i,th_id,elQ;							/* elQ == Global index of the current boundary element */
	register double e1[3],e2[3],e3[3];					/* Orthonormal companion basis associated with element elQ */
	register double xeQ[mx_nd_elt][3];					/* xeQ[k][0:2] == Cartesian coordinates of the k-th vertex of the triangle elQ */
	double aQ,bQ,cQ;								/* Geometric parameters of element elQ: bQ = base, aQ = height, */
	/* cQ = relative position of the 3rd node of element elQ        */
	register double r1[3],xi,zt,eth,et,ets,q[3],qs[3],rho[3],chi[3],gamma[3],gam;
	register double d[3];		//added by xzf
	register double vij[2][3],aQs,bmc,theta[3],alpha1, alpha2,alpha3,Rh[3],Rhs[3],ch11,ch12,ch22,ch23,ch33,ch31;
	double kQ3; /* kQ3 == Slope of Edge3 of element elQ */
	register double shc[mx_nd_elt],sncf1,sncf2,sncf3,cscf1,cscf2,cscf3,qetx2[3];
	double p11,p12,p22,p23,p33,p31,x3,z3,rh1s,rh2s,rh3s,p11s,p12s,p33s,p11rh1,p31rh1,omega;
	double cr11,cr12,zn1,zd1,cr22,cr23,zn2,zd2,cr31,cr33,zn3,zd3,p12rh2,p22rh2,p23rh3,p33rh3,theta0,ThetGam;

	//added by xzf
	register double prho1, prho2, prho3, qd[3], dsin[3], dcos[3], delta, deltas, deltac;
	double omes, omec, i3, i3et, i3e, i3d, i5ets, i5e, i5d;

	/*** Note: Below, the Projection of the collocation point ipts[i][0:2] means the Projection of the collocation point ipts[i][0:2] on
	*** Note: the plane of element elQ in the direction e3. e3 is the unit normal to element elQ */


	elQ = blockDim.x * blockIdx.x + threadIdx.x;

	if(elQ == 0)
		for (k=0; k<nField; k++) 
		{
			fldvl[k] = 0.0; /* Initialize the potential fldvl[k] at all interior points */
			g[k*3+0] = g[k*3+1] = g[k*3+2] = 0.0; /* Initialize the gradients g[k] at all interior points */
		}

		__syncthreads();

		//for (elQ=0; elQ<nElm; elQ++) /* Loop over all boundary elements */
		if(elQ < nElm)
		{
#if 1
			for (k=0; k<mx_nd_elt; k++) /* Loop over all vertices of element elQ */
				for (i=0; i<3; i++) /* Loop over all components */
					xeQ[k][i] = pBpt[(pElm[elQ*3+k] - 1)*3 + i];
			//xeQ[k][i] = bnd[belt[elQ][k]][i]; /* xeQ[k][0:2] == Cartesian coordinates of the k-th vertex on the triangle elQ */
			/* belt[elQ][k] == Global index of the k-th node on the boundary element elQ   */

			for (k=0; k<3; k++) /* Loop over all components */
			{
				vij[0][k] = xeQ[1][k] - xeQ[0][k]; /* y2-y1 <== Side 1 (Edge1) of element elQ */
				vij[1][k] = xeQ[2][k] - xeQ[0][k]; /* y3-y1 <== Side 3 (Edge3) of element elQ */
			}

			ortho_comp_basis(vij[0], vij[1],e1,e2,e3); /* Deternime the orthonormal companion basis associated with element elQ */

			bQ = vij[0][0]*e1[0] + vij[0][1]*e1[1] + vij[0][2]*e1[2]; /* bQ>0.0 always */
			aQ = vij[1][0]*e2[0] + vij[1][1]*e2[1] + vij[1][2]*e2[2]; /* aQ>0.0 always */
			cQ = vij[1][0]*e1[0] + vij[1][1]*e1[1] + vij[1][2]*e1[2];

			bmc = bQ-cQ; aQs = aQ*aQ;

			theta[0] = acosf(cQ/sqrt(cQ*cQ + aQs)); theta[1] = acosf(bmc/sqrt(bmc*bmc + aQs)); theta[2] = pi - (theta[0]+theta[1]);

			alpha1 = 0.0;	//added by xzf
			alpha2 = pi-theta[1]; alpha3 = pi+theta[0]; /* alpha1 = 0.0 */

			/*       kQ2 = -bmc/aQ;  kQ2 == Slope of Side 2 of element elQ */
			kQ3 = cQ/aQ; /* kQ3 == Slope of Side 3 of element elQ */

			cscf1 = cosf(alpha1); sncf1 = sinf(alpha1); cscf2 = cosf(alpha2); sncf2 = sinf(alpha2); cscf3 = cosf(alpha3); sncf3 = sinf(alpha3);
#else 
			aQ = pAQ[elQ], bQ = pBQ[elQ], cQ = pCQ[elQ];
			kQ3 = cQ/aQ;
			for (k=0; k<3; k++)
			{
				xeQ[0][k] = pXeQ[elQ*3+k];
				e1[k] = pE1[elQ*3+k];
				e2[k] = pE2[elQ*3+k];
				e3[k] = pE3[elQ*3+k];
			}
			theta[0] = pTheta[elQ*3+0], theta[1] = pTheta[elQ*3+1], theta[2] = pTheta[elQ*3+2];
			cscf1 = pCscf[elQ*3+0], sncf1 = pSncf[elQ*3+0];
			cscf2 = pCscf[elQ*3+1], sncf2 = pSncf[elQ*3+1];
			cscf3 = pCscf[elQ*3+2], sncf3 = pSncf[elQ*3+2];
#endif

			for (i=0; i<nField; i++) /* Loop over all interior points */
			{
				for (k=0; k<3; k++) /* Loop over all components */
					r1[k] = pFpt[i*3+k] - xeQ[0][k]; /* r1 = ipts[i][0:2]-y1 */

				//r1(x,y1),r1*el, r1*e2, r1*e3 (xzf)
				xi = r1[0]*e1[0] + r1[1]*e1[1] + r1[2]*e1[2];
				zt = r1[0]*e2[0] + r1[1]*e2[1] + r1[2]*e2[2];
				eth = r1[0]*e3[0] + r1[1]*e3[1] + r1[2]*e3[2];

				//the coordinates of three vertices of elQ, y1:(p11,q11), y2:(p12,q12), y3(p13,q13) and q11 = q12 (xzf)
				p11 = -xi; p12 = bQ-xi; q[0] = -zt; et = -eth;

				if (fabs(et) < gm_eps) et = 0.0;	 //the point ipt[i] is on the plane defined by elQ

				//here, p13 = x3, q13 = z3 (xzf)
				x3 = cQ+p11; z3 = aQ+q[0];
				p22 = p12*cscf2 + q[0]*sncf2; p23 = x3*cscf2 + z3*sncf2; q[1] = q[0]*cscf2 - p12*sncf2;
				p31 = p11*cscf3 + q[0]*sncf3; p33 = x3*cscf3 + z3*sncf3; q[2] = q[0]*cscf3 - p11*sncf3;

				ets = et*et; p11s = p11*p11; p12s = p12*p12; p33s = p33*p33; qs[0] = q[0]*q[0]; qs[1] = q[1]*q[1]; qs[2] = q[2]*q[2];

				//Rhs[0~2], the square of distance between the original point and three vertices (xzf)
				//rh1s, rh2s, rh3s, the square of distance between the source point and three vertices (xzf)
				Rhs[0] = p11s + qs[0]; Rhs[1] = p12s + qs[0]; Rhs[2] = p33s + qs[2]; /* Rhs(3) = x3*x3 + z3*z3 */
				Rh[0] = sqrtf(Rhs[0]); Rh[1] = sqrtf(Rhs[1]); Rh[2] = sqrtf(Rhs[2]);
				rh1s = Rhs[0] + ets; rh2s = Rhs[1] + ets; rh3s = Rhs[2] + ets;
				rho[0] = sqrtf(rh1s); rho[1] = sqrtf(rh2s); rho[2] = sqrtf(rh3s);

				shc[2] = -q[0]/aQ; shc[1] = (kQ3*q[0] - p11)/bQ; 
				shc[0] = 1.0 - (shc[2]+shc[1]); /* shc[0] = 1.0 + (p11 - kQ2*q[0])/bQ */

				if ( fabs(shc[0]) < gm_eps ) { shc[0] = 0.0; q[1] = 0.0; }

				if ( fabs(shc[1]) < gm_eps ) { shc[1] = 0.0; q[2] = 0.0; }
				else if ( fabs(shc[1]-1.0) < gm_eps ) shc[1] = 1.0;

				if ( fabs(shc[2]) < gm_eps ) { shc[2] = 0.0; q[0] = 0.0; }
				else if ( fabs(shc[2]-1.0) < gm_eps ) shc[2] = 1.0;

				theta0 = pt_tria_tst(shc[0],shc[1],shc[2],theta,&th_id); /* Point in triangle test */

				if ( th_id == 4 ) /* if the Projection of the collocation point ipts[i][0:2] is at vertex y1 */
				{
					Rh[0] = 0.0;
					if (et == 0.0) rho[0] = 0.0; /* ipts[i][0:2] = y1 */
				}
				else if ( th_id == 5 ) /* if the Projection of the collocation point ipts[i][0:2] is at vertex y2 */
				{
					Rh[1] = 0.0;
					if (et == 0.0) rho[1] = 0.0; /* ipts[i][0:2] = y2 */
				}
				else if ( th_id == 6 ) /* if the Projection of the collocation point ipts[i][0:2] is at vertex y3 */
				{
					Rh[2] = 0.0;
					if (et == 0.0) rho[2] = 0.0; /* ipts[i][0:2] = y3 */
				}

				qetx2[0] = 2.0*q[0]*et; qetx2[1] = 2.0*q[1]*et; qetx2[2] = 2.0*q[2]*et;

				p11rh1 = p11*rho[0]; p12rh2 = p12*rho[1];
				p22rh2 = p22*rho[1]; p23rh3 = p23*rho[2];
				p33rh3 = p33*rho[2]; p31rh1 = p31*rho[0];

				//added by xzf
				d[0] = qs[0]+ets; d[1] = qs[1]+ets; d[2] = qs[2]+ets;

				if ( rho[0] == 0.0 ) /* if The collocation point ipts[i][0:2] is at vertex y1, i.e., ipts[i][0:2] = y1 */
				{
					gam = 0.0;
					omega = q[1]*log( (p22+rho[1])/(p23+rho[2]) ); /* chi[1] = log( (p22+rho[1])/(p23+rho[2]) ) */
				}
				else if ( rho[1] == 0.0 ) /* if The collocation point ipts[i][0:2] is at vertex y2, i.e., ipts[i][0:2] = y2 */
				{
					gam = 0.0;
					omega = q[2]*log( (p33+rho[2])/(p31+rho[0]) ); /* chi[2] = log( (p33+rho[2])/(p31+rho[0]) ) */
				}
				else if ( rho[2] == 0.0 ) /* if The collocation point ipts[i][0:2] is at vertex y3, i.e., ipts[i][0:2] = y3 */
				{
					gam = 0.0;
					omega = q[0]*log( (p11+rho[0])/(p12+rho[1]) ); /* chi[0] = log( (p11+rho[0])/(p12+rho[1]) ) */
				}
				else
				{
					if ( Rh[0] == 0.0 ) /* if the Projection of the collocation point ipts[i][0:2] is at vertex y1 */
					{
						cr22 = qs[1]*rh2s-p22*p22*ets; cr23 = qs[1]*rh3s-p23*p23*ets;
						zn2 = qetx2[1]*(p23rh3*cr22 - p22rh2*cr23); zd2 = cr22*cr23 + qetx2[1]*qetx2[1]*p22rh2*p23rh3;
						gam = gammai(p22,p23,qetx2[1],cr22,cr23,zn2,zd2); /* gam = gamma[1]; gamma[0] = -gamma[2] = +pi or -pi */
					}
					else if ( Rh[1] == 0.0 ) /* if the Projection of the collocation point ipts[i][0:2] is at vertex y2 */
					{
						cr33 = qs[2]*rh3s - p33s*ets; cr31 = qs[2]*rh1s - p31*p31*ets;
						zn3 = qetx2[2]*(p31rh1*cr33 - p33rh3*cr31); zd3 = cr31*cr33 + qetx2[2]*qetx2[2]*p31rh1*p33rh3;
						gam = gammai(p33,p31,qetx2[2],cr33,cr31,zn3,zd3); /* gam = gamma[2]; gamma[0] = -gamma[1] = +pi or -pi */
					}
					else if ( Rh[2] == 0.0 ) /* if the Projection of the collocation point ipts[i][0:2] is at vertex y3 */
					{
						cr11 = qs[0]*rh1s - p11s*ets; cr12 = qs[0]*rh2s - p12s*ets;
						zn1 = qetx2[0]*(p12rh2*cr11 - p11rh1*cr12); zd1 = cr11*cr12 + qetx2[0]*qetx2[0]*p11rh1*p12rh2;
						gam = gammai(p11,p12,qetx2[0],cr11,cr12,zn1,zd1); /* gam = gamma[0]; gamma[1] = -gamma[2] = +pi or -pi */
					}
					else if ( q[0] == 0.0 ) /* if the Projection of the collocation point ipts[i][0:2] is on the line parallel to edge [y1,y2] */
					{
						cr22 = qs[1]*rh2s-p22*p22*ets; cr23 = qs[1]*rh3s-p23*p23*ets;
						zn2 = qetx2[1]*(p23rh3*cr22 - p22rh2*cr23); zd2 = cr22*cr23 + qetx2[1]*qetx2[1]*p22rh2*p23rh3;

						cr33 = qs[2]*rh3s - p33s*ets; cr31 = qs[2]*rh1s - p31*p31*ets;
						zn3 = qetx2[2]*(p31rh1*cr33 - p33rh3*cr31); zd3 = cr31*cr33 + qetx2[2]*qetx2[2]*p31rh1*p33rh3;

						gamma[1] = gammai(p22,p23,qetx2[1],cr22,cr23,zn2,zd2);
						gamma[2] = gammai(p33,p31,qetx2[2],cr33,cr31,zn3,zd3);
						gam = gamma[1]+gamma[2];
					}
					else if ( q[1] == 0.0 ) /* if the Projection of the collocation point ipts[i][0:2] is on the line parallel to edge [y2,y3] */
					{
						cr11 = qs[0]*rh1s - p11s*ets; cr12 = qs[0]*rh2s - p12s*ets;
						zn1 = qetx2[0]*(p12rh2*cr11 - p11rh1*cr12); zd1 = cr11*cr12 + qetx2[0]*qetx2[0]*p11rh1*p12rh2;

						cr33 = qs[2]*rh3s - p33s*ets; cr31 = qs[2]*rh1s - p31*p31*ets;
						zn3 = qetx2[2]*(p31rh1*cr33 - p33rh3*cr31); zd3 = cr31*cr33 + qetx2[2]*qetx2[2]*p31rh1*p33rh3;

						gamma[0] = gammai(p11,p12,qetx2[0],cr11,cr12,zn1,zd1);
						gamma[2] = gammai(p33,p31,qetx2[2],cr33,cr31,zn3,zd3);
						gam = gamma[0]+gamma[2];
					}
					else if ( q[2] == 0.0 ) /* if the Projection of the collocation point ipts[i][0:2] is on the line parallel to edge [y3,y1] */
					{
						cr11 = qs[0]*rh1s - p11s*ets; cr12 = qs[0]*rh2s - p12s*ets;
						zn1 = qetx2[0]*(p12rh2*cr11 - p11rh1*cr12); zd1 = cr11*cr12 + qetx2[0]*qetx2[0]*p11rh1*p12rh2;

						cr22 = qs[1]*rh2s-p22*p22*ets; cr23 = qs[1]*rh3s-p23*p23*ets;
						zn2 = qetx2[1]*(p23rh3*cr22 - p22rh2*cr23); zd2 = cr22*cr23 + qetx2[1]*qetx2[1]*p22rh2*p23rh3;

						gamma[0] = gammai(p11,p12,qetx2[0],cr11,cr12,zn1,zd1);
						gamma[1] = gammai(p22,p23,qetx2[1],cr22,cr23,zn2,zd2);
						gam = gamma[0]+gamma[1];
					}
					else
					{
						// alpha - beta = arctan[(x-y)/(1+xy)]   (xzf)
						cr11 = qs[0]*rh1s - p11s*ets; cr12 = qs[0]*rh2s - p12s*ets;
						zn1 = qetx2[0]*(p12rh2*cr11 - p11rh1*cr12); zd1 = cr11*cr12 + qetx2[0]*qetx2[0]*p11rh1*p12rh2;

						cr22 = qs[1]*rh2s-p22*p22*ets; cr23 = qs[1]*rh3s-p23*p23*ets;
						zn2 = qetx2[1]*(p23rh3*cr22 - p22rh2*cr23); zd2 = cr22*cr23 + qetx2[1]*qetx2[1]*p22rh2*p23rh3;

						cr33 = qs[2]*rh3s - p33s*ets; cr31 = qs[2]*rh1s - p31*p31*ets;
						zn3 = qetx2[2]*(p31rh1*cr33 - p33rh3*cr31); zd3 = cr31*cr33 + qetx2[2]*qetx2[2]*p31rh1*p33rh3;

						gamma[0] = gammai(p11,p12,qetx2[0],cr11,cr12,zn1,zd1);
						gamma[1] = gammai(p22,p23,qetx2[1],cr22,cr23,zn2,zd2);
						gamma[2] = gammai(p33,p31,qetx2[2],cr33,cr31,zn3,zd3);
						gam = gamma[0]+gamma[1]+gamma[2];
					}


					ch11 = p11+rho[0]; ch12 = p12+rho[1];
					ch22 = p22+rho[1]; ch23 = p23+rho[2];
					ch33 = p33+rho[2]; ch31 = p31+rho[0];

					if ( (ch11 >= f_eps) && (ch12 >= f_eps) )
						chi[0] = log(ch11/ch12);
					else
						chi[0] = asym_chi(p11,p12,qs[0] + ets);

					if ( (ch22 >= f_eps) && (ch23 >= f_eps) )
						chi[1] = log(ch22/ch23);
					else
						chi[1] = asym_chi(p22,p23,qs[1] + ets);

					if ( (ch33 >= f_eps) && (ch31 >= f_eps) )
						chi[2] = log(ch33/ch31);
					else
						chi[2] = asym_chi(p33,p31,qs[2] + ets);

					omega = q[0]*chi[0] + q[1]*chi[1] + q[2]*chi[2];

					//added by xzf
					// 		   prho1 = p11/rho[0]-p12/rho[1];
					// 		   prho2 = p22/rho[1]-p23/rho[2];
					// 		   prho3 = p33/rho[2]-p31/rho[0];
					// 
					// 		   qd[0] = q[0]/d[0];
					// 		   qd[1] = q[1]/d[1];
					// 		   qd[2] = q[2]/d[2];
					// 
					// 		   dsin[0] = sncf1/d[0];
					// 		   dsin[1] = sncf2/d[1];
					// 		   dsin[2] = sncf3/d[2];
					// 
					// 		   dcos[0] = cscf1/d[0];
					// 		   dcos[1] = cscf2/d[1];
					// 		   dcos[2] = cscf3/d[2];
					// 
					// 		   delta = qd[0]*prho1 + qd[1]*prho2 + qd[2]*prho3;
					// 		   deltas = dsin[0]*prho1 + dsin[1]*prho2 + dsin[2]*prho3;
					// 		   deltac = dcos[0]*prho1 + dcos[1]*prho2 + dcos[2]*prho3;
					// 
					// 		   omes = sncf1*chi[0] + sncf2*chi[1] + sncf3*chi[2];
					// 		   omec = cscf1*chi[0] + cscf2*chi[1] + cscf3*chi[2];

				} /* if The collocation point ipts[i][0:2] is at vertex y1, i.e., ipts[i][0:2] = y1 */

				if ( et >= 0.0 )
					ThetGam = 0.5*gam + theta0;
				else
					ThetGam = 0.5*gam - theta0;

				double flx, phi;
				if (pBc[elQ*2+0] == 1)
				{
					phi = pBc[elQ*2+1];
					flx = pU[elQ];
				}
				else
				{
					phi = pU[elQ];
					flx = pBc[elQ*2+1];
				}

				double vfldvl = (omega - et*ThetGam)*flx + ThetGam*phi;
				//fldvl[elQ] = vfldvl;

				atomicAdd(&fldvl[i], vfldvl);

				//fldvl[i] = fldvl[i] + (omega - et*ThetGam)*flx + ThetGam*phi;
			} /* Loop over all interior points */

		} /* Loop over all boundary elements */

		//     for (k=0; k<nField; k++)
		// 	{
		// 		fldvl[k] *= fac; /* Correct the potential at all interior points by the factor fac */
		// 
		// 		g[k*3+0] *= fac;
		// 		g[k*3+1] *= fac;
		// 		g[k*3+2] *= fac;
		// 	}
#endif
}

PotentialBEM3D::PotentialBEM3D(void)
{
	initatan2();
}


PotentialBEM3D::~PotentialBEM3D(void)
{
}

void PotentialBEM3D::initatan2()
{
	SIZE                 = /*1024*/4096;
	STRETCH            = pi;
	// Output will swing from -STRETCH to STRETCH (default: Math.PI)
	// Useful to change to 1 if you would normally do "atan2(y, x) / Math.PI"

	// Inverse of SIZE
	EZIS            = -SIZE;
	ATAN2_TABLE_PPY    = new double[SIZE + 1];
	ATAN2_TABLE_PPX    = new double[SIZE + 1];
	ATAN2_TABLE_PNY    = new double[SIZE + 1];
	ATAN2_TABLE_PNX    = new double[SIZE + 1];
	ATAN2_TABLE_NPY    = new double[SIZE + 1];
	ATAN2_TABLE_NPX    = new double[SIZE + 1];
	ATAN2_TABLE_NNY    = new double[SIZE + 1];
	ATAN2_TABLE_NNX    = new double[SIZE + 1];

	for (int i = 0; i <= SIZE; i++)
	{
		double f = (double)i / SIZE;
		ATAN2_TABLE_PPY[i] = (double)(atan(f) * STRETCH / pi);
		ATAN2_TABLE_PPX[i] = STRETCH * 0.5f - ATAN2_TABLE_PPY[i];
		ATAN2_TABLE_PNY[i] = -ATAN2_TABLE_PPY[i];
		ATAN2_TABLE_PNX[i] = ATAN2_TABLE_PPY[i] - STRETCH * 0.5f;
		ATAN2_TABLE_NPY[i] = STRETCH - ATAN2_TABLE_PPY[i];
		ATAN2_TABLE_NPX[i] = ATAN2_TABLE_PPY[i] + STRETCH * 0.5f;
		ATAN2_TABLE_NNY[i] = ATAN2_TABLE_PPY[i] - STRETCH;
		ATAN2_TABLE_NNX[i] = -STRETCH * 0.5f - ATAN2_TABLE_PPY[i];
	}
}

void PotentialBEM3D::SurfEvalBEM(int nElm, int nNod, int* pElm, double* pPt, 
	double* pBc, double* pNorm, double* pU)
{
	int i, tmp;
	FILE *fout = NULL, *fin = NULL;
	char str[128];

	//for testing
	int nfield = 0/*51*/;
	double *fpnt = NULL, *fldv, *g;

	fldv = new double[nfield];
	fpnt = new double[nfield*3];
	g = new double[nfield*3];

	fout = fopen("input.dat", "w");
	if (!fout)
	{
		spdlog::info("Error: cann't open file!\n");
		return;
	}

	fprintf(fout, "BLMesh Version 1.0\n");
	fprintf(fout, "\t%d \t%d \t%d\n", nElm, nNod, nfield);
	fprintf(fout, "\t1.D0 \t0.D0\n");
	fprintf(fout, "$ Nodes:\n");

	for (i=0; i<nNod; i++)
		fprintf(fout, "\t%d \t%e \t%e \t%e\n", i+1, pPt[i*3+0], pPt[i*3+1], pPt[i*3+2]);

	fprintf(fout, "$ Elements and Boundary Conditions:\n");

	for (i=0; i<nElm; i++)
	{
		int bct = pBc[2*i+0];
		fprintf(fout, "\t%10d %10d %10d %10d %10d \t%25.15e\n", i+1, pElm[i*3+0], pElm[i*3+1], pElm[i*3+2], bct, pBc[2*i+1]);
	}

	//for testing (box)
	if(nfield > 0)
	{
		fprintf(fout, "$ Field Points Inside Domain\n");
		double itv = 1.0/(nfield-1);
		for (i=0; i<nfield; i++)
		{
			fpnt[i*3+0] = itv*(i);
			fpnt[i*3+1] = 0.5;
			fpnt[i*3+2] = 0.5;

			fprintf(fout, "\t%10d \t%e \t%e \t%e\n", i+1, fpnt[i*3+0], fpnt[i*3+1], fpnt[i*3+2]);
		}
	}

	fprintf(fout, " $ End of the File\n");

	fclose(fout);
	fout = NULL;

	system("3D_Potential_CHBIE_FMM.exe");

	fin = fopen("output.dat", "r");
	if (!fin)
	{
		spdlog::info("Error: cann't open file!\n");
		return;
	}
	fscanf(fin, "%s %s %s %s %s\n", str, str, str, str, str);
	fscanf(fin, "%s %s %s %s\n", str, str, str, str);
	for (i=0; i<nElm; i++)
	{
		double fv, gv;
		fscanf(fin, "%d %lf %lf\n", &tmp, &fv, &gv);

		double av = pBc[i*2+1];
		if (pBc[i*2+0] == 1)
		{
			pU[i] = gv;
		}
		else
		{
			pU[i] = fv;
		}
	}

	fclose(fin);
	fin = NULL;

	//for testing (box)
	if(nfield > 0)
	{
		FieldEvalBEM(nfield, fpnt, fldv, g, nElm, nNod, pPt, pBc, pElm, pNorm, pU);
		fout = fopen("filed.out", "w");
		for (i=0; i<nfield; i++)
		{
			BLVector vec(0.0, 0.0, 0.0);
			vec.x = g[3*i+0];
			vec.y = g[3*i+1];
			vec.z = g[3*i+2];

			vec.normalize();

			fprintf(fout, "%d \t%e (%lf %lf %lf)\n", i+1, fldv[i], vec.x, vec.y, vec.z);
		}
		fclose(fout);
		fout = NULL;
	}
}

void PotentialBEM3D::SurfEvalBEM(int nElm, double* pBc, double* pU, char* filename)
{
	FILE *fin = NULL;
	char str[32];
	int i, tmp;

	fin = fopen(filename, "r");
	if (!fin)
	{
		spdlog::info("Error: cann't open file!\n");
		return;
	}
	fscanf(fin, "%s %s %s %s %s\n", str, str, str, str, str);
	fscanf(fin, "%s %s %s %s\n", str, str, str, str);
	for (i=0; i<nElm; i++)
	{
		double fv, gv;
		fscanf(fin, "%d %lf %lf\n", &tmp, &fv, &gv);

		double av = pBc[i*2+1];
		if (pBc[i*2+0] == 1)
		{
			pU[i] = gv;
		}
		else
		{
			pU[i] = fv;
		}
	}

	fclose(fin);
	fin = NULL;
}

void PotentialBEM3D::SurfElmInfo(int nElm, double* pBpt, int* pElm)
{
	int elQ, k, mx_nd_elt = 3, i;
	double aQ,bQ,cQ,bmc, xeQ[3][3], vij[3][3], aQs, theta[3], e1[3], e2[3], e3[3];
	double alpha1, alpha2, alpha3, kQ3, cscf1, sncf1, cscf2, sncf2, cscf3, sncf3;

	//allocate memeory
	m_nElm = nElm;
	m_pAQ = new double[nElm]; m_pBQ = new double[nElm]; m_pCQ = new double[nElm];
	m_pXeQ = new double[nElm*3];
	m_pE1 = new double[nElm*3];
	m_pE2 = new double[nElm*3];
	m_pE3 = new double[nElm*3];
	m_pTheta = new double[nElm*3];
	m_pCscf = new double[nElm*3];
	m_pSncf = new double[nElm*3];

	for (elQ=0; elQ<nElm; elQ++)
	{
		for (k=0; k<mx_nd_elt; k++) /* Loop over all vertices of element elQ */
			for (i=0; i<3; i++) /* Loop over all components */
				xeQ[k][i] = pBpt[(pElm[elQ*3+k] - 1)*3 + i];
		//xeQ[k][i] = bnd[belt[elQ][k]][i]; /* xeQ[k][0:2] == Cartesian coordinates of the k-th vertex on the triangle elQ */
		/* belt[elQ][k] == Global index of the k-th node on the boundary element elQ   */

		for (k=0; k<3; k++) /* Loop over all components */
		{
			vij[0][k] = xeQ[1][k] - xeQ[0][k]; /* y2-y1 <== Side 1 (Edge1) of element elQ */
			vij[1][k] = xeQ[2][k] - xeQ[0][k]; /* y3-y1 <== Side 3 (Edge3) of element elQ */
		}

		ortho_comp_basis(vij[0], vij[1],e1,e2,e3); /* Deternime the orthonormal companion basis associated with element elQ */

		bQ = vij[0][0]*e1[0] + vij[0][1]*e1[1] + vij[0][2]*e1[2]; /* bQ>0.0 always */
		aQ = vij[1][0]*e2[0] + vij[1][1]*e2[1] + vij[1][2]*e2[2]; /* aQ>0.0 always */
		cQ = vij[1][0]*e1[0] + vij[1][1]*e1[1] + vij[1][2]*e1[2];

		bmc = bQ-cQ; aQs = aQ*aQ;

		theta[0] = acos(cQ/sqrt(cQ*cQ + aQs)); theta[1] = acos(bmc/sqrt(bmc*bmc + aQs)); theta[2] = pi - (theta[0]+theta[1]);

		alpha1 = 0.0;	//added by xzf
		alpha2 = pi-theta[1]; alpha3 = pi+theta[0]; /* alpha1 = 0.0 */

		/*       kQ2 = -bmc/aQ;  kQ2 == Slope of Side 2 of element elQ */
		kQ3 = cQ/aQ; /* kQ3 == Slope of Side 3 of element elQ */

		cscf1 = cos(alpha1); sncf1 = sin(alpha1); cscf2 = cos(alpha2); sncf2 = sin(alpha2); cscf3 = cos(alpha3); sncf3 = sin(alpha3);

		m_pAQ[elQ] = aQ; m_pBQ[elQ] = bQ; m_pCQ[elQ] = cQ;
		for (k=0; k<3; k++)
		{
			m_pXeQ[elQ*3+k] = xeQ[0][k];
			m_pE1[elQ*3+k] = e1[k];
			m_pE2[elQ*3+k] = e2[k];
			m_pE3[elQ*3+k] = e3[k];

			m_pTheta[elQ*3+k] = theta[k];
		}
		m_pCscf[elQ*3+0] = cscf1; m_pSncf[elQ*3+0] = sncf1;
		m_pCscf[elQ*3+1] = cscf2; m_pSncf[elQ*3+1] = sncf2;
		m_pCscf[elQ*3+2] = cscf3; m_pSncf[elQ*3+2] = sncf3;
	}
}

double *pFpt_d = NULL, *pBpt_d, *pBc_d, *pNorm_d, *pU_d;
double *pAQ_d, *pBQ_d, *pCQ_d, *pXeQ_d, *pE1_d, *pE2_d, *pE3_d, *pTheta_d, *pCscf_d, *pSncf_d;
int *pElm_d;
double *fldvl_d, *g_d, *fval, *gval;

void PotentialBEM3D::FieldEvalBEM(int nField, double* pFpt, double* fldvl, double* g, int nElm,
	int nBpt, double* pBpt, double* pBc, int* pElm, double* pNorm, double* pU)
{
	// 	double *pFpt_d = NULL, *pBpt_d, *pBc_d, *pNorm_d, *pU_d;
	// 	double *pAQ_d, *pBQ_d, *pCQ_d, *pXeQ_d, *pE1_d, *pE2_d, *pE3_d, *pTheta_d, *pCscf_d, *pSncf_d;
	// 	int *pElm_d;
	// 	double *fldvl_d, *g_d, *fval, *gval;
	int i;
	int s1, s2, s3, sp1, sp3, sbp3;
	s1 = sizeof(double)*nElm;
	s2 = s1*2;
	s3 = s1*3;
	sp1 = sizeof(double)*nField;
	sp3 = sp1*3;
	sbp3 = sizeof(double)*nBpt*3;
	double ppi = pi;
	//spdlog::info("pi: %.15lf\n", ppi);
	double t1, t2, t3, t4, t5, t6;
	t1 = clock();
	for(i=0; i<nField; i++)
		fldvl[i] = 0.0;

	if(!pFpt_d)
	{
		fval = (double*)malloc(sizeof(double)*nField);
		gval = (double*)malloc(sizeof(double)*3);

		cudaMalloc(&pFpt_d, sp3);
		cudaMalloc(&fldvl_d, sizeof(double)*nField);
		cudaMalloc(&g_d, sizeof(double)*3);
		cudaMalloc(&pBpt_d, sbp3);

		cudaMalloc(&pBc_d, s2);
		cudaMalloc(&pElm_d, sizeof(int)*nElm*3);
		cudaMalloc(&pNorm_d, s3);
		cudaMalloc(&pU_d, s1);

	}

	cudaMemcpy(pBpt_d, pBpt, sbp3, cudaMemcpyHostToDevice);
	cudaMemcpy(pElm_d, pElm, sizeof(int)*nElm*3, cudaMemcpyHostToDevice);
	cudaMemcpy(pFpt_d, pFpt, sp3, cudaMemcpyHostToDevice);
	cudaMemcpy(pNorm_d, pNorm, s3, cudaMemcpyHostToDevice);

	cudaMemcpy(pBc_d, pBc, s2, cudaMemcpyHostToDevice);
	cudaMemcpy(pU_d, pU, s1, cudaMemcpyHostToDevice);
	cudaThreadSynchronize();
	t2 = clock();
	// Invoke kernel 
	int threadsPerBlock = (256 < nElm) ? 256 : nElm; 
	int blocksPerGrid = (nElm + threadsPerBlock - 1) / threadsPerBlock;

	PotentialBEMGPU<<<blocksPerGrid,threadsPerBlock>>>(nField, pFpt_d, fldvl_d, g_d, nElm, pBpt_d, pBc_d, pElm_d, 
		pNorm_d, pU_d/*, pAQ_d, pBQ_d, pCQ_d, pXeQ_d, pE1_d, pE2_d, pE3_d, pTheta_d, pCscf_d, pSncf_d*/);
	cudaThreadSynchronize();
	t3 = clock();

	cudaMemcpy(fval, fldvl_d, sizeof(double)*nField, cudaMemcpyDeviceToHost);
	//cudaMemcpy(gval, g_d, sizeof(double)*3, cudaMemcpyDeviceToHost);
	cudaThreadSynchronize();
	t4=clock();
	for (int k=0; k<nField; k++)
	{
		fldvl[k] = fval[k];

		fldvl[k] *= fac; /* Correct the potential at all interior points by the factor fac */
	}
	// 	spdlog::info("Init time: {}s\n", (t2-t1)/CLOCKS_PER_SEC);
	// 	spdlog::info("Compt time: {}s\n", (t3-t2)/CLOCKS_PER_SEC);
	// 	spdlog::info("Copy time: {}s\n", (t4-t3)/CLOCKS_PER_SEC);
	//free
#if 0
	{
		free(fval);
		free(gval);

		cudaFree(pFpt_d);
		cudaFree(fldvl_d);
		cudaFree(g_d);
		cudaFree(pBpt_d);

		cudaFree(pBc_d);
		cudaFree(pElm_d);
		cudaFree(pNorm_d);
		cudaFree(pU_d);

		cudaFree(pAQ_d);
		cudaFree(pBQ_d);
		cudaFree(pCQ_d);
		cudaFree(pXeQ_d);
		cudaFree(pE1_d);
		cudaFree(pE2_d);
		cudaFree(pE3_d);
		cudaFree(pTheta_d);
		cudaFree(pCscf_d);
		cudaFree(pSncf_d);
	}
#endif
}

void PotentialBEM3D::ReleaseMemory()
{
	// 	free(fval);
	// 	free(gval);
	// 
	// 	cudaFree(pFpt_d);
	// 	cudaFree(fldvl_d);
	// 	cudaFree(g_d);
	// 	cudaFree(pBpt_d);
	// 
	// 	cudaFree(pBc_d);
	// 	cudaFree(pElm_d);
	// 	cudaFree(pNorm_d);
	// 	cudaFree(pU_d);
	// 	
	// 	cudaFree(pAQ_d);
	// 	cudaFree(pBQ_d);
	// 	cudaFree(pCQ_d);
	// 	cudaFree(pXeQ_d);
	// 	cudaFree(pE1_d);
	// 	cudaFree(pE2_d);
	// 	cudaFree(pE3_d);
	// 	cudaFree(pTheta_d);
	// 	cudaFree(pCscf_d);
	// 	cudaFree(pSncf_d);
}
#endif
