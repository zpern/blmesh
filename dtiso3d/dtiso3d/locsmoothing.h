#ifndef __loc_smoothing_h__
#define __loc_smoothing_h__

#pragma warning(disable:4786)

#include "iso3d_define.h"
#include "spr.h"
#include "optdihangle.h"
#include <algorithm>
#include <limits>

class LocSmoother
{
public:
	LocSmoother(SPRImpl* aSprimpl);

	/* 初始化函数 */
	int setMesh_LocSmoothing(divide *divi, REAL *quals);
	int setMesh_LocSmoothing(INTEGER elems[][4], int numOfElems, REAL *quals);
	int setCentNode_LocSmoothing(int iCentNod);
	int getMeshQuality_LocSmoothing(double *elemQuality, int *size);

	void free_LocSmoothing();

	/* returns the basis B of S union M. S is a set of points known
	   to be in the basis */
	void findbasis(REAL S[][3],
		REAL M[][3],
		REAL B[][3],
		int sizeS,
		int sizeM,
		int *sizeB);

	/* finds the point on the convex hull of P nearest
	   the origin */
	void minconvexhullpoint(REAL P[][3],
		int sizeP,
		REAL nearest[]);

	/* given two values a and b and their gradients, compute the
	   gradient of their product grad(a*b) */
	void gradproduct(REAL a,
		REAL b,
		REAL grada[3],
		REAL gradb[3],
		REAL prod[3]);

	/* given two values top and bottom and their gradients, compute the
	   gradient of their quotient grad(top / bottom) */
	void gradquotient(REAL top,
		REAL bot,
		REAL gradtop[3],
		REAL gradbot[3],
		REAL quot[3]);

	/* compute Z, a quantity associated with circumradius computation
	   TODO this code is lifted from Jonathan's tetcircumcenter computation
	   in primitives.c */
	REAL getZ(REAL *tetorg,
		REAL *tetdest,
		REAL *tetfapex,
		REAL *tettapex);

	/* compute the (square) of the minimum sine
	   of all the dihedral angles in the tet defined
	   by the four vertices (vtx1, vtx2, vtx3, vtx4)
	*/
	REAL minsine(APOINT vtx1, APOINT vtx2, APOINT vtx3, APOINT vtx4);

	REAL radtodeg(REAL inangle);

	/* compute the minimum or maximum angle of the tet defined
	  by the four vertices (vtx1, vtx2, vtx3, vtx4)
   */
	REAL minmaxangle(APOINT vtx1, APOINT vtx2, APOINT vtx3, APOINT vtx4, bool max);

	/* warp the sine of the dihedral angle to penalize obtuse angles more than acute */
	REAL warpsine(REAL sine);

	/* compute the (square) of the minimum sine
	   of all the dihedral angles in the tet defined
	   by the four vertices (vtx1, vtx2, vtx3, vtx4)
	*/
	REAL warpedminsine(APOINT vtx1, APOINT vtx2, APOINT vtx3, APOINT vtx4);

	/* compute the (square) of the minimum sine
	   of all the dihedral angles in the tet defined
	   by the four vertices (vtx1, vtx2, vtx3, vtx4)
	*/
	REAL minsineandedgeratio(APOINT vtx1, APOINT vtx2, APOINT vtx3, APOINT vtx4);

	/* compute the mean of the sines
	   of all the dihedral angles in the tet defined
	   by the four vertices (vtx1, vtx2, vtx3, vtx4)
	*/
	REAL meansine(APOINT vtx1, APOINT vtx2, APOINT vtx3, APOINT vtx4);

	/* the inradius to circumradius ratio */
	REAL radiusratio(APOINT vtx1, APOINT vtx2, APOINT vtx3, APOINT vtx4);

	/* compute the ratio of the tet volume to the cube of
	   the rms edge length */
	REAL vlrms3ratio(APOINT vtx1, APOINT vtx2, APOINT vtx3, APOINT vtx4);

	REAL tetquality(APOINT vtx1, APOINT vtx2, APOINT vtx3, APOINT vtx4, int qualmeasure);

	/* 网格质量优化代码，以下函数将网格优化定义在一个局部网格上 */
	void getactiveset(Sphere sph,
		int nSph,
		REAL quals[],
		REAL qualgrads[][3],
		REAL activegrads[][3],
		int *numactive,
		REAL worstqual,
		int qualmeasure);

	REAL getinitialalpha(INTEGER iNod,
		Sphere sph,
		int nSph,
		REAL quals[],
		REAL qualgrads[][3],
		REAL d[3],
		REAL r,
		REAL worstqual,
		int qualmeasure);

	void nonsmoothlinesearch(INTEGER iNod,
		Sphere sph,
		int nSph,
		REAL d[],
		REAL inworstqual,
		REAL *ouworstqual,
		REAL *alpha,
		REAL r,
		int qualmeasure);

	void getoptinfo(INTEGER iNod,
		INTEGER iElem,
		REAL *qual,
		REAL qualgrad[][3],
		REAL *volume,
		REAL volumegrad[3],
		int qualmeasure);

	bool combinedSmoothing(INTEGER iNod,
		Sphere sph,
		int nSph,
		bool bSmartLap,
		REAL inworstqual,
		REAL qualthreshold,
		REAL *outworstqual,
		REAL *outquals,
		int qualmeasure);

public:
	Elem *g_pLocSmoothingElems;
	MBLNode *g_pLocSmoothingNodes;
	TetraElemQual *g_pLocSmoothingQuals;
	int g_nLocSmoothingElems;
	int g_nLocSmoothingNodes;
	int g_nLocSmoothingCentNode;

	SPRImpl* sprimpl;
};


/* 辅助函数 */
inline int compareFloatValues(const void *arg1, const void *arg2) {
	float val1 = *((float*)arg1);
	float val2 = *((float*)arg2);

	if (val1 < val2)		 return -1;
	else if (val1 > val2) return  1;
	else				  	 return  0;
}

/* 注意qualArray1和qualArray2中的值都是增序排列的 */
inline int compareQualArray(float qualArray1[], float qualArray2[], int size1, int size2, float epsilon)
{
	int i = 0, size = size1 < size2 ? size1 : size2;
	assert(size1 > 0 && size2 > 0);

	while (fabs(qualArray1[i] - qualArray2[i]) <= epsilon && i < size)
		i++;

	return i < size ? (qualArray1[i] < qualArray2[i] ? -1 : 1) : (size1 == size2 ? 0 : (size1 < size2 ? 1 : -1));
}

inline int compareDoubleValues(const void *arg1, const void *arg2)
{
	double val1 = *((double*)arg1);
	double val2 = *((double*)arg2);

	if (val1 < val2)		 return -1;
	else if (val1 > val2) return  1;
	else				  	 return  0;
}

/* 注意qualArray1和qualArray2中的值都是增序排列的 */
inline int compareQualArray(double qualArray1[], double qualArray2[], int size1, int size2, double epsilon)
{
	int i = 0, size = size1 < size2 ? size1 : size2;
	assert(size1 > 0 && size2 > 0);

	while (fabs(qualArray1[i] - qualArray2[i]) <= epsilon && i < size)
		i++;

	return i < size ? (qualArray1[i] < qualArray2[i] ? -1 : 1) : (size1 == size2 ? 0 : (size1 < size2 ? 1 : -1));
}

#endif