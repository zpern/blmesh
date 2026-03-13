// Octree.cpp: implementation of the Octree class.
//

#include <spdlog/spdlog.h> 
 #include "../include/Octree.h"
#include <queue>
#include <vector>
#include <algorithm>
//#include "geom_func.h"
#ifndef __APPLE__
#include <omp.h>
#endif
#include <assert.h>
//#include "boost/multiprecision/float128.hpp"
//for test (need to be deleted later)
#include "../include/common.h"
////
using namespace OCT;

namespace OCT{
//int Octree::time_machine = 0;tri_inter_with_a_line

bool tri_inter_with_a_line_test(double* A, double* B, double* C, double* P, double* Q);
void printpls(MBLNode* node, int *elm, int m, int n)
{
	int i, idx0, idx1, pidx;
	FILE *fout = NULL;

	fout = fopen("ti.pls", "w");

	idx0 = m;
	idx1 = n;

	fprintf(fout, "%d %d 0 0 0 0\n", 2, 6);
	for (i = 0; i < 3; i++)
	{
		pidx = elm[4 * idx0 + i];
		fprintf(fout, "%d %20.15lf %20.15lf %20.15lf \n", i + 1, node[pidx].coord[0], node[pidx].coord[1], node[pidx].coord[2]);
	}

	for (i = 0; i < 3; i++)
	{
		pidx = elm[4 * idx1 + i];
		fprintf(fout, "%d %20.15lf %20.15lf %20.15lf \n", i + 4, node[pidx].coord[0], node[pidx].coord[1], node[pidx].coord[2]);
	}

	fprintf(fout, "1 %d %d %d %d\n", 1, 2, 3, 1);
	fprintf(fout, "2 %d %d %d %d\n", 4, 5, 6, 1);

	fclose(fout);
}


void printpls(double *A, double *B, double*C, double*O, double *P, double*Q,std::string filename="tii")
{
	int i, idx0, idx1, pidx;
	FILE *fout = NULL;

	fout = fopen(filename.data(), "w");

	fprintf(fout, "%d %d 0 0 0 0\n", 2, 6);

	fprintf(fout, "%d %20.15lf %20.15lf %20.15lf \n", 1, A[0], A[1], A[2]);
	fprintf(fout, "%d %20.15lf %20.15lf %20.15lf \n", 2, B[0], B[1], B[2]);
	fprintf(fout, "%d %20.15lf %20.15lf %20.15lf \n", 3, C[0], C[1], C[2]);


	fprintf(fout, "%d %20.15lf %20.15lf %20.15lf \n", 4, O[0], O[1], O[2]);
	fprintf(fout, "%d %20.15lf %20.15lf %20.15lf \n", 5, P[0], P[1], P[2]);
	fprintf(fout, "%d %20.15lf %20.15lf %20.15lf \n", 6, Q[0], Q[1], Q[2]);

	fprintf(fout, "1 %d %d %d %d\n", 1, 2, 3, 1);
	fprintf(fout, "2 %d %d %d %d\n", 4, 5, 6, 1);

	fclose(fout);
}

//******************************************************************************/

using namespace std;
#pragma region predicates
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef CPU86
#include <double.h>
#endif /* CPU86 **/
#ifdef LINUX
#include <fpu_control.h>
#endif /* LINUX **/
#define INEXACT                          /* Nothing **/
/* #define INEXACT volatile **/

/* #define REAL double **/                      /* double or double **/
#define REALPRINT doubleprint
#define REALRAND doublerand
#define NARROWRAND narrowdoublerand
#define UNIFORMRAND uniformdoublerand

											   /* Which of the following two methods of finding the absolute values is      **/
											   /*   fastest is compiler-dependent.  A few compilers can inline and optimize **/
											   /*   the fabs() call; but most will incur the overhead of a function call,   **/
											   /*   which is disastrously slow.  A faster way on IEEE machines might be to  **/
											   /*   mask the appropriate bit, but that's difficult to do in C.              **/

#define Absolute(a)  ((a) >= 0.0 ? (a) : -(a))
											   /* #define Absolute(a)  fabs(a) **/

											   /* Many of the operations are broken up into two pieces, a main part that    **/
											   /*   performs an approximate operation, and a "tail" that computes the       **/
											   /*   roundoff error of that operation.                                       **/
											   /*                                                                           **/
											   /* The operations Fast_Two_Sum(), Fast_Two_Diff(), Two_Sum(), Two_Diff(),    **/
											   /*   Split(), and Two_Product() are all implemented as described in the      **/
											   /*   reference.  Each of these macros requires certain variables to be       **/
											   /*   defined in the calling routine.  The variables `bvirt', `c', `abig',    **/
											   /*   `_i', `_j', `_k', `_l', `_m', and `_n' are declared `INEXACT' because   **/
											   /*   they store the result of an operation that may incur roundoff error.    **/
											   /*   The input parameter `x' (or the highest numbered `x_' parameter) must   **/
											   /*   also be declared `INEXACT'.                                             **/

#define Fast_Two_Sum_Tail(a, b, x, y) \
  bvirt = x - a; \
  y = b - bvirt

#define Fast_Two_Sum(a, b, x, y) \
  x = (REAL) (a + b); \
  Fast_Two_Sum_Tail(a, b, x, y)

#define Fast_Two_Diff_Tail(a, b, x, y) \
  bvirt = a - x; \
  y = bvirt - b

#define Fast_Two_Diff(a, b, x, y) \
  x = (REAL) (a - b); \
  Fast_Two_Diff_Tail(a, b, x, y)

#define Two_Sum_Tail(a, b, x, y) \
  bvirt = (REAL) (x - a); \
  avirt = x - bvirt; \
  bround = b - bvirt; \
  around = a - avirt; \
  y = around + bround

#define Two_Sum(a, b, x, y) \
  x = (REAL) (a + b); \
  Two_Sum_Tail(a, b, x, y)

#define Two_Diff_Tail(a, b, x, y) \
  bvirt = (REAL) (a - x); \
  avirt = x + bvirt; \
  bround = bvirt - b; \
  around = a - avirt; \
  y = around + bround

#define Two_Diff(a, b, x, y) \
  x = (REAL) (a - b); \
  Two_Diff_Tail(a, b, x, y)

#define Split(a, ahi, alo) \
  c = (REAL) (splitter * a); \
  abig = (REAL) (c - a); \
  ahi = c - abig; \
  alo = a - ahi

#define Two_Product_Tail(a, b, x, y) \
  Split(a, ahi, alo); \
  Split(b, bhi, blo); \
  err1 = x - (ahi * bhi); \
  err2 = err1 - (alo * bhi); \
  err3 = err2 - (ahi * blo); \
  y = (alo * blo) - err3

#define Two_Product(a, b, x, y) \
  x = (REAL) (a * b); \
  Two_Product_Tail(a, b, x, y)

											   /* Two_Product_Presplit() is Two_Product() where one of the inputs has       **/
											   /*   already been split.  Avoids redundant splitting.                        **/

#define Two_Product_Presplit(a, b, bhi, blo, x, y) \
  x = (REAL) (a * b); \
  Split(a, ahi, alo); \
  err1 = x - (ahi * bhi); \
  err2 = err1 - (alo * bhi); \
  err3 = err2 - (ahi * blo); \
  y = (alo * blo) - err3

											   /* Two_Product_2Presplit() is Two_Product() where both of the inputs have    **/
											   /*   already been split.  Avoids redundant splitting.                        **/

#define Two_Product_2Presplit(a, ahi, alo, b, bhi, blo, x, y) \
  x = (REAL) (a * b); \
  err1 = x - (ahi * bhi); \
  err2 = err1 - (alo * bhi); \
  err3 = err2 - (ahi * blo); \
  y = (alo * blo) - err3

											   /* Square() can be done more quickly than Two_Product().                     **/

#define Square_Tail(a, x, y) \
  Split(a, ahi, alo); \
  err1 = x - (ahi * ahi); \
  err3 = err1 - ((ahi + ahi) * alo); \
  y = (alo * alo) - err3

#define Square(a, x, y) \
  x = (REAL) (a * a); \
  Square_Tail(a, x, y)

											   /* Macros for summing expansions of various fixed lengths.  These are all    **/
											   /*   unrolled versions of Expansion_Sum().                                   **/

#define Two_One_Sum(a1, a0, b, x2, x1, x0) \
  Two_Sum(a0, b , _i, x0); \
  Two_Sum(a1, _i, x2, x1)

#define Two_One_Diff(a1, a0, b, x2, x1, x0) \
  Two_Diff(a0, b , _i, x0); \
  Two_Sum( a1, _i, x2, x1)

#define Two_Two_Sum(a1, a0, b1, b0, x3, x2, x1, x0) \
  Two_One_Sum(a1, a0, b0, _j, _0, x0); \
  Two_One_Sum(_j, _0, b1, x3, x2, x1)

#define Two_Two_Diff(a1, a0, b1, b0, x3, x2, x1, x0) \
  Two_One_Diff(a1, a0, b0, _j, _0, x0); \
  Two_One_Diff(_j, _0, b1, x3, x2, x1)

#define Four_One_Sum(a3, a2, a1, a0, b, x4, x3, x2, x1, x0) \
  Two_One_Sum(a1, a0, b , _j, x1, x0); \
  Two_One_Sum(a3, a2, _j, x4, x3, x2)

#define Four_Two_Sum(a3, a2, a1, a0, b1, b0, x5, x4, x3, x2, x1, x0) \
  Four_One_Sum(a3, a2, a1, a0, b0, _k, _2, _1, _0, x0); \
  Four_One_Sum(_k, _2, _1, _0, b1, x5, x4, x3, x2, x1)

#define Four_Four_Sum(a3, a2, a1, a0, b4, b3, b1, b0, x7, x6, x5, x4, x3, x2, \
                      x1, x0) \
  Four_Two_Sum(a3, a2, a1, a0, b1, b0, _l, _2, _1, _0, x1, x0); \
  Four_Two_Sum(_l, _2, _1, _0, b4, b3, x7, x6, x5, x4, x3, x2)

#define Eight_One_Sum(a7, a6, a5, a4, a3, a2, a1, a0, b, x8, x7, x6, x5, x4, \
                      x3, x2, x1, x0) \
  Four_One_Sum(a3, a2, a1, a0, b , _j, x3, x2, x1, x0); \
  Four_One_Sum(a7, a6, a5, a4, _j, x8, x7, x6, x5, x4)

#define Eight_Two_Sum(a7, a6, a5, a4, a3, a2, a1, a0, b1, b0, x9, x8, x7, \
                      x6, x5, x4, x3, x2, x1, x0) \
  Eight_One_Sum(a7, a6, a5, a4, a3, a2, a1, a0, b0, _k, _6, _5, _4, _3, _2, \
                _1, _0, x0); \
  Eight_One_Sum(_k, _6, _5, _4, _3, _2, _1, _0, b1, x9, x8, x7, x6, x5, x4, \
                x3, x2, x1)

#define Eight_Four_Sum(a7, a6, a5, a4, a3, a2, a1, a0, b4, b3, b1, b0, x11, \
                       x10, x9, x8, x7, x6, x5, x4, x3, x2, x1, x0) \
  Eight_Two_Sum(a7, a6, a5, a4, a3, a2, a1, a0, b1, b0, _l, _6, _5, _4, _3, \
                _2, _1, _0, x1, x0); \
  Eight_Two_Sum(_l, _6, _5, _4, _3, _2, _1, _0, b4, b3, x11, x10, x9, x8, \
                x7, x6, x5, x4, x3, x2)

											   /* Macros for multiplying expansions of various fixed lengths.               **/

#define Two_One_Product(a1, a0, b, x3, x2, x1, x0) \
  Split(b, bhi, blo); \
  Two_Product_Presplit(a0, b, bhi, blo, _i, x0); \
  Two_Product_Presplit(a1, b, bhi, blo, _j, _0); \
  Two_Sum(_i, _0, _k, x1); \
  Fast_Two_Sum(_j, _k, x3, x2)

#define Four_One_Product(a3, a2, a1, a0, b, x7, x6, x5, x4, x3, x2, x1, x0) \
  Split(b, bhi, blo); \
  Two_Product_Presplit(a0, b, bhi, blo, _i, x0); \
  Two_Product_Presplit(a1, b, bhi, blo, _j, _0); \
  Two_Sum(_i, _0, _k, x1); \
  Fast_Two_Sum(_j, _k, _i, x2); \
  Two_Product_Presplit(a2, b, bhi, blo, _j, _0); \
  Two_Sum(_i, _0, _k, x3); \
  Fast_Two_Sum(_j, _k, _i, x4); \
  Two_Product_Presplit(a3, b, bhi, blo, _j, _0); \
  Two_Sum(_i, _0, _k, x5); \
  Fast_Two_Sum(_j, _k, x7, x6)

#define Two_Two_Product(a1, a0, b1, b0, x7, x6, x5, x4, x3, x2, x1, x0) \
  Split(a0, a0hi, a0lo); \
  Split(b0, bhi, blo); \
  Two_Product_2Presplit(a0, a0hi, a0lo, b0, bhi, blo, _i, x0); \
  Split(a1, a1hi, a1lo); \
  Two_Product_2Presplit(a1, a1hi, a1lo, b0, bhi, blo, _j, _0); \
  Two_Sum(_i, _0, _k, _1); \
  Fast_Two_Sum(_j, _k, _l, _2); \
  Split(b1, bhi, blo); \
  Two_Product_2Presplit(a0, a0hi, a0lo, b1, bhi, blo, _i, _0); \
  Two_Sum(_1, _0, _k, x1); \
  Two_Sum(_2, _k, _j, _1); \
  Two_Sum(_l, _j, _m, _2); \
  Two_Product_2Presplit(a1, a1hi, a1lo, b1, bhi, blo, _j, _0); \
  Two_Sum(_i, _0, _n, _0); \
  Two_Sum(_1, _0, _i, x2); \
  Two_Sum(_2, _i, _k, _1); \
  Two_Sum(_m, _k, _l, _2); \
  Two_Sum(_j, _n, _k, _0); \
  Two_Sum(_1, _0, _j, x3); \
  Two_Sum(_2, _j, _i, _1); \
  Two_Sum(_l, _i, _m, _2); \
  Two_Sum(_1, _k, _i, x4); \
  Two_Sum(_2, _i, _k, x5); \
  Two_Sum(_m, _k, x7, x6)

											   /* An expansion of length two can be squared more quickly than finding the   **/
											   /*   product of two different expansions of length two, and the result is    **/
											   /*   guaranteed to have no more than six (rather than eight) components.     **/

#define Two_Square(a1, a0, x5, x4, x3, x2, x1, x0) \
  Square(a0, _j, x0); \
  _0 = a0 + a0; \
  Two_Product(a1, _0, _k, _1); \
  Two_One_Sum(_k, _1, _j, _l, _2, x1); \
  Square(a1, _j, _1); \
  Two_Two_Sum(_j, _1, _l, _2, x5, x4, x3, x2)

#define REAL double
											   /* splitter = 2^ceiling(p / 2) + 1.  Used to split floats in half.           **/
static REAL splitter;
static REAL epsilon;         /* = 2^(-p).  Used to estimate roundoff errors. **/
							 /* A set of coefficients used to calculate maximum roundoff errors.          **/
static REAL resulterrbound;
static REAL ccwerrboundA, ccwerrboundB, ccwerrboundC;
static REAL o3derrboundA, o3derrboundB, o3derrboundC;
static REAL iccerrboundA, iccerrboundB, iccerrboundC;
static REAL isperrboundA, isperrboundB, isperrboundC;
REAL exactinit()
{
	REAL half;
	REAL check, lastcheck;
	int every_other;
#ifdef LINUX
	int cword;
#endif /* LINUX **/

#ifdef CPU86
#ifdef SINGLE
	_control87(_PC_24, _MCW_PC); /* Set FPU control word for single precision. **/
#else /* not SINGLE **/
	_control87(_PC_53, _MCW_PC); /* Set FPU control word for double precision. **/
#endif /* not SINGLE **/
#endif /* CPU86 **/
#ifdef LINUX
#ifdef SINGLE
								 /*  cword = 4223; **/
	cword = 4210;                 /* set FPU control word for single precision **/
#else /* not SINGLE **/
								 /*  cword = 4735; **/
	cword = 4722;                 /* set FPU control word for double precision **/
#endif /* not SINGLE **/
	_FPU_SETCW(cword);
#endif /* LINUX **/

	every_other = 1;
	half = 0.5;
	epsilon = 1.0;
	splitter = 1.0;
	check = 1.0;
	/* Repeatedly divide `epsilon' by two until it is too small to add to    **/
	/*   one without causing roundoff.  (Also check if the sum is equal to   **/
	/*   the previous sum, for machines that round up instead of using exact **/
	/*   rounding.  Not that this library will work on such machines anyway. **/
	do {
		lastcheck = check;
		epsilon *= half;
		if (every_other) {
			splitter *= 2.0;
		}
		every_other = !every_other;
		check = 1.0 + epsilon;
	} while ((check != 1.0) && (check != lastcheck));
	splitter += 1.0;

	/* Error bounds for orientation and incircle tests. **/
	resulterrbound = (3.0 + 8.0 * epsilon) * epsilon;
	ccwerrboundA = (3.0 + 16.0 * epsilon) * epsilon;
	ccwerrboundB = (2.0 + 12.0 * epsilon) * epsilon;
	ccwerrboundC = (9.0 + 64.0 * epsilon) * epsilon * epsilon;
	o3derrboundA = (7.0 + 56.0 * epsilon) * epsilon;
	o3derrboundB = (3.0 + 28.0 * epsilon) * epsilon;
	o3derrboundC = (26.0 + 288.0 * epsilon) * epsilon * epsilon;
	iccerrboundA = (10.0 + 96.0 * epsilon) * epsilon;
	iccerrboundB = (4.0 + 48.0 * epsilon) * epsilon;
	iccerrboundC = (44.0 + 576.0 * epsilon) * epsilon * epsilon;
	isperrboundA = (16.0 + 224.0 * epsilon) * epsilon;
	isperrboundB = (5.0 + 72.0 * epsilon) * epsilon;
	isperrboundC = (71.0 + 1408.0 * epsilon) * epsilon * epsilon;
	return epsilon; /* Added by H. Si 30 Juli, 2004. **/
}

/******************************************************************************/
/*                                                                           **/
/*  grow_expansion()   Add a scalar to an expansion.                         **/
/*                                                                           **/
/*  Sets h = e + b.  See the long version of my paper for details.           **/
/*                                                                           **/
/*  Maintains the nonoverlapping property.  If round-to-even is used (as     **/
/*  with IEEE 754), maintains the strongly nonoverlapping and nonadjacent    **/
/*  properties as well.  (That is, if e has one of these properties, so      **/
/*  will h.)                                                                 **/
/*                                                                           **/
/******************************************************************************/

int grow_expansion(int elen, REAL *e, REAL b, REAL *h)
/* e and h can be the same. **/
{
	REAL Q;
	INEXACT REAL Qnew;
	int eindex;
	REAL enow;
	INEXACT REAL bvirt;
	REAL avirt, bround, around;

	Q = b;
	for (eindex = 0; eindex < elen; eindex++) {
		enow = e[eindex];
		Two_Sum(Q, enow, Qnew, h[eindex]);
		Q = Qnew;
	}
	h[eindex] = Q;
	return eindex + 1;
}

/******************************************************************************/
/*                                                                           **/
/*  grow_expansion_zeroelim()   Add a scalar to an expansion, eliminating    **/
/*                              zero components from the output expansion.   **/
/*                                                                           **/
/*  Sets h = e + b.  See the long version of my paper for details.           **/
/*                                                                           **/
/*  Maintains the nonoverlapping property.  If round-to-even is used (as     **/
/*  with IEEE 754), maintains the strongly nonoverlapping and nonadjacent    **/
/*  properties as well.  (That is, if e has one of these properties, so      **/
/*  will h.)                                                                 **/
/*                                                                           **/
/******************************************************************************/

int grow_expansion_zeroelim(int elen, REAL *e, REAL b, REAL *h)
/* e and h can be the same. **/
{
	REAL Q, hh;
	INEXACT REAL Qnew;
	int eindex, hindex;
	REAL enow;
	INEXACT REAL bvirt;
	REAL avirt, bround, around;

	hindex = 0;
	Q = b;
	for (eindex = 0; eindex < elen; eindex++) {
		enow = e[eindex];
		Two_Sum(Q, enow, Qnew, hh);
		Q = Qnew;
		if (hh != 0.0) {
			h[hindex++] = hh;
		}
	}
	if ((Q != 0.0) || (hindex == 0)) {
		h[hindex++] = Q;
	}
	return hindex;
}

/******************************************************************************/
/*                                                                           **/
/*  expansion_sum()   Sum two expansions.                                    **/
/*                                                                           **/
/*  Sets h = e + f.  See the long version of my paper for details.           **/
/*                                                                           **/
/*  Maintains the nonoverlapping property.  If round-to-even is used (as     **/
/*  with IEEE 754), maintains the nonadjacent property as well.  (That is,   **/
/*  if e has one of these properties, so will h.)  Does NOT maintain the     **/
/*  strongly nonoverlapping property.                                        **/
/*                                                                           **/
/******************************************************************************/

int expansion_sum(int elen, REAL *e, int flen, REAL *f, REAL *h)
/* e and h can be the same, but f and h cannot. **/
{
	REAL Q;
	INEXACT REAL Qnew;
	int findex, hindex, hlast;
	REAL hnow;
	INEXACT REAL bvirt;
	REAL avirt, bround, around;

	Q = f[0];
	for (hindex = 0; hindex < elen; hindex++) {
		hnow = e[hindex];
		Two_Sum(Q, hnow, Qnew, h[hindex]);
		Q = Qnew;
	}
	h[hindex] = Q;
	hlast = hindex;
	for (findex = 1; findex < flen; findex++) {
		Q = f[findex];
		for (hindex = findex; hindex <= hlast; hindex++) {
			hnow = h[hindex];
			Two_Sum(Q, hnow, Qnew, h[hindex]);
			Q = Qnew;
		}
		h[++hlast] = Q;
	}
	return hlast + 1;
}

/******************************************************************************/
/*                                                                           **/
/*  expansion_sum_zeroelim1()   Sum two expansions, eliminating zero         **/
/*                              components from the output expansion.        **/
/*                                                                           **/
/*  Sets h = e + f.  See the long version of my paper for details.           **/
/*                                                                           **/
/*  Maintains the nonoverlapping property.  If round-to-even is used (as     **/
/*  with IEEE 754), maintains the nonadjacent property as well.  (That is,   **/
/*  if e has one of these properties, so will h.)  Does NOT maintain the     **/
/*  strongly nonoverlapping property.                                        **/
/*                                                                           **/
/******************************************************************************/

int expansion_sum_zeroelim1(int elen, REAL *e, int flen, REAL *f, REAL *h)
/* e and h can be the same, but f and h cannot. **/
{
	REAL Q;
	INEXACT REAL Qnew;
	int index, findex, hindex, hlast;
	REAL hnow;
	INEXACT REAL bvirt;
	REAL avirt, bround, around;

	Q = f[0];
	for (hindex = 0; hindex < elen; hindex++) {
		hnow = e[hindex];
		Two_Sum(Q, hnow, Qnew, h[hindex]);
		Q = Qnew;
	}
	h[hindex] = Q;
	hlast = hindex;
	for (findex = 1; findex < flen; findex++) {
		Q = f[findex];
		for (hindex = findex; hindex <= hlast; hindex++) {
			hnow = h[hindex];
			Two_Sum(Q, hnow, Qnew, h[hindex]);
			Q = Qnew;
		}
		h[++hlast] = Q;
	}
	hindex = -1;
	for (index = 0; index <= hlast; index++) {
		hnow = h[index];
		if (hnow != 0.0) {
			h[++hindex] = hnow;
		}
	}
	if (hindex == -1) {
		return 1;
	}
	else {
		return hindex + 1;
	}
}

/******************************************************************************/
/*                                                                           **/
/*  expansion_sum_zeroelim2()   Sum two expansions, eliminating zero         **/
/*                              components from the output expansion.        **/
/*                                                                           **/
/*  Sets h = e + f.  See the long version of my paper for details.           **/
/*                                                                           **/
/*  Maintains the nonoverlapping property.  If round-to-even is used (as     **/
/*  with IEEE 754), maintains the nonadjacent property as well.  (That is,   **/
/*  if e has one of these properties, so will h.)  Does NOT maintain the     **/
/*  strongly nonoverlapping property.                                        **/
/*                                                                           **/
/******************************************************************************/

int expansion_sum_zeroelim2(int elen, REAL *e, int flen, REAL *f, REAL *h)
/* e and h can be the same, but f and h cannot. **/
{
	REAL Q, hh;
	INEXACT REAL Qnew;
	int eindex, findex, hindex, hlast;
	REAL enow;
	INEXACT REAL bvirt;
	REAL avirt, bround, around;

	hindex = 0;
	Q = f[0];
	for (eindex = 0; eindex < elen; eindex++) {
		enow = e[eindex];
		Two_Sum(Q, enow, Qnew, hh);
		Q = Qnew;
		if (hh != 0.0) {
			h[hindex++] = hh;
		}
	}
	h[hindex] = Q;
	hlast = hindex;
	for (findex = 1; findex < flen; findex++) {
		hindex = 0;
		Q = f[findex];
		for (eindex = 0; eindex <= hlast; eindex++) {
			enow = h[eindex];
			Two_Sum(Q, enow, Qnew, hh);
			Q = Qnew;
			if (hh != 0) {
				h[hindex++] = hh;
			}
		}
		h[hindex] = Q;
		hlast = hindex;
	}
	return hlast + 1;
}

/******************************************************************************/
/*                                                                           **/
/*  fast_expansion_sum()   Sum two expansions.                               **/
/*                                                                           **/
/*  Sets h = e + f.  See the long version of my paper for details.           **/
/*                                                                           **/
/*  If round-to-even is used (as with IEEE 754), maintains the strongly      **/
/*  nonoverlapping property.  (That is, if e is strongly nonoverlapping, h   **/
/*  will be also.)  Does NOT maintain the nonoverlapping or nonadjacent      **/
/*  properties.                                                              **/
/*                                                                           **/
/******************************************************************************/

int fast_expansion_sum(int elen, REAL *e, int flen, REAL *f, REAL *h)
/* h cannot be e or f. **/
{
	REAL Q;
	INEXACT REAL Qnew;
	INEXACT REAL bvirt;
	REAL avirt, bround, around;
	int eindex, findex, hindex;
	REAL enow, fnow;

	enow = e[0];
	fnow = f[0];
	eindex = findex = 0;
	if ((fnow > enow) == (fnow > -enow)) {
		Q = enow;
		enow = e[++eindex];
	}
	else {
		Q = fnow;
		fnow = f[++findex];
	}
	hindex = 0;
	if ((eindex < elen) && (findex < flen)) {
		if ((fnow > enow) == (fnow > -enow)) {
			Fast_Two_Sum(enow, Q, Qnew, h[0]);
			enow = e[++eindex];
		}
		else {
			Fast_Two_Sum(fnow, Q, Qnew, h[0]);
			fnow = f[++findex];
		}
		Q = Qnew;
		hindex = 1;
		while ((eindex < elen) && (findex < flen)) {
			if ((fnow > enow) == (fnow > -enow)) {
				Two_Sum(Q, enow, Qnew, h[hindex]);
				enow = e[++eindex];
			}
			else {
				Two_Sum(Q, fnow, Qnew, h[hindex]);
				fnow = f[++findex];
			}
			Q = Qnew;
			hindex++;
		}
	}
	while (eindex < elen) {
		Two_Sum(Q, enow, Qnew, h[hindex]);
		enow = e[++eindex];
		Q = Qnew;
		hindex++;
	}
	while (findex < flen) {
		Two_Sum(Q, fnow, Qnew, h[hindex]);
		fnow = f[++findex];
		Q = Qnew;
		hindex++;
	}
	h[hindex] = Q;
	return hindex + 1;
}

/******************************************************************************/
/*                                                                           **/
/*  fast_expansion_sum_zeroelim()   Sum two expansions, eliminating zero     **/
/*                                  components from the output expansion.    **/
/*                                                                           **/
/*  Sets h = e + f.  See the long version of my paper for details.           **/
/*                                                                           **/
/*  If round-to-even is used (as with IEEE 754), maintains the strongly      **/
/*  nonoverlapping property.  (That is, if e is strongly nonoverlapping, h   **/
/*  will be also.)  Does NOT maintain the nonoverlapping or nonadjacent      **/
/*  properties.                                                              **/
/*                                                                           **/
/******************************************************************************/

int fast_expansion_sum_zeroelim(int elen, REAL *e, int flen, REAL *f, REAL *h)
/* h cannot be e or f. **/
{
	REAL Q;
	INEXACT REAL Qnew;
	INEXACT REAL hh;
	INEXACT REAL bvirt;
	REAL avirt, bround, around;
	int eindex, findex, hindex;
	REAL enow, fnow;

	enow = e[0];
	fnow = f[0];
	eindex = findex = 0;
	if ((fnow > enow) == (fnow > -enow)) {
		Q = enow;
		enow = e[++eindex];
	}
	else {
		Q = fnow;
		fnow = f[++findex];
	}
	hindex = 0;
	if ((eindex < elen) && (findex < flen)) {
		if ((fnow > enow) == (fnow > -enow)) {
			Fast_Two_Sum(enow, Q, Qnew, hh);
			enow = e[++eindex];
		}
		else {
			Fast_Two_Sum(fnow, Q, Qnew, hh);
			fnow = f[++findex];
		}
		Q = Qnew;
		if (hh != 0.0) {
			h[hindex++] = hh;
		}
		while ((eindex < elen) && (findex < flen)) {
			if ((fnow > enow) == (fnow > -enow)) {
				Two_Sum(Q, enow, Qnew, hh);
				enow = e[++eindex];
			}
			else {
				Two_Sum(Q, fnow, Qnew, hh);
				fnow = f[++findex];
			}
			Q = Qnew;
			if (hh != 0.0) {
				h[hindex++] = hh;
			}
		}
	}
	while (eindex < elen) {
		Two_Sum(Q, enow, Qnew, hh);
		enow = e[++eindex];
		Q = Qnew;
		if (hh != 0.0) {
			h[hindex++] = hh;
		}
	}
	while (findex < flen) {
		Two_Sum(Q, fnow, Qnew, hh);
		fnow = f[++findex];
		Q = Qnew;
		if (hh != 0.0) {
			h[hindex++] = hh;
		}
	}
	if ((Q != 0.0) || (hindex == 0)) {
		h[hindex++] = Q;
	}
	return hindex;
}

/******************************************************************************/
/*                                                                           **/
/*  linear_expansion_sum()   Sum two expansions.                             **/
/*                                                                           **/
/*  Sets h = e + f.  See either version of my paper for details.             **/
/*                                                                           **/
/*  Maintains the nonoverlapping property.  (That is, if e is                **/
/*  nonoverlapping, h will be also.)                                         **/
/*                                                                           **/
/******************************************************************************/

int linear_expansion_sum(int elen, REAL *e, int flen, REAL *f, REAL *h)
/* h cannot be e or f. **/
{
	REAL Q, q;
	INEXACT REAL Qnew;
	INEXACT REAL R;
	INEXACT REAL bvirt;
	REAL avirt, bround, around;
	int eindex, findex, hindex;
	REAL enow, fnow;
	REAL g0;

	enow = e[0];
	fnow = f[0];
	eindex = findex = 0;
	if ((fnow > enow) == (fnow > -enow)) {
		g0 = enow;
		enow = e[++eindex];
	}
	else {
		g0 = fnow;
		fnow = f[++findex];
	}
	if ((eindex < elen) && ((findex >= flen)
		|| ((fnow > enow) == (fnow > -enow)))) {
		Fast_Two_Sum(enow, g0, Qnew, q);
		enow = e[++eindex];
	}
	else {
		Fast_Two_Sum(fnow, g0, Qnew, q);
		fnow = f[++findex];
	}
	Q = Qnew;
	for (hindex = 0; hindex < elen + flen - 2; hindex++) {
		if ((eindex < elen) && ((findex >= flen)
			|| ((fnow > enow) == (fnow > -enow)))) {
			Fast_Two_Sum(enow, q, R, h[hindex]);
			enow = e[++eindex];
		}
		else {
			Fast_Two_Sum(fnow, q, R, h[hindex]);
			fnow = f[++findex];
		}
		Two_Sum(Q, R, Qnew, q);
		Q = Qnew;
	}
	h[hindex] = q;
	h[hindex + 1] = Q;
	return hindex + 2;
}

/******************************************************************************/
/*                                                                           **/
/*  linear_expansion_sum_zeroelim()   Sum two expansions, eliminating zero   **/
/*                                    components from the output expansion.  **/
/*                                                                           **/
/*  Sets h = e + f.  See either version of my paper for details.             **/
/*                                                                           **/
/*  Maintains the nonoverlapping property.  (That is, if e is                **/
/*  nonoverlapping, h will be also.)                                         **/
/*                                                                           **/
/******************************************************************************/

int linear_expansion_sum_zeroelim(int elen, REAL *e, int flen, REAL *f,
	REAL *h)
	/* h cannot be e or f. **/
{
	REAL Q, q, hh;
	INEXACT REAL Qnew;
	INEXACT REAL R;
	INEXACT REAL bvirt;
	REAL avirt, bround, around;
	int eindex, findex, hindex;
	int count;
	REAL enow, fnow;
	REAL g0;

	enow = e[0];
	fnow = f[0];
	eindex = findex = 0;
	hindex = 0;
	if ((fnow > enow) == (fnow > -enow)) {
		g0 = enow;
		enow = e[++eindex];
	}
	else {
		g0 = fnow;
		fnow = f[++findex];
	}
	if ((eindex < elen) && ((findex >= flen)
		|| ((fnow > enow) == (fnow > -enow)))) {
		Fast_Two_Sum(enow, g0, Qnew, q);
		enow = e[++eindex];
	}
	else {
		Fast_Two_Sum(fnow, g0, Qnew, q);
		fnow = f[++findex];
	}
	Q = Qnew;
	for (count = 2; count < elen + flen; count++) {
		if ((eindex < elen) && ((findex >= flen)
			|| ((fnow > enow) == (fnow > -enow)))) {
			Fast_Two_Sum(enow, q, R, hh);
			enow = e[++eindex];
		}
		else {
			Fast_Two_Sum(fnow, q, R, hh);
			fnow = f[++findex];
		}
		Two_Sum(Q, R, Qnew, q);
		Q = Qnew;
		if (hh != 0) {
			h[hindex++] = hh;
		}
	}
	if (q != 0) {
		h[hindex++] = q;
	}
	if ((Q != 0.0) || (hindex == 0)) {
		h[hindex++] = Q;
	}
	return hindex;
}

/******************************************************************************/
/*                                                                           **/
/*  scale_expansion()   Multiply an expansion by a scalar.                   **/
/*                                                                           **/
/*  Sets h = be.  See either version of my paper for details.                **/
/*                                                                           **/
/*  Maintains the nonoverlapping property.  If round-to-even is used (as     **/
/*  with IEEE 754), maintains the strongly nonoverlapping and nonadjacent    **/
/*  properties as well.  (That is, if e has one of these properties, so      **/
/*  will h.)                                                                 **/
/*                                                                           **/
/******************************************************************************/

int scale_expansion(int elen, REAL *e, REAL b, REAL *h)
/* e and h cannot be the same. **/
{
	INEXACT REAL Q;
	INEXACT REAL sum;
	INEXACT REAL product1;
	REAL product0;
	int eindex, hindex;
	REAL enow;
	INEXACT REAL bvirt;
	REAL avirt, bround, around;
	INEXACT REAL c;
	INEXACT REAL abig;
	REAL ahi, alo, bhi, blo;
	REAL err1, err2, err3;

	Split(b, bhi, blo);
	Two_Product_Presplit(e[0], b, bhi, blo, Q, h[0]);
	hindex = 1;
	for (eindex = 1; eindex < elen; eindex++) {
		enow = e[eindex];
		Two_Product_Presplit(enow, b, bhi, blo, product1, product0);
		Two_Sum(Q, product0, sum, h[hindex]);
		hindex++;
		Two_Sum(product1, sum, Q, h[hindex]);
		hindex++;
	}
	h[hindex] = Q;
	return elen + elen;
}

/******************************************************************************/
/*                                                                           **/
/*  scale_expansion_zeroelim()   Multiply an expansion by a scalar,          **/
/*                               eliminating zero components from the        **/
/*                               output expansion.                           **/
/*                                                                           **/
/*  Sets h = be.  See either version of my paper for details.                **/
/*                                                                           **/
/*  Maintains the nonoverlapping property.  If round-to-even is used (as     **/
/*  with IEEE 754), maintains the strongly nonoverlapping and nonadjacent    **/
/*  properties as well.  (That is, if e has one of these properties, so      **/
/*  will h.)                                                                 **/
/*                                                                           **/
/******************************************************************************/

int scale_expansion_zeroelim(int elen, REAL *e, REAL b, REAL *h)
/* e and h cannot be the same. **/
{
	INEXACT REAL Q, sum;
	REAL hh;
	INEXACT REAL product1;
	REAL product0;
	int eindex, hindex;
	REAL enow;
	INEXACT REAL bvirt;
	REAL avirt, bround, around;
	INEXACT REAL c;
	INEXACT REAL abig;
	REAL ahi, alo, bhi, blo;
	REAL err1, err2, err3;

	Split(b, bhi, blo);
	Two_Product_Presplit(e[0], b, bhi, blo, Q, hh);
	hindex = 0;
	if (hh != 0) {
		h[hindex++] = hh;
	}
	for (eindex = 1; eindex < elen; eindex++) {
		enow = e[eindex];
		Two_Product_Presplit(enow, b, bhi, blo, product1, product0);
		Two_Sum(Q, product0, sum, hh);
		if (hh != 0) {
			h[hindex++] = hh;
		}
		Fast_Two_Sum(product1, sum, Q, hh);
		if (hh != 0) {
			h[hindex++] = hh;
		}
	}
	if ((Q != 0.0) || (hindex == 0)) {
		h[hindex++] = Q;
	}
	return hindex;
}

/******************************************************************************/
/*                                                                           **/
/*  compression()   compression an expansion.                                      **/
/*                                                                           **/
/*  See the long version of my paper for details.                            **/
/*                                                                           **/
/*  Maintains the nonoverlapping property.  If round-to-even is used (as     **/
/*  with IEEE 754), then any nonoverlapping expansion is converted to a      **/
/*  nonadjacent expansion.                                                   **/
/*                                                                           **/
/******************************************************************************/

int compression(int elen, REAL *e, REAL *h)
/* e and h may be the same. **/
{
	REAL Q, q;
	INEXACT REAL Qnew;
	int eindex, hindex;
	INEXACT REAL bvirt;
	REAL enow, hnow;
	int top, bottom;

	bottom = elen - 1;
	Q = e[bottom];
	for (eindex = elen - 2; eindex >= 0; eindex--) {
		enow = e[eindex];
		Fast_Two_Sum(Q, enow, Qnew, q);
		if (q != 0) {
			h[bottom--] = Qnew;
			Q = q;
		}
		else {
			Q = Qnew;
		}
	}
	top = 0;
	for (hindex = bottom + 1; hindex < elen; hindex++) {
		hnow = h[hindex];
		Fast_Two_Sum(hnow, Q, Qnew, q);
		if (q != 0) {
			h[top++] = q;
		}
		Q = Qnew;
	}
	h[top] = Q;
	return top + 1;
}

/******************************************************************************/
/*                                                                           **/
/*  estimate()   Produce a one-word estimate of an expansion's value.        **/
/*                                                                           **/
/*  See either version of my paper for details.                              **/
/*                                                                           **/
/******************************************************************************/

REAL estimate(int elen, REAL *e)
{
	REAL Q;
	int eindex;

	Q = e[0];
	for (eindex = 1; eindex < elen; eindex++) {
		Q += e[eindex];
	}
	return Q;
}

REAL Two_Divide(REAL a, REAL b)
{
	INEXACT REAL bvirt;
	REAL avirt, bround, around;
	INEXACT REAL c;
	INEXACT REAL abig;
	REAL ahi, alo, bhi, blo;
	REAL err1, err2, err3;
	INEXACT REAL _i, _j, _k;
	REAL _0, _1, _2;


	REAL x3, x2, x1, x0, avir1, avir2;
	REAL result1, result2, result3, result4;
	REAL temp;
	x3 = a / b;
	Two_Product(x3, b, avir1, avir2);
	Two_One_Sum(a, (-avir1), (-avir2), x2, x1, x0);
	x2 /= b;
	x1 /= b;
	x0 /= b;
	Two_Two_Sum(x3, x2, x1, x0, result1, result2, result3, result4);
	return result1;
}

/* Add by Wu Bin **/
/* For Test **/
REAL fixedSplitPoint(REAL s1, REAL s2, REAL pnt1, REAL pnt2)
{
	INEXACT REAL bvirt;
	REAL avirt, bround, around;
	INEXACT REAL c;
	INEXACT REAL abig;
	REAL ahi, alo, bhi, blo;
	REAL err1, err2, err3;
	INEXACT REAL _i, _j, _k;
	REAL _0, _1, _2;

	REAL result = 0.0;
	REAL wu[30];
	Two_Sum(s1, s2, wu[0], wu[1]);
	REAL bili = Two_Divide(s1, wu[0]);
	//bili = s1/(wu[0]+wu[1]);
	Two_Sum(1.0, -bili, wu[2], wu[3]);
	Two_Product(bili, pnt2, wu[4], wu[5]);
	Two_Product(wu[2], pnt1, wu[6], wu[7]);
	Two_Product(wu[3], pnt1, wu[8], wu[9]);
	Four_Two_Sum(wu[4], wu[5], wu[6], wu[7], wu[8], wu[9], wu[10], wu[11], wu[12], wu[13], wu[14], wu[15]);

	//result = wu_e/(wu_i+wu_j) + wu_f/(wu_i+wu_j) + wu_g/(wu_i+wu_j) + wu_h/(wu_i+wu_j);
	return result = wu[10] + wu[11];
}
/* End of add **/

/******************************************************************************/
/*                                                                           **/
/*  orient2dfast()   Approximate 2D orientation test.  Nonrobust.            **/
/*  orient2dexact()   Exact 2D orientation test.  Robust.                    **/
/*  orient2dslow()   Another exact 2D orientation test.  Robust.             **/
/*  orient2d()   Adaptive exact 2D orientation test.  Robust.                **/
/*                                                                           **/
/*               Return a positive value if the points pa, pb, and pc occur  **/
/*               in counterclockwise order; a negative value if they occur   **/
/*               in clockwise order; and zero if they are collinear.  The    **/
/*               result is also a rough approximation of twice the signed    **/
/*               area of the triangle defined by the three points.           **/
/*                                                                           **/
/*  Only the first and last routine should be used; the middle two are for   **/
/*  timings.                                                                 **/
/*                                                                           **/
/*  The last three use exact arithmetic to ensure a correct answer.  The     **/
/*  result returned is the determinant of a matrix.  In orient2d() only,     **/
/*  this determinant is computed adaptively, in the sense that exact         **/
/*  arithmetic is used only to the degree it is needed to ensure that the    **/
/*  returned value has the correct sign.  Hence, orient2d() is usually quite **/
/*  fast, but will run more slowly when the input points are collinear or    **/
/*  nearly so.                                                               **/
/*                                                                           **/
/******************************************************************************/

REAL orient2dfast(REAL *pa, REAL *pb, REAL *pc)
{
	REAL acx, bcx, acy, bcy;

	acx = pa[0] - pc[0];
	bcx = pb[0] - pc[0];
	acy = pa[1] - pc[1];
	bcy = pb[1] - pc[1];
	return acx * bcy - acy * bcx;
}

REAL orient2dexact(REAL *pa, REAL *pb, REAL *pc)
{
	INEXACT REAL axby1, axcy1, bxcy1, bxay1, cxay1, cxby1;
	REAL axby0, axcy0, bxcy0, bxay0, cxay0, cxby0;
	REAL aterms[4], bterms[4], cterms[4];
	INEXACT REAL aterms3, bterms3, cterms3;
	REAL v[8], w[12];
	int vlength, wlength;

	INEXACT REAL bvirt;
	REAL avirt, bround, around;
	INEXACT REAL c;
	INEXACT REAL abig;
	REAL ahi, alo, bhi, blo;
	REAL err1, err2, err3;
	INEXACT REAL _i, _j;
	REAL _0;

	Two_Product(pa[0], pb[1], axby1, axby0);
	Two_Product(pa[0], pc[1], axcy1, axcy0);
	Two_Two_Diff(axby1, axby0, axcy1, axcy0,
		aterms3, aterms[2], aterms[1], aterms[0]);
	aterms[3] = aterms3;

	Two_Product(pb[0], pc[1], bxcy1, bxcy0);
	Two_Product(pb[0], pa[1], bxay1, bxay0);
	Two_Two_Diff(bxcy1, bxcy0, bxay1, bxay0,
		bterms3, bterms[2], bterms[1], bterms[0]);
	bterms[3] = bterms3;

	Two_Product(pc[0], pa[1], cxay1, cxay0);
	Two_Product(pc[0], pb[1], cxby1, cxby0);
	Two_Two_Diff(cxay1, cxay0, cxby1, cxby0,
		cterms3, cterms[2], cterms[1], cterms[0]);
	cterms[3] = cterms3;

	vlength = fast_expansion_sum_zeroelim(4, aterms, 4, bterms, v);
	wlength = fast_expansion_sum_zeroelim(vlength, v, 4, cterms, w);

	return w[wlength - 1];
}

REAL orient2dslow(REAL *pa, REAL *pb, REAL *pc)
{
	INEXACT REAL acx, acy, bcx, bcy;
	REAL acxtail, acytail;
	REAL bcxtail, bcytail;
	REAL negate, negatetail;
	REAL axby[8], bxay[8];
	INEXACT REAL axby7, bxay7;
	REAL deter[16];
	int deterlen;

	INEXACT REAL bvirt;
	REAL avirt, bround, around;
	INEXACT REAL c;
	INEXACT REAL abig;
	REAL a0hi, a0lo, a1hi, a1lo, bhi, blo;
	REAL err1, err2, err3;
	INEXACT REAL _i, _j, _k, _l, _m, _n;
	REAL _0, _1, _2;

	Two_Diff(pa[0], pc[0], acx, acxtail);
	Two_Diff(pa[1], pc[1], acy, acytail);
	Two_Diff(pb[0], pc[0], bcx, bcxtail);
	Two_Diff(pb[1], pc[1], bcy, bcytail);

	Two_Two_Product(acx, acxtail, bcy, bcytail,
		axby7, axby[6], axby[5], axby[4],
		axby[3], axby[2], axby[1], axby[0]);
	axby[7] = axby7;
	negate = -acy;
	negatetail = -acytail;
	Two_Two_Product(bcx, bcxtail, negate, negatetail,
		bxay7, bxay[6], bxay[5], bxay[4],
		bxay[3], bxay[2], bxay[1], bxay[0]);
	bxay[7] = bxay7;

	deterlen = fast_expansion_sum_zeroelim(8, axby, 8, bxay, deter);

	return deter[deterlen - 1];
}

REAL orient2dadapt(REAL *pa, REAL *pb, REAL *pc, REAL detsum)
{
	INEXACT REAL acx, acy, bcx, bcy;
	REAL acxtail, acytail, bcxtail, bcytail;
	INEXACT REAL detleft, detright;
	REAL detlefttail, detrighttail;
	REAL det, errbound;
	REAL B[4], C1[8], C2[12], D[16];
	INEXACT REAL B3;
	int C1length, C2length, Dlength;
	REAL u[4];
	INEXACT REAL u3;
	INEXACT REAL s1, t1;
	REAL s0, t0;

	INEXACT REAL bvirt;
	REAL avirt, bround, around;
	INEXACT REAL c;
	INEXACT REAL abig;
	REAL ahi, alo, bhi, blo;
	REAL err1, err2, err3;
	INEXACT REAL _i, _j;
	REAL _0;

	acx = (REAL)(pa[0] - pc[0]);
	bcx = (REAL)(pb[0] - pc[0]);
	acy = (REAL)(pa[1] - pc[1]);
	bcy = (REAL)(pb[1] - pc[1]);

	Two_Product(acx, bcy, detleft, detlefttail);
	Two_Product(acy, bcx, detright, detrighttail);

	Two_Two_Diff(detleft, detlefttail, detright, detrighttail,
		B3, B[2], B[1], B[0]);
	B[3] = B3;

	det = estimate(4, B);
	errbound = ccwerrboundB * detsum;
	if ((det >= errbound) || (-det >= errbound)) {
		return det;
	}

	Two_Diff_Tail(pa[0], pc[0], acx, acxtail);
	Two_Diff_Tail(pb[0], pc[0], bcx, bcxtail);
	Two_Diff_Tail(pa[1], pc[1], acy, acytail);
	Two_Diff_Tail(pb[1], pc[1], bcy, bcytail);

	if ((acxtail == 0.0) && (acytail == 0.0)
		&& (bcxtail == 0.0) && (bcytail == 0.0)) {
		return det;
	}

	errbound = ccwerrboundC * detsum + resulterrbound * Absolute(det);
	det += (acx * bcytail + bcy * acxtail)
		- (acy * bcxtail + bcx * acytail);
	if ((det >= errbound) || (-det >= errbound)) {
		return det;
	}

	Two_Product(acxtail, bcy, s1, s0);
	Two_Product(acytail, bcx, t1, t0);
	Two_Two_Diff(s1, s0, t1, t0, u3, u[2], u[1], u[0]);
	u[3] = u3;
	C1length = fast_expansion_sum_zeroelim(4, B, 4, u, C1);

	Two_Product(acx, bcytail, s1, s0);
	Two_Product(acy, bcxtail, t1, t0);
	Two_Two_Diff(s1, s0, t1, t0, u3, u[2], u[1], u[0]);
	u[3] = u3;
	C2length = fast_expansion_sum_zeroelim(C1length, C1, 4, u, C2);

	Two_Product(acxtail, bcytail, s1, s0);
	Two_Product(acytail, bcxtail, t1, t0);
	Two_Two_Diff(s1, s0, t1, t0, u3, u[2], u[1], u[0]);
	u[3] = u3;
	Dlength = fast_expansion_sum_zeroelim(C2length, C2, 4, u, D);

	return(D[Dlength - 1]);
}

REAL orient2d(REAL *pa, REAL *pb, REAL *pc)
{
	REAL detleft, detright, det;
	REAL detsum, errbound;

	detleft = (pa[0] - pc[0]) * (pb[1] - pc[1]);
	detright = (pa[1] - pc[1]) * (pb[0] - pc[0]);
	det = detleft - detright;

	if (detleft > 0.0) {
		if (detright <= 0.0) {
			return det;
		}
		else {
			detsum = detleft + detright;
		}
	}
	else if (detleft < 0.0) {
		if (detright >= 0.0) {
			return det;
		}
		else {
			detsum = -detleft - detright;
		}
	}
	else {
		return det;
	}

	errbound = ccwerrboundA * detsum;
	if ((det >= errbound) || (-det >= errbound)) {
		return det;
	}

	return orient2dadapt(pa, pb, pc, detsum);
}

/******************************************************************************/
/*                                                                           **/
/*  orient3dfast()   Approximate 3D orientation test.  Nonrobust.            **/
/*  orient3dexact()   Exact 3D orientation test.  Robust.                    **/
/*  orient3dslow()   Another exact 3D orientation test.  Robust.             **/
/*  orient3d()   Adaptive exact 3D orientation test.  Robust.                **/
/*                                                                           **/
/*               Return a positive value if the point pd lies below the      **/
/*               plane passing through pa, pb, and pc; "below" is defined so **/
/*               that pa, pb, and pc appear in counterclockwise order when   **/
/*               viewed from above the plane.  Returns a negative value if   **/
/*               pd lies above the plane.  Returns zero if the points are    **/
/*               coplanar.  The result is also a rough approximation of six  **/
/*               times the signed volume of the tetrahedron defined by the   **/
/*               four points.                                                **/
/*                                                                           **/
/*  Only the first and last routine should be used; the middle two are for   **/
/*  timings.                                                                 **/
/*                                                                           **/
/*  The last three use exact arithmetic to ensure a correct answer.  The     **/
/*  result returned is the determinant of a matrix.  In orient3d() only,     **/
/*  this determinant is computed adaptively, in the sense that exact         **/
/*  arithmetic is used only to the degree it is needed to ensure that the    **/
/*  returned value has the correct sign.  Hence, orient3d() is usually quite **/
/*  fast, but will run more slowly when the input points are coplanar or     **/
/*  nearly so.                                                               **/
/*                                                                           **/
/******************************************************************************/

REAL orient3dfast(REAL *pa, REAL *pb, REAL *pc, REAL *pd)
{
	REAL adx, bdx, cdx;
	REAL ady, bdy, cdy;
	REAL adz, bdz, cdz;

	adx = pa[0] - pd[0];
	bdx = pb[0] - pd[0];
	cdx = pc[0] - pd[0];
	ady = pa[1] - pd[1];
	bdy = pb[1] - pd[1];
	cdy = pc[1] - pd[1];
	adz = pa[2] - pd[2];
	bdz = pb[2] - pd[2];
	cdz = pc[2] - pd[2];

	return adx * (bdy * cdz - bdz * cdy)
		+ bdx * (cdy * adz - cdz * ady)
		+ cdx * (ady * bdz - adz * bdy);
}

REAL orient3dexact(REAL *pa, REAL *pb, REAL *pc, REAL *pd)

{
	INEXACT REAL axby1, bxcy1, cxdy1, dxay1, axcy1, bxdy1;
	INEXACT REAL bxay1, cxby1, dxcy1, axdy1, cxay1, dxby1;
	REAL axby0, bxcy0, cxdy0, dxay0, axcy0, bxdy0;
	REAL bxay0, cxby0, dxcy0, axdy0, cxay0, dxby0;
	REAL ab[4], bc[4], cd[4], da[4], ac[4], bd[4];
	REAL temp8[8];
	int templen;
	REAL abc[12], bcd[12], cda[12], dab[12];
	int abclen, bcdlen, cdalen, dablen;
	REAL adet[24], bdet[24], cdet[24], ddet[24];
	int alen, blen, clen, dlen;
	REAL abdet[48], cddet[48];
	int ablen, cdlen;
	REAL deter[96];
	int deterlen;
	int i;

	INEXACT REAL bvirt;
	REAL avirt, bround, around;
	INEXACT REAL c;
	INEXACT REAL abig;
	REAL ahi, alo, bhi, blo;
	REAL err1, err2, err3;
	INEXACT REAL _i, _j;
	REAL _0;

	Two_Product(pa[0], pb[1], axby1, axby0);
	Two_Product(pb[0], pa[1], bxay1, bxay0);
	Two_Two_Diff(axby1, axby0, bxay1, bxay0, ab[3], ab[2], ab[1], ab[0]);

	Two_Product(pb[0], pc[1], bxcy1, bxcy0);
	Two_Product(pc[0], pb[1], cxby1, cxby0);
	Two_Two_Diff(bxcy1, bxcy0, cxby1, cxby0, bc[3], bc[2], bc[1], bc[0]);

	Two_Product(pc[0], pd[1], cxdy1, cxdy0);
	Two_Product(pd[0], pc[1], dxcy1, dxcy0);
	Two_Two_Diff(cxdy1, cxdy0, dxcy1, dxcy0, cd[3], cd[2], cd[1], cd[0]);

	Two_Product(pd[0], pa[1], dxay1, dxay0);
	Two_Product(pa[0], pd[1], axdy1, axdy0);
	Two_Two_Diff(dxay1, dxay0, axdy1, axdy0, da[3], da[2], da[1], da[0]);

	Two_Product(pa[0], pc[1], axcy1, axcy0);
	Two_Product(pc[0], pa[1], cxay1, cxay0);
	Two_Two_Diff(axcy1, axcy0, cxay1, cxay0, ac[3], ac[2], ac[1], ac[0]);

	Two_Product(pb[0], pd[1], bxdy1, bxdy0);
	Two_Product(pd[0], pb[1], dxby1, dxby0);
	Two_Two_Diff(bxdy1, bxdy0, dxby1, dxby0, bd[3], bd[2], bd[1], bd[0]);

	templen = fast_expansion_sum_zeroelim(4, cd, 4, da, temp8);
	cdalen = fast_expansion_sum_zeroelim(templen, temp8, 4, ac, cda);
	templen = fast_expansion_sum_zeroelim(4, da, 4, ab, temp8);
	dablen = fast_expansion_sum_zeroelim(templen, temp8, 4, bd, dab);
	for (i = 0; i < 4; i++) {
		bd[i] = -bd[i];
		ac[i] = -ac[i];
	}
	templen = fast_expansion_sum_zeroelim(4, ab, 4, bc, temp8);
	abclen = fast_expansion_sum_zeroelim(templen, temp8, 4, ac, abc);
	templen = fast_expansion_sum_zeroelim(4, bc, 4, cd, temp8);
	bcdlen = fast_expansion_sum_zeroelim(templen, temp8, 4, bd, bcd);

	alen = scale_expansion_zeroelim(bcdlen, bcd, pa[2], adet);
	blen = scale_expansion_zeroelim(cdalen, cda, -pb[2], bdet);
	clen = scale_expansion_zeroelim(dablen, dab, pc[2], cdet);
	dlen = scale_expansion_zeroelim(abclen, abc, -pd[2], ddet);

	ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, abdet);
	cdlen = fast_expansion_sum_zeroelim(clen, cdet, dlen, ddet, cddet);
	deterlen = fast_expansion_sum_zeroelim(ablen, abdet, cdlen, cddet, deter);

	return deter[deterlen - 1];
}

REAL orient3dslow(REAL *pa, REAL *pb, REAL *pc, REAL *pd)
{
	INEXACT REAL adx, ady, adz, bdx, bdy, bdz, cdx, cdy, cdz;
	REAL adxtail, adytail, adztail;
	REAL bdxtail, bdytail, bdztail;
	REAL cdxtail, cdytail, cdztail;
	REAL negate, negatetail;
	INEXACT REAL axby7, bxcy7, axcy7, bxay7, cxby7, cxay7;
	REAL axby[8], bxcy[8], axcy[8], bxay[8], cxby[8], cxay[8];
	REAL temp16[16], temp32[32], temp32t[32];
	int temp16len, temp32len, temp32tlen;
	REAL adet[64], bdet[64], cdet[64];
	int alen, blen, clen;
	REAL abdet[128];
	int ablen;
	REAL deter[192];
	int deterlen;

	INEXACT REAL bvirt;
	REAL avirt, bround, around;
	INEXACT REAL c;
	INEXACT REAL abig;
	REAL a0hi, a0lo, a1hi, a1lo, bhi, blo;
	REAL err1, err2, err3;
	INEXACT REAL _i, _j, _k, _l, _m, _n;
	REAL _0, _1, _2;

	Two_Diff(pa[0], pd[0], adx, adxtail);
	Two_Diff(pa[1], pd[1], ady, adytail);
	Two_Diff(pa[2], pd[2], adz, adztail);
	Two_Diff(pb[0], pd[0], bdx, bdxtail);
	Two_Diff(pb[1], pd[1], bdy, bdytail);
	Two_Diff(pb[2], pd[2], bdz, bdztail);
	Two_Diff(pc[0], pd[0], cdx, cdxtail);
	Two_Diff(pc[1], pd[1], cdy, cdytail);
	Two_Diff(pc[2], pd[2], cdz, cdztail);

	Two_Two_Product(adx, adxtail, bdy, bdytail,
		axby7, axby[6], axby[5], axby[4],
		axby[3], axby[2], axby[1], axby[0]);
	axby[7] = axby7;
	negate = -ady;
	negatetail = -adytail;
	Two_Two_Product(bdx, bdxtail, negate, negatetail,
		bxay7, bxay[6], bxay[5], bxay[4],
		bxay[3], bxay[2], bxay[1], bxay[0]);
	bxay[7] = bxay7;
	Two_Two_Product(bdx, bdxtail, cdy, cdytail,
		bxcy7, bxcy[6], bxcy[5], bxcy[4],
		bxcy[3], bxcy[2], bxcy[1], bxcy[0]);
	bxcy[7] = bxcy7;
	negate = -bdy;
	negatetail = -bdytail;
	Two_Two_Product(cdx, cdxtail, negate, negatetail,
		cxby7, cxby[6], cxby[5], cxby[4],
		cxby[3], cxby[2], cxby[1], cxby[0]);
	cxby[7] = cxby7;
	Two_Two_Product(cdx, cdxtail, ady, adytail,
		cxay7, cxay[6], cxay[5], cxay[4],
		cxay[3], cxay[2], cxay[1], cxay[0]);
	cxay[7] = cxay7;
	negate = -cdy;
	negatetail = -cdytail;
	Two_Two_Product(adx, adxtail, negate, negatetail,
		axcy7, axcy[6], axcy[5], axcy[4],
		axcy[3], axcy[2], axcy[1], axcy[0]);
	axcy[7] = axcy7;

	temp16len = fast_expansion_sum_zeroelim(8, bxcy, 8, cxby, temp16);
	temp32len = scale_expansion_zeroelim(temp16len, temp16, adz, temp32);
	temp32tlen = scale_expansion_zeroelim(temp16len, temp16, adztail, temp32t);
	alen = fast_expansion_sum_zeroelim(temp32len, temp32, temp32tlen, temp32t,
		adet);

	temp16len = fast_expansion_sum_zeroelim(8, cxay, 8, axcy, temp16);
	temp32len = scale_expansion_zeroelim(temp16len, temp16, bdz, temp32);
	temp32tlen = scale_expansion_zeroelim(temp16len, temp16, bdztail, temp32t);
	blen = fast_expansion_sum_zeroelim(temp32len, temp32, temp32tlen, temp32t,
		bdet);

	temp16len = fast_expansion_sum_zeroelim(8, axby, 8, bxay, temp16);
	temp32len = scale_expansion_zeroelim(temp16len, temp16, cdz, temp32);
	temp32tlen = scale_expansion_zeroelim(temp16len, temp16, cdztail, temp32t);
	clen = fast_expansion_sum_zeroelim(temp32len, temp32, temp32tlen, temp32t,
		cdet);

	ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, abdet);
	deterlen = fast_expansion_sum_zeroelim(ablen, abdet, clen, cdet, deter);

	return deter[deterlen - 1];
}

REAL orient3dadapt(REAL *pa, REAL *pb, REAL *pc, REAL *pd, REAL permanent)
{
	INEXACT REAL adx, bdx, cdx, ady, bdy, cdy, adz, bdz, cdz;
	REAL det, errbound;

	INEXACT REAL bdxcdy1, cdxbdy1, cdxady1, adxcdy1, adxbdy1, bdxady1;
	REAL bdxcdy0, cdxbdy0, cdxady0, adxcdy0, adxbdy0, bdxady0;
	REAL bc[4], ca[4], ab[4];
	INEXACT REAL bc3, ca3, ab3;
	REAL adet[8], bdet[8], cdet[8];
	int alen, blen, clen;
	REAL abdet[16];
	int ablen;
	REAL *finnow, *finother, *finswap;
	REAL fin1[192], fin2[192];
	int finlength;

	REAL adxtail, bdxtail, cdxtail;
	REAL adytail, bdytail, cdytail;
	REAL adztail, bdztail, cdztail;
	INEXACT REAL at_blarge, at_clarge;
	INEXACT REAL bt_clarge, bt_alarge;
	INEXACT REAL ct_alarge, ct_blarge;
	REAL at_b[4], at_c[4], bt_c[4], bt_a[4], ct_a[4], ct_b[4];
	int at_blen, at_clen, bt_clen, bt_alen, ct_alen, ct_blen;
	INEXACT REAL bdxt_cdy1, cdxt_bdy1, cdxt_ady1;
	INEXACT REAL adxt_cdy1, adxt_bdy1, bdxt_ady1;
	REAL bdxt_cdy0, cdxt_bdy0, cdxt_ady0;
	REAL adxt_cdy0, adxt_bdy0, bdxt_ady0;
	INEXACT REAL bdyt_cdx1, cdyt_bdx1, cdyt_adx1;
	INEXACT REAL adyt_cdx1, adyt_bdx1, bdyt_adx1;
	REAL bdyt_cdx0, cdyt_bdx0, cdyt_adx0;
	REAL adyt_cdx0, adyt_bdx0, bdyt_adx0;
	REAL bct[8], cat[8], abt[8];
	int bctlen, catlen, abtlen;
	INEXACT REAL bdxt_cdyt1, cdxt_bdyt1, cdxt_adyt1;
	INEXACT REAL adxt_cdyt1, adxt_bdyt1, bdxt_adyt1;
	REAL bdxt_cdyt0, cdxt_bdyt0, cdxt_adyt0;
	REAL adxt_cdyt0, adxt_bdyt0, bdxt_adyt0;
	REAL u[4], v[12], w[16];
	INEXACT REAL u3;
	int vlength, wlength;
	REAL negate;

	INEXACT REAL bvirt;
	REAL avirt, bround, around;
	INEXACT REAL c;
	INEXACT REAL abig;
	REAL ahi, alo, bhi, blo;
	REAL err1, err2, err3;
	INEXACT REAL _i, _j, _k;
	REAL _0;

	adx = (REAL)(pa[0] - pd[0]);
	bdx = (REAL)(pb[0] - pd[0]);
	cdx = (REAL)(pc[0] - pd[0]);
	ady = (REAL)(pa[1] - pd[1]);
	bdy = (REAL)(pb[1] - pd[1]);
	cdy = (REAL)(pc[1] - pd[1]);
	adz = (REAL)(pa[2] - pd[2]);
	bdz = (REAL)(pb[2] - pd[2]);
	cdz = (REAL)(pc[2] - pd[2]);

	Two_Product(bdx, cdy, bdxcdy1, bdxcdy0);
	Two_Product(cdx, bdy, cdxbdy1, cdxbdy0);
	Two_Two_Diff(bdxcdy1, bdxcdy0, cdxbdy1, cdxbdy0, bc3, bc[2], bc[1], bc[0]);
	bc[3] = bc3;
	alen = scale_expansion_zeroelim(4, bc, adz, adet);

	Two_Product(cdx, ady, cdxady1, cdxady0);
	Two_Product(adx, cdy, adxcdy1, adxcdy0);
	Two_Two_Diff(cdxady1, cdxady0, adxcdy1, adxcdy0, ca3, ca[2], ca[1], ca[0]);
	ca[3] = ca3;
	blen = scale_expansion_zeroelim(4, ca, bdz, bdet);

	Two_Product(adx, bdy, adxbdy1, adxbdy0);
	Two_Product(bdx, ady, bdxady1, bdxady0);
	Two_Two_Diff(adxbdy1, adxbdy0, bdxady1, bdxady0, ab3, ab[2], ab[1], ab[0]);
	ab[3] = ab3;
	clen = scale_expansion_zeroelim(4, ab, cdz, cdet);

	ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, abdet);
	finlength = fast_expansion_sum_zeroelim(ablen, abdet, clen, cdet, fin1);

	det = estimate(finlength, fin1);
	errbound = o3derrboundB * permanent;
	if ((det >= errbound) || (-det >= errbound)) {
		return det;
	}

	Two_Diff_Tail(pa[0], pd[0], adx, adxtail);
	Two_Diff_Tail(pb[0], pd[0], bdx, bdxtail);
	Two_Diff_Tail(pc[0], pd[0], cdx, cdxtail);
	Two_Diff_Tail(pa[1], pd[1], ady, adytail);
	Two_Diff_Tail(pb[1], pd[1], bdy, bdytail);
	Two_Diff_Tail(pc[1], pd[1], cdy, cdytail);
	Two_Diff_Tail(pa[2], pd[2], adz, adztail);
	Two_Diff_Tail(pb[2], pd[2], bdz, bdztail);
	Two_Diff_Tail(pc[2], pd[2], cdz, cdztail);

	if ((adxtail == 0.0) && (bdxtail == 0.0) && (cdxtail == 0.0)
		&& (adytail == 0.0) && (bdytail == 0.0) && (cdytail == 0.0)
		&& (adztail == 0.0) && (bdztail == 0.0) && (cdztail == 0.0)) {
		return det;
	}

	errbound = o3derrboundC * permanent + resulterrbound * Absolute(det);
	det += (adz * ((bdx * cdytail + cdy * bdxtail)
		- (bdy * cdxtail + cdx * bdytail))
		+ adztail * (bdx * cdy - bdy * cdx))
		+ (bdz * ((cdx * adytail + ady * cdxtail)
			- (cdy * adxtail + adx * cdytail))
			+ bdztail * (cdx * ady - cdy * adx))
		+ (cdz * ((adx * bdytail + bdy * adxtail)
			- (ady * bdxtail + bdx * adytail))
			+ cdztail * (adx * bdy - ady * bdx));
	if ((det >= errbound) || (-det >= errbound)) {
		return det;
	}

	finnow = fin1;
	finother = fin2;

	if (adxtail == 0.0) {
		if (adytail == 0.0) {
			at_b[0] = 0.0;
			at_blen = 1;
			at_c[0] = 0.0;
			at_clen = 1;
		}
		else {
			negate = -adytail;
			Two_Product(negate, bdx, at_blarge, at_b[0]);
			at_b[1] = at_blarge;
			at_blen = 2;
			Two_Product(adytail, cdx, at_clarge, at_c[0]);
			at_c[1] = at_clarge;
			at_clen = 2;
		}
	}
	else {
		if (adytail == 0.0) {
			Two_Product(adxtail, bdy, at_blarge, at_b[0]);
			at_b[1] = at_blarge;
			at_blen = 2;
			negate = -adxtail;
			Two_Product(negate, cdy, at_clarge, at_c[0]);
			at_c[1] = at_clarge;
			at_clen = 2;
		}
		else {
			Two_Product(adxtail, bdy, adxt_bdy1, adxt_bdy0);
			Two_Product(adytail, bdx, adyt_bdx1, adyt_bdx0);
			Two_Two_Diff(adxt_bdy1, adxt_bdy0, adyt_bdx1, adyt_bdx0,
				at_blarge, at_b[2], at_b[1], at_b[0]);
			at_b[3] = at_blarge;
			at_blen = 4;
			Two_Product(adytail, cdx, adyt_cdx1, adyt_cdx0);
			Two_Product(adxtail, cdy, adxt_cdy1, adxt_cdy0);
			Two_Two_Diff(adyt_cdx1, adyt_cdx0, adxt_cdy1, adxt_cdy0,
				at_clarge, at_c[2], at_c[1], at_c[0]);
			at_c[3] = at_clarge;
			at_clen = 4;
		}
	}
	if (bdxtail == 0.0) {
		if (bdytail == 0.0) {
			bt_c[0] = 0.0;
			bt_clen = 1;
			bt_a[0] = 0.0;
			bt_alen = 1;
		}
		else {
			negate = -bdytail;
			Two_Product(negate, cdx, bt_clarge, bt_c[0]);
			bt_c[1] = bt_clarge;
			bt_clen = 2;
			Two_Product(bdytail, adx, bt_alarge, bt_a[0]);
			bt_a[1] = bt_alarge;
			bt_alen = 2;
		}
	}
	else {
		if (bdytail == 0.0) {
			Two_Product(bdxtail, cdy, bt_clarge, bt_c[0]);
			bt_c[1] = bt_clarge;
			bt_clen = 2;
			negate = -bdxtail;
			Two_Product(negate, ady, bt_alarge, bt_a[0]);
			bt_a[1] = bt_alarge;
			bt_alen = 2;
		}
		else {
			Two_Product(bdxtail, cdy, bdxt_cdy1, bdxt_cdy0);
			Two_Product(bdytail, cdx, bdyt_cdx1, bdyt_cdx0);
			Two_Two_Diff(bdxt_cdy1, bdxt_cdy0, bdyt_cdx1, bdyt_cdx0,
				bt_clarge, bt_c[2], bt_c[1], bt_c[0]);
			bt_c[3] = bt_clarge;
			bt_clen = 4;
			Two_Product(bdytail, adx, bdyt_adx1, bdyt_adx0);
			Two_Product(bdxtail, ady, bdxt_ady1, bdxt_ady0);
			Two_Two_Diff(bdyt_adx1, bdyt_adx0, bdxt_ady1, bdxt_ady0,
				bt_alarge, bt_a[2], bt_a[1], bt_a[0]);
			bt_a[3] = bt_alarge;
			bt_alen = 4;
		}
	}
	if (cdxtail == 0.0) {
		if (cdytail == 0.0) {
			ct_a[0] = 0.0;
			ct_alen = 1;
			ct_b[0] = 0.0;
			ct_blen = 1;
		}
		else {
			negate = -cdytail;
			Two_Product(negate, adx, ct_alarge, ct_a[0]);
			ct_a[1] = ct_alarge;
			ct_alen = 2;
			Two_Product(cdytail, bdx, ct_blarge, ct_b[0]);
			ct_b[1] = ct_blarge;
			ct_blen = 2;
		}
	}
	else {
		if (cdytail == 0.0) {
			Two_Product(cdxtail, ady, ct_alarge, ct_a[0]);
			ct_a[1] = ct_alarge;
			ct_alen = 2;
			negate = -cdxtail;
			Two_Product(negate, bdy, ct_blarge, ct_b[0]);
			ct_b[1] = ct_blarge;
			ct_blen = 2;
		}
		else {
			Two_Product(cdxtail, ady, cdxt_ady1, cdxt_ady0);
			Two_Product(cdytail, adx, cdyt_adx1, cdyt_adx0);
			Two_Two_Diff(cdxt_ady1, cdxt_ady0, cdyt_adx1, cdyt_adx0,
				ct_alarge, ct_a[2], ct_a[1], ct_a[0]);
			ct_a[3] = ct_alarge;
			ct_alen = 4;
			Two_Product(cdytail, bdx, cdyt_bdx1, cdyt_bdx0);
			Two_Product(cdxtail, bdy, cdxt_bdy1, cdxt_bdy0);
			Two_Two_Diff(cdyt_bdx1, cdyt_bdx0, cdxt_bdy1, cdxt_bdy0,
				ct_blarge, ct_b[2], ct_b[1], ct_b[0]);
			ct_b[3] = ct_blarge;
			ct_blen = 4;
		}
	}

	bctlen = fast_expansion_sum_zeroelim(bt_clen, bt_c, ct_blen, ct_b, bct);
	wlength = scale_expansion_zeroelim(bctlen, bct, adz, w);
	finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w,
		finother);
	finswap = finnow; finnow = finother; finother = finswap;

	catlen = fast_expansion_sum_zeroelim(ct_alen, ct_a, at_clen, at_c, cat);
	wlength = scale_expansion_zeroelim(catlen, cat, bdz, w);
	finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w,
		finother);
	finswap = finnow; finnow = finother; finother = finswap;

	abtlen = fast_expansion_sum_zeroelim(at_blen, at_b, bt_alen, bt_a, abt);
	wlength = scale_expansion_zeroelim(abtlen, abt, cdz, w);
	finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w,
		finother);
	finswap = finnow; finnow = finother; finother = finswap;

	if (adztail != 0.0) {
		vlength = scale_expansion_zeroelim(4, bc, adztail, v);
		finlength = fast_expansion_sum_zeroelim(finlength, finnow, vlength, v,
			finother);
		finswap = finnow; finnow = finother; finother = finswap;
	}
	if (bdztail != 0.0) {
		vlength = scale_expansion_zeroelim(4, ca, bdztail, v);
		finlength = fast_expansion_sum_zeroelim(finlength, finnow, vlength, v,
			finother);
		finswap = finnow; finnow = finother; finother = finswap;
	}
	if (cdztail != 0.0) {
		vlength = scale_expansion_zeroelim(4, ab, cdztail, v);
		finlength = fast_expansion_sum_zeroelim(finlength, finnow, vlength, v,
			finother);
		finswap = finnow; finnow = finother; finother = finswap;
	}

	if (adxtail != 0.0) {
		if (bdytail != 0.0) {
			Two_Product(adxtail, bdytail, adxt_bdyt1, adxt_bdyt0);
			Two_One_Product(adxt_bdyt1, adxt_bdyt0, cdz, u3, u[2], u[1], u[0]);
			u[3] = u3;
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
				finother);
			finswap = finnow; finnow = finother; finother = finswap;
			if (cdztail != 0.0) {
				Two_One_Product(adxt_bdyt1, adxt_bdyt0, cdztail, u3, u[2], u[1], u[0]);
				u[3] = u3;
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
					finother);
				finswap = finnow; finnow = finother; finother = finswap;
			}
		}
		if (cdytail != 0.0) {
			negate = -adxtail;
			Two_Product(negate, cdytail, adxt_cdyt1, adxt_cdyt0);
			Two_One_Product(adxt_cdyt1, adxt_cdyt0, bdz, u3, u[2], u[1], u[0]);
			u[3] = u3;
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
				finother);
			finswap = finnow; finnow = finother; finother = finswap;
			if (bdztail != 0.0) {
				Two_One_Product(adxt_cdyt1, adxt_cdyt0, bdztail, u3, u[2], u[1], u[0]);
				u[3] = u3;
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
					finother);
				finswap = finnow; finnow = finother; finother = finswap;
			}
		}
	}
	if (bdxtail != 0.0) {
		if (cdytail != 0.0) {
			Two_Product(bdxtail, cdytail, bdxt_cdyt1, bdxt_cdyt0);
			Two_One_Product(bdxt_cdyt1, bdxt_cdyt0, adz, u3, u[2], u[1], u[0]);
			u[3] = u3;
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
				finother);
			finswap = finnow; finnow = finother; finother = finswap;
			if (adztail != 0.0) {
				Two_One_Product(bdxt_cdyt1, bdxt_cdyt0, adztail, u3, u[2], u[1], u[0]);
				u[3] = u3;
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
					finother);
				finswap = finnow; finnow = finother; finother = finswap;
			}
		}
		if (adytail != 0.0) {
			negate = -bdxtail;
			Two_Product(negate, adytail, bdxt_adyt1, bdxt_adyt0);
			Two_One_Product(bdxt_adyt1, bdxt_adyt0, cdz, u3, u[2], u[1], u[0]);
			u[3] = u3;
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
				finother);
			finswap = finnow; finnow = finother; finother = finswap;
			if (cdztail != 0.0) {
				Two_One_Product(bdxt_adyt1, bdxt_adyt0, cdztail, u3, u[2], u[1], u[0]);
				u[3] = u3;
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
					finother);
				finswap = finnow; finnow = finother; finother = finswap;
			}
		}
	}
	if (cdxtail != 0.0) {
		if (adytail != 0.0) {
			Two_Product(cdxtail, adytail, cdxt_adyt1, cdxt_adyt0);
			Two_One_Product(cdxt_adyt1, cdxt_adyt0, bdz, u3, u[2], u[1], u[0]);
			u[3] = u3;
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
				finother);
			finswap = finnow; finnow = finother; finother = finswap;
			if (bdztail != 0.0) {
				Two_One_Product(cdxt_adyt1, cdxt_adyt0, bdztail, u3, u[2], u[1], u[0]);
				u[3] = u3;
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
					finother);
				finswap = finnow; finnow = finother; finother = finswap;
			}
		}
		if (bdytail != 0.0) {
			negate = -cdxtail;
			Two_Product(negate, bdytail, cdxt_bdyt1, cdxt_bdyt0);
			Two_One_Product(cdxt_bdyt1, cdxt_bdyt0, adz, u3, u[2], u[1], u[0]);
			u[3] = u3;
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
				finother);
			finswap = finnow; finnow = finother; finother = finswap;
			if (adztail != 0.0) {
				Two_One_Product(cdxt_bdyt1, cdxt_bdyt0, adztail, u3, u[2], u[1], u[0]);
				u[3] = u3;
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
					finother);
				finswap = finnow; finnow = finother; finother = finswap;
			}
		}
	}

	if (adztail != 0.0) {
		wlength = scale_expansion_zeroelim(bctlen, bct, adztail, w);
		finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w,
			finother);
		finswap = finnow; finnow = finother; finother = finswap;
	}
	if (bdztail != 0.0) {
		wlength = scale_expansion_zeroelim(catlen, cat, bdztail, w);
		finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w,
			finother);
		finswap = finnow; finnow = finother; finother = finswap;
	}
	if (cdztail != 0.0) {
		wlength = scale_expansion_zeroelim(abtlen, abt, cdztail, w);
		finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w,
			finother);
		finswap = finnow; finnow = finother; finother = finswap;
	}

	return finnow[finlength - 1];
}

inline REAL orient3d(REAL *&pa, REAL *&pb, REAL *&pc, REAL *&pd)
{
	REAL adx = pa[0] - pd[0];
	REAL bdx = pb[0] - pd[0];
	REAL cdx = pc[0] - pd[0];
	REAL ady = pa[1] - pd[1];
	REAL bdy = pb[1] - pd[1];
	REAL cdy = pc[1] - pd[1];
	REAL adz = pa[2] - pd[2];
	REAL bdz = pb[2] - pd[2];
	REAL cdz = pc[2] - pd[2];


	REAL det= adx * (bdy * cdz - bdz * cdy)
		+ bdx * (cdy * adz - cdz * ady)
		+ cdx * (ady * bdz - adz * bdy);
	constexpr static double bound = 1e-12;//这个值经过测试
	if (det > bound&& -det<bound )
		return det;

	REAL bdxcdy = bdx * cdy;
	REAL cdxbdy = cdx * bdy;

	REAL cdxady = cdx * ady;
	REAL adxcdy = adx * cdy;

	REAL adxbdy = adx * bdy;
	REAL bdxady = bdx * ady;


	REAL permanent = (Absolute(bdxcdy) + Absolute(cdxbdy)) * Absolute(adz)
		+ (Absolute(cdxady) + Absolute(adxcdy)) * Absolute(bdz)
		+ (Absolute(adxbdy) + Absolute(bdxady)) * Absolute(cdz);
	REAL  errbound = o3derrboundA * permanent;
	if ((det > errbound) || (-det > errbound)) {
		return det;
	}

	return orient3dadapt(pa, pb, pc, pd, permanent);
}


/******************************************************************************/
/*                                                                           **/
/*  incirclefast()   Approximate 2D incircle test.  Nonrobust.               **/
/*  incircleexact()   Exact 2D incircle test.  Robust.                       **/
/*  incircleslow()   Another exact 2D incircle test.  Robust.                **/
/*  incircle()   Adaptive exact 2D incircle test.  Robust.                   **/
/*                                                                           **/
/*               Return a positive value if the point pd lies inside the     **/
/*               circle passing through pa, pb, and pc; a negative value if  **/
/*               it lies outside; and zero if the four points are cocircular.**/
/*               The points pa, pb, and pc must be in counterclockwise       **/
/*               order, or the sign of the result will be reversed.          **/
/*                                                                           **/
/*  Only the first and last routine should be used; the middle two are for   **/
/*  timings.                                                                 **/
/*                                                                           **/
/*  The last three use exact arithmetic to ensure a correct answer.  The     **/
/*  result returned is the determinant of a matrix.  In incircle() only,     **/
/*  this determinant is computed adaptively, in the sense that exact         **/
/*  arithmetic is used only to the degree it is needed to ensure that the    **/
/*  returned value has the correct sign.  Hence, incircle() is usually quite **/
/*  fast, but will run more slowly when the input points are cocircular or   **/
/*  nearly so.                                                               **/
/*                                                                           **/
/******************************************************************************/

REAL incirclefast(REAL *pa, REAL *pb, REAL *pc, REAL *pd)
{
	REAL adx, ady, bdx, bdy, cdx, cdy;
	REAL abdet, bcdet, cadet;
	REAL alift, blift, clift;

	adx = pa[0] - pd[0];
	ady = pa[1] - pd[1];
	bdx = pb[0] - pd[0];
	bdy = pb[1] - pd[1];
	cdx = pc[0] - pd[0];
	cdy = pc[1] - pd[1];

	abdet = adx * bdy - bdx * ady;
	bcdet = bdx * cdy - cdx * bdy;
	cadet = cdx * ady - adx * cdy;
	alift = adx * adx + ady * ady;
	blift = bdx * bdx + bdy * bdy;
	clift = cdx * cdx + cdy * cdy;

	return alift * bcdet + blift * cadet + clift * abdet;
}

REAL incircleexact(REAL *pa, REAL *pb, REAL *pc, REAL *pd)
{
	INEXACT REAL axby1, bxcy1, cxdy1, dxay1, axcy1, bxdy1;
	INEXACT REAL bxay1, cxby1, dxcy1, axdy1, cxay1, dxby1;
	REAL axby0, bxcy0, cxdy0, dxay0, axcy0, bxdy0;
	REAL bxay0, cxby0, dxcy0, axdy0, cxay0, dxby0;
	REAL ab[4], bc[4], cd[4], da[4], ac[4], bd[4];
	REAL temp8[8];
	int templen;
	REAL abc[12], bcd[12], cda[12], dab[12];
	int abclen, bcdlen, cdalen, dablen;
	REAL det24x[24], det24y[24], det48x[48], det48y[48];
	int xlen, ylen;
	REAL adet[96], bdet[96], cdet[96], ddet[96];
	int alen, blen, clen, dlen;
	REAL abdet[192], cddet[192];
	int ablen, cdlen;
	REAL deter[384];
	int deterlen;
	int i;

	INEXACT REAL bvirt;
	REAL avirt, bround, around;
	INEXACT REAL c;
	INEXACT REAL abig;
	REAL ahi, alo, bhi, blo;
	REAL err1, err2, err3;
	INEXACT REAL _i, _j;
	REAL _0;

	Two_Product(pa[0], pb[1], axby1, axby0);
	Two_Product(pb[0], pa[1], bxay1, bxay0);
	Two_Two_Diff(axby1, axby0, bxay1, bxay0, ab[3], ab[2], ab[1], ab[0]);

	Two_Product(pb[0], pc[1], bxcy1, bxcy0);
	Two_Product(pc[0], pb[1], cxby1, cxby0);
	Two_Two_Diff(bxcy1, bxcy0, cxby1, cxby0, bc[3], bc[2], bc[1], bc[0]);

	Two_Product(pc[0], pd[1], cxdy1, cxdy0);
	Two_Product(pd[0], pc[1], dxcy1, dxcy0);
	Two_Two_Diff(cxdy1, cxdy0, dxcy1, dxcy0, cd[3], cd[2], cd[1], cd[0]);

	Two_Product(pd[0], pa[1], dxay1, dxay0);
	Two_Product(pa[0], pd[1], axdy1, axdy0);
	Two_Two_Diff(dxay1, dxay0, axdy1, axdy0, da[3], da[2], da[1], da[0]);

	Two_Product(pa[0], pc[1], axcy1, axcy0);
	Two_Product(pc[0], pa[1], cxay1, cxay0);
	Two_Two_Diff(axcy1, axcy0, cxay1, cxay0, ac[3], ac[2], ac[1], ac[0]);

	Two_Product(pb[0], pd[1], bxdy1, bxdy0);
	Two_Product(pd[0], pb[1], dxby1, dxby0);
	Two_Two_Diff(bxdy1, bxdy0, dxby1, dxby0, bd[3], bd[2], bd[1], bd[0]);

	templen = fast_expansion_sum_zeroelim(4, cd, 4, da, temp8);
	cdalen = fast_expansion_sum_zeroelim(templen, temp8, 4, ac, cda);
	templen = fast_expansion_sum_zeroelim(4, da, 4, ab, temp8);
	dablen = fast_expansion_sum_zeroelim(templen, temp8, 4, bd, dab);
	for (i = 0; i < 4; i++) {
		bd[i] = -bd[i];
		ac[i] = -ac[i];
	}
	templen = fast_expansion_sum_zeroelim(4, ab, 4, bc, temp8);
	abclen = fast_expansion_sum_zeroelim(templen, temp8, 4, ac, abc);
	templen = fast_expansion_sum_zeroelim(4, bc, 4, cd, temp8);
	bcdlen = fast_expansion_sum_zeroelim(templen, temp8, 4, bd, bcd);

	xlen = scale_expansion_zeroelim(bcdlen, bcd, pa[0], det24x);
	xlen = scale_expansion_zeroelim(xlen, det24x, pa[0], det48x);
	ylen = scale_expansion_zeroelim(bcdlen, bcd, pa[1], det24y);
	ylen = scale_expansion_zeroelim(ylen, det24y, pa[1], det48y);
	alen = fast_expansion_sum_zeroelim(xlen, det48x, ylen, det48y, adet);

	xlen = scale_expansion_zeroelim(cdalen, cda, pb[0], det24x);
	xlen = scale_expansion_zeroelim(xlen, det24x, -pb[0], det48x);
	ylen = scale_expansion_zeroelim(cdalen, cda, pb[1], det24y);
	ylen = scale_expansion_zeroelim(ylen, det24y, -pb[1], det48y);
	blen = fast_expansion_sum_zeroelim(xlen, det48x, ylen, det48y, bdet);

	xlen = scale_expansion_zeroelim(dablen, dab, pc[0], det24x);
	xlen = scale_expansion_zeroelim(xlen, det24x, pc[0], det48x);
	ylen = scale_expansion_zeroelim(dablen, dab, pc[1], det24y);
	ylen = scale_expansion_zeroelim(ylen, det24y, pc[1], det48y);
	clen = fast_expansion_sum_zeroelim(xlen, det48x, ylen, det48y, cdet);

	xlen = scale_expansion_zeroelim(abclen, abc, pd[0], det24x);
	xlen = scale_expansion_zeroelim(xlen, det24x, -pd[0], det48x);
	ylen = scale_expansion_zeroelim(abclen, abc, pd[1], det24y);
	ylen = scale_expansion_zeroelim(ylen, det24y, -pd[1], det48y);
	dlen = fast_expansion_sum_zeroelim(xlen, det48x, ylen, det48y, ddet);

	ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, abdet);
	cdlen = fast_expansion_sum_zeroelim(clen, cdet, dlen, ddet, cddet);
	deterlen = fast_expansion_sum_zeroelim(ablen, abdet, cdlen, cddet, deter);

	return deter[deterlen - 1];
}

REAL incircleslow(REAL *pa, REAL *pb, REAL *pc, REAL *pd)
{
	INEXACT REAL adx, bdx, cdx, ady, bdy, cdy;
	REAL adxtail, bdxtail, cdxtail;
	REAL adytail, bdytail, cdytail;
	REAL negate, negatetail;
	INEXACT REAL axby7, bxcy7, axcy7, bxay7, cxby7, cxay7;
	REAL axby[8], bxcy[8], axcy[8], bxay[8], cxby[8], cxay[8];
	REAL temp16[16];
	int temp16len;
	REAL detx[32], detxx[64], detxt[32], detxxt[64], detxtxt[64];
	int xlen, xxlen, xtlen, xxtlen, xtxtlen;
	REAL x1[128], x2[192];
	int x1len, x2len;
	REAL dety[32], detyy[64], detyt[32], detyyt[64], detytyt[64];
	int ylen, yylen, ytlen, yytlen, ytytlen;
	REAL y1[128], y2[192];
	int y1len, y2len;
	REAL adet[384], bdet[384], cdet[384], abdet[768], deter[1152];
	int alen, blen, clen, ablen, deterlen;
	int i;

	INEXACT REAL bvirt;
	REAL avirt, bround, around;
	INEXACT REAL c;
	INEXACT REAL abig;
	REAL a0hi, a0lo, a1hi, a1lo, bhi, blo;
	REAL err1, err2, err3;
	INEXACT REAL _i, _j, _k, _l, _m, _n;
	REAL _0, _1, _2;

	Two_Diff(pa[0], pd[0], adx, adxtail);
	Two_Diff(pa[1], pd[1], ady, adytail);
	Two_Diff(pb[0], pd[0], bdx, bdxtail);
	Two_Diff(pb[1], pd[1], bdy, bdytail);
	Two_Diff(pc[0], pd[0], cdx, cdxtail);
	Two_Diff(pc[1], pd[1], cdy, cdytail);

	Two_Two_Product(adx, adxtail, bdy, bdytail,
		axby7, axby[6], axby[5], axby[4],
		axby[3], axby[2], axby[1], axby[0]);
	axby[7] = axby7;
	negate = -ady;
	negatetail = -adytail;
	Two_Two_Product(bdx, bdxtail, negate, negatetail,
		bxay7, bxay[6], bxay[5], bxay[4],
		bxay[3], bxay[2], bxay[1], bxay[0]);
	bxay[7] = bxay7;
	Two_Two_Product(bdx, bdxtail, cdy, cdytail,
		bxcy7, bxcy[6], bxcy[5], bxcy[4],
		bxcy[3], bxcy[2], bxcy[1], bxcy[0]);
	bxcy[7] = bxcy7;
	negate = -bdy;
	negatetail = -bdytail;
	Two_Two_Product(cdx, cdxtail, negate, negatetail,
		cxby7, cxby[6], cxby[5], cxby[4],
		cxby[3], cxby[2], cxby[1], cxby[0]);
	cxby[7] = cxby7;
	Two_Two_Product(cdx, cdxtail, ady, adytail,
		cxay7, cxay[6], cxay[5], cxay[4],
		cxay[3], cxay[2], cxay[1], cxay[0]);
	cxay[7] = cxay7;
	negate = -cdy;
	negatetail = -cdytail;
	Two_Two_Product(adx, adxtail, negate, negatetail,
		axcy7, axcy[6], axcy[5], axcy[4],
		axcy[3], axcy[2], axcy[1], axcy[0]);
	axcy[7] = axcy7;


	temp16len = fast_expansion_sum_zeroelim(8, bxcy, 8, cxby, temp16);

	xlen = scale_expansion_zeroelim(temp16len, temp16, adx, detx);
	xxlen = scale_expansion_zeroelim(xlen, detx, adx, detxx);
	xtlen = scale_expansion_zeroelim(temp16len, temp16, adxtail, detxt);
	xxtlen = scale_expansion_zeroelim(xtlen, detxt, adx, detxxt);
	for (i = 0; i < xxtlen; i++) {
		detxxt[i] *= 2.0;
	}
	xtxtlen = scale_expansion_zeroelim(xtlen, detxt, adxtail, detxtxt);
	x1len = fast_expansion_sum_zeroelim(xxlen, detxx, xxtlen, detxxt, x1);
	x2len = fast_expansion_sum_zeroelim(x1len, x1, xtxtlen, detxtxt, x2);

	ylen = scale_expansion_zeroelim(temp16len, temp16, ady, dety);
	yylen = scale_expansion_zeroelim(ylen, dety, ady, detyy);
	ytlen = scale_expansion_zeroelim(temp16len, temp16, adytail, detyt);
	yytlen = scale_expansion_zeroelim(ytlen, detyt, ady, detyyt);
	for (i = 0; i < yytlen; i++) {
		detyyt[i] *= 2.0;
	}
	ytytlen = scale_expansion_zeroelim(ytlen, detyt, adytail, detytyt);
	y1len = fast_expansion_sum_zeroelim(yylen, detyy, yytlen, detyyt, y1);
	y2len = fast_expansion_sum_zeroelim(y1len, y1, ytytlen, detytyt, y2);

	alen = fast_expansion_sum_zeroelim(x2len, x2, y2len, y2, adet);


	temp16len = fast_expansion_sum_zeroelim(8, cxay, 8, axcy, temp16);

	xlen = scale_expansion_zeroelim(temp16len, temp16, bdx, detx);
	xxlen = scale_expansion_zeroelim(xlen, detx, bdx, detxx);
	xtlen = scale_expansion_zeroelim(temp16len, temp16, bdxtail, detxt);
	xxtlen = scale_expansion_zeroelim(xtlen, detxt, bdx, detxxt);
	for (i = 0; i < xxtlen; i++) {
		detxxt[i] *= 2.0;
	}
	xtxtlen = scale_expansion_zeroelim(xtlen, detxt, bdxtail, detxtxt);
	x1len = fast_expansion_sum_zeroelim(xxlen, detxx, xxtlen, detxxt, x1);
	x2len = fast_expansion_sum_zeroelim(x1len, x1, xtxtlen, detxtxt, x2);

	ylen = scale_expansion_zeroelim(temp16len, temp16, bdy, dety);
	yylen = scale_expansion_zeroelim(ylen, dety, bdy, detyy);
	ytlen = scale_expansion_zeroelim(temp16len, temp16, bdytail, detyt);
	yytlen = scale_expansion_zeroelim(ytlen, detyt, bdy, detyyt);
	for (i = 0; i < yytlen; i++) {
		detyyt[i] *= 2.0;
	}
	ytytlen = scale_expansion_zeroelim(ytlen, detyt, bdytail, detytyt);
	y1len = fast_expansion_sum_zeroelim(yylen, detyy, yytlen, detyyt, y1);
	y2len = fast_expansion_sum_zeroelim(y1len, y1, ytytlen, detytyt, y2);

	blen = fast_expansion_sum_zeroelim(x2len, x2, y2len, y2, bdet);


	temp16len = fast_expansion_sum_zeroelim(8, axby, 8, bxay, temp16);

	xlen = scale_expansion_zeroelim(temp16len, temp16, cdx, detx);
	xxlen = scale_expansion_zeroelim(xlen, detx, cdx, detxx);
	xtlen = scale_expansion_zeroelim(temp16len, temp16, cdxtail, detxt);
	xxtlen = scale_expansion_zeroelim(xtlen, detxt, cdx, detxxt);
	for (i = 0; i < xxtlen; i++) {
		detxxt[i] *= 2.0;
	}
	xtxtlen = scale_expansion_zeroelim(xtlen, detxt, cdxtail, detxtxt);
	x1len = fast_expansion_sum_zeroelim(xxlen, detxx, xxtlen, detxxt, x1);
	x2len = fast_expansion_sum_zeroelim(x1len, x1, xtxtlen, detxtxt, x2);

	ylen = scale_expansion_zeroelim(temp16len, temp16, cdy, dety);
	yylen = scale_expansion_zeroelim(ylen, dety, cdy, detyy);
	ytlen = scale_expansion_zeroelim(temp16len, temp16, cdytail, detyt);
	yytlen = scale_expansion_zeroelim(ytlen, detyt, cdy, detyyt);
	for (i = 0; i < yytlen; i++) {
		detyyt[i] *= 2.0;
	}
	ytytlen = scale_expansion_zeroelim(ytlen, detyt, cdytail, detytyt);
	y1len = fast_expansion_sum_zeroelim(yylen, detyy, yytlen, detyyt, y1);
	y2len = fast_expansion_sum_zeroelim(y1len, y1, ytytlen, detytyt, y2);

	clen = fast_expansion_sum_zeroelim(x2len, x2, y2len, y2, cdet);

	ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, abdet);
	deterlen = fast_expansion_sum_zeroelim(ablen, abdet, clen, cdet, deter);

	return deter[deterlen - 1];
}

REAL incircleadapt(REAL *pa, REAL *pb, REAL *pc, REAL *pd, REAL permanent)
{
	INEXACT REAL adx, bdx, cdx, ady, bdy, cdy;
	REAL det, errbound;

	INEXACT REAL bdxcdy1, cdxbdy1, cdxady1, adxcdy1, adxbdy1, bdxady1;
	REAL bdxcdy0, cdxbdy0, cdxady0, adxcdy0, adxbdy0, bdxady0;
	REAL bc[4], ca[4], ab[4];
	INEXACT REAL bc3, ca3, ab3;
	REAL axbc[8], axxbc[16], aybc[8], ayybc[16], adet[32];
	int axbclen, axxbclen, aybclen, ayybclen, alen;
	REAL bxca[8], bxxca[16], byca[8], byyca[16], bdet[32];
	int bxcalen, bxxcalen, bycalen, byycalen, blen;
	REAL cxab[8], cxxab[16], cyab[8], cyyab[16], cdet[32];
	int cxablen, cxxablen, cyablen, cyyablen, clen;
	REAL abdet[64];
	int ablen;
	REAL fin1[1152], fin2[1152];
	REAL *finnow, *finother, *finswap;
	int finlength;

	REAL adxtail, bdxtail, cdxtail, adytail, bdytail, cdytail;
	INEXACT REAL adxadx1, adyady1, bdxbdx1, bdybdy1, cdxcdx1, cdycdy1;
	REAL adxadx0, adyady0, bdxbdx0, bdybdy0, cdxcdx0, cdycdy0;
	REAL aa[4], bb[4], cc[4];
	INEXACT REAL aa3, bb3, cc3;
	INEXACT REAL ti1, tj1;
	REAL ti0, tj0;
	REAL u[4], v[4];
	INEXACT REAL u3, v3;
	REAL temp8[8], temp16a[16], temp16b[16], temp16c[16];
	REAL temp32a[32], temp32b[32], temp48[48], temp64[64];
	int temp8len, temp16alen, temp16blen, temp16clen;
	int temp32alen, temp32blen, temp48len, temp64len;
	REAL axtbb[8], axtcc[8], aytbb[8], aytcc[8];
	int axtbblen, axtcclen, aytbblen, aytcclen;
	REAL bxtaa[8], bxtcc[8], bytaa[8], bytcc[8];
	int bxtaalen, bxtcclen, bytaalen, bytcclen;
	REAL cxtaa[8], cxtbb[8], cytaa[8], cytbb[8];
	int cxtaalen, cxtbblen, cytaalen, cytbblen;
	REAL axtbc[8], aytbc[8], bxtca[8], bytca[8], cxtab[8], cytab[8];
	int axtbclen, aytbclen, bxtcalen, bytcalen, cxtablen, cytablen;
	REAL axtbct[16], aytbct[16], bxtcat[16], bytcat[16], cxtabt[16], cytabt[16];
	int axtbctlen, aytbctlen, bxtcatlen, bytcatlen, cxtabtlen, cytabtlen;
	REAL axtbctt[8], aytbctt[8], bxtcatt[8];
	REAL bytcatt[8], cxtabtt[8], cytabtt[8];
	int axtbcttlen, aytbcttlen, bxtcattlen, bytcattlen, cxtabttlen, cytabttlen;
	REAL abt[8], bct[8], cat[8];
	int abtlen, bctlen, catlen;
	REAL abtt[4], bctt[4], catt[4];
	int abttlen, bcttlen, cattlen;
	INEXACT REAL abtt3, bctt3, catt3;
	REAL negate;

	INEXACT REAL bvirt;
	REAL avirt, bround, around;
	INEXACT REAL c;
	INEXACT REAL abig;
	REAL ahi, alo, bhi, blo;
	REAL err1, err2, err3;
	INEXACT REAL _i, _j;
	REAL _0;

	adx = (REAL)(pa[0] - pd[0]);
	bdx = (REAL)(pb[0] - pd[0]);
	cdx = (REAL)(pc[0] - pd[0]);
	ady = (REAL)(pa[1] - pd[1]);
	bdy = (REAL)(pb[1] - pd[1]);
	cdy = (REAL)(pc[1] - pd[1]);

	Two_Product(bdx, cdy, bdxcdy1, bdxcdy0);
	Two_Product(cdx, bdy, cdxbdy1, cdxbdy0);
	Two_Two_Diff(bdxcdy1, bdxcdy0, cdxbdy1, cdxbdy0, bc3, bc[2], bc[1], bc[0]);
	bc[3] = bc3;
	axbclen = scale_expansion_zeroelim(4, bc, adx, axbc);
	axxbclen = scale_expansion_zeroelim(axbclen, axbc, adx, axxbc);
	aybclen = scale_expansion_zeroelim(4, bc, ady, aybc);
	ayybclen = scale_expansion_zeroelim(aybclen, aybc, ady, ayybc);
	alen = fast_expansion_sum_zeroelim(axxbclen, axxbc, ayybclen, ayybc, adet);

	Two_Product(cdx, ady, cdxady1, cdxady0);
	Two_Product(adx, cdy, adxcdy1, adxcdy0);
	Two_Two_Diff(cdxady1, cdxady0, adxcdy1, adxcdy0, ca3, ca[2], ca[1], ca[0]);
	ca[3] = ca3;
	bxcalen = scale_expansion_zeroelim(4, ca, bdx, bxca);
	bxxcalen = scale_expansion_zeroelim(bxcalen, bxca, bdx, bxxca);
	bycalen = scale_expansion_zeroelim(4, ca, bdy, byca);
	byycalen = scale_expansion_zeroelim(bycalen, byca, bdy, byyca);
	blen = fast_expansion_sum_zeroelim(bxxcalen, bxxca, byycalen, byyca, bdet);

	Two_Product(adx, bdy, adxbdy1, adxbdy0);
	Two_Product(bdx, ady, bdxady1, bdxady0);
	Two_Two_Diff(adxbdy1, adxbdy0, bdxady1, bdxady0, ab3, ab[2], ab[1], ab[0]);
	ab[3] = ab3;
	cxablen = scale_expansion_zeroelim(4, ab, cdx, cxab);
	cxxablen = scale_expansion_zeroelim(cxablen, cxab, cdx, cxxab);
	cyablen = scale_expansion_zeroelim(4, ab, cdy, cyab);
	cyyablen = scale_expansion_zeroelim(cyablen, cyab, cdy, cyyab);
	clen = fast_expansion_sum_zeroelim(cxxablen, cxxab, cyyablen, cyyab, cdet);

	ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, abdet);
	finlength = fast_expansion_sum_zeroelim(ablen, abdet, clen, cdet, fin1);

	det = estimate(finlength, fin1);
	errbound = iccerrboundB * permanent;
	if ((det >= errbound) || (-det >= errbound)) {
		return det;
	}

	Two_Diff_Tail(pa[0], pd[0], adx, adxtail);
	Two_Diff_Tail(pa[1], pd[1], ady, adytail);
	Two_Diff_Tail(pb[0], pd[0], bdx, bdxtail);
	Two_Diff_Tail(pb[1], pd[1], bdy, bdytail);
	Two_Diff_Tail(pc[0], pd[0], cdx, cdxtail);
	Two_Diff_Tail(pc[1], pd[1], cdy, cdytail);
	if ((adxtail == 0.0) && (bdxtail == 0.0) && (cdxtail == 0.0)
		&& (adytail == 0.0) && (bdytail == 0.0) && (cdytail == 0.0)) {
		return det;
	}

	errbound = iccerrboundC * permanent + resulterrbound * Absolute(det);
	det += ((adx * adx + ady * ady) * ((bdx * cdytail + cdy * bdxtail)
		- (bdy * cdxtail + cdx * bdytail))
		+ 2.0 * (adx * adxtail + ady * adytail) * (bdx * cdy - bdy * cdx))
		+ ((bdx * bdx + bdy * bdy) * ((cdx * adytail + ady * cdxtail)
			- (cdy * adxtail + adx * cdytail))
			+ 2.0 * (bdx * bdxtail + bdy * bdytail) * (cdx * ady - cdy * adx))
		+ ((cdx * cdx + cdy * cdy) * ((adx * bdytail + bdy * adxtail)
			- (ady * bdxtail + bdx * adytail))
			+ 2.0 * (cdx * cdxtail + cdy * cdytail) * (adx * bdy - ady * bdx));
	if ((det >= errbound) || (-det >= errbound)) {
		return det;
	}

	finnow = fin1;
	finother = fin2;

	if ((bdxtail != 0.0) || (bdytail != 0.0)
		|| (cdxtail != 0.0) || (cdytail != 0.0)) {
		Square(adx, adxadx1, adxadx0);
		Square(ady, adyady1, adyady0);
		Two_Two_Sum(adxadx1, adxadx0, adyady1, adyady0, aa3, aa[2], aa[1], aa[0]);
		aa[3] = aa3;
	}
	if ((cdxtail != 0.0) || (cdytail != 0.0)
		|| (adxtail != 0.0) || (adytail != 0.0)) {
		Square(bdx, bdxbdx1, bdxbdx0);
		Square(bdy, bdybdy1, bdybdy0);
		Two_Two_Sum(bdxbdx1, bdxbdx0, bdybdy1, bdybdy0, bb3, bb[2], bb[1], bb[0]);
		bb[3] = bb3;
	}
	if ((adxtail != 0.0) || (adytail != 0.0)
		|| (bdxtail != 0.0) || (bdytail != 0.0)) {
		Square(cdx, cdxcdx1, cdxcdx0);
		Square(cdy, cdycdy1, cdycdy0);
		Two_Two_Sum(cdxcdx1, cdxcdx0, cdycdy1, cdycdy0, cc3, cc[2], cc[1], cc[0]);
		cc[3] = cc3;
	}

	if (adxtail != 0.0) {
		axtbclen = scale_expansion_zeroelim(4, bc, adxtail, axtbc);
		temp16alen = scale_expansion_zeroelim(axtbclen, axtbc, 2.0 * adx,
			temp16a);

		axtcclen = scale_expansion_zeroelim(4, cc, adxtail, axtcc);
		temp16blen = scale_expansion_zeroelim(axtcclen, axtcc, bdy, temp16b);

		axtbblen = scale_expansion_zeroelim(4, bb, adxtail, axtbb);
		temp16clen = scale_expansion_zeroelim(axtbblen, axtbb, -cdy, temp16c);

		temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
			temp16blen, temp16b, temp32a);
		temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c,
			temp32alen, temp32a, temp48);
		finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
			temp48, finother);
		finswap = finnow; finnow = finother; finother = finswap;
	}
	if (adytail != 0.0) {
		aytbclen = scale_expansion_zeroelim(4, bc, adytail, aytbc);
		temp16alen = scale_expansion_zeroelim(aytbclen, aytbc, 2.0 * ady,
			temp16a);

		aytbblen = scale_expansion_zeroelim(4, bb, adytail, aytbb);
		temp16blen = scale_expansion_zeroelim(aytbblen, aytbb, cdx, temp16b);

		aytcclen = scale_expansion_zeroelim(4, cc, adytail, aytcc);
		temp16clen = scale_expansion_zeroelim(aytcclen, aytcc, -bdx, temp16c);

		temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
			temp16blen, temp16b, temp32a);
		temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c,
			temp32alen, temp32a, temp48);
		finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
			temp48, finother);
		finswap = finnow; finnow = finother; finother = finswap;
	}
	if (bdxtail != 0.0) {
		bxtcalen = scale_expansion_zeroelim(4, ca, bdxtail, bxtca);
		temp16alen = scale_expansion_zeroelim(bxtcalen, bxtca, 2.0 * bdx,
			temp16a);

		bxtaalen = scale_expansion_zeroelim(4, aa, bdxtail, bxtaa);
		temp16blen = scale_expansion_zeroelim(bxtaalen, bxtaa, cdy, temp16b);

		bxtcclen = scale_expansion_zeroelim(4, cc, bdxtail, bxtcc);
		temp16clen = scale_expansion_zeroelim(bxtcclen, bxtcc, -ady, temp16c);

		temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
			temp16blen, temp16b, temp32a);
		temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c,
			temp32alen, temp32a, temp48);
		finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
			temp48, finother);
		finswap = finnow; finnow = finother; finother = finswap;
	}
	if (bdytail != 0.0) {
		bytcalen = scale_expansion_zeroelim(4, ca, bdytail, bytca);
		temp16alen = scale_expansion_zeroelim(bytcalen, bytca, 2.0 * bdy,
			temp16a);

		bytcclen = scale_expansion_zeroelim(4, cc, bdytail, bytcc);
		temp16blen = scale_expansion_zeroelim(bytcclen, bytcc, adx, temp16b);

		bytaalen = scale_expansion_zeroelim(4, aa, bdytail, bytaa);
		temp16clen = scale_expansion_zeroelim(bytaalen, bytaa, -cdx, temp16c);

		temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
			temp16blen, temp16b, temp32a);
		temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c,
			temp32alen, temp32a, temp48);
		finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
			temp48, finother);
		finswap = finnow; finnow = finother; finother = finswap;
	}
	if (cdxtail != 0.0) {
		cxtablen = scale_expansion_zeroelim(4, ab, cdxtail, cxtab);
		temp16alen = scale_expansion_zeroelim(cxtablen, cxtab, 2.0 * cdx,
			temp16a);

		cxtbblen = scale_expansion_zeroelim(4, bb, cdxtail, cxtbb);
		temp16blen = scale_expansion_zeroelim(cxtbblen, cxtbb, ady, temp16b);

		cxtaalen = scale_expansion_zeroelim(4, aa, cdxtail, cxtaa);
		temp16clen = scale_expansion_zeroelim(cxtaalen, cxtaa, -bdy, temp16c);

		temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
			temp16blen, temp16b, temp32a);
		temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c,
			temp32alen, temp32a, temp48);
		finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
			temp48, finother);
		finswap = finnow; finnow = finother; finother = finswap;
	}
	if (cdytail != 0.0) {
		cytablen = scale_expansion_zeroelim(4, ab, cdytail, cytab);
		temp16alen = scale_expansion_zeroelim(cytablen, cytab, 2.0 * cdy,
			temp16a);

		cytaalen = scale_expansion_zeroelim(4, aa, cdytail, cytaa);
		temp16blen = scale_expansion_zeroelim(cytaalen, cytaa, bdx, temp16b);

		cytbblen = scale_expansion_zeroelim(4, bb, cdytail, cytbb);
		temp16clen = scale_expansion_zeroelim(cytbblen, cytbb, -adx, temp16c);

		temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
			temp16blen, temp16b, temp32a);
		temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c,
			temp32alen, temp32a, temp48);
		finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
			temp48, finother);
		finswap = finnow; finnow = finother; finother = finswap;
	}

	if ((adxtail != 0.0) || (adytail != 0.0)) {
		if ((bdxtail != 0.0) || (bdytail != 0.0)
			|| (cdxtail != 0.0) || (cdytail != 0.0)) {
			Two_Product(bdxtail, cdy, ti1, ti0);
			Two_Product(bdx, cdytail, tj1, tj0);
			Two_Two_Sum(ti1, ti0, tj1, tj0, u3, u[2], u[1], u[0]);
			u[3] = u3;
			negate = -bdy;
			Two_Product(cdxtail, negate, ti1, ti0);
			negate = -bdytail;
			Two_Product(cdx, negate, tj1, tj0);
			Two_Two_Sum(ti1, ti0, tj1, tj0, v3, v[2], v[1], v[0]);
			v[3] = v3;
			bctlen = fast_expansion_sum_zeroelim(4, u, 4, v, bct);

			Two_Product(bdxtail, cdytail, ti1, ti0);
			Two_Product(cdxtail, bdytail, tj1, tj0);
			Two_Two_Diff(ti1, ti0, tj1, tj0, bctt3, bctt[2], bctt[1], bctt[0]);
			bctt[3] = bctt3;
			bcttlen = 4;
		}
		else {
			bct[0] = 0.0;
			bctlen = 1;
			bctt[0] = 0.0;
			bcttlen = 1;
		}

		if (adxtail != 0.0) {
			temp16alen = scale_expansion_zeroelim(axtbclen, axtbc, adxtail, temp16a);
			axtbctlen = scale_expansion_zeroelim(bctlen, bct, adxtail, axtbct);
			temp32alen = scale_expansion_zeroelim(axtbctlen, axtbct, 2.0 * adx,
				temp32a);
			temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a,
				temp32alen, temp32a, temp48);
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
				temp48, finother);
			finswap = finnow; finnow = finother; finother = finswap;
			if (bdytail != 0.0) {
				temp8len = scale_expansion_zeroelim(4, cc, adxtail, temp8);
				temp16alen = scale_expansion_zeroelim(temp8len, temp8, bdytail,
					temp16a);
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen,
					temp16a, finother);
				finswap = finnow; finnow = finother; finother = finswap;
			}
			if (cdytail != 0.0) {
				temp8len = scale_expansion_zeroelim(4, bb, -adxtail, temp8);
				temp16alen = scale_expansion_zeroelim(temp8len, temp8, cdytail,
					temp16a);
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen,
					temp16a, finother);
				finswap = finnow; finnow = finother; finother = finswap;
			}

			temp32alen = scale_expansion_zeroelim(axtbctlen, axtbct, adxtail,
				temp32a);
			axtbcttlen = scale_expansion_zeroelim(bcttlen, bctt, adxtail, axtbctt);
			temp16alen = scale_expansion_zeroelim(axtbcttlen, axtbctt, 2.0 * adx,
				temp16a);
			temp16blen = scale_expansion_zeroelim(axtbcttlen, axtbctt, adxtail,
				temp16b);
			temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
				temp16blen, temp16b, temp32b);
			temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a,
				temp32blen, temp32b, temp64);
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len,
				temp64, finother);
			finswap = finnow; finnow = finother; finother = finswap;
		}
		if (adytail != 0.0) {
			temp16alen = scale_expansion_zeroelim(aytbclen, aytbc, adytail, temp16a);
			aytbctlen = scale_expansion_zeroelim(bctlen, bct, adytail, aytbct);
			temp32alen = scale_expansion_zeroelim(aytbctlen, aytbct, 2.0 * ady,
				temp32a);
			temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a,
				temp32alen, temp32a, temp48);
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
				temp48, finother);
			finswap = finnow; finnow = finother; finother = finswap;


			temp32alen = scale_expansion_zeroelim(aytbctlen, aytbct, adytail,
				temp32a);
			aytbcttlen = scale_expansion_zeroelim(bcttlen, bctt, adytail, aytbctt);
			temp16alen = scale_expansion_zeroelim(aytbcttlen, aytbctt, 2.0 * ady,
				temp16a);
			temp16blen = scale_expansion_zeroelim(aytbcttlen, aytbctt, adytail,
				temp16b);
			temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
				temp16blen, temp16b, temp32b);
			temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a,
				temp32blen, temp32b, temp64);
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len,
				temp64, finother);
			finswap = finnow; finnow = finother; finother = finswap;
		}
	}
	if ((bdxtail != 0.0) || (bdytail != 0.0)) {
		if ((cdxtail != 0.0) || (cdytail != 0.0)
			|| (adxtail != 0.0) || (adytail != 0.0)) {
			Two_Product(cdxtail, ady, ti1, ti0);
			Two_Product(cdx, adytail, tj1, tj0);
			Two_Two_Sum(ti1, ti0, tj1, tj0, u3, u[2], u[1], u[0]);
			u[3] = u3;
			negate = -cdy;
			Two_Product(adxtail, negate, ti1, ti0);
			negate = -cdytail;
			Two_Product(adx, negate, tj1, tj0);
			Two_Two_Sum(ti1, ti0, tj1, tj0, v3, v[2], v[1], v[0]);
			v[3] = v3;
			catlen = fast_expansion_sum_zeroelim(4, u, 4, v, cat);

			Two_Product(cdxtail, adytail, ti1, ti0);
			Two_Product(adxtail, cdytail, tj1, tj0);
			Two_Two_Diff(ti1, ti0, tj1, tj0, catt3, catt[2], catt[1], catt[0]);
			catt[3] = catt3;
			cattlen = 4;
		}
		else {
			cat[0] = 0.0;
			catlen = 1;
			catt[0] = 0.0;
			cattlen = 1;
		}

		if (bdxtail != 0.0) {
			temp16alen = scale_expansion_zeroelim(bxtcalen, bxtca, bdxtail, temp16a);
			bxtcatlen = scale_expansion_zeroelim(catlen, cat, bdxtail, bxtcat);
			temp32alen = scale_expansion_zeroelim(bxtcatlen, bxtcat, 2.0 * bdx,
				temp32a);
			temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a,
				temp32alen, temp32a, temp48);
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
				temp48, finother);
			finswap = finnow; finnow = finother; finother = finswap;
			if (cdytail != 0.0) {
				temp8len = scale_expansion_zeroelim(4, aa, bdxtail, temp8);
				temp16alen = scale_expansion_zeroelim(temp8len, temp8, cdytail,
					temp16a);
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen,
					temp16a, finother);
				finswap = finnow; finnow = finother; finother = finswap;
			}
			if (adytail != 0.0) {
				temp8len = scale_expansion_zeroelim(4, cc, -bdxtail, temp8);
				temp16alen = scale_expansion_zeroelim(temp8len, temp8, adytail,
					temp16a);
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen,
					temp16a, finother);
				finswap = finnow; finnow = finother; finother = finswap;
			}

			temp32alen = scale_expansion_zeroelim(bxtcatlen, bxtcat, bdxtail,
				temp32a);
			bxtcattlen = scale_expansion_zeroelim(cattlen, catt, bdxtail, bxtcatt);
			temp16alen = scale_expansion_zeroelim(bxtcattlen, bxtcatt, 2.0 * bdx,
				temp16a);
			temp16blen = scale_expansion_zeroelim(bxtcattlen, bxtcatt, bdxtail,
				temp16b);
			temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
				temp16blen, temp16b, temp32b);
			temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a,
				temp32blen, temp32b, temp64);
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len,
				temp64, finother);
			finswap = finnow; finnow = finother; finother = finswap;
		}
		if (bdytail != 0.0) {
			temp16alen = scale_expansion_zeroelim(bytcalen, bytca, bdytail, temp16a);
			bytcatlen = scale_expansion_zeroelim(catlen, cat, bdytail, bytcat);
			temp32alen = scale_expansion_zeroelim(bytcatlen, bytcat, 2.0 * bdy,
				temp32a);
			temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a,
				temp32alen, temp32a, temp48);
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
				temp48, finother);
			finswap = finnow; finnow = finother; finother = finswap;


			temp32alen = scale_expansion_zeroelim(bytcatlen, bytcat, bdytail,
				temp32a);
			bytcattlen = scale_expansion_zeroelim(cattlen, catt, bdytail, bytcatt);
			temp16alen = scale_expansion_zeroelim(bytcattlen, bytcatt, 2.0 * bdy,
				temp16a);
			temp16blen = scale_expansion_zeroelim(bytcattlen, bytcatt, bdytail,
				temp16b);
			temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
				temp16blen, temp16b, temp32b);
			temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a,
				temp32blen, temp32b, temp64);
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len,
				temp64, finother);
			finswap = finnow; finnow = finother; finother = finswap;
		}
	}
	if ((cdxtail != 0.0) || (cdytail != 0.0)) {
		if ((adxtail != 0.0) || (adytail != 0.0)
			|| (bdxtail != 0.0) || (bdytail != 0.0)) {
			Two_Product(adxtail, bdy, ti1, ti0);
			Two_Product(adx, bdytail, tj1, tj0);
			Two_Two_Sum(ti1, ti0, tj1, tj0, u3, u[2], u[1], u[0]);
			u[3] = u3;
			negate = -ady;
			Two_Product(bdxtail, negate, ti1, ti0);
			negate = -adytail;
			Two_Product(bdx, negate, tj1, tj0);
			Two_Two_Sum(ti1, ti0, tj1, tj0, v3, v[2], v[1], v[0]);
			v[3] = v3;
			abtlen = fast_expansion_sum_zeroelim(4, u, 4, v, abt);

			Two_Product(adxtail, bdytail, ti1, ti0);
			Two_Product(bdxtail, adytail, tj1, tj0);
			Two_Two_Diff(ti1, ti0, tj1, tj0, abtt3, abtt[2], abtt[1], abtt[0]);
			abtt[3] = abtt3;
			abttlen = 4;
		}
		else {
			abt[0] = 0.0;
			abtlen = 1;
			abtt[0] = 0.0;
			abttlen = 1;
		}

		if (cdxtail != 0.0) {
			temp16alen = scale_expansion_zeroelim(cxtablen, cxtab, cdxtail, temp16a);
			cxtabtlen = scale_expansion_zeroelim(abtlen, abt, cdxtail, cxtabt);
			temp32alen = scale_expansion_zeroelim(cxtabtlen, cxtabt, 2.0 * cdx,
				temp32a);
			temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a,
				temp32alen, temp32a, temp48);
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
				temp48, finother);
			finswap = finnow; finnow = finother; finother = finswap;
			if (adytail != 0.0) {
				temp8len = scale_expansion_zeroelim(4, bb, cdxtail, temp8);
				temp16alen = scale_expansion_zeroelim(temp8len, temp8, adytail,
					temp16a);
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen,
					temp16a, finother);
				finswap = finnow; finnow = finother; finother = finswap;
			}
			if (bdytail != 0.0) {
				temp8len = scale_expansion_zeroelim(4, aa, -cdxtail, temp8);
				temp16alen = scale_expansion_zeroelim(temp8len, temp8, bdytail,
					temp16a);
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen,
					temp16a, finother);
				finswap = finnow; finnow = finother; finother = finswap;
			}

			temp32alen = scale_expansion_zeroelim(cxtabtlen, cxtabt, cdxtail,
				temp32a);
			cxtabttlen = scale_expansion_zeroelim(abttlen, abtt, cdxtail, cxtabtt);
			temp16alen = scale_expansion_zeroelim(cxtabttlen, cxtabtt, 2.0 * cdx,
				temp16a);
			temp16blen = scale_expansion_zeroelim(cxtabttlen, cxtabtt, cdxtail,
				temp16b);
			temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
				temp16blen, temp16b, temp32b);
			temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a,
				temp32blen, temp32b, temp64);
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len,
				temp64, finother);
			finswap = finnow; finnow = finother; finother = finswap;
		}
		if (cdytail != 0.0) {
			temp16alen = scale_expansion_zeroelim(cytablen, cytab, cdytail, temp16a);
			cytabtlen = scale_expansion_zeroelim(abtlen, abt, cdytail, cytabt);
			temp32alen = scale_expansion_zeroelim(cytabtlen, cytabt, 2.0 * cdy,
				temp32a);
			temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a,
				temp32alen, temp32a, temp48);
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
				temp48, finother);
			finswap = finnow; finnow = finother; finother = finswap;


			temp32alen = scale_expansion_zeroelim(cytabtlen, cytabt, cdytail,
				temp32a);
			cytabttlen = scale_expansion_zeroelim(abttlen, abtt, cdytail, cytabtt);
			temp16alen = scale_expansion_zeroelim(cytabttlen, cytabtt, 2.0 * cdy,
				temp16a);
			temp16blen = scale_expansion_zeroelim(cytabttlen, cytabtt, cdytail,
				temp16b);
			temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
				temp16blen, temp16b, temp32b);
			temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a,
				temp32blen, temp32b, temp64);
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len,
				temp64, finother);
			finswap = finnow; finnow = finother; finother = finswap;
		}
	}

	return finnow[finlength - 1];
}

REAL incircle(REAL *pa, REAL *pb, REAL *pc, REAL *pd)
{
	REAL adx, bdx, cdx, ady, bdy, cdy;
	REAL bdxcdy, cdxbdy, cdxady, adxcdy, adxbdy, bdxady;
	REAL alift, blift, clift;
	REAL det;
	REAL permanent, errbound;

	adx = pa[0] - pd[0];
	bdx = pb[0] - pd[0];
	cdx = pc[0] - pd[0];
	ady = pa[1] - pd[1];
	bdy = pb[1] - pd[1];
	cdy = pc[1] - pd[1];

	bdxcdy = bdx * cdy;
	cdxbdy = cdx * bdy;
	alift = adx * adx + ady * ady;

	cdxady = cdx * ady;
	adxcdy = adx * cdy;
	blift = bdx * bdx + bdy * bdy;

	adxbdy = adx * bdy;
	bdxady = bdx * ady;
	clift = cdx * cdx + cdy * cdy;

	det = alift * (bdxcdy - cdxbdy)
		+ blift * (cdxady - adxcdy)
		+ clift * (adxbdy - bdxady);

	permanent = (Absolute(bdxcdy) + Absolute(cdxbdy)) * alift
		+ (Absolute(cdxady) + Absolute(adxcdy)) * blift
		+ (Absolute(adxbdy) + Absolute(bdxady)) * clift;
	errbound = iccerrboundA * permanent;
	if ((det > errbound) || (-det > errbound)) {
		return det;
	}

	return incircleadapt(pa, pb, pc, pd, permanent);
}

/******************************************************************************/
/*                                                                           **/
/*  inspherefast()   Approximate 3D insphere test.  Nonrobust.               **/
/*  insphereexact()   Exact 3D insphere test.  Robust.                       **/
/*  insphereslow()   Another exact 3D insphere test.  Robust.                **/
/*  insphere()   Adaptive exact 3D insphere test.  Robust.                   **/
/*                                                                           **/
/*               Return a positive value if the point pe lies inside the     **/
/*               sphere passing through pa, pb, pc, and pd; a negative value **/
/*               if it lies outside; and zero if the five points are         **/
/*               cospherical.  The points pa, pb, pc, and pd must be ordered **/
/*               so that they have a positive orientation (as defined by     **/
/*               orient3d()), or the sign of the result will be reversed.    **/
/*                                                                           **/
/*  Only the first and last routine should be used; the middle two are for   **/
/*  timings.                                                                 **/
/*                                                                           **/
/*  The last three use exact arithmetic to ensure a correct answer.  The     **/
/*  result returned is the determinant of a matrix.  In insphere() only,     **/
/*  this determinant is computed adaptively, in the sense that exact         **/
/*  arithmetic is used only to the degree it is needed to ensure that the    **/
/*  returned value has the correct sign.  Hence, insphere() is usually quite **/
/*  fast, but will run more slowly when the input points are cospherical or  **/
/*  nearly so.                                                               **/
/*                                                                           **/
/******************************************************************************/

REAL inspherefast(REAL *pa, REAL *pb, REAL *pc, REAL *pd, REAL *pe)
{
	REAL aex, bex, cex, dex;
	REAL aey, bey, cey, dey;
	REAL aez, bez, cez, dez;
	REAL alift, blift, clift, dlift;
	REAL ab, bc, cd, da, ac, bd;
	REAL abc, bcd, cda, dab;

	aex = pa[0] - pe[0];
	bex = pb[0] - pe[0];
	cex = pc[0] - pe[0];
	dex = pd[0] - pe[0];
	aey = pa[1] - pe[1];
	bey = pb[1] - pe[1];
	cey = pc[1] - pe[1];
	dey = pd[1] - pe[1];
	aez = pa[2] - pe[2];
	bez = pb[2] - pe[2];
	cez = pc[2] - pe[2];
	dez = pd[2] - pe[2];

	ab = aex * bey - bex * aey;
	bc = bex * cey - cex * bey;
	cd = cex * dey - dex * cey;
	da = dex * aey - aex * dey;

	ac = aex * cey - cex * aey;
	bd = bex * dey - dex * bey;

	abc = aez * bc - bez * ac + cez * ab;
	bcd = bez * cd - cez * bd + dez * bc;
	cda = cez * da + dez * ac + aez * cd;
	dab = dez * ab + aez * bd + bez * da;

	alift = aex * aex + aey * aey + aez * aez;
	blift = bex * bex + bey * bey + bez * bez;
	clift = cex * cex + cey * cey + cez * cez;
	dlift = dex * dex + dey * dey + dez * dez;

	return (dlift * abc - clift * dab) + (blift * cda - alift * bcd);
}

REAL insphereexact(REAL *pa, REAL *pb, REAL *pc, REAL *pd, REAL *pe)
{
	INEXACT REAL axby1, bxcy1, cxdy1, dxey1, exay1;
	INEXACT REAL bxay1, cxby1, dxcy1, exdy1, axey1;
	INEXACT REAL axcy1, bxdy1, cxey1, dxay1, exby1;
	INEXACT REAL cxay1, dxby1, excy1, axdy1, bxey1;
	REAL axby0, bxcy0, cxdy0, dxey0, exay0;
	REAL bxay0, cxby0, dxcy0, exdy0, axey0;
	REAL axcy0, bxdy0, cxey0, dxay0, exby0;
	REAL cxay0, dxby0, excy0, axdy0, bxey0;
	REAL ab[4], bc[4], cd[4], de[4], ea[4];
	REAL ac[4], bd[4], ce[4], da[4], eb[4];
	REAL temp8a[8], temp8b[8], temp16[16];
	int temp8alen, temp8blen, temp16len;
	REAL abc[24], bcd[24], cde[24], dea[24], eab[24];
	REAL abd[24], bce[24], cda[24], deb[24], eac[24];
	int abclen, bcdlen, cdelen, dealen, eablen;
	int abdlen, bcelen, cdalen, deblen, eaclen;
	REAL temp48a[48], temp48b[48];
	int temp48alen, temp48blen;
	REAL abcd[96], bcde[96], cdea[96], deab[96], eabc[96];
	int abcdlen, bcdelen, cdealen, deablen, eabclen;
	REAL temp192[192];
	REAL det384x[384], det384y[384], det384z[384];
	int xlen, ylen, zlen;
	REAL detxy[768];
	int xylen;
	REAL adet[1152], bdet[1152], cdet[1152], ddet[1152], edet[1152];
	int alen, blen, clen, dlen, elen;
	REAL abdet[2304], cddet[2304], cdedet[3456];
	int ablen, cdlen;
	REAL deter[5760];
	int deterlen;
	int i;

	INEXACT REAL bvirt;
	REAL avirt, bround, around;
	INEXACT REAL c;
	INEXACT REAL abig;
	REAL ahi, alo, bhi, blo;
	REAL err1, err2, err3;
	INEXACT REAL _i, _j;
	REAL _0;

	Two_Product(pa[0], pb[1], axby1, axby0);
	Two_Product(pb[0], pa[1], bxay1, bxay0);
	Two_Two_Diff(axby1, axby0, bxay1, bxay0, ab[3], ab[2], ab[1], ab[0]);

	Two_Product(pb[0], pc[1], bxcy1, bxcy0);
	Two_Product(pc[0], pb[1], cxby1, cxby0);
	Two_Two_Diff(bxcy1, bxcy0, cxby1, cxby0, bc[3], bc[2], bc[1], bc[0]);

	Two_Product(pc[0], pd[1], cxdy1, cxdy0);
	Two_Product(pd[0], pc[1], dxcy1, dxcy0);
	Two_Two_Diff(cxdy1, cxdy0, dxcy1, dxcy0, cd[3], cd[2], cd[1], cd[0]);

	Two_Product(pd[0], pe[1], dxey1, dxey0);
	Two_Product(pe[0], pd[1], exdy1, exdy0);
	Two_Two_Diff(dxey1, dxey0, exdy1, exdy0, de[3], de[2], de[1], de[0]);

	Two_Product(pe[0], pa[1], exay1, exay0);
	Two_Product(pa[0], pe[1], axey1, axey0);
	Two_Two_Diff(exay1, exay0, axey1, axey0, ea[3], ea[2], ea[1], ea[0]);

	Two_Product(pa[0], pc[1], axcy1, axcy0);
	Two_Product(pc[0], pa[1], cxay1, cxay0);
	Two_Two_Diff(axcy1, axcy0, cxay1, cxay0, ac[3], ac[2], ac[1], ac[0]);

	Two_Product(pb[0], pd[1], bxdy1, bxdy0);
	Two_Product(pd[0], pb[1], dxby1, dxby0);
	Two_Two_Diff(bxdy1, bxdy0, dxby1, dxby0, bd[3], bd[2], bd[1], bd[0]);

	Two_Product(pc[0], pe[1], cxey1, cxey0);
	Two_Product(pe[0], pc[1], excy1, excy0);
	Two_Two_Diff(cxey1, cxey0, excy1, excy0, ce[3], ce[2], ce[1], ce[0]);

	Two_Product(pd[0], pa[1], dxay1, dxay0);
	Two_Product(pa[0], pd[1], axdy1, axdy0);
	Two_Two_Diff(dxay1, dxay0, axdy1, axdy0, da[3], da[2], da[1], da[0]);

	Two_Product(pe[0], pb[1], exby1, exby0);
	Two_Product(pb[0], pe[1], bxey1, bxey0);
	Two_Two_Diff(exby1, exby0, bxey1, bxey0, eb[3], eb[2], eb[1], eb[0]);

	temp8alen = scale_expansion_zeroelim(4, bc, pa[2], temp8a);
	temp8blen = scale_expansion_zeroelim(4, ac, -pb[2], temp8b);
	temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b,
		temp16);
	temp8alen = scale_expansion_zeroelim(4, ab, pc[2], temp8a);
	abclen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16,
		abc);

	temp8alen = scale_expansion_zeroelim(4, cd, pb[2], temp8a);
	temp8blen = scale_expansion_zeroelim(4, bd, -pc[2], temp8b);
	temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b,
		temp16);
	temp8alen = scale_expansion_zeroelim(4, bc, pd[2], temp8a);
	bcdlen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16,
		bcd);

	temp8alen = scale_expansion_zeroelim(4, de, pc[2], temp8a);
	temp8blen = scale_expansion_zeroelim(4, ce, -pd[2], temp8b);
	temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b,
		temp16);
	temp8alen = scale_expansion_zeroelim(4, cd, pe[2], temp8a);
	cdelen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16,
		cde);

	temp8alen = scale_expansion_zeroelim(4, ea, pd[2], temp8a);
	temp8blen = scale_expansion_zeroelim(4, da, -pe[2], temp8b);
	temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b,
		temp16);
	temp8alen = scale_expansion_zeroelim(4, de, pa[2], temp8a);
	dealen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16,
		dea);

	temp8alen = scale_expansion_zeroelim(4, ab, pe[2], temp8a);
	temp8blen = scale_expansion_zeroelim(4, eb, -pa[2], temp8b);
	temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b,
		temp16);
	temp8alen = scale_expansion_zeroelim(4, ea, pb[2], temp8a);
	eablen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16,
		eab);

	temp8alen = scale_expansion_zeroelim(4, bd, pa[2], temp8a);
	temp8blen = scale_expansion_zeroelim(4, da, pb[2], temp8b);
	temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b,
		temp16);
	temp8alen = scale_expansion_zeroelim(4, ab, pd[2], temp8a);
	abdlen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16,
		abd);

	temp8alen = scale_expansion_zeroelim(4, ce, pb[2], temp8a);
	temp8blen = scale_expansion_zeroelim(4, eb, pc[2], temp8b);
	temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b,
		temp16);
	temp8alen = scale_expansion_zeroelim(4, bc, pe[2], temp8a);
	bcelen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16,
		bce);

	temp8alen = scale_expansion_zeroelim(4, da, pc[2], temp8a);
	temp8blen = scale_expansion_zeroelim(4, ac, pd[2], temp8b);
	temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b,
		temp16);
	temp8alen = scale_expansion_zeroelim(4, cd, pa[2], temp8a);
	cdalen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16,
		cda);

	temp8alen = scale_expansion_zeroelim(4, eb, pd[2], temp8a);
	temp8blen = scale_expansion_zeroelim(4, bd, pe[2], temp8b);
	temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b,
		temp16);
	temp8alen = scale_expansion_zeroelim(4, de, pb[2], temp8a);
	deblen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16,
		deb);

	temp8alen = scale_expansion_zeroelim(4, ac, pe[2], temp8a);
	temp8blen = scale_expansion_zeroelim(4, ce, pa[2], temp8b);
	temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b,
		temp16);
	temp8alen = scale_expansion_zeroelim(4, ea, pc[2], temp8a);
	eaclen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16,
		eac);

	temp48alen = fast_expansion_sum_zeroelim(cdelen, cde, bcelen, bce, temp48a);
	temp48blen = fast_expansion_sum_zeroelim(deblen, deb, bcdlen, bcd, temp48b);
	for (i = 0; i < temp48blen; i++) {
		temp48b[i] = -temp48b[i];
	}
	bcdelen = fast_expansion_sum_zeroelim(temp48alen, temp48a,
		temp48blen, temp48b, bcde);
	xlen = scale_expansion_zeroelim(bcdelen, bcde, pa[0], temp192);
	xlen = scale_expansion_zeroelim(xlen, temp192, pa[0], det384x);
	ylen = scale_expansion_zeroelim(bcdelen, bcde, pa[1], temp192);
	ylen = scale_expansion_zeroelim(ylen, temp192, pa[1], det384y);
	zlen = scale_expansion_zeroelim(bcdelen, bcde, pa[2], temp192);
	zlen = scale_expansion_zeroelim(zlen, temp192, pa[2], det384z);
	xylen = fast_expansion_sum_zeroelim(xlen, det384x, ylen, det384y, detxy);
	alen = fast_expansion_sum_zeroelim(xylen, detxy, zlen, det384z, adet);

	temp48alen = fast_expansion_sum_zeroelim(dealen, dea, cdalen, cda, temp48a);
	temp48blen = fast_expansion_sum_zeroelim(eaclen, eac, cdelen, cde, temp48b);
	for (i = 0; i < temp48blen; i++) {
		temp48b[i] = -temp48b[i];
	}
	cdealen = fast_expansion_sum_zeroelim(temp48alen, temp48a,
		temp48blen, temp48b, cdea);
	xlen = scale_expansion_zeroelim(cdealen, cdea, pb[0], temp192);
	xlen = scale_expansion_zeroelim(xlen, temp192, pb[0], det384x);
	ylen = scale_expansion_zeroelim(cdealen, cdea, pb[1], temp192);
	ylen = scale_expansion_zeroelim(ylen, temp192, pb[1], det384y);
	zlen = scale_expansion_zeroelim(cdealen, cdea, pb[2], temp192);
	zlen = scale_expansion_zeroelim(zlen, temp192, pb[2], det384z);
	xylen = fast_expansion_sum_zeroelim(xlen, det384x, ylen, det384y, detxy);
	blen = fast_expansion_sum_zeroelim(xylen, detxy, zlen, det384z, bdet);

	temp48alen = fast_expansion_sum_zeroelim(eablen, eab, deblen, deb, temp48a);
	temp48blen = fast_expansion_sum_zeroelim(abdlen, abd, dealen, dea, temp48b);
	for (i = 0; i < temp48blen; i++) {
		temp48b[i] = -temp48b[i];
	}
	deablen = fast_expansion_sum_zeroelim(temp48alen, temp48a,
		temp48blen, temp48b, deab);
	xlen = scale_expansion_zeroelim(deablen, deab, pc[0], temp192);
	xlen = scale_expansion_zeroelim(xlen, temp192, pc[0], det384x);
	ylen = scale_expansion_zeroelim(deablen, deab, pc[1], temp192);
	ylen = scale_expansion_zeroelim(ylen, temp192, pc[1], det384y);
	zlen = scale_expansion_zeroelim(deablen, deab, pc[2], temp192);
	zlen = scale_expansion_zeroelim(zlen, temp192, pc[2], det384z);
	xylen = fast_expansion_sum_zeroelim(xlen, det384x, ylen, det384y, detxy);
	clen = fast_expansion_sum_zeroelim(xylen, detxy, zlen, det384z, cdet);

	temp48alen = fast_expansion_sum_zeroelim(abclen, abc, eaclen, eac, temp48a);
	temp48blen = fast_expansion_sum_zeroelim(bcelen, bce, eablen, eab, temp48b);
	for (i = 0; i < temp48blen; i++) {
		temp48b[i] = -temp48b[i];
	}
	eabclen = fast_expansion_sum_zeroelim(temp48alen, temp48a,
		temp48blen, temp48b, eabc);
	xlen = scale_expansion_zeroelim(eabclen, eabc, pd[0], temp192);
	xlen = scale_expansion_zeroelim(xlen, temp192, pd[0], det384x);
	ylen = scale_expansion_zeroelim(eabclen, eabc, pd[1], temp192);
	ylen = scale_expansion_zeroelim(ylen, temp192, pd[1], det384y);
	zlen = scale_expansion_zeroelim(eabclen, eabc, pd[2], temp192);
	zlen = scale_expansion_zeroelim(zlen, temp192, pd[2], det384z);
	xylen = fast_expansion_sum_zeroelim(xlen, det384x, ylen, det384y, detxy);
	dlen = fast_expansion_sum_zeroelim(xylen, detxy, zlen, det384z, ddet);

	temp48alen = fast_expansion_sum_zeroelim(bcdlen, bcd, abdlen, abd, temp48a);
	temp48blen = fast_expansion_sum_zeroelim(cdalen, cda, abclen, abc, temp48b);
	for (i = 0; i < temp48blen; i++) {
		temp48b[i] = -temp48b[i];
	}
	abcdlen = fast_expansion_sum_zeroelim(temp48alen, temp48a,
		temp48blen, temp48b, abcd);
	xlen = scale_expansion_zeroelim(abcdlen, abcd, pe[0], temp192);
	xlen = scale_expansion_zeroelim(xlen, temp192, pe[0], det384x);
	ylen = scale_expansion_zeroelim(abcdlen, abcd, pe[1], temp192);
	ylen = scale_expansion_zeroelim(ylen, temp192, pe[1], det384y);
	zlen = scale_expansion_zeroelim(abcdlen, abcd, pe[2], temp192);
	zlen = scale_expansion_zeroelim(zlen, temp192, pe[2], det384z);
	xylen = fast_expansion_sum_zeroelim(xlen, det384x, ylen, det384y, detxy);
	elen = fast_expansion_sum_zeroelim(xylen, detxy, zlen, det384z, edet);

	ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, abdet);
	cdlen = fast_expansion_sum_zeroelim(clen, cdet, dlen, ddet, cddet);
	cdelen = fast_expansion_sum_zeroelim(cdlen, cddet, elen, edet, cdedet);
	deterlen = fast_expansion_sum_zeroelim(ablen, abdet, cdelen, cdedet, deter);

	return deter[deterlen - 1];
}

REAL insphereslow(REAL *pa, REAL *pb, REAL *pc, REAL *pd, REAL *pe)
{
	INEXACT REAL aex, bex, cex, dex, aey, bey, cey, dey, aez, bez, cez, dez;
	REAL aextail, bextail, cextail, dextail;
	REAL aeytail, beytail, ceytail, deytail;
	REAL aeztail, beztail, ceztail, deztail;
	REAL negate, negatetail;
	INEXACT REAL axby7, bxcy7, cxdy7, dxay7, axcy7, bxdy7;
	INEXACT REAL bxay7, cxby7, dxcy7, axdy7, cxay7, dxby7;
	REAL axby[8], bxcy[8], cxdy[8], dxay[8], axcy[8], bxdy[8];
	REAL bxay[8], cxby[8], dxcy[8], axdy[8], cxay[8], dxby[8];
	REAL ab[16], bc[16], cd[16], da[16], ac[16], bd[16];
	int ablen, bclen, cdlen, dalen, aclen, bdlen;
	REAL temp32a[32], temp32b[32], temp64a[64], temp64b[64], temp64c[64];
	int temp32alen, temp32blen, temp64alen, temp64blen, temp64clen;
	REAL temp128[128], temp192[192];
	int temp128len, temp192len;
	REAL detx[384], detxx[768], detxt[384], detxxt[768], detxtxt[768];
	int xlen, xxlen, xtlen, xxtlen, xtxtlen;
	REAL x1[1536], x2[2304];
	int x1len, x2len;
	REAL dety[384], detyy[768], detyt[384], detyyt[768], detytyt[768];
	int ylen, yylen, ytlen, yytlen, ytytlen;
	REAL y1[1536], y2[2304];
	int y1len, y2len;
	REAL detz[384], detzz[768], detzt[384], detzzt[768], detztzt[768];
	int zlen, zzlen, ztlen, zztlen, ztztlen;
	REAL z1[1536], z2[2304];
	int z1len, z2len;
	REAL detxy[4608];
	int xylen;
	REAL adet[6912], bdet[6912], cdet[6912], ddet[6912];
	int alen, blen, clen, dlen;
	REAL abdet[13824], cddet[13824], deter[27648];
	int deterlen;
	int i;

	INEXACT REAL bvirt;
	REAL avirt, bround, around;
	INEXACT REAL c;
	INEXACT REAL abig;
	REAL a0hi, a0lo, a1hi, a1lo, bhi, blo;
	REAL err1, err2, err3;
	INEXACT REAL _i, _j, _k, _l, _m, _n;
	REAL _0, _1, _2;

	Two_Diff(pa[0], pe[0], aex, aextail);
	Two_Diff(pa[1], pe[1], aey, aeytail);
	Two_Diff(pa[2], pe[2], aez, aeztail);
	Two_Diff(pb[0], pe[0], bex, bextail);
	Two_Diff(pb[1], pe[1], bey, beytail);
	Two_Diff(pb[2], pe[2], bez, beztail);
	Two_Diff(pc[0], pe[0], cex, cextail);
	Two_Diff(pc[1], pe[1], cey, ceytail);
	Two_Diff(pc[2], pe[2], cez, ceztail);
	Two_Diff(pd[0], pe[0], dex, dextail);
	Two_Diff(pd[1], pe[1], dey, deytail);
	Two_Diff(pd[2], pe[2], dez, deztail);

	Two_Two_Product(aex, aextail, bey, beytail,
		axby7, axby[6], axby[5], axby[4],
		axby[3], axby[2], axby[1], axby[0]);
	axby[7] = axby7;
	negate = -aey;
	negatetail = -aeytail;
	Two_Two_Product(bex, bextail, negate, negatetail,
		bxay7, bxay[6], bxay[5], bxay[4],
		bxay[3], bxay[2], bxay[1], bxay[0]);
	bxay[7] = bxay7;
	ablen = fast_expansion_sum_zeroelim(8, axby, 8, bxay, ab);
	Two_Two_Product(bex, bextail, cey, ceytail,
		bxcy7, bxcy[6], bxcy[5], bxcy[4],
		bxcy[3], bxcy[2], bxcy[1], bxcy[0]);
	bxcy[7] = bxcy7;
	negate = -bey;
	negatetail = -beytail;
	Two_Two_Product(cex, cextail, negate, negatetail,
		cxby7, cxby[6], cxby[5], cxby[4],
		cxby[3], cxby[2], cxby[1], cxby[0]);
	cxby[7] = cxby7;
	bclen = fast_expansion_sum_zeroelim(8, bxcy, 8, cxby, bc);
	Two_Two_Product(cex, cextail, dey, deytail,
		cxdy7, cxdy[6], cxdy[5], cxdy[4],
		cxdy[3], cxdy[2], cxdy[1], cxdy[0]);
	cxdy[7] = cxdy7;
	negate = -cey;
	negatetail = -ceytail;
	Two_Two_Product(dex, dextail, negate, negatetail,
		dxcy7, dxcy[6], dxcy[5], dxcy[4],
		dxcy[3], dxcy[2], dxcy[1], dxcy[0]);
	dxcy[7] = dxcy7;
	cdlen = fast_expansion_sum_zeroelim(8, cxdy, 8, dxcy, cd);
	Two_Two_Product(dex, dextail, aey, aeytail,
		dxay7, dxay[6], dxay[5], dxay[4],
		dxay[3], dxay[2], dxay[1], dxay[0]);
	dxay[7] = dxay7;
	negate = -dey;
	negatetail = -deytail;
	Two_Two_Product(aex, aextail, negate, negatetail,
		axdy7, axdy[6], axdy[5], axdy[4],
		axdy[3], axdy[2], axdy[1], axdy[0]);
	axdy[7] = axdy7;
	dalen = fast_expansion_sum_zeroelim(8, dxay, 8, axdy, da);
	Two_Two_Product(aex, aextail, cey, ceytail,
		axcy7, axcy[6], axcy[5], axcy[4],
		axcy[3], axcy[2], axcy[1], axcy[0]);
	axcy[7] = axcy7;
	negate = -aey;
	negatetail = -aeytail;
	Two_Two_Product(cex, cextail, negate, negatetail,
		cxay7, cxay[6], cxay[5], cxay[4],
		cxay[3], cxay[2], cxay[1], cxay[0]);
	cxay[7] = cxay7;
	aclen = fast_expansion_sum_zeroelim(8, axcy, 8, cxay, ac);
	Two_Two_Product(bex, bextail, dey, deytail,
		bxdy7, bxdy[6], bxdy[5], bxdy[4],
		bxdy[3], bxdy[2], bxdy[1], bxdy[0]);
	bxdy[7] = bxdy7;
	negate = -bey;
	negatetail = -beytail;
	Two_Two_Product(dex, dextail, negate, negatetail,
		dxby7, dxby[6], dxby[5], dxby[4],
		dxby[3], dxby[2], dxby[1], dxby[0]);
	dxby[7] = dxby7;
	bdlen = fast_expansion_sum_zeroelim(8, bxdy, 8, dxby, bd);

	temp32alen = scale_expansion_zeroelim(cdlen, cd, -bez, temp32a);
	temp32blen = scale_expansion_zeroelim(cdlen, cd, -beztail, temp32b);
	temp64alen = fast_expansion_sum_zeroelim(temp32alen, temp32a,
		temp32blen, temp32b, temp64a);
	temp32alen = scale_expansion_zeroelim(bdlen, bd, cez, temp32a);
	temp32blen = scale_expansion_zeroelim(bdlen, bd, ceztail, temp32b);
	temp64blen = fast_expansion_sum_zeroelim(temp32alen, temp32a,
		temp32blen, temp32b, temp64b);
	temp32alen = scale_expansion_zeroelim(bclen, bc, -dez, temp32a);
	temp32blen = scale_expansion_zeroelim(bclen, bc, -deztail, temp32b);
	temp64clen = fast_expansion_sum_zeroelim(temp32alen, temp32a,
		temp32blen, temp32b, temp64c);
	temp128len = fast_expansion_sum_zeroelim(temp64alen, temp64a,
		temp64blen, temp64b, temp128);
	temp192len = fast_expansion_sum_zeroelim(temp64clen, temp64c,
		temp128len, temp128, temp192);
	xlen = scale_expansion_zeroelim(temp192len, temp192, aex, detx);
	xxlen = scale_expansion_zeroelim(xlen, detx, aex, detxx);
	xtlen = scale_expansion_zeroelim(temp192len, temp192, aextail, detxt);
	xxtlen = scale_expansion_zeroelim(xtlen, detxt, aex, detxxt);
	for (i = 0; i < xxtlen; i++) {
		detxxt[i] *= 2.0;
	}
	xtxtlen = scale_expansion_zeroelim(xtlen, detxt, aextail, detxtxt);
	x1len = fast_expansion_sum_zeroelim(xxlen, detxx, xxtlen, detxxt, x1);
	x2len = fast_expansion_sum_zeroelim(x1len, x1, xtxtlen, detxtxt, x2);
	ylen = scale_expansion_zeroelim(temp192len, temp192, aey, dety);
	yylen = scale_expansion_zeroelim(ylen, dety, aey, detyy);
	ytlen = scale_expansion_zeroelim(temp192len, temp192, aeytail, detyt);
	yytlen = scale_expansion_zeroelim(ytlen, detyt, aey, detyyt);
	for (i = 0; i < yytlen; i++) {
		detyyt[i] *= 2.0;
	}
	ytytlen = scale_expansion_zeroelim(ytlen, detyt, aeytail, detytyt);
	y1len = fast_expansion_sum_zeroelim(yylen, detyy, yytlen, detyyt, y1);
	y2len = fast_expansion_sum_zeroelim(y1len, y1, ytytlen, detytyt, y2);
	zlen = scale_expansion_zeroelim(temp192len, temp192, aez, detz);
	zzlen = scale_expansion_zeroelim(zlen, detz, aez, detzz);
	ztlen = scale_expansion_zeroelim(temp192len, temp192, aeztail, detzt);
	zztlen = scale_expansion_zeroelim(ztlen, detzt, aez, detzzt);
	for (i = 0; i < zztlen; i++) {
		detzzt[i] *= 2.0;
	}
	ztztlen = scale_expansion_zeroelim(ztlen, detzt, aeztail, detztzt);
	z1len = fast_expansion_sum_zeroelim(zzlen, detzz, zztlen, detzzt, z1);
	z2len = fast_expansion_sum_zeroelim(z1len, z1, ztztlen, detztzt, z2);
	xylen = fast_expansion_sum_zeroelim(x2len, x2, y2len, y2, detxy);
	alen = fast_expansion_sum_zeroelim(z2len, z2, xylen, detxy, adet);

	temp32alen = scale_expansion_zeroelim(dalen, da, cez, temp32a);
	temp32blen = scale_expansion_zeroelim(dalen, da, ceztail, temp32b);
	temp64alen = fast_expansion_sum_zeroelim(temp32alen, temp32a,
		temp32blen, temp32b, temp64a);
	temp32alen = scale_expansion_zeroelim(aclen, ac, dez, temp32a);
	temp32blen = scale_expansion_zeroelim(aclen, ac, deztail, temp32b);
	temp64blen = fast_expansion_sum_zeroelim(temp32alen, temp32a,
		temp32blen, temp32b, temp64b);
	temp32alen = scale_expansion_zeroelim(cdlen, cd, aez, temp32a);
	temp32blen = scale_expansion_zeroelim(cdlen, cd, aeztail, temp32b);
	temp64clen = fast_expansion_sum_zeroelim(temp32alen, temp32a,
		temp32blen, temp32b, temp64c);
	temp128len = fast_expansion_sum_zeroelim(temp64alen, temp64a,
		temp64blen, temp64b, temp128);
	temp192len = fast_expansion_sum_zeroelim(temp64clen, temp64c,
		temp128len, temp128, temp192);
	xlen = scale_expansion_zeroelim(temp192len, temp192, bex, detx);
	xxlen = scale_expansion_zeroelim(xlen, detx, bex, detxx);
	xtlen = scale_expansion_zeroelim(temp192len, temp192, bextail, detxt);
	xxtlen = scale_expansion_zeroelim(xtlen, detxt, bex, detxxt);
	for (i = 0; i < xxtlen; i++) {
		detxxt[i] *= 2.0;
	}
	xtxtlen = scale_expansion_zeroelim(xtlen, detxt, bextail, detxtxt);
	x1len = fast_expansion_sum_zeroelim(xxlen, detxx, xxtlen, detxxt, x1);
	x2len = fast_expansion_sum_zeroelim(x1len, x1, xtxtlen, detxtxt, x2);
	ylen = scale_expansion_zeroelim(temp192len, temp192, bey, dety);
	yylen = scale_expansion_zeroelim(ylen, dety, bey, detyy);
	ytlen = scale_expansion_zeroelim(temp192len, temp192, beytail, detyt);
	yytlen = scale_expansion_zeroelim(ytlen, detyt, bey, detyyt);
	for (i = 0; i < yytlen; i++) {
		detyyt[i] *= 2.0;
	}
	ytytlen = scale_expansion_zeroelim(ytlen, detyt, beytail, detytyt);
	y1len = fast_expansion_sum_zeroelim(yylen, detyy, yytlen, detyyt, y1);
	y2len = fast_expansion_sum_zeroelim(y1len, y1, ytytlen, detytyt, y2);
	zlen = scale_expansion_zeroelim(temp192len, temp192, bez, detz);
	zzlen = scale_expansion_zeroelim(zlen, detz, bez, detzz);
	ztlen = scale_expansion_zeroelim(temp192len, temp192, beztail, detzt);
	zztlen = scale_expansion_zeroelim(ztlen, detzt, bez, detzzt);
	for (i = 0; i < zztlen; i++) {
		detzzt[i] *= 2.0;
	}
	ztztlen = scale_expansion_zeroelim(ztlen, detzt, beztail, detztzt);
	z1len = fast_expansion_sum_zeroelim(zzlen, detzz, zztlen, detzzt, z1);
	z2len = fast_expansion_sum_zeroelim(z1len, z1, ztztlen, detztzt, z2);
	xylen = fast_expansion_sum_zeroelim(x2len, x2, y2len, y2, detxy);
	blen = fast_expansion_sum_zeroelim(z2len, z2, xylen, detxy, bdet);

	temp32alen = scale_expansion_zeroelim(ablen, ab, -dez, temp32a);
	temp32blen = scale_expansion_zeroelim(ablen, ab, -deztail, temp32b);
	temp64alen = fast_expansion_sum_zeroelim(temp32alen, temp32a,
		temp32blen, temp32b, temp64a);
	temp32alen = scale_expansion_zeroelim(bdlen, bd, -aez, temp32a);
	temp32blen = scale_expansion_zeroelim(bdlen, bd, -aeztail, temp32b);
	temp64blen = fast_expansion_sum_zeroelim(temp32alen, temp32a,
		temp32blen, temp32b, temp64b);
	temp32alen = scale_expansion_zeroelim(dalen, da, -bez, temp32a);
	temp32blen = scale_expansion_zeroelim(dalen, da, -beztail, temp32b);
	temp64clen = fast_expansion_sum_zeroelim(temp32alen, temp32a,
		temp32blen, temp32b, temp64c);
	temp128len = fast_expansion_sum_zeroelim(temp64alen, temp64a,
		temp64blen, temp64b, temp128);
	temp192len = fast_expansion_sum_zeroelim(temp64clen, temp64c,
		temp128len, temp128, temp192);
	xlen = scale_expansion_zeroelim(temp192len, temp192, cex, detx);
	xxlen = scale_expansion_zeroelim(xlen, detx, cex, detxx);
	xtlen = scale_expansion_zeroelim(temp192len, temp192, cextail, detxt);
	xxtlen = scale_expansion_zeroelim(xtlen, detxt, cex, detxxt);
	for (i = 0; i < xxtlen; i++) {
		detxxt[i] *= 2.0;
	}
	xtxtlen = scale_expansion_zeroelim(xtlen, detxt, cextail, detxtxt);
	x1len = fast_expansion_sum_zeroelim(xxlen, detxx, xxtlen, detxxt, x1);
	x2len = fast_expansion_sum_zeroelim(x1len, x1, xtxtlen, detxtxt, x2);
	ylen = scale_expansion_zeroelim(temp192len, temp192, cey, dety);
	yylen = scale_expansion_zeroelim(ylen, dety, cey, detyy);
	ytlen = scale_expansion_zeroelim(temp192len, temp192, ceytail, detyt);
	yytlen = scale_expansion_zeroelim(ytlen, detyt, cey, detyyt);
	for (i = 0; i < yytlen; i++) {
		detyyt[i] *= 2.0;
	}
	ytytlen = scale_expansion_zeroelim(ytlen, detyt, ceytail, detytyt);
	y1len = fast_expansion_sum_zeroelim(yylen, detyy, yytlen, detyyt, y1);
	y2len = fast_expansion_sum_zeroelim(y1len, y1, ytytlen, detytyt, y2);
	zlen = scale_expansion_zeroelim(temp192len, temp192, cez, detz);
	zzlen = scale_expansion_zeroelim(zlen, detz, cez, detzz);
	ztlen = scale_expansion_zeroelim(temp192len, temp192, ceztail, detzt);
	zztlen = scale_expansion_zeroelim(ztlen, detzt, cez, detzzt);
	for (i = 0; i < zztlen; i++) {
		detzzt[i] *= 2.0;
	}
	ztztlen = scale_expansion_zeroelim(ztlen, detzt, ceztail, detztzt);
	z1len = fast_expansion_sum_zeroelim(zzlen, detzz, zztlen, detzzt, z1);
	z2len = fast_expansion_sum_zeroelim(z1len, z1, ztztlen, detztzt, z2);
	xylen = fast_expansion_sum_zeroelim(x2len, x2, y2len, y2, detxy);
	clen = fast_expansion_sum_zeroelim(z2len, z2, xylen, detxy, cdet);

	temp32alen = scale_expansion_zeroelim(bclen, bc, aez, temp32a);
	temp32blen = scale_expansion_zeroelim(bclen, bc, aeztail, temp32b);
	temp64alen = fast_expansion_sum_zeroelim(temp32alen, temp32a,
		temp32blen, temp32b, temp64a);
	temp32alen = scale_expansion_zeroelim(aclen, ac, -bez, temp32a);
	temp32blen = scale_expansion_zeroelim(aclen, ac, -beztail, temp32b);
	temp64blen = fast_expansion_sum_zeroelim(temp32alen, temp32a,
		temp32blen, temp32b, temp64b);
	temp32alen = scale_expansion_zeroelim(ablen, ab, cez, temp32a);
	temp32blen = scale_expansion_zeroelim(ablen, ab, ceztail, temp32b);
	temp64clen = fast_expansion_sum_zeroelim(temp32alen, temp32a,
		temp32blen, temp32b, temp64c);
	temp128len = fast_expansion_sum_zeroelim(temp64alen, temp64a,
		temp64blen, temp64b, temp128);
	temp192len = fast_expansion_sum_zeroelim(temp64clen, temp64c,
		temp128len, temp128, temp192);
	xlen = scale_expansion_zeroelim(temp192len, temp192, dex, detx);
	xxlen = scale_expansion_zeroelim(xlen, detx, dex, detxx);
	xtlen = scale_expansion_zeroelim(temp192len, temp192, dextail, detxt);
	xxtlen = scale_expansion_zeroelim(xtlen, detxt, dex, detxxt);
	for (i = 0; i < xxtlen; i++) {
		detxxt[i] *= 2.0;
	}
	xtxtlen = scale_expansion_zeroelim(xtlen, detxt, dextail, detxtxt);
	x1len = fast_expansion_sum_zeroelim(xxlen, detxx, xxtlen, detxxt, x1);
	x2len = fast_expansion_sum_zeroelim(x1len, x1, xtxtlen, detxtxt, x2);
	ylen = scale_expansion_zeroelim(temp192len, temp192, dey, dety);
	yylen = scale_expansion_zeroelim(ylen, dety, dey, detyy);
	ytlen = scale_expansion_zeroelim(temp192len, temp192, deytail, detyt);
	yytlen = scale_expansion_zeroelim(ytlen, detyt, dey, detyyt);
	for (i = 0; i < yytlen; i++) {
		detyyt[i] *= 2.0;
	}
	ytytlen = scale_expansion_zeroelim(ytlen, detyt, deytail, detytyt);
	y1len = fast_expansion_sum_zeroelim(yylen, detyy, yytlen, detyyt, y1);
	y2len = fast_expansion_sum_zeroelim(y1len, y1, ytytlen, detytyt, y2);
	zlen = scale_expansion_zeroelim(temp192len, temp192, dez, detz);
	zzlen = scale_expansion_zeroelim(zlen, detz, dez, detzz);
	ztlen = scale_expansion_zeroelim(temp192len, temp192, deztail, detzt);
	zztlen = scale_expansion_zeroelim(ztlen, detzt, dez, detzzt);
	for (i = 0; i < zztlen; i++) {
		detzzt[i] *= 2.0;
	}
	ztztlen = scale_expansion_zeroelim(ztlen, detzt, deztail, detztzt);
	z1len = fast_expansion_sum_zeroelim(zzlen, detzz, zztlen, detzzt, z1);
	z2len = fast_expansion_sum_zeroelim(z1len, z1, ztztlen, detztzt, z2);
	xylen = fast_expansion_sum_zeroelim(x2len, x2, y2len, y2, detxy);
	dlen = fast_expansion_sum_zeroelim(z2len, z2, xylen, detxy, ddet);

	ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, abdet);
	cdlen = fast_expansion_sum_zeroelim(clen, cdet, dlen, ddet, cddet);
	deterlen = fast_expansion_sum_zeroelim(ablen, abdet, cdlen, cddet, deter);

	return deter[deterlen - 1];
}

REAL insphereadapt(REAL *pa, REAL *pb, REAL *pc, REAL *pd, REAL *pe,
	REAL permanent)
{
	INEXACT REAL aex, bex, cex, dex, aey, bey, cey, dey, aez, bez, cez, dez;
	REAL det, errbound;

	INEXACT REAL aexbey1, bexaey1, bexcey1, cexbey1;
	INEXACT REAL cexdey1, dexcey1, dexaey1, aexdey1;
	INEXACT REAL aexcey1, cexaey1, bexdey1, dexbey1;
	REAL aexbey0, bexaey0, bexcey0, cexbey0;
	REAL cexdey0, dexcey0, dexaey0, aexdey0;
	REAL aexcey0, cexaey0, bexdey0, dexbey0;
	REAL ab[4], bc[4], cd[4], da[4], ac[4], bd[4];
	INEXACT REAL ab3, bc3, cd3, da3, ac3, bd3;
	REAL abeps, bceps, cdeps, daeps, aceps, bdeps;
	REAL temp8a[8], temp8b[8], temp8c[8], temp16[16], temp24[24], temp48[48];
	int temp8alen, temp8blen, temp8clen, temp16len, temp24len, temp48len;
	REAL xdet[96], ydet[96], zdet[96], xydet[192];
	int xlen, ylen, zlen, xylen;
	REAL adet[288], bdet[288], cdet[288], ddet[288];
	int alen, blen, clen, dlen;
	REAL abdet[576], cddet[576];
	int ablen, cdlen;
	REAL fin1[1152];
	int finlength;

	REAL aextail, bextail, cextail, dextail;
	REAL aeytail, beytail, ceytail, deytail;
	REAL aeztail, beztail, ceztail, deztail;

	INEXACT REAL bvirt;
	REAL avirt, bround, around;
	INEXACT REAL c;
	INEXACT REAL abig;
	REAL ahi, alo, bhi, blo;
	REAL err1, err2, err3;
	INEXACT REAL _i, _j;
	REAL _0;

	aex = (REAL)(pa[0] - pe[0]);
	bex = (REAL)(pb[0] - pe[0]);
	cex = (REAL)(pc[0] - pe[0]);
	dex = (REAL)(pd[0] - pe[0]);
	aey = (REAL)(pa[1] - pe[1]);
	bey = (REAL)(pb[1] - pe[1]);
	cey = (REAL)(pc[1] - pe[1]);
	dey = (REAL)(pd[1] - pe[1]);
	aez = (REAL)(pa[2] - pe[2]);
	bez = (REAL)(pb[2] - pe[2]);
	cez = (REAL)(pc[2] - pe[2]);
	dez = (REAL)(pd[2] - pe[2]);

	Two_Product(aex, bey, aexbey1, aexbey0);
	Two_Product(bex, aey, bexaey1, bexaey0);
	Two_Two_Diff(aexbey1, aexbey0, bexaey1, bexaey0, ab3, ab[2], ab[1], ab[0]);
	ab[3] = ab3;

	Two_Product(bex, cey, bexcey1, bexcey0);
	Two_Product(cex, bey, cexbey1, cexbey0);
	Two_Two_Diff(bexcey1, bexcey0, cexbey1, cexbey0, bc3, bc[2], bc[1], bc[0]);
	bc[3] = bc3;

	Two_Product(cex, dey, cexdey1, cexdey0);
	Two_Product(dex, cey, dexcey1, dexcey0);
	Two_Two_Diff(cexdey1, cexdey0, dexcey1, dexcey0, cd3, cd[2], cd[1], cd[0]);
	cd[3] = cd3;

	Two_Product(dex, aey, dexaey1, dexaey0);
	Two_Product(aex, dey, aexdey1, aexdey0);
	Two_Two_Diff(dexaey1, dexaey0, aexdey1, aexdey0, da3, da[2], da[1], da[0]);
	da[3] = da3;

	Two_Product(aex, cey, aexcey1, aexcey0);
	Two_Product(cex, aey, cexaey1, cexaey0);
	Two_Two_Diff(aexcey1, aexcey0, cexaey1, cexaey0, ac3, ac[2], ac[1], ac[0]);
	ac[3] = ac3;

	Two_Product(bex, dey, bexdey1, bexdey0);
	Two_Product(dex, bey, dexbey1, dexbey0);
	Two_Two_Diff(bexdey1, bexdey0, dexbey1, dexbey0, bd3, bd[2], bd[1], bd[0]);
	bd[3] = bd3;

	temp8alen = scale_expansion_zeroelim(4, cd, bez, temp8a);
	temp8blen = scale_expansion_zeroelim(4, bd, -cez, temp8b);
	temp8clen = scale_expansion_zeroelim(4, bc, dez, temp8c);
	temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a,
		temp8blen, temp8b, temp16);
	temp24len = fast_expansion_sum_zeroelim(temp8clen, temp8c,
		temp16len, temp16, temp24);
	temp48len = scale_expansion_zeroelim(temp24len, temp24, aex, temp48);
	xlen = scale_expansion_zeroelim(temp48len, temp48, -aex, xdet);
	temp48len = scale_expansion_zeroelim(temp24len, temp24, aey, temp48);
	ylen = scale_expansion_zeroelim(temp48len, temp48, -aey, ydet);
	temp48len = scale_expansion_zeroelim(temp24len, temp24, aez, temp48);
	zlen = scale_expansion_zeroelim(temp48len, temp48, -aez, zdet);
	xylen = fast_expansion_sum_zeroelim(xlen, xdet, ylen, ydet, xydet);
	alen = fast_expansion_sum_zeroelim(xylen, xydet, zlen, zdet, adet);

	temp8alen = scale_expansion_zeroelim(4, da, cez, temp8a);
	temp8blen = scale_expansion_zeroelim(4, ac, dez, temp8b);
	temp8clen = scale_expansion_zeroelim(4, cd, aez, temp8c);
	temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a,
		temp8blen, temp8b, temp16);
	temp24len = fast_expansion_sum_zeroelim(temp8clen, temp8c,
		temp16len, temp16, temp24);
	temp48len = scale_expansion_zeroelim(temp24len, temp24, bex, temp48);
	xlen = scale_expansion_zeroelim(temp48len, temp48, bex, xdet);
	temp48len = scale_expansion_zeroelim(temp24len, temp24, bey, temp48);
	ylen = scale_expansion_zeroelim(temp48len, temp48, bey, ydet);
	temp48len = scale_expansion_zeroelim(temp24len, temp24, bez, temp48);
	zlen = scale_expansion_zeroelim(temp48len, temp48, bez, zdet);
	xylen = fast_expansion_sum_zeroelim(xlen, xdet, ylen, ydet, xydet);
	blen = fast_expansion_sum_zeroelim(xylen, xydet, zlen, zdet, bdet);

	temp8alen = scale_expansion_zeroelim(4, ab, dez, temp8a);
	temp8blen = scale_expansion_zeroelim(4, bd, aez, temp8b);
	temp8clen = scale_expansion_zeroelim(4, da, bez, temp8c);
	temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a,
		temp8blen, temp8b, temp16);
	temp24len = fast_expansion_sum_zeroelim(temp8clen, temp8c,
		temp16len, temp16, temp24);
	temp48len = scale_expansion_zeroelim(temp24len, temp24, cex, temp48);
	xlen = scale_expansion_zeroelim(temp48len, temp48, -cex, xdet);
	temp48len = scale_expansion_zeroelim(temp24len, temp24, cey, temp48);
	ylen = scale_expansion_zeroelim(temp48len, temp48, -cey, ydet);
	temp48len = scale_expansion_zeroelim(temp24len, temp24, cez, temp48);
	zlen = scale_expansion_zeroelim(temp48len, temp48, -cez, zdet);
	xylen = fast_expansion_sum_zeroelim(xlen, xdet, ylen, ydet, xydet);
	clen = fast_expansion_sum_zeroelim(xylen, xydet, zlen, zdet, cdet);

	temp8alen = scale_expansion_zeroelim(4, bc, aez, temp8a);
	temp8blen = scale_expansion_zeroelim(4, ac, -bez, temp8b);
	temp8clen = scale_expansion_zeroelim(4, ab, cez, temp8c);
	temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a,
		temp8blen, temp8b, temp16);
	temp24len = fast_expansion_sum_zeroelim(temp8clen, temp8c,
		temp16len, temp16, temp24);
	temp48len = scale_expansion_zeroelim(temp24len, temp24, dex, temp48);
	xlen = scale_expansion_zeroelim(temp48len, temp48, dex, xdet);
	temp48len = scale_expansion_zeroelim(temp24len, temp24, dey, temp48);
	ylen = scale_expansion_zeroelim(temp48len, temp48, dey, ydet);
	temp48len = scale_expansion_zeroelim(temp24len, temp24, dez, temp48);
	zlen = scale_expansion_zeroelim(temp48len, temp48, dez, zdet);
	xylen = fast_expansion_sum_zeroelim(xlen, xdet, ylen, ydet, xydet);
	dlen = fast_expansion_sum_zeroelim(xylen, xydet, zlen, zdet, ddet);

	ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, abdet);
	cdlen = fast_expansion_sum_zeroelim(clen, cdet, dlen, ddet, cddet);
	finlength = fast_expansion_sum_zeroelim(ablen, abdet, cdlen, cddet, fin1);

	det = estimate(finlength, fin1);
	errbound = isperrboundB * permanent;
	if ((det >= errbound) || (-det >= errbound)) {
		return det;
	}

	Two_Diff_Tail(pa[0], pe[0], aex, aextail);
	Two_Diff_Tail(pa[1], pe[1], aey, aeytail);
	Two_Diff_Tail(pa[2], pe[2], aez, aeztail);
	Two_Diff_Tail(pb[0], pe[0], bex, bextail);
	Two_Diff_Tail(pb[1], pe[1], bey, beytail);
	Two_Diff_Tail(pb[2], pe[2], bez, beztail);
	Two_Diff_Tail(pc[0], pe[0], cex, cextail);
	Two_Diff_Tail(pc[1], pe[1], cey, ceytail);
	Two_Diff_Tail(pc[2], pe[2], cez, ceztail);
	Two_Diff_Tail(pd[0], pe[0], dex, dextail);
	Two_Diff_Tail(pd[1], pe[1], dey, deytail);
	Two_Diff_Tail(pd[2], pe[2], dez, deztail);
	if ((aextail == 0.0) && (aeytail == 0.0) && (aeztail == 0.0)
		&& (bextail == 0.0) && (beytail == 0.0) && (beztail == 0.0)
		&& (cextail == 0.0) && (ceytail == 0.0) && (ceztail == 0.0)
		&& (dextail == 0.0) && (deytail == 0.0) && (deztail == 0.0)) {
		return det;
	}

	errbound = isperrboundC * permanent + resulterrbound * Absolute(det);
	abeps = (aex * beytail + bey * aextail)
		- (aey * bextail + bex * aeytail);
	bceps = (bex * ceytail + cey * bextail)
		- (bey * cextail + cex * beytail);
	cdeps = (cex * deytail + dey * cextail)
		- (cey * dextail + dex * ceytail);
	daeps = (dex * aeytail + aey * dextail)
		- (dey * aextail + aex * deytail);
	aceps = (aex * ceytail + cey * aextail)
		- (aey * cextail + cex * aeytail);
	bdeps = (bex * deytail + dey * bextail)
		- (bey * dextail + dex * beytail);
	det += (((bex * bex + bey * bey + bez * bez)
		* ((cez * daeps + dez * aceps + aez * cdeps)
			+ (ceztail * da3 + deztail * ac3 + aeztail * cd3))
		+ (dex * dex + dey * dey + dez * dez)
		* ((aez * bceps - bez * aceps + cez * abeps)
			+ (aeztail * bc3 - beztail * ac3 + ceztail * ab3)))
		- ((aex * aex + aey * aey + aez * aez)
			* ((bez * cdeps - cez * bdeps + dez * bceps)
				+ (beztail * cd3 - ceztail * bd3 + deztail * bc3))
			+ (cex * cex + cey * cey + cez * cez)
			* ((dez * abeps + aez * bdeps + bez * daeps)
				+ (deztail * ab3 + aeztail * bd3 + beztail * da3))))
		+ 2.0 * (((bex * bextail + bey * beytail + bez * beztail)
			* (cez * da3 + dez * ac3 + aez * cd3)
			+ (dex * dextail + dey * deytail + dez * deztail)
			* (aez * bc3 - bez * ac3 + cez * ab3))
			- ((aex * aextail + aey * aeytail + aez * aeztail)
				* (bez * cd3 - cez * bd3 + dez * bc3)
				+ (cex * cextail + cey * ceytail + cez * ceztail)
				* (dez * ab3 + aez * bd3 + bez * da3)));
	if ((det >= errbound) || (-det >= errbound)) {
		return det;
	}

	return insphereexact(pa, pb, pc, pd, pe);
}

REAL insphere(REAL *pa, REAL *pb, REAL *pc, REAL *pd, REAL *pe)
{
	REAL aex, bex, cex, dex;
	REAL aey, bey, cey, dey;
	REAL aez, bez, cez, dez;
	REAL aexbey, bexaey, bexcey, cexbey, cexdey, dexcey, dexaey, aexdey;
	REAL aexcey, cexaey, bexdey, dexbey;
	REAL alift, blift, clift, dlift;
	REAL ab, bc, cd, da, ac, bd;
	REAL abc, bcd, cda, dab;
	REAL aezplus, bezplus, cezplus, dezplus;
	REAL aexbeyplus, bexaeyplus, bexceyplus, cexbeyplus;
	REAL cexdeyplus, dexceyplus, dexaeyplus, aexdeyplus;
	REAL aexceyplus, cexaeyplus, bexdeyplus, dexbeyplus;
	REAL det;
	REAL permanent, errbound;

	aex = pa[0] - pe[0];
	bex = pb[0] - pe[0];
	cex = pc[0] - pe[0];
	dex = pd[0] - pe[0];
	aey = pa[1] - pe[1];
	bey = pb[1] - pe[1];
	cey = pc[1] - pe[1];
	dey = pd[1] - pe[1];
	aez = pa[2] - pe[2];
	bez = pb[2] - pe[2];
	cez = pc[2] - pe[2];
	dez = pd[2] - pe[2];

	aexbey = aex * bey;
	bexaey = bex * aey;
	ab = aexbey - bexaey;
	bexcey = bex * cey;
	cexbey = cex * bey;
	bc = bexcey - cexbey;
	cexdey = cex * dey;
	dexcey = dex * cey;
	cd = cexdey - dexcey;
	dexaey = dex * aey;
	aexdey = aex * dey;
	da = dexaey - aexdey;

	aexcey = aex * cey;
	cexaey = cex * aey;
	ac = aexcey - cexaey;
	bexdey = bex * dey;
	dexbey = dex * bey;
	bd = bexdey - dexbey;

	abc = aez * bc - bez * ac + cez * ab;
	bcd = bez * cd - cez * bd + dez * bc;
	cda = cez * da + dez * ac + aez * cd;
	dab = dez * ab + aez * bd + bez * da;

	alift = aex * aex + aey * aey + aez * aez;
	blift = bex * bex + bey * bey + bez * bez;
	clift = cex * cex + cey * cey + cez * cez;
	dlift = dex * dex + dey * dey + dez * dez;

	det = (dlift * abc - clift * dab) + (blift * cda - alift * bcd);

	aezplus = Absolute(aez);
	bezplus = Absolute(bez);
	cezplus = Absolute(cez);
	dezplus = Absolute(dez);
	aexbeyplus = Absolute(aexbey);
	bexaeyplus = Absolute(bexaey);
	bexceyplus = Absolute(bexcey);
	cexbeyplus = Absolute(cexbey);
	cexdeyplus = Absolute(cexdey);
	dexceyplus = Absolute(dexcey);
	dexaeyplus = Absolute(dexaey);
	aexdeyplus = Absolute(aexdey);
	aexceyplus = Absolute(aexcey);
	cexaeyplus = Absolute(cexaey);
	bexdeyplus = Absolute(bexdey);
	dexbeyplus = Absolute(dexbey);
	permanent = ((cexdeyplus + dexceyplus) * bezplus
		+ (dexbeyplus + bexdeyplus) * cezplus
		+ (bexceyplus + cexbeyplus) * dezplus)
		* alift
		+ ((dexaeyplus + aexdeyplus) * cezplus
			+ (aexceyplus + cexaeyplus) * dezplus
			+ (cexdeyplus + dexceyplus) * aezplus)
		* blift
		+ ((aexbeyplus + bexaeyplus) * dezplus
			+ (bexdeyplus + dexbeyplus) * aezplus
			+ (dexaeyplus + aexdeyplus) * bezplus)
		* clift
		+ ((bexceyplus + cexbeyplus) * aezplus
			+ (cexaeyplus + aexceyplus) * bezplus
			+ (aexbeyplus + bexaeyplus) * cezplus)
		* dlift;
	errbound = isperrboundA * permanent;
	if ((det > errbound) || (-det > errbound)) {
		return det;
	}

	return insphereadapt(pa, pb, pc, pd, pe, permanent);
}


#ifdef __cplusplus
extern "C" {
#endif
	REAL orient2d_c(REAL *pa, REAL *pb, REAL *pc)
	{
		return orient2d(pa, pb, pc);
	}
#ifdef __cplusplus
}
#endif

#pragma endregion


/* this edge to edge test is based on Franlin Antonio's gem:
   "Faster Line Segment Intersection", in Graphics Gems III,
   pp. 199-202 **/
#define EDGE_EDGE_TEST(V0,U0,U1)                      \
  Bx=U0[i0]-U1[i0];                                   \
  By=U0[i1]-U1[i1];                                   \
  Cx=V0[i0]-U0[i0];                                   \
  Cy=V0[i1]-U0[i1];                                   \
  f=Ay*Bx-Ax*By;                                      \
  d=By*Cx-Bx*Cy;                                      \
  if((f>0 && d>=0 && d<=f) || (f<0 && d<=0 && d>=f))  \
  {                                                   \
    e=Ax*Cy-Ay*Cx;                                    \
    if(f>0)                                           \
    {                                                 \
      if(e>=0 && e<=f) return 1;                      \
    }                                                 \
    else                                              \
    {                                                 \
      if(e<=0 && e>=f) return 1;                      \
    }                                                 \
  }

#define EDGE_AGAINST_TRI_EDGES(V0,V1,U0,U1,U2) \
{                                              \
  double Ax,Ay,Bx,By,Cx,Cy,e,d,f;               \
  Ax=V1[i0]-V0[i0];                            \
  Ay=V1[i1]-V0[i1];                            \
  /* test edge U0,U1 against V0,V1 **/          \
  EDGE_EDGE_TEST(V0,U0,U1);                    \
  /* test edge U1,U2 against V0,V1 **/          \
  EDGE_EDGE_TEST(V0,U1,U2);                    \
  /* test edge U2,U1 against V0,V1 **/          \
  EDGE_EDGE_TEST(V0,U2,U0);                    \
}

#define POINT_IN_TRI(V0,U0,U1,U2)           \
{                                           \
  double a,b,c,d0,d1,d2;                     \
  /* is T1 completly inside T2? **/          \
  /* check if V0 is inside tri(U0,U1,U2) **/ \
  a=U1[i1]-U0[i1];                          \
  b=-(U1[i0]-U0[i0]);                       \
  c=-a*U0[i0]-b*U0[i1];                     \
  d0=a*V0[i0]+b*V0[i1]+c;                   \
                                            \
  a=U2[i1]-U1[i1];                          \
  b=-(U2[i0]-U1[i0]);                       \
  c=-a*U1[i0]-b*U1[i1];                     \
  d1=a*V0[i0]+b*V0[i1]+c;                   \
                                            \
  a=U0[i1]-U2[i1];                          \
  b=-(U0[i0]-U2[i0]);                       \
  c=-a*U2[i0]-b*U2[i1];                     \
  d2=a*V0[i0]+b*V0[i1]+c;                   \
  if(d0*d1>0.0)                             \
  {                                         \
    if(d0*d2>0.0) return 1;                 \
  }                                         \
}


void printpls(double* coord, int *elm, int m, int n)
{
	int i, idx0, idx1, pidx;
	FILE *fout = NULL;

	fout = fopen("ti.pls", "w");

	idx0 = m;
	idx1 = n;

	fprintf(fout, "%d %d 0 0 0 0\n", 2, 6);
	for (i = 0; i < 3; i++)
	{
		pidx = elm[4 * idx0 + i];
		fprintf(fout, "%d %20.15lf %20.15lf %20.15lf \n", i + 1, coord[pidx * 3 + 0], coord[pidx * 3 + 1], coord[pidx * 3 + 2]);
	}

	for (i = 0; i < 3; i++)
	{
		pidx = elm[4 * idx1 + i];
		fprintf(fout, "%d %20.15lf %20.15lf %20.15lf \n", i + 4, coord[pidx * 3 + 0], coord[pidx * 3 + 1], coord[pidx * 3 + 2]);
	}

	fprintf(fout, "1 %d %d %d %d\n", 1, 2, 3, 1);
	fprintf(fout, "2 %d %d %d %d\n", 4, 5, 6, 1);

	fclose(fout);
}

//title
#define REAL double
typedef REAL *point;

// Labels that signify the result of triangle-triangle intersection test.
enum interresult {
	DISJOINT = 0, INTERSECT, SHAREVERT, SHAREEDGE, SHAREFACE,
	TOUCHEDGE, TOUCHFACE, ACROSSVERT, ACROSSEDGE, ACROSSFACE,
	COLLISIONFACE, ACROSSSEG, ACROSSSUB
};
#define SETVECTOR3(V, a0, a1, a2) (V)[0] = (a0); (V)[1] = (a1); (V)[2] = (a2)
//start
void cross(REAL* v1, REAL* v2, REAL* n)
{
	n[0] = v1[1] * v2[2] - v2[1] * v1[2];
	n[1] = -(v1[0] * v2[2] - v2[0] * v1[2]);
	n[2] = v1[0] * v2[1] - v2[0] * v1[1];
}
REAL dot(REAL* v1, REAL* v2)
{
	return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}
void facenormal(double * pa, double * pb, double * pc, REAL *n, int pivot,
	REAL* lav)
{
	REAL v1[3], v2[3], v3[3], *pv1, *pv2;
	REAL L1, L2, L3;

	v1[0] = pb[0] - pa[0];  // edge vector v1: a->b
	v1[1] = pb[1] - pa[1];
	v1[2] = pb[2] - pa[2];
	v2[0] = pa[0] - pc[0];  // edge vector v2: c->a
	v2[1] = pa[1] - pc[1];
	v2[2] = pa[2] - pc[2];

	// Default, normal is calculated by: v1 x (-v2) (see Fig. fnormal).
	if (pivot > 0) {
		// Choose edge vectors by Burdakov's algorithm.
		v3[0] = pc[0] - pb[0];  // edge vector v3: b->c
		v3[1] = pc[1] - pb[1];
		v3[2] = pc[2] - pb[2];
		L1 = dot(v1, v1);
		L2 = dot(v2, v2);
		L3 = dot(v3, v3);
		// Sort the three edge lengths.
		if (L1 < L2) {
			if (L2 < L3) {
				pv1 = v1; pv2 = v2; // n = v1 x (-v2).
			}
			else {
				pv1 = v3; pv2 = v1; // n = v3 x (-v1).
			}
		}
		else {
			if (L1 < L3) {
				pv1 = v1; pv2 = v2; // n = v1 x (-v2).
			}
			else {
				pv1 = v2; pv2 = v3; // n = v2 x (-v3).
			}
		}
		if (lav) {
			// return the average edge length.
			*lav = (sqrt(L1) + sqrt(L2) + sqrt(L3)) / 3.0;
		}
	}
	else {
		pv1 = v1; pv2 = v2; // n = v1 x (-v2).
	}

	// Calculate the face normal.
	cross(pv1, pv2, n);
	// Inverse the direction;
	n[0] = -n[0];
	n[1] = -n[1];
	n[2] = -n[2];
}
REAL distance(REAL* p1, REAL* p2)
{
	return sqrt((p2[0] - p1[0]) * (p2[0] - p1[0]) +
		(p2[1] - p1[1]) * (p2[1] - p1[1]) +
		(p2[2] - p1[2]) * (p2[2] - p1[2]));
}
int tri_edge_2d(double* A, double* B, double* C, double* P, double* Q,
	double* R, int level, int *types, int *pos)

{
	point U[3], V[3];  // The permuted vectors of points.
	int pu[3], pv[3];  // The original positions of points.
	REAL abovept[3];
	REAL sA, sB, sC;
	REAL s1, s2, s3, s4;
	int z1;

	if (R == NULL) {
		// Calculate a lift point.
		if (1) {
			REAL n[3], len;
			// Calculate a lift point, saved in dummypoint.
			facenormal(A, B, C, n, 1, NULL);
			len = sqrt(dot(n, n));
			if (len != 0) {
				n[0] /= len;
				n[1] /= len;
				n[2] /= len;
				len = distance(A, B);
				len += distance(B, C);
				len += distance(C, A);
				len /= 3.0;
				R = abovept; //dummypoint;
				R[0] = A[0] + len * n[0];
				R[1] = A[1] + len * n[1];
				R[2] = A[2] + len * n[2];
			}
			else {
				// The triangle [A,B,C] is (nearly) degenerate, i.e., it is (close)
				//   to a line.  We need a line-line intersection test.
				//assert(0);
				// !!! A non-save return value.!!!
				return 0;  // DISJOINT
			}
		}
	}

	// Test A's, B's, and C's orientations wrt plane PQR.
	sA = orient3d(P, Q, R, A);
	sB = orient3d(P, Q, R, B);
	sC = orient3d(P, Q, R, C);


	if (sA < 0) {
		if (sB < 0) {
			if (sC < 0) { // (---).
				return 0;
			}
			else {
				if (sC > 0) { // (--+).
							  // All points are in the right positions.
					SETVECTOR3(U, A, B, C);  // I3
					SETVECTOR3(V, P, Q, R);  // I2
					SETVECTOR3(pu, 0, 1, 2);
					SETVECTOR3(pv, 0, 1, 2);
					z1 = 0;
				}
				else { // (--0).
					SETVECTOR3(U, A, B, C);  // I3
					SETVECTOR3(V, P, Q, R);  // I2
					SETVECTOR3(pu, 0, 1, 2);
					SETVECTOR3(pv, 0, 1, 2);
					z1 = 1;
				}
			}
		}
		else {
			if (sB > 0) {
				if (sC < 0) { // (-+-).
					SETVECTOR3(U, C, A, B);  // PT = ST
					SETVECTOR3(V, P, Q, R);  // I2
					SETVECTOR3(pu, 2, 0, 1);
					SETVECTOR3(pv, 0, 1, 2);
					z1 = 0;
				}
				else {
					if (sC > 0) { // (-++).
						SETVECTOR3(U, B, C, A);  // PT = ST x ST
						SETVECTOR3(V, Q, P, R);  // PL = SL
						SETVECTOR3(pu, 1, 2, 0);
						SETVECTOR3(pv, 1, 0, 2);
						z1 = 0;
					}
					else { // (-+0).
						SETVECTOR3(U, C, A, B);  // PT = ST
						SETVECTOR3(V, P, Q, R);  // I2
						SETVECTOR3(pu, 2, 0, 1);
						SETVECTOR3(pv, 0, 1, 2);
						z1 = 2;
					}
				}
			}
			else {
				if (sC < 0) { // (-0-).
					SETVECTOR3(U, C, A, B);  // PT = ST
					SETVECTOR3(V, P, Q, R);  // I2
					SETVECTOR3(pu, 2, 0, 1);
					SETVECTOR3(pv, 0, 1, 2);
					z1 = 1;
				}
				else {
					if (sC > 0) { // (-0+).
						SETVECTOR3(U, B, C, A);  // PT = ST x ST
						SETVECTOR3(V, Q, P, R);  // PL = SL
						SETVECTOR3(pu, 1, 2, 0);
						SETVECTOR3(pv, 1, 0, 2);
						z1 = 2;
					}
					else { // (-00).
						SETVECTOR3(U, B, C, A);  // PT = ST x ST
						SETVECTOR3(V, Q, P, R);  // PL = SL
						SETVECTOR3(pu, 1, 2, 0);
						SETVECTOR3(pv, 1, 0, 2);
						z1 = 3;
					}
				}
			}
		}
	}
	else {
		if (sA > 0) {
			if (sB < 0) {
				if (sC < 0) { // (+--).
					SETVECTOR3(U, B, C, A);  // PT = ST x ST
					SETVECTOR3(V, P, Q, R);  // I2
					SETVECTOR3(pu, 1, 2, 0);
					SETVECTOR3(pv, 0, 1, 2);
					z1 = 0;
				}
				else {
					if (sC > 0) { // (+-+).
						SETVECTOR3(U, C, A, B);  // PT = ST
						SETVECTOR3(V, Q, P, R);  // PL = SL
						SETVECTOR3(pu, 2, 0, 1);
						SETVECTOR3(pv, 1, 0, 2);
						z1 = 0;
					}
					else { // (+-0).
						SETVECTOR3(U, C, A, B);  // PT = ST
						SETVECTOR3(V, Q, P, R);  // PL = SL
						SETVECTOR3(pu, 2, 0, 1);
						SETVECTOR3(pv, 1, 0, 2);
						z1 = 2;
					}
				}
			}
			else {
				if (sB > 0) {
					if (sC < 0) { // (++-).
						SETVECTOR3(U, A, B, C);  // I3
						SETVECTOR3(V, Q, P, R);  // PL = SL
						SETVECTOR3(pu, 0, 1, 2);
						SETVECTOR3(pv, 1, 0, 2);
						z1 = 0;
					}
					else {
						if (sC > 0) { // (+++).
							return 0;
						}
						else { // (++0).
							SETVECTOR3(U, A, B, C);  // I3
							SETVECTOR3(V, Q, P, R);  // PL = SL
							SETVECTOR3(pu, 0, 1, 2);
							SETVECTOR3(pv, 1, 0, 2);
							z1 = 1;
						}
					}
				}
				else { // (+0#)
					if (sC < 0) { // (+0-).
						SETVECTOR3(U, B, C, A);  // PT = ST x ST
						SETVECTOR3(V, P, Q, R);  // I2
						SETVECTOR3(pu, 1, 2, 0);
						SETVECTOR3(pv, 0, 1, 2);
						z1 = 2;
					}
					else {
						if (sC > 0) { // (+0+).
							SETVECTOR3(U, C, A, B);  // PT = ST
							SETVECTOR3(V, Q, P, R);  // PL = SL
							SETVECTOR3(pu, 2, 0, 1);
							SETVECTOR3(pv, 1, 0, 2);
							z1 = 1;
						}
						else { // (+00).
							SETVECTOR3(U, B, C, A);  // PT = ST x ST
							SETVECTOR3(V, P, Q, R);  // I2
							SETVECTOR3(pu, 1, 2, 0);
							SETVECTOR3(pv, 0, 1, 2);
							z1 = 3;
						}
					}
				}
			}
		}
		else {
			if (sB < 0) {
				if (sC < 0) { // (0--).
					SETVECTOR3(U, B, C, A);  // PT = ST x ST
					SETVECTOR3(V, P, Q, R);  // I2
					SETVECTOR3(pu, 1, 2, 0);
					SETVECTOR3(pv, 0, 1, 2);
					z1 = 1;
				}
				else {
					if (sC > 0) { // (0-+).
						SETVECTOR3(U, A, B, C);  // I3
						SETVECTOR3(V, P, Q, R);  // I2
						SETVECTOR3(pu, 0, 1, 2);
						SETVECTOR3(pv, 0, 1, 2);
						z1 = 2;
					}
					else { // (0-0).
						SETVECTOR3(U, C, A, B);  // PT = ST
						SETVECTOR3(V, Q, P, R);  // PL = SL
						SETVECTOR3(pu, 2, 0, 1);
						SETVECTOR3(pv, 1, 0, 2);
						z1 = 3;
					}
				}
			}
			else {
				if (sB > 0) {
					if (sC < 0) { // (0+-).
						SETVECTOR3(U, A, B, C);  // I3
						SETVECTOR3(V, Q, P, R);  // PL = SL
						SETVECTOR3(pu, 0, 1, 2);
						SETVECTOR3(pv, 1, 0, 2);
						z1 = 2;
					}
					else {
						if (sC > 0) { // (0++).
							SETVECTOR3(U, B, C, A);  // PT = ST x ST
							SETVECTOR3(V, Q, P, R);  // PL = SL
							SETVECTOR3(pu, 1, 2, 0);
							SETVECTOR3(pv, 1, 0, 2);
							z1 = 1;
						}
						else { // (0+0).
							SETVECTOR3(U, C, A, B);  // PT = ST
							SETVECTOR3(V, P, Q, R);  // I2
							SETVECTOR3(pu, 2, 0, 1);
							SETVECTOR3(pv, 0, 1, 2);
							z1 = 3;
						}
					}
				}
				else { // (00#)
					if (sC < 0) { // (00-).
						SETVECTOR3(U, A, B, C);  // I3
						SETVECTOR3(V, Q, P, R);  // PL = SL
						SETVECTOR3(pu, 0, 1, 2);
						SETVECTOR3(pv, 1, 0, 2);
						z1 = 3;
					}
					else {
						if (sC > 0) { // (00+).
							SETVECTOR3(U, A, B, C);  // I3
							SETVECTOR3(V, P, Q, R);  // I2
							SETVECTOR3(pu, 0, 1, 2);
							SETVECTOR3(pv, 0, 1, 2);
							z1 = 3;
						}
						else { // (000)
							   // Not possible unless ABC is degenerate.
							   // Avoiding compiler warnings.
							SETVECTOR3(U, A, B, C);  // I3
							SETVECTOR3(V, P, Q, R);  // I2
							SETVECTOR3(pu, 0, 1, 2);
							SETVECTOR3(pv, 0, 1, 2);
							z1 = 4;
						}
					}
				}
			}
		}
	}

	s1 = orient3d(U[0], U[2], R, V[1]);  // A, C, R, Q
	s2 = orient3d(U[1], U[2], R, V[0]);  // B, C, R, P

	if (s1 > 0) {
		return 0;
	}
	if (s2 < 0) {
		return 0;
	}

	if (level == 0) {
		return 1;  // They are intersected.
	}
	//```assert(z1 != 4); // SELF_CHECK

	if (z1 == 1) {
		if (s1 == 0) {  // (0###)
						// C = Q.
			types[0] = (int)SHAREVERT;
			pos[0] = pu[2]; // C
			pos[1] = pv[1]; // Q
			types[1] = (int)DISJOINT;
		}
		else {
			if (s2 == 0) { // (#0##)
						   // C = P.
				types[0] = (int)SHAREVERT;
				pos[0] = pu[2]; // C
				pos[1] = pv[0]; // P
				types[1] = (int)DISJOINT;
			}
			else { // (-+##)
				   // C in [P, Q].
				types[0] = (int)ACROSSVERT;
				pos[0] = pu[2]; // C
				pos[1] = pv[0]; // [P, Q]
				types[1] = (int)DISJOINT;
			}
		}
		return 4;
	}

	s3 = orient3d(U[0], U[2], R, V[0]);  // A, C, R, P
	s4 = orient3d(U[1], U[2], R, V[1]);  // B, C, R, Q

	if (z1 == 0) {  // (tritri-03)
		if (s1 < 0) {
			if (s3 > 0) {
				assert(s2 > 0); // SELF_CHECK
				if (s4 > 0) {
					// [P, Q] overlaps [k, l] (-+++).
					types[0] = (int)ACROSSEDGE;
					pos[0] = pu[2]; // [C, A]
					pos[1] = pv[0]; // [P, Q]
					types[1] = (int)TOUCHFACE;
					pos[2] = 3;     // [A, B, C]
					pos[3] = pv[1]; // Q
				}
				else {
					if (s4 == 0) {
						// Q = l, [P, Q] contains [k, l] (-++0).
						types[0] = (int)ACROSSEDGE;
						pos[0] = pu[2]; // [C, A]
						pos[1] = pv[0]; // [P, Q]
						types[1] = (int)TOUCHEDGE;
						pos[2] = pu[1]; // [B, C]
						pos[3] = pv[1]; // Q
					}
					else { // s4 < 0
						   // [P, Q] contains [k, l] (-++-).
						types[0] = (int)ACROSSEDGE;
						pos[0] = pu[2]; // [C, A]
						pos[1] = pv[0]; // [P, Q]
						types[1] = (int)ACROSSEDGE;
						pos[2] = pu[1]; // [B, C]
						pos[3] = pv[0]; // [P, Q]
					}
				}
			}
			else {
				if (s3 == 0) {
					assert(s2 > 0); // SELF_CHECK
					if (s4 > 0) {
						// P = k, [P, Q] in [k, l] (-+0+).
						types[0] = (int)TOUCHEDGE;
						pos[0] = pu[2]; // [C, A]
						pos[1] = pv[0]; // P
						types[1] = (int)TOUCHFACE;
						pos[2] = 3;     // [A, B, C]
						pos[3] = pv[1]; // Q
					}
					else {
						if (s4 == 0) {
							// [P, Q] = [k, l] (-+00).
							types[0] = (int)TOUCHEDGE;
							pos[0] = pu[2]; // [C, A]
							pos[1] = pv[0]; // P
							types[1] = (int)TOUCHEDGE;
							pos[2] = pu[1]; // [B, C]
							pos[3] = pv[1]; // Q
						}
						else {
							// P = k, [P, Q] contains [k, l] (-+0-).
							types[0] = (int)TOUCHEDGE;
							pos[0] = pu[2]; // [C, A]
							pos[1] = pv[0]; // P
							types[1] = (int)ACROSSEDGE;
							pos[2] = pu[1]; // [B, C]
							pos[3] = pv[0]; // [P, Q]
						}
					}
				}
				else { // s3 < 0
					if (s2 > 0) {
						if (s4 > 0) {
							// [P, Q] in [k, l] (-+-+).
							types[0] = (int)TOUCHFACE;
							pos[0] = 3;     // [A, B, C]
							pos[1] = pv[0]; // P
							types[1] = (int)TOUCHFACE;
							pos[2] = 3;     // [A, B, C]
							pos[3] = pv[1]; // Q
						}
						else {
							if (s4 == 0) {
								// Q = l, [P, Q] in [k, l] (-+-0).
								types[0] = (int)TOUCHFACE;
								pos[0] = 3;     // [A, B, C]
								pos[1] = pv[0]; // P
								types[1] = (int)TOUCHEDGE;
								pos[2] = pu[1]; // [B, C]
								pos[3] = pv[1]; // Q
							}
							else { // s4 < 0
								   // [P, Q] overlaps [k, l] (-+--).
								types[0] = (int)TOUCHFACE;
								pos[0] = 3;     // [A, B, C]
								pos[1] = pv[0]; // P
								types[1] = (int)ACROSSEDGE;
								pos[2] = pu[1]; // [B, C]
								pos[3] = pv[0]; // [P, Q]
							}
						}
					}
					else { // s2 == 0
						   // P = l (#0##).
						types[0] = (int)TOUCHEDGE;
						pos[0] = pu[1]; // [B, C]
						pos[1] = pv[0]; // P
						types[1] = (int)DISJOINT;
					}
				}
			}
		}
		else { // s1 == 0
			   // Q = k (0####)
			types[0] = (int)TOUCHEDGE;
			pos[0] = pu[2]; // [C, A]
			pos[1] = pv[1]; // Q
			types[1] = (int)DISJOINT;
		}
	}
	else if (z1 == 2) {  // (tritri-23)
		if (s1 < 0) {
			if (s3 > 0) {
				assert(s2 > 0); // SELF_CHECK
				if (s4 > 0) {
					// [P, Q] overlaps [A, l] (-+++).
					types[0] = (int)ACROSSVERT;
					pos[0] = pu[0]; // A
					pos[1] = pv[0]; // [P, Q]
					types[1] = (int)TOUCHFACE;
					pos[2] = 3;     // [A, B, C]
					pos[3] = pv[1]; // Q
				}
				else {
					if (s4 == 0) {
						// Q = l, [P, Q] contains [A, l] (-++0).
						types[0] = (int)ACROSSVERT;
						pos[0] = pu[0]; // A
						pos[1] = pv[0]; // [P, Q]
						types[1] = (int)TOUCHEDGE;
						pos[2] = pu[1]; // [B, C]
						pos[3] = pv[1]; // Q
					}
					else { // s4 < 0
						   // [P, Q] contains [A, l] (-++-).
						types[0] = (int)ACROSSVERT;
						pos[0] = pu[0]; // A
						pos[1] = pv[0]; // [P, Q]
						types[1] = (int)ACROSSEDGE;
						pos[2] = pu[1]; // [B, C]
						pos[3] = pv[0]; // [P, Q]
					}
				}
			}
			else {
				if (s3 == 0) {
					assert(s2 > 0); // SELF_CHECK
					if (s4 > 0) {
						// P = A, [P, Q] in [A, l] (-+0+).
						types[0] = (int)SHAREVERT;
						pos[0] = pu[0]; // A
						pos[1] = pv[0]; // P
						types[1] = (int)TOUCHFACE;
						pos[2] = 3;     // [A, B, C]
						pos[3] = pv[1]; // Q
					}
					else {
						if (s4 == 0) {
							// [P, Q] = [A, l] (-+00).
							types[0] = (int)SHAREVERT;
							pos[0] = pu[0]; // A
							pos[1] = pv[0]; // P
							types[1] = (int)TOUCHEDGE;
							pos[2] = pu[1]; // [B, C]
							pos[3] = pv[1]; // Q
						}
						else { // s4 < 0
							   // Q = l, [P, Q] in [A, l] (-+0-).
							types[0] = (int)SHAREVERT;
							pos[0] = pu[0]; // A
							pos[1] = pv[0]; // P
							types[1] = (int)ACROSSEDGE;
							pos[2] = pu[1]; // [B, C]
							pos[3] = pv[0]; // [P, Q]
						}
					}
				}
				else { // s3 < 0
					if (s2 > 0) {
						if (s4 > 0) {
							// [P, Q] in [A, l] (-+-+).
							types[0] = (int)TOUCHFACE;
							pos[0] = 3;     // [A, B, C]
							pos[1] = pv[0]; // P
							types[0] = (int)TOUCHFACE;
							pos[0] = 3;     // [A, B, C]
							pos[1] = pv[1]; // Q
						}
						else {
							if (s4 == 0) {
								// Q = l, [P, Q] in [A, l] (-+-0).
								types[0] = (int)TOUCHFACE;
								pos[0] = 3;     // [A, B, C]
								pos[1] = pv[0]; // P
								types[0] = (int)TOUCHEDGE;
								pos[0] = pu[1]; // [B, C]
								pos[1] = pv[1]; // Q
							}
							else { // s4 < 0
								   // [P, Q] overlaps [A, l] (-+--).
								types[0] = (int)TOUCHFACE;
								pos[0] = 3;     // [A, B, C]
								pos[1] = pv[0]; // P
								types[0] = (int)ACROSSEDGE;
								pos[0] = pu[1]; // [B, C]
								pos[1] = pv[0]; // [P, Q]
							}
						}
					}
					else { // s2 == 0
						   // P = l (#0##).
						types[0] = (int)TOUCHEDGE;
						pos[0] = pu[1]; // [B, C]
						pos[1] = pv[0]; // P
						types[1] = (int)DISJOINT;
					}
				}
			}
		}
		else { // s1 == 0
			   // Q = A (0###).
			types[0] = (int)SHAREVERT;
			pos[0] = pu[0]; // A
			pos[1] = pv[1]; // Q
			types[1] = (int)DISJOINT;
		}
	}
	else if (z1 == 3) {  // (tritri-33)
		if (s1 < 0) {
			if (s3 > 0) {
				assert(s2 > 0); // SELF_CHECK
				if (s4 > 0) {
					// [P, Q] overlaps [A, B] (-+++).
					types[0] = (int)ACROSSVERT;
					pos[0] = pu[0]; // A
					pos[1] = pv[0]; // [P, Q]
					types[1] = (int)TOUCHEDGE;
					pos[2] = pu[0]; // [A, B]
					pos[3] = pv[1]; // Q
				}
				else {
					if (s4 == 0) {
						// Q = B, [P, Q] contains [A, B] (-++0).
						types[0] = (int)ACROSSVERT;
						pos[0] = pu[0]; // A
						pos[1] = pv[0]; // [P, Q]
						types[1] = (int)SHAREVERT;
						pos[2] = pu[1]; // B
						pos[3] = pv[1]; // Q
					}
					else { // s4 < 0
						   // [P, Q] contains [A, B] (-++-).
						types[0] = (int)ACROSSVERT;
						pos[0] = pu[0]; // A
						pos[1] = pv[0]; // [P, Q]
						types[1] = (int)ACROSSVERT;
						pos[2] = pu[1]; // B
						pos[3] = pv[0]; // [P, Q]
					}
				}
			}
			else {
				if (s3 == 0) {
					assert(s2 > 0); // SELF_CHECK
					if (s4 > 0) {
						// P = A, [P, Q] in [A, B] (-+0+).
						types[0] = (int)SHAREVERT;
						pos[0] = pu[0]; // A
						pos[1] = pv[0]; // P
						types[1] = (int)TOUCHEDGE;
						pos[2] = pu[0]; // [A, B]
						pos[3] = pv[1]; // Q
					}
					else {
						if (s4 == 0) {
							// [P, Q] = [A, B] (-+00).
							types[0] = (int)SHAREEDGE;
							pos[0] = pu[0]; // [A, B]
							pos[1] = pv[0]; // [P, Q]
							types[1] = (int)DISJOINT;
						}
						else { // s4 < 0
							   // P= A, [P, Q] in [A, B] (-+0-).
							types[0] = (int)SHAREVERT;
							pos[0] = pu[0]; // A
							pos[1] = pv[0]; // P
							types[1] = (int)ACROSSVERT;
							pos[2] = pu[1]; // B
							pos[3] = pv[0]; // [P, Q]
						}
					}
				}
				else { // s3 < 0
					if (s2 > 0) {
						if (s4 > 0) {
							// [P, Q] in [A, B] (-+-+).
							types[0] = (int)TOUCHEDGE;
							pos[0] = pu[0]; // [A, B]
							pos[1] = pv[0]; // P
							types[1] = (int)TOUCHEDGE;
							pos[2] = pu[0]; // [A, B]
							pos[3] = pv[1]; // Q
						}
						else {
							if (s4 == 0) {
								// Q = B, [P, Q] in [A, B] (-+-0).
								types[0] = (int)TOUCHEDGE;
								pos[0] = pu[0]; // [A, B]
								pos[1] = pv[0]; // P
								types[1] = (int)SHAREVERT;
								pos[2] = pu[1]; // B
								pos[3] = pv[1]; // Q
							}
							else { // s4 < 0
								   // [P, Q] overlaps [A, B] (-+--).
								types[0] = (int)TOUCHEDGE;
								pos[0] = pu[0]; // [A, B]
								pos[1] = pv[0]; // P
								types[1] = (int)ACROSSVERT;
								pos[2] = pu[1]; // B
								pos[3] = pv[0]; // [P, Q]
							}
						}
					}
					else { // s2 == 0
						   // P = B (#0##).
						types[0] = (int)SHAREVERT;
						pos[0] = pu[1]; // B
						pos[1] = pv[0]; // P
						types[1] = (int)DISJOINT;
					}
				}
			}
		}
		else { // s1 == 0
			   // Q = A (0###).
			types[0] = (int)SHAREVERT;
			pos[0] = pu[0]; // A
			pos[1] = pv[1]; // Q
			types[1] = (int)DISJOINT;
		}
	}

	return 4;
}
int tri_edge_tail(double *A, double *B, double *C, double *P, double *Q, double *R, double sP, double sQ, int level, int *types, int *pos)

{
	point U[3], V[3]; //, Ptmp;
	int pu[3], pv[3]; //, itmp;
	REAL s1, s2, s3;
	int z1;


	if (sP < 0) {
		if (sQ < 0) { // (--) disjoint
			return 0;
		}
		else {
			if (sQ > 0) { // (-+)
				SETVECTOR3(U, A, B, C);
				SETVECTOR3(V, P, Q, R);
				SETVECTOR3(pu, 0, 1, 2);
				SETVECTOR3(pv, 0, 1, 2);
				z1 = 0;
			}
			else { // (-0)
				SETVECTOR3(U, A, B, C);
				SETVECTOR3(V, P, Q, R);
				SETVECTOR3(pu, 0, 1, 2);
				SETVECTOR3(pv, 0, 1, 2);
				z1 = 1;
			}
		}
	}
	else {
		if (sP > 0) { // (+-)
			if (sQ < 0) {
				SETVECTOR3(U, A, B, C);
				SETVECTOR3(V, Q, P, R);  // P and Q are flipped.
				SETVECTOR3(pu, 0, 1, 2);
				SETVECTOR3(pv, 1, 0, 2);
				z1 = 0;
			}
			else {
				if (sQ > 0) { // (++) disjoint
					return 0;
				}
				else { // (+0)
					SETVECTOR3(U, B, A, C); // A and B are flipped.
					SETVECTOR3(V, P, Q, R);
					SETVECTOR3(pu, 1, 0, 2);
					SETVECTOR3(pv, 0, 1, 2);
					z1 = 1;
				}
			}
		}
		else { // sP == 0
			if (sQ < 0) { // (0-)
				SETVECTOR3(U, A, B, C);
				SETVECTOR3(V, Q, P, R);  // P and Q are flipped.
				SETVECTOR3(pu, 0, 1, 2);
				SETVECTOR3(pv, 1, 0, 2);
				z1 = 1;
			}
			else {
				if (sQ > 0) { // (0+)
					SETVECTOR3(U, B, A, C);  // A and B are flipped.
					SETVECTOR3(V, Q, P, R);  // P and Q are flipped.
					SETVECTOR3(pu, 1, 0, 2);
					SETVECTOR3(pv, 1, 0, 2);
					z1 = 1;
				}
				else { // (00)
					   // A, B, C, P, and Q are coplanar.
					z1 = 2;
				}
			}
		}
	}

	if (z1 == 2) {
		return tri_edge_2d(A, B, C, P, Q, R, level, types, pos);
	}

	s1 = orient3d(U[0], U[1], V[0], V[1]);
	if (s1 < 0) {
		return 0;
	}

	s2 = orient3d(U[1], U[2], V[0], V[1]);
	if (s2 < 0) {
		return 0;
	}

	s3 = orient3d(U[2], U[0], V[0], V[1]);
	if (s3 < 0) {
		return 0;
	}

	if (level == 0) {
		return 1;  // The are intersected.
	}

	types[1] = (int)DISJOINT; // No second intersection point.

	if (z1 == 0) {
		if (s1 > 0) {
			if (s2 > 0) {
				if (s3 > 0) { // (+++)
							  // [P, Q] passes interior of [A, B, C].
					types[0] = (int)ACROSSFACE;
					pos[0] = 3;  // interior of [A, B, C]
					pos[1] = 0;  // [P, Q]
				}
				else { // s3 == 0 (++0)
					   // [P, Q] intersects [C, A].
					types[0] = (int)ACROSSEDGE;
					pos[0] = pu[2];  // [C, A]
					pos[1] = 0;  // [P, Q]
				}
			}
			else { // s2 == 0
				if (s3 > 0) { // (+0+)
							  // [P, Q] intersects [B, C].
					types[0] = (int)ACROSSEDGE;
					pos[0] = pu[1];  // [B, C]
					pos[1] = 0;  // [P, Q]
				}
				else { // s3 == 0 (+00)
					   // [P, Q] passes C.
					types[0] = (int)ACROSSVERT;
					pos[0] = pu[2];  // C
					pos[1] = 0;  // [P, Q]
				}
			}
		}
		else { // s1 == 0
			if (s2 > 0) {
				if (s3 > 0) { // (0++)
							  // [P, Q] intersects [A, B].
					types[0] = (int)ACROSSEDGE;
					pos[0] = pu[0];  // [A, B]
					pos[1] = 0;  // [P, Q]
				}
				else { // s3 == 0 (0+0)
					   // [P, Q] passes A.
					types[0] = (int)ACROSSVERT;
					pos[0] = pu[0];  // A
					pos[1] = 0;  // [P, Q]
				}
			}
			else { // s2 == 0
				if (s3 > 0) { // (00+)
							  // [P, Q] passes B.
					types[0] = (int)ACROSSVERT;
					pos[0] = pu[1];  // B
					pos[1] = 0;  // [P, Q]
				}
				else { // s3 == 0 (000)
					   // Impossible.
					//assert(0);
				}
			}
		}
	}
	else { // z1 == 1
		if (s1 > 0) {
			if (s2 > 0) {
				if (s3 > 0) { // (+++)
							  // Q lies in [A, B, C].
					types[0] = (int)TOUCHFACE;
					pos[0] = 0; // [A, B, C]
					pos[1] = pv[1]; // Q
				}
				else { // s3 == 0 (++0)
					   // Q lies on [C, A].
					types[0] = (int)TOUCHEDGE;
					pos[0] = pu[2]; // [C, A]
					pos[1] = pv[1]; // Q
				}
			}
			else { // s2 == 0
				if (s3 > 0) { // (+0+)
							  // Q lies on [B, C].
					types[0] = (int)TOUCHEDGE;
					pos[0] = pu[1]; // [B, C]
					pos[1] = pv[1]; // Q
				}
				else { // s3 == 0 (+00)
					   // Q = C.
					types[0] = (int)SHAREVERT;
					pos[0] = pu[2]; // C
					pos[1] = pv[1]; // Q
				}
			}
		}
		else { // s1 == 0
			if (s2 > 0) {
				if (s3 > 0) { // (0++)
							  // Q lies on [A, B].
					types[0] = (int)TOUCHEDGE;
					pos[0] = pu[0]; // [A, B]
					pos[1] = pv[1]; // Q
				}
				else { // s3 == 0 (0+0)
					   // Q = A.
					types[0] = (int)SHAREVERT;
					pos[0] = pu[0]; // A
					pos[1] = pv[1]; // Q
				}
			}
			else { // s2 == 0
				if (s3 > 0) { // (00+)
							  // Q = B.
					types[0] = (int)SHAREVERT;
					pos[0] = pu[1]; // B
					pos[1] = pv[1]; // Q
				}
				else { // s3 == 0 (000)
					   // Impossible.
					//assert(0);
				}
			}
		}
	}

	// T and E intersect in a single point.
	return 2;
}
int tri_edge_inter_tail(double* A, double* B, double* C, double* P,
	double* Q, REAL s_p, REAL s_q)
{
	int types[2], pos[4];
	double inter_pos[3];
	int ni;  // =0, 2, 4

	//ni = tri_edge_tail(A, B, C, P, Q, NULL, s_p, s_q, 1, types, pos);

	double face[3][3];
	double line[2][3];
	for (int i = 0; i < 3; i++) {
		face[0][i] = A[i];
		face[1][i] = B[i];
		face[2][i] = C[i];

		line[0][i] = P[i];
		line[1][i] = Q[i];
	}
	//double linep[2][3], double facep[3][3], int *intTyp, int *intCod, double intPnt[3], bool bEpsilon)
	int intCod;
	//ni=GEOM_FUNC::lin_tri_intersect3d_check(line,face);
	ni = GEOM_FUNC::lin_tri_intersect3d(line, face, types, &intCod, inter_pos, false);
	return ni;
	//ni=tri_edge_inter_tail();

	// if(ni==2){
	//     /*check**/
	//     std::vector<std::pair<double*,double*>> check_pair;
	//     check_pair.push_back(std::make_pair(A,P));
	//     check_pair.push_back(std::make_pair(B,P));
	//     check_pair.push_back(std::make_pair(C,P));
	//     check_pair.push_back(std::make_pair(A,Q));
	//     check_pair.push_back(std::make_pair(B,Q));
	//     check_pair.push_back(std::make_pair(C,Q));
	//     bool true_value=false;
	//     for(auto j:check_pair){
	//         bool same=true;
	//         for(int i=0;i<3;i++){
	//             if(j.first[i]!=j.second[i]){
	//                 same=false;
	//             }
	//         }
	//         true_value=true_value||same;
	//     }
	//     if(!true_value){
	//         ni=0;
	//     }

	// }
	if (ni > 0) {
		if (ni == 2) {
			// Get the intersection type.
			if (types[0] == (int)SHAREVERT) {
				return (int)SHAREVERT;
			}
			else {
				return (int)INTERSECT;
			}
		}
		else if (ni == 4) {
			// There may be two intersections.
			if (types[0] == (int)SHAREVERT) {
				if (types[1] == (int)DISJOINT) {
					return (int)SHAREVERT;
				}
				else {
					assert(types[1] != (int)SHAREVERT);
					return (int)INTERSECT;
				}
			}
			else {
				if (types[0] == (int)SHAREEDGE) {
					return (int)SHAREEDGE;
				}
				else {
					return (int)INTERSECT;
				}
			}
		}
		else {
			assert(0);
		}
	}

	return (int)DISJOINT;

}
//#define newmethod
// A and O share
int share_one_node_tri_tri_inter(double *A, double *B, double*C, double*O, double *P, double*Q, int& i, int& j)
{
	REAL *P1=O, *P2=Q;
	if (j == 0)
	{
		P1 = P; 
	}
	else if (j == 2)
	{
		 P2 = P;
	}

	REAL s_p = orient3d(A, B, C, P1);
	REAL s_q = orient3d(A, B, C, P2);
	REAL *P3, *P4;
	if (s_p*s_q > 0.0) {
		return 0;
	}
	if (i == 0)
	{
		P3 = B; P4 = C;
	}
	else if (i == 1)
	{
		P3 = A; P4 = C;
	}
	else if (i == 2)
	{
		P3 = A; P4 = B;
	}

	REAL s_b = orient3d(O, P, Q, P3);
	REAL s_c = orient3d(O, P, Q, P4);
	if (s_b*s_c > 0.0) {
		return 0;
	}

		int abcop, opqab;
#ifndef newmethod
		abcop = tri_edge_inter_tail(A, B, C, P1, P2, s_p, s_q);
#else
		abcop = tri_inter_with_a_line_test(A, B, C, P1, P2);
#endif
		if (abcop) {
			//cout << "1#1212";
			return (int)INTERSECT;
		}
#ifndef newmethod
		opqab = tri_edge_inter_tail(A, B, C, P1, P2, s_p, s_q);
#else
		opqab = tri_inter_with_a_line_test(O, P, Q, P3, P4);
#endif
		if (opqab) {
			return (int)INTERSECT;
		}

		return 0; // DISJOINT;
	return 0; // DISJOINT;
}
bool pnt_project_in_tri(REAL* pa, REAL* pb, REAL* pc, REAL* pd) {
	

	REAL adx = pa[0] - pd[0];
	REAL bdx = pb[0] - pd[0];
	REAL cdx = pc[0] - pd[0];
	REAL ady = pa[1] - pd[1];
	REAL bdy = pb[1] - pd[1];
	REAL cdy = pc[1] - pd[1];
	REAL adz = pa[2] - pd[2];
	REAL bdz = pb[2] - pd[2];
	REAL cdz = pc[2] - pd[2];

	REAL pa_d[3];
	pa_d[0] = pa[0]+ bdy * cdz - bdz * cdy;
	pa_d[1] = pa[1]+ cdy * adz - cdz * ady;
	pa_d[2] = pa[2]+ ady * bdz - adz * bdy;
	REAL *pe=pa_d;
	REAL o_1 = orient3d(pa, pe, pb, pc);
	REAL o_2 = orient3d(pa, pe, pc, pd);
	REAL o_3 = orient3d(pa, pe, pd, pb);
	//return true;
	return o_1*o_2>=0&&o_2*o_3>=0&&o_1*o_3>=0;

}
REAL tri_area(REAL* A, REAL* B, REAL* C) {
	double area = -1;
	double dis;//三角形的高
	double side[3];//存储三条边的长度;

	side[0] = sqrt(pow(A[0] - B[0], 2) + pow(A[1] - B[1], 2) + pow(A[2] - B[2], 2));
	side[1] = sqrt(pow(A[0] - C[0], 2) + pow(A[1] - C[1], 2) + pow(A[2] - C[2], 2));
	side[2] = sqrt(pow(C[0] - B[0], 2) + pow(C[1] - B[1], 2) + pow(C[2] - B[2], 2));


	//利用海伦公式。s=sqr(p*(p-a)(p-b)(p-c)); 
	double p = (side[0] + side[1] + side[2]) / 2; //半周长;
	area = sqrt(p * (p - side[0]) * (p - side[1]) * (p - side[2]));

	return area;

}
int tri_tri_inter(REAL* A, REAL* B, REAL* C, REAL* O, REAL* P, REAL* Q)
{
	REAL s_o, s_p, s_q;
	REAL s_a, s_b, s_c;

	s_o = orient3d(A, B, C, O);
	s_p = orient3d(A, B, C, P);
	s_q = orient3d(A, B, C, Q);
	//cf.min_volumn_eps =1e-5;
	//if (abs(s_o) < 1) {


	//	bool b = pnt_project_in_tri(O, A, B, C);
	//	if (b) {
	//		if (abs(s_o) / tri_area(A, B, C) < cf.min_volumn_eps) {
	//			cout << abs(s_o) / tri_area(A, B, C) << "#" << endl;
	//			cout << b << endl;
	//		}
	//	}

	//}

	//if (abs(s_p) < 1) {


	//	bool b = pnt_project_in_tri(P, A, B, C);
	//	if (b) {
	//		if (abs(s_p) / tri_area(A, B, C) < cf.min_volumn_eps) {
	//			cout << abs(s_p) / tri_area(A, B, C) << "#" << endl;
	//			cout << b << endl;
	//		}
	//	}
	//}

	//if (abs(s_q) < 1) {


	//	bool b = pnt_project_in_tri(Q, A, B, C);
	//	if (b) {
	//		if (abs(s_q) / tri_area(A, B, C) < cf.min_volumn_eps) {
	//			cout << abs(s_q) / tri_area(A, B, C) << "#" << endl;
	//			cout << b << endl;
	//		}
	//	}
	//}


	if ((s_o * s_p > 0.0) && (s_o * s_q > 0.0)) {
		// o, p, q are all in the same halfspace of ABC.

		return 0; // DISJOINT;
	}

	s_a = orient3d(O, P, Q, A);
	s_b = orient3d(O, P, Q, B);
	s_c = orient3d(O, P, Q, C);
	if ((s_a * s_b > 0.0) && (s_a * s_c > 0.0)) {
		// a, b, c are all in the same halfspace of OPQ.
		return 0; // DISJOINT;
	}

	int abcop, abcpq, abcqo;
	int shareedge = 0;

	abcop = tri_edge_inter_tail(A, B, C, O, P, s_o, s_p);
	if (abcop == (int)INTERSECT) {
		return (int)INTERSECT;
	}
	else if (abcop == (int)SHAREEDGE) {
		shareedge++;
	}
	abcpq = tri_edge_inter_tail(A, B, C, P, Q, s_p, s_q);
	if (abcpq == (int)INTERSECT) {
		//cout << "2";
		return (int)INTERSECT;
	}
	else if (abcpq == (int)SHAREEDGE) {
		shareedge++;
	}
	abcqo = tri_edge_inter_tail(A, B, C, Q, O, s_q, s_o);
	if (abcqo == (int)INTERSECT) {
		//cout << "3";
		return (int)INTERSECT;
	}
	else if (abcqo == (int)SHAREEDGE) {
		//cout << "4";
		shareedge++;
	}
	if (shareedge == 3) {
		// opq are coincident with abc.
		//cout << "5";
		return (int)SHAREFACE;
	}

	// Continue to detect whether opq and abc are intersecting or not.
	int opqab, opqbc, opqca;

	opqab = tri_edge_inter_tail(O, P, Q, A, B, s_a, s_b);
	if (opqab == (int)INTERSECT) {
		//cout << "6";
		return (int)INTERSECT;
	}
	opqbc = tri_edge_inter_tail(O, P, Q, B, C, s_b, s_c);
	if (opqbc == (int)INTERSECT) {
		///cout << "7";
		return (int)INTERSECT;
	}
	opqca = tri_edge_inter_tail(O, P, Q, C, A, s_c, s_a);
	if (opqca == (int)INTERSECT) {
		//cout << "8";
		return (int)INTERSECT;
	}

	// At this point, two triangles are not intersecting and not coincident.
	//   They may be share an edge, or share a vertex, or disjoint.
	if (abcop == (int)SHAREEDGE) {
		//cout << "9";
		// op is coincident with an edge of abc.
		return (int)SHAREEDGE;
	}
	if (abcpq == (int)SHAREEDGE) {
		//cout << "10";
		// pq is coincident with an edge of abc.
		return (int)SHAREEDGE;
	}
	if (abcqo == (int)SHAREEDGE) {
		//cout << "11";
		// qo is coincident with an edge of abc.
		return (int)SHAREEDGE;
	}

	// They may share a vertex or disjoint.
	if (abcop == (int)SHAREVERT) {
		//cout << "12";
		return (int)SHAREVERT;
	}
	if (abcpq == (int)SHAREVERT) {
		//cout << "13";
		// q is the coincident vertex.
		return (int)SHAREVERT;
	}

	// They are disjoint.
	return (int)DISJOINT;
}




//////////////////**********************************************************/


#define TRI_NODE_NUM 3

bool tri_inter_with_a_line_test(double* A, double* B, double* C, double* P, double* Q) {

	BLVector start(P[0], P[1], P[2]);
	BLVector end(Q[0], Q[1], Q[2]);

	BLVector pos[3];
	pos[0] = BLVector(A[0], A[1], A[2]);
	pos[1] = BLVector(B[0], B[1], B[2]);
	pos[2] = BLVector(C[0], C[1], C[2]);
	BLVector e1 = pos[1] - pos[0];
	BLVector e2 = pos[2] - pos[0];
	BLVector q = BLVector::crossProduct(end - start, e2);
	double a = BLVector::dotProduct(q, e1);
	if (abs(a) < FLT_EPSILON) {
		return false;
	}
	double f = double(1.0) / a;
	BLVector s = start - pos[0];
	double u = f * BLVector::dotProduct(q, s);
	if (u < -FLT_EPSILON) {
		return false;
	}
	BLVector rp = BLVector::crossProduct(s, e1);
	double v = f * BLVector::dotProduct(rp, end - start);
	if (v < -FLT_EPSILON || u + v>1 + FLT_EPSILON) {
		return false;
	}
	double tt = f * BLVector::dotProduct(rp, e2);
	if (tt < -FLT_EPSILON || tt>1 + FLT_EPSILON) {
		return false;
	}
	//BLVector interPos = start + tt * (end-start);
	return true;


}


double tri_inter_with_a_line(std::array<int, 3>& triIdx1, MBLNode *pNods, BLVector start, BLVector end) {

	BLVector pos[3];
	for (int k = 0; k < 3; k++) {
		pos[k] = BLVector(pNods[triIdx1[k]].coord);
	}
	BLVector e1 = pos[1] - pos[0];
	BLVector e2 = pos[2] - pos[0];
	BLVector q = BLVector::crossProduct(end - start, e2);
	double a = BLVector::dotProduct(q, e1);
	if (abs(a) < double(1e-7)) {
		return 1.0;
	}
	double f = double(1.0) / a;
	BLVector s = start - pos[0];
	double u = f * BLVector::dotProduct(q, s);
	if (u<double(0.0)) {
		return 1.0;
	}
	BLVector rp = BLVector::crossProduct(s, e1);
	double v = f * BLVector::dotProduct(rp, end - start);
	if (v < -FLT_EPSILON || u + v>1 + FLT_EPSILON) {
		return 1.0;
	}
	double tt = f * BLVector::dotProduct(rp, e2);
	if (tt < 0) {
		return 1.0;
	}
	//BLVector interPos = start + tt * (end-start);
	return min(1.0, tt);


}
/*
* @note: 这个函数的特点是不会有包围盒求交
**/

int  Octree::tri_overlap_test_no_box(std::array<int, 3>& triIdx1, std::array<int, 3>& triIdx2, MBLNode *&pNods,int &share_node_num,int &share1,int& share2) {


	


}
/*
* @note：因为这个函数的效率直接决定了整个程序的效率，因此我对其进行了深度优化，特别是前面的部分，请勿轻易修改，除非你已经看懂--yhf
**/

int  Octree::tri_overlap_test(int triIdx1, int triIdx2, MBLNode *pNods, int *ele)
{
	int share_node_num = 0, i, j;
	//int k, ptIdx1, ptIdx2;
	//double p1[TRI_NODE_NUM][3], p2[TRI_NODE_NUM][3];
	//
	//double min1[3] = { std::numeric_limits<double>::max(),std::numeric_limits<double>::max(),std::numeric_limits<double>::max() }, max1[3] = { -std::numeric_limits<double>::max(),-std::numeric_limits<double>::max(),-std::numeric_limits<double>::max() }, min2[3] = { std::numeric_limits<double>::max(),std::numeric_limits<double>::max(),std::numeric_limits<double>::max() }, max2[3] = { -std::numeric_limits<double>::max(),-std::numeric_limits<double>::max(),-std::numeric_limits<double>::max() };
	int share1, share2;
	//for (i=0; i<TRI_NODE_NUM; i++)
	//{
	//	ptIdx1 = ele[triIdx1*(TRI_NODE_NUM+1) + i];
	//	ptIdx2 = ele[triIdx2*(TRI_NODE_NUM+1) + i];
	//	for (j=0; j<3; j++)
	//	{			
	//		p1[i][j] = pNods[ptIdx1].coord[j];
	//		p2[i][j] = pNods[ptIdx2].coord[j];
	//		if(p1[i][j] < min1[j])
	//			min1[j] = p1[i][j];
	//		if(p1[i][j] > max1[j])
	//			max1[j] = p1[i][j];
	//		if(p2[i][j] < min2[j])
	//			min2[j] = p2[i][j];
	//		if(p2[i][j] > max2[j])
	//			max2[j] = p2[i][j];
	//	}
	//}
	//
	//if( (min1[0] >= max2[0] || min1[1] >= max2[1] || min1[2] >= max2[2]) ||
	//	(max1[0] <= min2[0] || max1[1] <= min2[1] || max1[2] <= min2[2]))
	//	return false;
	//
	double* pt[6];
	pt[0] = pNods[ele[triIdx1*(TRI_NODE_NUM) + 0]].coord;
	pt[1] = pNods[ele[triIdx1*(TRI_NODE_NUM) + 1]].coord;
	pt[2] = pNods[ele[triIdx1*(TRI_NODE_NUM) + 2]].coord;

	pt[3] = pNods[ele[triIdx2*(TRI_NODE_NUM) + 0]].coord;
	pt[4] = pNods[ele[triIdx2*(TRI_NODE_NUM) + 1]].coord;
	pt[5] = pNods[ele[triIdx2*(TRI_NODE_NUM) + 2]].coord;
	double max_1, min_1,max_2,min_2;//不做初始化，减少时间
	for (i = 0; i < 3; i++) {

		if (pt[1][i] > pt[0][i]) {
			max_1 = pt[1][i]; min_1 = pt[0][i];
		}
		else {
			max_1 = pt[0][i]; min_1 = pt[1][i];
		}
		if (pt[2][i] > max_1) {
			max_1 = pt[2][i];
		}
		else if (pt[2][i] < min_1) {
			min_1 = pt[2][i];
		}


		if (pt[3][i] > pt[4][i]) {
			max_2 = pt[3][i]; min_2 = pt[4][i];
		}
		else {
			max_2 = pt[4][i]; min_2 = pt[3][i];
		}
		if (pt[5][i] > max_2) {
			max_2 = pt[5][i];
		}
		else if (pt[5][i] < min_2) {
			min_2 = pt[5][i];
		}

		if (max_2 > min_1&&max_1 > min_2)//这里不取等于是因为，等于要么共享两个点，要么共享一个点但不相交
			continue;
		return false;
	}

	for (i = 0; i < TRI_NODE_NUM; i++) {
		for (j = 0; j < 3; j++) {
			if (ele[triIdx1*(TRI_NODE_NUM) + i] == ele[triIdx2*(TRI_NODE_NUM) + j])
			{
				share_node_num++;
				share1 = i;
				share2 = j;
			}
		}
	}

	if (share_node_num > 1)
		return 0;
	if (!share_node_num)
	{
		return tri_tri_inter(pt[0], pt[1], pt[2], pt[3], pt[4], pt[5]);
	}
	return OCT::share_one_node_tri_tri_inter(pt[0], pt[1], pt[2], pt[3], pt[4], pt[5], share1, share2);
}

//bool tri_overlap_test_octree(TreeNode *node, Node *pNods, int *ele)
//{
//	int i, j, k;
//	TreeNode* pChild, *pFather;
//	bool bNoChild = true, ret = false;
//	if( node != NULL )
//	{
//		for (i=0; i<OC_SUBNODES; i++)
//		{
//			pChild = node->getChild(i);
//			
//			bNoChild = bNoChild && (pChild == NULL);
//		}
//
//		if ( !bNoChild )
//		{
//			for (i=0; i<OC_SUBNODES; i++)
//				ret |= tri_overlap_test_octree(node->getChild(i), pNods, ele);
//		}
//		else
//		{
//			for (auto i: node->getData())
//			{
//				for (auto j : node->getData())
//				{
//					//test overlap
//					if (j <= i)
//						continue;
//					if( tri_overlap_test(i, j, pNods, ele) )
//					{
//						printf("patch idx: %d %d\n",i, j);
//						ret = true;
//					}
//				}
//			}
//		}
//	}
//
//	return ret;
//}


TreeNode::TreeNode()
{
	layer = -1;
	vecData.clear();

	for (int i = 0; i < OC_SUBNODES; i++)
	{
		pChildNodes[i] = NULL;
	}

	pFather = NULL;
}

Eigen::Vector3d TreeNode::getPCA()
{
	Eigen::Vector3d ans;
	//

	//target_direction = BLVector(1, 0, 0);
	//return ans;



	//PCA
	Eigen::Matrix<double, 3, Eigen::Dynamic> M;
	M.resize(3, 3 * order_data.size());

	BLVector averagy(0, 0, 0);
	int count = 0;
	for (auto id_t : order_data) {
		auto i = id_t.id;
		BLVector centor(0, 0, 0);
		BLVector node[3];
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				centor[j] += (pOctreeAgent)->getcoord(i)[k][j];
				node[k][j] = (pOctreeAgent)->getcoord(i)[k][j];
			}
		}

		centor = centor / 3;

		averagy = averagy + centor;
		for (int k = 0; k < 3; k++) {
			for (int j = 0; j < 3; j++) {
				M(j, count) = centor[j] / 3 + 2.0 / 3 * node[k][j];
			}
			count++;
		}

	}
	averagy = averagy / order_data.size();
	for (int i = 0; i < M.cols(); i++) {
		for (int j = 0; j < 3; j++) {
			M(j, i) = M(j, i) - averagy[j];
		}
	}

	auto cov = (1.0 / count)*(M * M.transpose());



	//cout << "cove=" << endl;
	//cout << cov << endl;

	Eigen::EigenSolver<Eigen::Matrix3d> es(cov);
	auto D = es.eigenvalues();
	auto V = es.eigenvectors();

	//cout << D;

	Eigen::MatrixXd evalsReal;//注意这里定义的MatrixXd里没有c
	evalsReal = D.real();//获取特征值实数部分
	Eigen::MatrixXf::Index evalsMax;
	evalsReal.rowwise().sum().maxCoeff(&evalsMax);//得到最大特征值的位置

	static double d1 = 0, d2 = 0, d3 = 0;
	//d1 += D(0, 0).real();
	//d2 += D(1, 1).real();
	//d3 += D(2, 2).real();

	//if (rand() % 10 == 0) {
	//	cout << D << D(0, 0).real() << " " << D(1, 1).real() <<" "<<D(2,2).real()<< endl;
	//	cout << "================================" << endl;
	//	cout << V << endl;
	//	cout << d1 << " " << d2 << " " << d3 << endl;
	//}
	







	double sum = 0;
	//target_direction = BLVector(1, 0, 0);
	for (int j = 0; j < 3; j++) {

		target_direction[j] = V(j, evalsMax).real();
		ans[j] = V(j, evalsMax).real();
		//sum += V(0, j).real()*V(0, j).real();
	}

	//target_direction = BLVector(1.0, 0.0, 0.0);
	if (target_direction[0] < 0)
		target_direction = -1 * target_direction;
	//cout << target_direction;


	return ans;
}


TreeNode::TreeNode(const OCCUBE& cube, const std::vector<int>& data, TreeNode* father, int layer, int fnormal, OctreeAgent * pOctree):target_direction(1.0,0,0)
{
	num_obj = data.size();
	this->ocCube.lower = cube.lower;
	this->ocCube.upper = cube.upper;


	most_dispersed_direction = fnormal;
	pOctreeAgent = pOctree;
	
	order_data.reserve(data.size() * 1.5);





	for (auto i : data) {
		Tri t;
		//vecData.insert(i);
		auto coord = pOctreeAgent->getcoord(i);
		for (int j = 0; j < 3;j++) {
			for (int k = 0; k < 3; k++) {
				t.box[k] = std::max(t.box[k], coord[j][k]);
				t.box[3+k] = std::min(t.box[3+k], coord[j][k]);
			}
		}
		t.id = i;


		order_data.push_back(t);
	}
	//if (order_data.size() < 200&& order_data.size()>5)
	getPCA();

	for (auto &t : order_data) {
		auto coord = pOctreeAgent->getcoord(t.id);
		t.min_value = std::min(target_direction*coord[0], target_direction*coord[1]);
		t.max_value = std::max(target_direction*coord[0], target_direction*coord[1]);
		t.min_value = std::min(t.min_value, target_direction*coord[2]);
		t.max_value = std::max(t.max_value, target_direction*coord[2]);
	}

	sort(order_data.begin(), order_data.end());
	double cur_max = std::numeric_limits<double>::lowest();;
	for (auto i = order_data.begin(); i != order_data.end(); i++) {
		cur_max = std::max(cur_max, i->max_value);
		i->max_value = cur_max;
	}

	this->pFather = father;
	this->layer = layer;
	for (int i = 0; i < OC_SUBNODES; i++)
	{
		pChildNodes[i] = NULL;
	}
	this->changed = true;
}

TreeNode::~TreeNode()
{
	layer = -1;
	//vecData.clear();

	for (int i = 0; i < OC_SUBNODES; i++)
	{
		if (pChildNodes[i])
		{
			delete pChildNodes[i];
			pChildNodes[i] = NULL;
		}
	}

	pFather = NULL;
}



TreeNode* TreeNode::getFather()
{
	return pFather;
}

TreeNode*& TreeNode::getChild(int i)
{
	return pChildNodes[i];
}

void TreeNode::appendData(vector<int>& data)
{
	for (auto i : data) {
		Tri t;
		pOctreeAgent->getBox(i, t.box);

		const auto& coord = pOctreeAgent->getcoord(i);
		t.min_value = std::min(target_direction*coord[0], target_direction*coord[1]);
		t.max_value = std::max(target_direction*coord[0], target_direction*coord[1]);
		t.min_value = std::min(t.min_value, target_direction*coord[2]);
		t.max_value = std::max(t.max_value, target_direction*coord[2]);

		t.id = i;
		order_data.push_back(t);
	}

	sort(order_data.begin(), order_data.end());
	double cur_max = std::numeric_limits<double>::lowest();;
	for (auto i = order_data.begin(); i != order_data.end(); i++) {
		cur_max = std::max(cur_max, i->max_value);
		i->max_value = cur_max;
	}

}

void TreeNode::appendData_unsorted(int i, std::array<double, 6>& box)
{
	static int p = 0;

	Tri t(box,i);

	const auto& coord = pOctreeAgent->getcoord(i);

	t.min_value = std::min(target_direction*coord[0], target_direction*coord[1]);
	t.max_value = std::max(target_direction*coord[0], target_direction*coord[1]);
	t.min_value = std::min(t.min_value, target_direction*coord[2]);
	t.max_value = std::max(t.max_value, target_direction*coord[2]);

#ifdef BRUTE_FORCE_CHECK
	auto back = order_data;
#endif
	if (memory_queue.empty()) {
		order_data.push_back(t);
	}
	else {
		order_data[memory_queue.back()]=t;
		memory_queue.pop_back();
	}



}
void TreeNode::appendData(int i,std::array<double, 6>& box)
{
	static int p = 0;
	
	Tri t;
	t.box = box;
	t.id = i;

	auto coord = pOctreeAgent->getcoord(i);



	t.min_value = std::min(target_direction*coord[0], target_direction*coord[1]);
	t.max_value = std::max(target_direction*coord[0], target_direction*coord[1]);
	t.min_value = std::min(t.min_value, target_direction*coord[2]);
	t.max_value = std::max(t.max_value, target_direction*coord[2]);

#ifdef BRUTE_FORCE_CHECK
	auto back = order_data;
#endif
	int start = 0;
	int end = order_data.size();
	int size = end;
	int mid;

	//找到第一个大于它的值，相当于std::lower_bound;
	while (start < end) {
		mid = (start + end) / 2;
		if (order_data[mid] < t) {
			start = mid + 1;
		}
		else {
			end = mid;
		}
	}

	int j=end;

		//avoid duplicated element
		if (j < end&&pOctreeAgent->getElm()[order_data[j].id][0] == pOctreeAgent->getElm()[i][0] && \
			pOctreeAgent->getElm()[order_data[j].id][1] == pOctreeAgent->getElm()[i][1] && pOctreeAgent->getElm()[order_data[j].id][2] == pOctreeAgent->getElm()[i][2]) {
			order_data[j] = t;
			return;
		}
		order_data.insert(order_data.begin() + j, t);
		//if(order_Data.size()>200)
		//cout << order_Data.size()<<" ";

		if (j)
			if (order_data[j - 1].max_value > t.max_value) {
				order_data[j].max_value = order_data[j - 1].max_value;
				t.max_value = order_data[j - 1].max_value;
			}
		for (j++; j < size + 1; j++) {//这里插入了一个，size增加了1
			if (order_data[j].max_value < t.max_value)
				order_data[j].max_value = t.max_value;
			else
				break;
		}
#ifdef BRUTE_FORCE_CHECK
	for (short i = order_data.size() - 1; i > 0; i--) {
		if (order_data[i].min_value < order_data[i - 1].min_value \
			|| order_data[i].max_value < order_data[i - 1].max_value) {
			cout << "not ordered data!" << endl;
			for (auto i : back) {
				cout << i.min_value << " " << i.max_value << "###";
			}
			cout << endl;
			for (auto i : order_data) {
				cout << i.min_value << " " << i.max_value << "###";
			}
			cout << endl << " id=" << i;
			throw std::exception("not ordered data!");
		}
	}
#endif



}
void TreeNode::rmData_unsorted(int i)
{
	//auto coord = pOctreeAgent->getcoord(i, most_dispersed_direction);
	//double max_value = -1e20;
	//double min_value = 1e20;
	//for (auto j : coord) {
	//	max_value = std::max(max_value, j);
	//	min_value = std::min(min_value, j);
	//}
#ifdef BRUTE_FORCE_CHECK
	auto back = order_data;
#endif
	int start = 0;
	int end = order_data.size();
	int size = end;
	int mid;
	//找到第一个大于它的值，相当于std::lower_bound;

	for (start=size-1; start >= 0; start--) {
		if (order_data[start].id == i) {
			order_data[start].valid = false;
			memory_queue.push_back(start);
			break;
		}
	}



#ifdef BRUTE_FORCE_CHECK
	for (short i = order_data.size()-1; i > 0; i--) {
		if (order_data[i].min_value < order_data[i - 1].min_value \
			|| order_data[i].max_value < order_data[i - 1].max_value) {
			cout << "not ordered data!" << endl;
			for (auto i : back) {
				cout << i.min_value << " " << i.max_value << "###";
			}
			cout << endl;
			for (auto i : order_data) {
				cout << i.min_value << " " << i.max_value << "###";
			}
			cout << endl << " id=" << i;
			throw std::exception("not ordered data!");
		}
	}
#endif

	//while (it != order_Data.end()) {
	//	if () {
	//		order_Data.erase(it);
	//		//order_Data.
	//		return;
	//	}
	//	it++;
	//}
	//vecData.erase(i);



}
void TreeNode::rmData(int i)
{
	//auto coord = pOctreeAgent->getcoord(i, most_dispersed_direction);
	//double max_value = -1e20;
	//double min_value = 1e20;
	//for (auto j : coord) {
	//	max_value = std::max(max_value, j);
	//	min_value = std::min(min_value, j);
	//}
#ifdef BRUTE_FORCE_CHECK
	auto back = order_data;
#endif
	int start = 0;
	int end = order_data.size();
	int size = end;
	int mid;
	//找到第一个大于它的值，相当于std::lower_bound;

	for (start=size-1; start >= 0; start--) {
		if (order_data[start].id == i) {
			
			order_data.erase(order_data.begin()+start);
			if (order_data.size()*1.0 / order_data.capacity() < 0.4)
				order_data.shrink_to_fit();//删除冗余内存
			break;
		}
	}



#ifdef BRUTE_FORCE_CHECK
	for (short i = order_data.size()-1; i > 0; i--) {
		if (order_data[i].min_value < order_data[i - 1].min_value \
			|| order_data[i].max_value < order_data[i - 1].max_value) {
			cout << "not ordered data!" << endl;
			for (auto i : back) {
				cout << i.min_value << " " << i.max_value << "###";
			}
			cout << endl;
			for (auto i : order_data) {
				cout << i.min_value << " " << i.max_value << "###";
			}
			cout << endl << " id=" << i;
			throw std::exception("not ordered data!");
		}
	}
#endif

	//while (it != order_Data.end()) {
	//	if () {
	//		order_Data.erase(it);
	//		//order_Data.
	//		return;
	//	}
	//	it++;
	//}
	//vecData.erase(i);



}

void TreeNode::sortElement()
{

	sort(order_data.begin(), order_data.end());
	double cur_max = std::numeric_limits<double>::lowest();;
	for (auto i = order_data.begin(); i != order_data.end(); i++) {
		cur_max = std::max(cur_max, i->max_value);
		i->max_value = cur_max;
	}

}

Octree::Octree() :pRoot(nullptr)

#ifdef _DEBUG
,transfer_time_cost_(0),
num_inter(0),
call_times_(0),
transfer_times_(0)
#endif
{

	depth = -1;
	sorted = true;
	special_id = -1;
}

Octree::~Octree()
{
	if(pRoot)
	destroyNode(pRoot);
	if (node_before) {
		node_before = nullptr;
	}
	if (pRoot) {
		pRoot = nullptr;
	}

	//delete pOctreeAgent;
}

int Octree::buildOcTree(const OCCUBE& cube, const std::vector<int>& vecData, int maxObjNum)
{
	OCT::exactinit();

	max_depth_ = depth;

	std::vector<pair<double, int>> lengths = {std::make_pair((cube.upper - cube.lower).x,0),std::make_pair((cube.upper - cube.lower).y,1) ,std::make_pair((cube.upper - cube.lower).z,2) };
	sort(lengths.begin(), lengths.end());
	last_length_ = lengths[0].second;
	second_length_ = lengths[1].second;
	most_length_ = 3 - last_length_ - second_length_;

	createBranch(&pRoot, cube, vecData, 0, maxObjNum, max_depth_);
	node_before = pRoot;
	//trim(pRoot);

	return 0;
}

void Octree::rebuild()
{

	std::set<int> datas;
	auto leafnodes_bak = leaf_nodes;
	for (auto i : leafnodes_bak) {
		for (auto j : i->order_data) {
			datas.insert(j.id);
		}
	}
	cout << "Rebuild the octree." << endl;
	vector<int> datas_vec;
	for (auto i : datas) {
		datas_vec.push_back(i);
	}
	auto cube = pRoot->getCube();
	destroyNode(pRoot);
	leaf_nodes.clear();
	pRoot = nullptr;
	buildOcTree(cube,datas_vec,100);

}

int Octree::createBranch(TreeNode** node, const OCCUBE& cube, const std::vector<int>& vecData, int layer, int maxObjNum, int maxDepth)
{

	bool isleaf = true;
	if (maxDepth != 0 && vecData.size() > maxObjNum)
		isleaf = false;

	if (!(*node))
	{
		if(isleaf)
			*node = new TreeNode(cube, vecData, NULL, layer, most_length_, pOctreeAgent);
		else
			*node = new TreeNode(cube, std::vector<int>(), NULL, layer, most_length_, pOctreeAgent);
		if (layer > depth)
			depth = layer;
	}

	if (!isleaf)
	{
		std::vector<int> subdata[OC_SUBNODES ];
		OCCUBE subcube[OC_SUBNODES];

		pOctreeAgent->divideData(cube, vecData, subdata);

		(*node)->vecData.clear();
		(*node)->order_data.shrink_to_fit();

		//divide cube to eight sub-cubes
		CUBECOORD coorcenter, len;
		coorcenter = (cube.lower + cube.upper) / 2;
		len = cube.upper-coorcenter;

		//subcube[UNE]
		subcube[UNE].lower = coorcenter;
		subcube[UNE].upper = coorcenter + len;

		//subcube[UNW]
		subcube[UNW].lower = coorcenter;
		subcube[UNW].lower.x -= len.x;
		subcube[UNW].upper = coorcenter;
		subcube[UNW].upper.y += len.y;
		subcube[UNW].upper.z += len.z;

		//subcube[USW]
		subcube[USW].lower = coorcenter;
		subcube[USW].lower.x -= len.x;
		subcube[USW].lower.y -= len.y;
		subcube[USW].upper = coorcenter;
		subcube[USW].upper.z += len.z;

		//subcube[USE]
		subcube[USE].lower = coorcenter;
		subcube[USE].lower.y -= len.y;
		subcube[USE].upper = coorcenter;
		subcube[USE].upper.x += len.x;
		subcube[USE].upper.z += len.z;

		//subcube[LNE]
		subcube[LNE].lower = coorcenter;
		subcube[LNE].lower.z -= len.z;
		subcube[LNE].upper = coorcenter;
		subcube[LNE].upper.x += len.x;
		subcube[LNE].upper.y += len.y;

		//subcube[LNW]
		subcube[LNW].lower = coorcenter;
		subcube[LNW].lower.x -= len.x;
		subcube[LNW].lower.z -= len.z;
		subcube[LNW].upper = coorcenter;
		subcube[LNW].upper.y += len.y;

		//subcube[LSW]
		subcube[LSW].lower = coorcenter - len;
		subcube[LSW].upper = coorcenter;

		//subcube[LSE]
		subcube[LSE].lower = coorcenter;
		subcube[LSE].lower.y -= len.y;
		subcube[LSE].lower.z -= len.z;
		subcube[LSE].upper = coorcenter;
		subcube[LSE].upper.x += len.x;


		createBranch(&(*node)->pChildNodes[UNE], subcube[UNE], subdata[UNE], layer + 1, maxObjNum, maxDepth - 1);
		createBranch(&(*node)->pChildNodes[UNW], subcube[UNW], subdata[UNW], layer + 1, maxObjNum, maxDepth - 1);
		createBranch(&(*node)->pChildNodes[USW], subcube[USW], subdata[USW], layer + 1, maxObjNum, maxDepth - 1);
		createBranch(&(*node)->pChildNodes[USE], subcube[USE], subdata[USE], layer + 1, maxObjNum, maxDepth - 1);

		createBranch(&(*node)->pChildNodes[LNE], subcube[LNE], subdata[LNE], layer + 1, maxObjNum, maxDepth - 1);
		createBranch(&(*node)->pChildNodes[LNW], subcube[LNW], subdata[LNW], layer + 1, maxObjNum, maxDepth - 1);
		createBranch(&(*node)->pChildNodes[LSW], subcube[LSW], subdata[LSW], layer + 1, maxObjNum, maxDepth - 1);
		createBranch(&(*node)->pChildNodes[LSE], subcube[LSE], subdata[LSE], layer + 1, maxObjNum, maxDepth - 1);

#if 0
		spdlog::info("data size: {} layer: {}\n", (*node)->vecData.size(), (*node)->layer);
#endif

		for (int i = 0; i < OC_SUBNODES; i++)
		{
			(*node)->pChildNodes[i]->pFather = *node;

#if 0
			spdlog::info("data size: {} layer: {}\n", (*node)->pChildNodes[i]->vecData.size(), (*node)->pChildNodes[i]->layer);
#endif
		}
	}
	else {// is leaf node
		leaf_nodes.push_back(*node);
		//cout << func_normal << endl;;
	}

	return 1;
}

void Octree::traverseOcTree()
{
	int i, j, size = 0;
	bool bnl = false;
	std::queue<TreeNode*> quNodes;
	if (!pRoot)
		spdlog::info("The tree is empty!\n");

	printf("Layer: %d, data_size: %d\n\n", pRoot->layer, pRoot->vecData.size());
	size += pRoot->vecData.size();


#if 0
	for (j = 0; j < pRoot->vecData.size(); j++)
	{
		spdlog::info("{} ", pRoot->vecData[j]);
	}
	spdlog::info("\n");
#endif

	quNodes.push(pRoot);

	while (quNodes.size() > 0)
	{
		TreeNode* node = quNodes.front();
		quNodes.pop();

		//while ()
		{
			for (i = 0; i < OC_SUBNODES; i++)
			{
				if (node->pChildNodes[i])
				{
					bnl = true;
					printf("Layer: %d, data_size: %d\n", node->pChildNodes[i]->layer, node->pChildNodes[i]->vecData.size());

					size += node->pChildNodes[i]->vecData.size();
#if 0
					for (j = 0; j < node->pChildNodes[i]->vecData.size(); j++)
					{
						spdlog::info("{} ", node->pChildNodes[i]->vecData[j]);
					}
					spdlog::info("\n");
#endif

					quNodes.push(node->pChildNodes[i]);
				}
			}
			if (bnl)
			{
				spdlog::info("\n");
				bnl = false;
			}
		}
	}

	spdlog::info("total size: {}\n", size);
}
void Octree::printNodeSize()
{
	vector<int> p(100, 0);
	printNodeSize(getRootNode(),p);
	cout << endl;
	for (int i = 0; i < 40; i++)
		cout <<" #"<< i << ": " << p[i];
	cout << endl;
}
void Octree::printNodeSize(TreeNode* r,vector<int> & p) {
	if (r->getChild(0)) {
		for (int i = 0; i < 8; i++) {
			printNodeSize(r->getChild(i),p);
		}
	}
	else {
		p[r->order_data.size() / 100]++;
	}
}
void Octree::printElement(std::string filename)
{
	//return;
	queue<TreeNode *> list;
	std::map<int, int> index_map;
	std::map<int, int> point_map;
	vector<int> points;
	vector<int> elems;
	list.push(pRoot);
	while (!list.empty()) {
		TreeNode* node = list.front();
		list.pop();
		if (node->pChildNodes[UNE]) {
			list.push(node->pChildNodes[UNE]);
			list.push(node->pChildNodes[UNW]);
			list.push(node->pChildNodes[USW]);
			list.push(node->pChildNodes[USE]);
			list.push(node->pChildNodes[LNE]);
			list.push(node->pChildNodes[LNW]);
			list.push(node->pChildNodes[LSW]);
			list.push(node->pChildNodes[LSE]);
		}
		for (auto i : node->order_data) {
			if (index_map.find(i.id) == index_map.end()) {
				index_map[i.id] = index_map.size();
				elems.push_back(i.id);
			}
		}
	}
	for (auto i : index_map) {
		for (int j = 0; j < 3; j++) {
			if (point_map.find(pOctreeAgent->getElm()[i.first][j]) == point_map.end()) {
				point_map[pOctreeAgent->getElm()[i.first][j]] = point_map.size();
				points.push_back(pOctreeAgent->getElm()[i.first][j]);
			}
		}
	}
	std::ofstream fout(filename);
	fout << index_map.size() << " " << point_map.size() << " 0 0 0 0" << std::endl;
	for (int i = 0; i< points.size();i++) {
		fout << i + 1 << " " << pOctreeAgent->getNod()[points[i]].coord[0] << " " << pOctreeAgent->getNod()[points[i]].coord[1] << " " << pOctreeAgent->getNod()[points[i]].coord[2] << std::endl;
	}
	for (int i = 0; i < elems.size(); i++) {
		fout << i + 1 << " " << point_map[pOctreeAgent->getElm()[elems[i]][0]]+1 << " " << point_map[pOctreeAgent->getElm()[elems[i]][1]]+1 << " " << point_map[pOctreeAgent->getElm()[elems[i]][2]]+1<<" "<<1 << std::endl;
	}
	fout.close();
}
bool Octree::TriIn(int id) {
	queue<TreeNode *> list;
	list.push(pRoot);
	while (!list.empty()) {
		TreeNode* node = list.front();
		list.pop();
		if (node->pChildNodes[UNE]) {
			list.push(node->pChildNodes[UNE]);
			list.push(node->pChildNodes[UNW]);
			list.push(node->pChildNodes[USW]);
			list.push(node->pChildNodes[USE]);
			list.push(node->pChildNodes[LNE]);
			list.push(node->pChildNodes[LNW]);
			list.push(node->pChildNodes[LSW]);
			list.push(node->pChildNodes[LSE]);
		}
		else
		{
			for (auto i : node->order_data) {
				if (i.id == id)
					return true;
			}
		}
	
	}
	return false;
}
int Octree::getTriNum()
{
	queue<TreeNode *> list;
	std::set<int> index_map;
	list.push(pRoot);
	while (!list.empty()) {
		TreeNode* node = list.front();
		list.pop();
		if (node->pChildNodes[UNE]) {
			list.push(node->pChildNodes[UNE]);
			list.push(node->pChildNodes[UNW]);
			list.push(node->pChildNodes[USW]);
			list.push(node->pChildNodes[USE]);
			list.push(node->pChildNodes[LNE]);
			list.push(node->pChildNodes[LNW]);
			list.push(node->pChildNodes[LSW]);
			list.push(node->pChildNodes[LSE]);
		}
		else {
			for (auto i : node->order_data) {
				index_map.insert(i.id);
			}
		}
	}
	return index_map.size();
}
void Octree::saveOCTreeVTK(std::string filename)
{
	//return;
	queue<TreeNode *> list;
	vector<TreeNode*> leaf_node;

	map<std::array<double, 6>, int> point2idx;
	vector<std::array<double, 6>> points;

	list.push(pRoot);
	while (!list.empty()) {
		TreeNode* node = list.front();
		list.pop();
		if (node->pChildNodes[UNE]) {
			list.push(node->pChildNodes[UNE]);
			list.push(node->pChildNodes[UNW]);
			list.push(node->pChildNodes[USW]);
			list.push(node->pChildNodes[USE]);
			list.push(node->pChildNodes[LNE]);
			list.push(node->pChildNodes[LNW]);
			list.push(node->pChildNodes[LSW]);
			list.push(node->pChildNodes[LSE]);
		}
		else
			leaf_node.push_back(node);

	}
	std::vector<std::array<std::array<double, 6>, 8>> cubes(leaf_node.size());
	for (int i = 0; i < leaf_node.size(); i++) {
		auto& corner = cubes[i];

		for (int j = 0; j < 8; j++) {
			if (j < 4)
				corner[j][2] = leaf_node[i]->getCube().lower[2];
			else
				corner[j][2] = leaf_node[i]->getCube().upper[2];

			if (j % 4 < 2)
				corner[j][1] = leaf_node[i]->getCube().lower[1];
			else
				corner[j][1] = leaf_node[i]->getCube().upper[1];

			if (j % 2 < 1)
				corner[j][0] = leaf_node[i]->getCube().lower[0];
			else
				corner[j][0] = leaf_node[i]->getCube().upper[0];

		}

		for (int j = 0; j < 8; j++) {
			if (point2idx.find(corner[j]) == point2idx.end()) {
				point2idx[corner[j]] = points.size();
				points.push_back(corner[j]);
			}
		}
	}

	std::ofstream fout(filename);
	fout << "# vtk DataFile Version 3.0" << endl;
	fout << "Really cool data" << endl;
	fout << "ASCII" << endl;
	fout << "DATASET UNSTRUCTURED_GRID" << endl;
	fout << "POINTS " << point2idx.size() << " double" << endl;
	for (int i = 0; i < points.size(); i++) {
		fout << points[i][0] << " " << points[i][1] << " " << points[i][2] << endl;
	}
	fout << "CELLS " << leaf_node.size() << " " << leaf_node.size() * 9 << endl;;
	for (int i = 0; i < leaf_node.size(); i++) {
		fout << 8;
		for (int j = 0; j < 8; j++) {
			fout << " " << point2idx[cubes[i][j]];
		}
		fout << endl;
	}
	fout << "CELL_TYPES " << leaf_node.size() << endl;;
	for (int i = 0; i < leaf_node.size(); i++) {
		fout << 11 << endl;
	}
	fout.close();
}
void Octree::saveBiggestOctantVectorVTK(std::string filename)
{
	//return;
	queue<TreeNode *> list;
	vector<TreeNode*> leaf_node;

	map<std::array<double, 6>, int> point2idx;
	vector<std::array<double, 6>> points;
	vector<std::array<int, 3>> tris;

	list.push(pRoot);
	while (!list.empty()) {
		TreeNode* node = list.front();
		list.pop();
		if (node->pChildNodes[UNE]) {
			list.push(node->pChildNodes[UNE]);
			list.push(node->pChildNodes[UNW]);
			list.push(node->pChildNodes[USW]);
			list.push(node->pChildNodes[USE]);
			list.push(node->pChildNodes[LNE]);
			list.push(node->pChildNodes[LNW]);
			list.push(node->pChildNodes[LSW]);
			list.push(node->pChildNodes[LSE]);
		}
		else
			leaf_node.push_back(node);

	}


	for (int i = 0; i < leaf_node.size(); i++) {
		if (leaf_node[i]->order_data.size() > 60) {
			leaf_node[0] = leaf_node[i];
			leaf_node.resize(1);
			break;
		}
	}
	for (int i = 0; i < leaf_node[0]->order_data.size(); i++) {
		std::array<double, 6> p;
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				p[k] = pOctreeAgent->getNod()[pOctreeAgent->getElm()[leaf_node[0]->order_data[i].id][j]].coord[k];
			}
			points.push_back(p);

		}
		tris.push_back(std::array<int, 3>{3*i,3*i+1,3*i+2});
	}

	std::vector<std::array<std::array<double, 6>, 8>> cubes(leaf_node.size());
	for (int i = 0; i < leaf_node.size(); i++) {
		auto& corner = cubes[i];

		for (int j = 0; j < 8; j++) {
			if (j < 4)
				corner[j][2] = leaf_node[i]->getCube().lower[2];
			else
				corner[j][2] = leaf_node[i]->getCube().upper[2];

			if (j % 4 < 2)
				corner[j][1] = leaf_node[i]->getCube().lower[1];
			else
				corner[j][1] = leaf_node[i]->getCube().upper[1];

			if (j % 2 < 1)
				corner[j][0] = leaf_node[i]->getCube().lower[0];
			else
				corner[j][0] = leaf_node[i]->getCube().upper[0];

		}

		for (int j = 0; j < 8; j++) {
			if (point2idx.find(corner[j]) == point2idx.end()) {
				point2idx[corner[j]] = points.size();
				points.push_back(corner[j]);
			}
		}
	}

	std::ofstream fout(filename);
	fout << "# vtk DataFile Version 3.0" << endl;
	fout << "Really cool data" << endl;
	fout << "ASCII" << endl;
	fout << "DATASET UNSTRUCTURED_GRID" << endl;
	fout << "POINTS " << points.size() + leaf_node.size() << " double" << endl;
	for (int i = 0; i < points.size(); i++) {
		fout << points[i][0] << " " << points[i][1] << " " << points[i][2] << endl;
	}
	for (int i = 0; i < leaf_node.size(); i++) {
		fout << (leaf_node[i]->getCube().lower + leaf_node[i]->getCube().upper)[0] / 2 << " " << (leaf_node[i]->getCube().lower + leaf_node[i]->getCube().upper)[1] / 2 << " " << (leaf_node[i]->getCube().lower + leaf_node[i]->getCube().upper)[2] / 2 << endl;
	}
	fout << "CELLS " << leaf_node.size()+ tris.size() << " " << leaf_node.size() * 9+ tris.size()*4 << endl;;
	for (int i = 0; i < leaf_node.size(); i++) {
		fout << 8;
		for (int j = 0; j < 8; j++) {
			fout << " " << point2idx[cubes[i][j]];
		}
		fout << endl;
	}
	for (int i = 0; i < tris.size(); i++) {
		fout << 3;
		for (int j = 0; j < 3; j++) {
			fout << " " << tris[i][j];
		}
		fout << endl;
	}
	fout << "CELL_TYPES " << leaf_node.size()+tris.size() << endl;;
	for (int i = 0; i < leaf_node.size(); i++) {
		fout << 11 << endl;
	}
	for (int i = 0; i < tris.size(); i++) {
		fout << 5 << endl;
	}


	/* find target triangles **/
	
	fout<<"CELL_DATA " << leaf_node.size() + tris.size() << endl;
	fout << "SCALARS surface int 1" << endl;

	fout << "LOOKUP_TABLE default" << endl;
	for (int i = 0; i < leaf_node.size(); i++) {
		fout << 6 << endl;
	}
	int target = 29;
	for (int i = 0; i < tris.size(); i++) {
		if (i == target)
			fout << 10 << endl;
		else {

			if (	leaf_node[0]->order_data[i].min_value< leaf_node[0]->order_data[target].max_value\
				&&  leaf_node[0]->order_data[target].min_value < leaf_node[0]->order_data[i].max_value) {
				fout << 4 << endl;
			}
			else

				fout << 3 << endl;
		}
	}


	fout << "POINT_DATA " << leaf_node.size()+ points.size() << endl;;
	fout << "NORMALS vectors double" << endl;
	for (int i = 0; i < leaf_node.size() + points.size()-1; i++) {
		fout << 0 << " " << 0 << " " << 0 << endl;
	}
	fout << leaf_node[0]->target_direction[0] << " " << leaf_node[0]->target_direction[1] << " " << leaf_node[0]->target_direction[2] << endl;
	fout.close();
}
void Octree::saveOCTreeVectorVTK(std::string filename)
{
	//return;
	queue<TreeNode *> list;
	vector<TreeNode*> leaf_node;

	map<std::array<double, 6>,int> point2idx;
	vector<std::array<double, 6>> points;

	list.push(pRoot);
	while (!list.empty()) {
		TreeNode* node = list.front();
		list.pop();
		if (node->pChildNodes[UNE]) {
			list.push(node->pChildNodes[UNE]);
			list.push(node->pChildNodes[UNW]);
			list.push(node->pChildNodes[USW]);
			list.push(node->pChildNodes[USE]);
			list.push(node->pChildNodes[LNE]);
			list.push(node->pChildNodes[LNW]);
			list.push(node->pChildNodes[LSW]);
			list.push(node->pChildNodes[LSE]);
		}
		else
			leaf_node.push_back(node);
		
	}
	std::vector<std::array<std::array<double, 6>, 8>> cubes(leaf_node.size());
	for (int i = 0; i < leaf_node.size();i++) {
		auto& corner= cubes[i];

		for (int j = 0; j < 8; j++) {
			if(j<4)
				corner[j][2] = leaf_node[i]->getCube().lower[2];
			else
				corner[j][2] = leaf_node[i]->getCube().upper[2];

			if (j%4 < 2)
				corner[j][1] = leaf_node[i]->getCube().lower[1];
			else
				corner[j][1] = leaf_node[i]->getCube().upper[1];

			if (j % 2 < 1)
				corner[j][0] = leaf_node[i]->getCube().lower[0];
			else
				corner[j][0] = leaf_node[i]->getCube().upper[0];

		}

		for (int j = 0; j < 8; j++) {
			if (point2idx.find(corner[j]) == point2idx.end()) {
				point2idx[corner[j]] = points.size();
				points.push_back(corner[j]);
			}
		}
	}

	std::ofstream fout(filename);
	fout << "# vtk DataFile Version 3.0" << endl;
	fout << "Really cool data" << endl;
	fout << "ASCII" << endl;
	fout << "DATASET UNSTRUCTURED_GRID" << endl;
	fout << "POINTS "<< points.size()+leaf_node.size()<<" double" << endl;
	for (int i = 0; i < points.size(); i++) {
		fout << points[i][0] << " " << points[i][1] << " " << points[i][2] << endl;
	}
	for (int i = 0; i < leaf_node.size(); i++) {
		fout << (leaf_node[i]->getCube().lower + leaf_node[i]->getCube().upper)[0] / 2 << " " << (leaf_node[i]->getCube().lower + leaf_node[i]->getCube().upper)[1] / 2 << " "<<(leaf_node[i]->getCube().lower + leaf_node[i]->getCube().upper)[2] / 2 << endl;
	}
	fout << "CELLS " << leaf_node.size() << " " << leaf_node.size() * 9 << endl;;
	for (int i = 0; i < leaf_node.size(); i++) {
		fout << 8 ;
		for (int j = 0; j < 8; j++) {
			fout << " " << point2idx[cubes[i][j]];
		}
		fout << endl;
	}
	fout << "CELL_TYPES " << leaf_node.size() << endl;;
	for (int i = 0; i < leaf_node.size(); i++) {
		fout << 11<<endl;
	}

	fout << "POINT_DATA "<< leaf_node.size()+ points.size() << endl;;
	fout << "NORMALS vectors double" <<  endl;
	for (int i = 0; i < points.size(); i++) {
		fout << 0 << " " << 0 << " " << 0 << endl;
	}
	for (int i = 0; i < leaf_node.size(); i++) {
		leaf_node[i]->target_direction.normalize();
		if(leaf_node[i]->target_direction[0]>0.999)
			fout << 0 << " " << 0 << " " << 0 << endl;
		else
			fout << leaf_node[i]->target_direction[0] << " " << leaf_node[i]->target_direction[1] << " " << leaf_node[i]->target_direction[2] << endl;
	}
	fout.close();
}

void Octree::printTransferCost()
{
#ifdef _DEBUG
	cout << endl << "=====================";
	cout << "timecost = " << transfer_time_cost_ << endl;
	cout << "total call times=" << call_times_ << endl;
	cout << "total transfer times=" << transfer_times_ << endl;
	cout << "averagy transfer times" << transfer_times_ * 1.0 / call_times_ << endl;
#endif
}


int Octree::trim(TreeNode* node)
{
	TreeNode* pChild;
	bool bNoChild = true;
	int i, j;
	if (node != NULL)
	{
		for (i = 0; i < OC_SUBNODES; i++)
		{
			pChild = node->pChildNodes[i];

			bNoChild = bNoChild && (pChild == NULL);
		}

		if (bNoChild && (node->vecData.size() == 0))
		{
			for (j = 0; j < OC_SUBNODES; j++)
			{
				if (node == node->pFather->pChildNodes[j])
					break;
			}

			node->pFather->pChildNodes[j] = NULL;
			delete node;
		}
		else if (!bNoChild)
		{
			for (j = 0; j < OC_SUBNODES; j++)
				trim(node->pChildNodes[j]);
		}
	}
	else
	{
		return 0;
	}
}

void Octree::insert(int data, TreeNode *node, const OCCUBE cube, std::array<double, 6>& box)
{
	TreeNode* pChild;

	if (data < 0)
	{
		spdlog::info("error data.\n");
		return;
	}
	if (box[0] < cube.upper.x&&box[1] < cube.upper.y&&box[2] < cube.upper.z&&\
		box[3] > cube.lower.x&&box[4] > cube.lower.y&&box[5] > cube.lower.z) {
		node_before = node;
	}
	if (node->pChildNodes[0])//如果有孩子节点
	{


		//divide cube to eight sub-cubes
		CUBECOORD coorcenter((cube.lower.x + cube.upper.x) / 2, (cube.lower.y + cube.upper.y) / 2, (cube.lower.z + cube.upper.z) / 2);
		CUBECOORD len(cube.upper.x - coorcenter.x, cube.upper.y - coorcenter.y, cube.upper.z - coorcenter.z);
		OCCUBE subcube(coorcenter, cube.upper);
		bool in[6] = { box[0] > coorcenter.x,  box[3] < coorcenter.x,box[1] > coorcenter.y, box[4] < coorcenter.y,box[2] > coorcenter.z, box[5] < coorcenter.z };


		if (in[0] && in[2] && in[4]) {
			//subcube[UNE]
			insert(data, node->pChildNodes[UNE], subcube, box);
		}

		if (in[1] && in[2] && in[4]) {
			//subcube[UNW]
			subcube.lower.x -= len.x;
			subcube.upper.x -= len.x;
			insert(data, node->pChildNodes[UNW], subcube, box);
		}

		if (in[1] && in[3] && in[4]) {
			//subcube[USW]
			subcube.lower = coorcenter;
			subcube.lower.x -= len.x;
			subcube.lower.y -= len.y;
			subcube.upper = coorcenter;
			subcube.upper.z += len.z;
			insert(data, node->pChildNodes[USW], subcube, box);
		}

		if (in[0] && in[3] && in[4]) {
			//subcube[USE]
			subcube.lower = coorcenter;
			subcube.lower.y -= len.y;
			subcube.upper = coorcenter;
			subcube.upper.x += len.x;
			subcube.upper.z += len.z;
			insert(data, node->pChildNodes[USE], subcube, box);
		}

		if (in[0] && in[2] && in[5]) {
			//subcube[LNE]
			subcube.lower = coorcenter;
			subcube.lower.z -= len.z;
			subcube.upper = coorcenter;
			subcube.upper.x += len.x;
			subcube.upper.y += len.y;
			insert(data, node->pChildNodes[LNE], subcube, box);
		}

		if (in[1] && in[2] && in[5]) {
			//subcube[LNW]
			subcube.lower = coorcenter;
			subcube.lower.x -= len.x;
			subcube.lower.z -= len.z;
			subcube.upper = coorcenter;
			subcube.upper.y += len.y;
			insert(data, node->pChildNodes[LNW], subcube, box);
		}


		if (in[1] && in[3] && in[5]) {
			//subcube[LSW]
			subcube.lower = cube.lower;
			subcube.upper = coorcenter;
			insert(data, node->pChildNodes[LSW], subcube, box);
		}



		if (in[0] && in[3] && in[5]) {

			//subcube[LSE]
			subcube.lower = coorcenter;
			subcube.lower.y -= len.y;
			subcube.lower.z -= len.z;
			subcube.upper = coorcenter;
			subcube.upper.x += len.x;

			insert(data, node->pChildNodes[LSE], subcube, box);
		}


	}
	else
	{//如果是叶子节点

		node->changed = true;
		if(sort_telerant)
			node->appendData_unsorted(data, box);
		else
			node->appendData(data, box);

		

		//if (node->order_Data.size() > 1000)
		//	cout << node->order_Data.size() << " ";
	//不采用动态树，因为维护动态的代价可能要比收益要高

		//if ((node->order_data.size() > 800&&node->layer<15)) {
		//   
		//    vector<int> vec_data1;
		//	//cout << node->vecData.size() << " "<< node->layer<<endl;
		//	for (auto i : node->order_data)
		//	{	
		//		vec_data1.push_back(i.id);
		//	}
		//	node->num_obj = vec_data1.size();
		//    createBranch(&node, node->getCube(), vec_data1, node->layer+1,100,1);
		//    node->vecData.clear();
		//	node->order_data.clear();
		//}


	}

}
void Octree::insert(vector<int>& data,TreeNode* node)
{

	if(node==nullptr)
		node=pRoot;
	OCCUBE cube = node->getCube();

	if (node->pChildNodes[0])//如果有孩子节点
	{
		vector<int> datas_in_cell[8];
		pOctreeAgent->divideData(cube, data, datas_in_cell);
		for (int j = 0; j < OC_SUBNODES; j++)
		{
			insert(datas_in_cell[j], node->pChildNodes[j]);
		}
	}
	else {
		node->changed = true;
		node->appendData(data);
	}

}
void Octree::insertPreProcess(int data)
{
#ifdef _DEBUG
	double t1 = clock();
#endif
	TreeNode* tmp = node_before;
	OCCUBE cube = tmp->getCube();
	if (data < 0 || !tmp)
	{
		//cout << "warning,nonpositive triangle index!" << endl;
		return;
	}
	


	std::array<double, 6> box = { pOctreeAgent->pNod[pOctreeAgent->pEle[data ][0]].coord[0],\
		pOctreeAgent->pNod[pOctreeAgent->pEle[data][0]].coord[1] , \
		pOctreeAgent->pNod[pOctreeAgent->pEle[data][0]].coord[2] , \
		pOctreeAgent->pNod[pOctreeAgent->pEle[data][0]].coord[0] , \
		pOctreeAgent->pNod[pOctreeAgent->pEle[data][0]].coord[1] , \
		pOctreeAgent->pNod[pOctreeAgent->pEle[data][0]].coord[2] };
	for (int i = 0; i < 3; i++) {
		box[i] = max(box[i], pOctreeAgent->pNod[pOctreeAgent->pEle[data][1]].coord[i]);
		box[3 + i] = min(box[3 + i], pOctreeAgent->pNod[pOctreeAgent->pEle[data][1]].coord[i]);
		box[i] = max(box[i], pOctreeAgent->pNod[pOctreeAgent->pEle[data][2]].coord[i]);
		box[3 + i] = min(box[3 + i], pOctreeAgent->pNod[pOctreeAgent->pEle[data][2]].coord[i]);
	}
#ifdef _DEBUG
	int path = 0;
#endif
	while (tmp != pRoot) {
		if (box[0] < cube.upper.x&&box[1] < cube.upper.y&&box[2] < cube.upper.z&&\
			box[3] > cube.lower.x&&box[4] > cube.lower.y&&box[5] > cube.lower.z)
			break;
		tmp = tmp->getFather();

		cube = tmp->getCube();
#ifdef _DEBUG
		path++;
#endif
	}
	//cout << "path=" << path << "layer=" << node_before->layer << endl;
#ifdef _DEBUG
	double t2 = clock();
	transfer_time_cost_ += (t2 - t1)/CLOCKS_PER_SEC;
	transfer_times_+= path;
	call_times_++;
#endif

	insert(data, tmp, cube , box);
}
void Octree::rmDataPreProcess(int data) {
#ifdef _DEBUG
	double t1 = clock();
#endif
	TreeNode* tmp = node_before;
	OCCUBE cube = tmp->getCube();
	if (data < 0 || !tmp)
	{
		return;
	}

	std::array<double, 6> box = { pOctreeAgent->pNod[pOctreeAgent->pEle[data][0]].coord[0],\
	pOctreeAgent->pNod[pOctreeAgent->pEle[data][0]].coord[1] , \
	pOctreeAgent->pNod[pOctreeAgent->pEle[data][0]].coord[2] , \
	pOctreeAgent->pNod[pOctreeAgent->pEle[data][0]].coord[0] , \
	pOctreeAgent->pNod[pOctreeAgent->pEle[data][0]].coord[1] , \
	pOctreeAgent->pNod[pOctreeAgent->pEle[data][0]].coord[2] };
	for (int i = 0; i < 3; i++) {
		box[i] = max(box[i], pOctreeAgent->pNod[pOctreeAgent->pEle[data][1]].coord[i]);
		box[3 + i] = min(box[3 + i], pOctreeAgent->pNod[pOctreeAgent->pEle[data][1]].coord[i]);
		box[i] = max(box[i], pOctreeAgent->pNod[pOctreeAgent->pEle[data][2]].coord[i]);
		box[3 + i] = min(box[3 + i], pOctreeAgent->pNod[pOctreeAgent->pEle[data][2]].coord[i]);
	}
	int path = 0;
	while (tmp != pRoot) {
		if (box[0] < cube.upper.x&&box[1] < cube.upper.y&&box[2] < cube.upper.z&&\
			box[3] > cube.lower.x&&box[4] > cube.lower.y&&box[5] > cube.lower.z)
			break;
		tmp = tmp->getFather();

		cube = tmp->getCube();
		path++;
	}
#ifdef _DEBUG
	double t2 = clock();
	transfer_time_cost_ += (t2 - t1) / CLOCKS_PER_SEC;
	transfer_times_ += path;
	call_times_++;
#endif
	rmData(data, tmp, cube,box);
}
void Octree::sortNode(TreeNode *node)
{
	if (node == nullptr)
		node = pRoot;
	if (node->pChildNodes[0]) {
		for (int i = 0; i < OC_SUBNODES; i++)
			sortNode(node->pChildNodes[i]);
	}
	else {
		node->sortElement();
	}
}
void Octree::rmData(int data, TreeNode *node, const OCCUBE cube, std::array<double, 6>& box)
{
	if (box[0] < cube.upper.x&&box[1] < cube.upper.y&&box[2] < cube.upper.z&&\
		box[3] > cube.lower.x&&box[4] > cube.lower.y&&box[5] > cube.lower.z) {
		node_before = node;
	}
	TreeNode* pChild;


	if (node->pChildNodes[0])
	{

		//divide cube to eight sub-cubes
		CUBECOORD coorcenter((cube.lower.x + cube.upper.x) / 2, (cube.lower.y + cube.upper.y) / 2, (cube.lower.z + cube.upper.z) / 2);
		CUBECOORD len(cube.upper.x - coorcenter.x, cube.upper.y - coorcenter.y, cube.upper.z - coorcenter.z);
		OCCUBE subcube(coorcenter, cube.upper);
		bool in[6] = { box[0] > coorcenter.x,  box[3] < coorcenter.x,box[1] > coorcenter.y, box[4] < coorcenter.y,box[2] > coorcenter.z, box[5] < coorcenter.z };


		if (in[0] && in[2] && in[4]) {
			//subcube[UNE]
			rmData(data, node->pChildNodes[UNE], subcube, box);
		}

		if (in[1] && in[2] && in[4]) {
			//subcube[UNW]
			subcube.lower.x -= len.x;
			subcube.upper.x -= len.x;
			rmData(data, node->pChildNodes[UNW], subcube, box);
		}

		if (in[1] && in[3] && in[4]) {
			//subcube[USW]
			subcube.lower = coorcenter;
			subcube.lower.x -= len.x;
			subcube.lower.y -= len.y;
			subcube.upper = coorcenter;
			subcube.upper.z += len.z;
			rmData(data, node->pChildNodes[USW], subcube, box);
		}

		if (in[0] && in[3] && in[4]) {
			//subcube[USE]
			subcube.lower = coorcenter;
			subcube.lower.y -= len.y;
			subcube.upper = coorcenter;
			subcube.upper.x += len.x;
			subcube.upper.z += len.z;
			rmData(data, node->pChildNodes[USE], subcube, box);
		}

		if (in[0] && in[2] && in[5]) {
			//subcube[LNE]
			subcube.lower = coorcenter;
			subcube.lower.z -= len.z;
			subcube.upper = coorcenter;
			subcube.upper.x += len.x;
			subcube.upper.y += len.y;
			rmData(data, node->pChildNodes[LNE], subcube, box);
		}

		if (in[1] && in[2] && in[5]) {
			//subcube[LNW]
			subcube.lower = coorcenter;
			subcube.lower.x -= len.x;
			subcube.lower.z -= len.z;
			subcube.upper = coorcenter;
			subcube.upper.y += len.y;
			rmData(data, node->pChildNodes[LNW], subcube, box);
		}


		if (in[1] && in[3] && in[5]) {

			//subcube[LSW]
			subcube.lower = cube.lower;
			subcube.upper = coorcenter;


			rmData(data, node->pChildNodes[LSW], subcube, box);
		}



		if (in[0] && in[3] && in[5]) {

			//subcube[LSE]
			subcube.lower = coorcenter;
			subcube.lower.y -= len.y;
			subcube.lower.z -= len.z;
			subcube.upper = coorcenter;
			subcube.upper.x += len.x;

			rmData(data, node->pChildNodes[LSE], subcube, box);
		}
	}
	else
	{

		if (sort_telerant)
			node->rmData_unsorted(data);
		else
			node->rmData(data);
		
	}
}
bool Octree::chckIntersectPreProcess(int data) {

	if (!node_before)
		node_before = pRoot;


	TreeNode* tmp = node_before;
#ifdef _DEBUG
	double t1 = clock();
#endif

	OCCUBE cube = tmp->getCube();
	if (data < 0 || !tmp)
	{
		return false;
	}

	std::array<double, 6> box = { pOctreeAgent->pNod[pOctreeAgent->pEle[data][0]].coord[0],\
	pOctreeAgent->pNod[pOctreeAgent->pEle[data][0]].coord[1] , \
	pOctreeAgent->pNod[pOctreeAgent->pEle[data][0]].coord[2] , \
	pOctreeAgent->pNod[pOctreeAgent->pEle[data][0]].coord[0] , \
	pOctreeAgent->pNod[pOctreeAgent->pEle[data][0]].coord[1] , \
	pOctreeAgent->pNod[pOctreeAgent->pEle[data][0]].coord[2] };


	for (int i = 0; i < 3; i++) {
		box[i] = max(box[i], pOctreeAgent->pNod[pOctreeAgent->pEle[data][1]].coord[i]);
		box[3 + i] = min(box[3 + i], pOctreeAgent->pNod[pOctreeAgent->pEle[data][1]].coord[i]);
		box[i] = max(box[i], pOctreeAgent->pNod[pOctreeAgent->pEle[data][2]].coord[i]);
		box[3 + i] = min(box[3 + i], pOctreeAgent->pNod[pOctreeAgent->pEle[data][2]].coord[i]);
	}
	int path = 0;
	while (tmp != pRoot) {
		if (box[0] < cube.upper.x&&box[1] < cube.upper.y&&box[2] < cube.upper.z&&\
			box[3] > cube.lower.x&&box[4] > cube.lower.y&&box[5] > cube.lower.z)
			break;
		tmp = tmp->getFather();
		cube = tmp->getCube();
		path++;
	}
#ifdef _DEBUG
	double t2 = clock();
	transfer_time_cost_ += (t2 - t1) / CLOCKS_PER_SEC;
	transfer_times_ += path;
	call_times_++;
#endif

	return chckIntersect(data, tmp, cube,box);
}
void Octree::destroyNode(TreeNode* p)
{

	for (int i = 0; i < OC_SUBNODES; i++) {
		if (p->getChild(i)) {
			destroyNode(p->getChild(i));
			p->getChild(i) = nullptr;
		}
	}
	delete p;
}
//这里cube不作引用
bool Octree::chckIntersect(int data, TreeNode *node, const OCCUBE cube,std::array<double,6>& box)
{

	if (box[0] < cube.upper.x&&box[1] < cube.upper.y&&box[2] < cube.upper.z&&\
		box[3] > cube.lower.x&&box[4] > cube.lower.y&&box[5] > cube.lower.z) {
		node_before = node;
	}
	
	if (node->pChildNodes[0])//如果不是子节点
	{
		

		//divide cube to eight sub-cubes
		CUBECOORD coorcenter((cube.lower.x + cube.upper.x) / 2, (cube.lower.y + cube.upper.y) / 2, (cube.lower.z + cube.upper.z) / 2);
		CUBECOORD len(cube.upper.x - coorcenter.x, cube.upper.y - coorcenter.y,cube.upper.z - coorcenter.z);
		OCCUBE subcube(coorcenter, cube.upper);
		bool in[6] = { box[0] > coorcenter.x,  box[3] < coorcenter.x,box[1] > coorcenter.y, box[4] < coorcenter.y,box[2] > coorcenter.z, box[5] < coorcenter.z };


		if (in[0] && in[2] && in[4]) {
			//subcube[UNE]
			if (chckIntersect(data, node->pChildNodes[UNE], subcube, box))
				return true;
		}


		if (in[1] && in[2] && in[4]) {
			//subcube[UNW]
			subcube.lower.x -= len.x;
			subcube.upper.x -= len.x;
			if (chckIntersect(data, node->pChildNodes[UNW], subcube, box))
				return true;
		}

		if (in[1] && in[3] && in[4]) {
			//subcube[USW]
			subcube.lower = coorcenter;
			subcube.lower.x -= len.x;
			subcube.lower.y -= len.y;
			subcube.upper = coorcenter;
			subcube.upper.z += len.z;
			if (chckIntersect(data, node->pChildNodes[USW], subcube, box))
				return true;
		}

		if (in[0] && in[3] && in[4]) {
			//subcube[USE]
			subcube.lower = coorcenter;
			subcube.lower.y -= len.y;
			subcube.upper = coorcenter;
			subcube.upper.x += len.x;
			subcube.upper.z += len.z;
			if (chckIntersect(data, node->pChildNodes[USE], subcube, box))
				return true;
		}

		if (in[0] && in[2] && in[5]) {
			//subcube[LNE]
			subcube.lower = coorcenter;
			subcube.lower.z -= len.z;
			subcube.upper = coorcenter;
			subcube.upper.x += len.x;
			subcube.upper.y += len.y;
			if (chckIntersect(data, node->pChildNodes[LNE], subcube, box))
				return true;
		}

		if (in[1] && in[2] && in[5]) {
			//subcube[LNW]
			subcube.lower = coorcenter;
			subcube.lower.x -= len.x;
			subcube.lower.z -= len.z;
			subcube.upper = coorcenter;
			subcube.upper.y += len.y;
			if (chckIntersect(data, node->pChildNodes[LNW], subcube, box))
				return true;
		}


		if (in[1] && in[3] && in[5]) {

			//subcube[LSW]
			subcube.lower = cube.lower;
			subcube.upper = coorcenter;


			if (chckIntersect(data, node->pChildNodes[LSW], subcube, box))
				return true;
		}



		if (in[0] && in[3] && in[5]) {

			//subcube[LSE]
			subcube.lower = coorcenter;
			subcube.lower.y -= len.y;
			subcube.lower.z -= len.z;
			subcube.upper = coorcenter;
			subcube.upper.x += len.x;

			if (chckIntersect(data, node->pChildNodes[LSE], subcube, box))
				return true;
		}
	}
	else//叶子节点
	{

		
		

		/*重写,这段代码做了很多深度优化，请勿随意修改，除非你能看懂--yhf2020**/
		const auto &cell_data = node->order_data;

		int start = 0;
		int end = cell_data.size();
		int mid;
		static std::size_t maxsize = 0;


		auto coord = pOctreeAgent->getcoord(data);

		double min_value = std::numeric_limits<double>::max();
		double max_value = std::numeric_limits<double>::lowest();
		for (int j = 0; j < 3; j++) {
			min_value = std::min(min_value, node->target_direction*coord[j]);
			max_value = std::max(max_value, node->target_direction*coord[j]);
		}
		//找到第一个大于它的值，相当于std::lower_bound;
		while (start < end) {
			mid = (start + end) / 2;
			if (cell_data[mid].min_value <= max_value) {
				start = mid + 1;
			}
			else {
				end = mid;
			}
		}
		int it_end = end;
		if (!end)
			return false;
		start = 0;
		while (start < end) {
			mid = (start + end) / 2;
			if (cell_data[mid].max_value >= min_value) {

				end = mid - 1;
			}
			else {
				if (start == mid) {
					break;

				}
				start = mid;
			}
		}

		const auto &triIdx1 = pOctreeAgent->pEle[data].indexs;

#ifdef BRUTE_FORCE_CHECK
		for (short i = cell_data.size() - 1; i > 0; i--) {
			if (cell_data[i].min_value < cell_data[i - 1].min_value \
				|| cell_data[i].max_value < cell_data[i - 1].max_value) {
				cout << "not ordered data!";
				throw std::exception("not ordered data!");
			}
		}
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				if (pOctreeAgent->pNod[triIdx1[i]].coord[j] > box[j] ||
					pOctreeAgent->pNod[triIdx1[i]].coord[j] < box[j + 3]) {
					cout << "box error!";
					throw std::exception("box error!");
				}
			}
		}


#endif
		//static int d = 0;
		//static int c= 0;
		//c+= (cell_data.size());
		//d+= (it_end - start) ;
		//if(rand()%10000==0)
		//	cout << d*100.0/c << "%" << endl;
		//for (short i = cell_data.size() - 1; i >= 0; i--) {
		//if (data == 77749116) {
		//	
		//	for (int i = 0; i < cell_data.size(); i++) {
		//		if (cell_data[i].box[1] > box[4] && box[1] > cell_data[i].box[4] && \
		//			cell_data[i].box[2] > box[5] && box[2] > cell_data[i].box[5] && \
		//			cell_data[i].box[0] > box[3] && box[0] > cell_data[i].box[3]) {
		//			const auto &pNods = pOctreeAgent->pNod;

		//			const auto &triIdx2 = pOctreeAgent->pEle[cell_data[i].id].indexs;

		//			short share_node_num = 0, k, j;
		//			int share1, share2;

		//			for (k = 0; k < 3; k++) {
		//				for (j = 0; j < 3; j++) {
		//					if (triIdx1[k] == triIdx2[j])
		//					{
		//						share_node_num++;
		//						share1 = k;
		//						share2 = j;
		//					}
		//				}
		//			}


		//			if (share_node_num > 0) {
		//				continue;
		//			}



		//			static int f = 0;
		//			if (f == 0)
		//				printElement("firsttryoctree.vtk");
		//			printpls(pNods[triIdx1[0]].coord, pNods[triIdx1[1]].coord, pNods[triIdx1[2]].coord, \
		//				pNods[triIdx2[0]].coord, pNods[triIdx2[1]].coord, pNods[triIdx2[2]].coord, std::to_string(f++) + ".pls");


		//			if (tri_tri_inter(pNods[triIdx1[0]].coord, pNods[triIdx1[1]].coord, pNods[triIdx1[2]].coord, \
		//				pNods[triIdx2[0]].coord, pNods[triIdx2[1]].coord, pNods[triIdx2[2]].coord)) {
		//				last_intersection_ = cell_data[i].id;
		//				exit(-99);
		//			}
		//		}
		//	}
		//}
#if ONED_FILTER ==1
		for (int i = it_end -1; i > start; i--) {
#else
		for (int i = cell_data.size()-1; i >=0; i--) {
#endif

			

			//if (false&&sorted) {
			//	if (time_data <= pOctreeAgent->pEle[cell_data[i].id].time_stamp)
			//		continue;
			//}
#if 0		
			if (
				
				cell_data[i].box[second_length_] >              box[second_length_ + 3] &&\
				             box[second_length_] > cell_data[i].box[second_length_ + 3] && \
				cell_data[i].box[last_length_] >                box[last_length_ + 3] &&\
				             box[last_length_] >   cell_data[i].box[last_length_ + 3] && \
				cell_data[i].box[most_length_] >                box[most_length_ + 3]) {//这里取5个比较而不是6个
				
#else
			if (cell_data[i].box[1] > box[4] && box[1] > cell_data[i].box[4] && \
				cell_data[i].box[2] > box[5] && box[2] > cell_data[i].box[5] && \
				cell_data[i].box[0] > box[3] && box[0] > cell_data[i].box[3]) {
#endif
				const auto &pNods = pOctreeAgent->pNod;
				
				const auto &triIdx2 = pOctreeAgent->pEle[cell_data[i].id].indexs;
				short share_node_num = 0, k, j;
				int share1, share2;

                for (k = 0; k < 3; k++) {
				for (j = 0; j < 3; j++) {
					if (triIdx1[k] == triIdx2[j])
					{
						share_node_num++;
						share1 = k;
						share2 = j;
					}
				}
				}


				if (share_node_num > 1) {
					continue;
				}


				
				if (!share_node_num)
				{
					if (tri_tri_inter(pNods[triIdx1[0]].coord, pNods[triIdx1[1]].coord, pNods[triIdx1[2]].coord, 
									  pNods[triIdx2[0]].coord, pNods[triIdx2[1]].coord, pNods[triIdx2[2]].coord)) {
						last_intersection_ = cell_data[i].id;
						return true;
					}

				}

				else if (OCT::share_one_node_tri_tri_inter(pNods[triIdx1[0]].coord, pNods[triIdx1[1]].coord, pNods[triIdx1[2]].coord, \
					pNods[triIdx2[0]].coord, pNods[triIdx2[1]].coord, pNods[triIdx2[2]].coord, share1, share2)) {
						last_intersection_ = cell_data[i].id;
						return true;
				}


			}
		}
		
	}

	return false;
}



double Octree::chckIntersectWithLine(BLVector start, BLVector end, TreeNode * node, const OCCUBE & cube)
{
	TreeNode* pChild;
	double ans = 1;
	bool bNoChild = true, ret = false;
	if (node != NULL && pOctreeAgent->LineInCube(start, end, cube))
	{
		if (pOctreeAgent->LineTotalInCube(start, end, cube)) {
			node_before = node;
		}
		for (int i = 0; i < OC_SUBNODES; i++)
		{


			pChild = node->pChildNodes[i];

			bNoChild = bNoChild && (pChild == NULL);
		}

		if (!bNoChild)//如果不是子节点
		{
			OCCUBE subcube[OC_SUBNODES];

			//divide cube to eight sub-cubes
			CUBECOORD coorcenter, len;
			coorcenter = (cube.lower + cube.upper) / 2;
			len = cube.upper-coorcenter;

			//subcube[UNE]
			subcube[UNE].lower = coorcenter;
			subcube[UNE].upper = coorcenter + len;

			//subcube[UNW]
			subcube[UNW].lower = coorcenter;
			subcube[UNW].lower.x -= len.x;
			subcube[UNW].upper = coorcenter;
			subcube[UNW].upper.y += len.y;
			subcube[UNW].upper.z += len.z;

			//subcube[USW]
			subcube[USW].lower = coorcenter;
			subcube[USW].lower.x -= len.x;
			subcube[USW].lower.y -= len.y;
			subcube[USW].upper = coorcenter;
			subcube[USW].upper.z += len.z;

			//subcube[USE]
			subcube[USE].lower = coorcenter;
			subcube[USE].lower.y -= len.y;
			subcube[USE].upper = coorcenter;
			subcube[USE].upper.x += len.x;
			subcube[USE].upper.z += len.z;

			//subcube[LNE]
			subcube[LNE].lower = coorcenter;
			subcube[LNE].lower.z -= len.z;
			subcube[LNE].upper = coorcenter;
			subcube[LNE].upper.x += len.x;
			subcube[LNE].upper.y += len.y;

			//subcube[LNW]
			subcube[LNW].lower = coorcenter;
			subcube[LNW].lower.x -= len.x;
			subcube[LNW].lower.z -= len.z;
			subcube[LNW].upper = coorcenter;
			subcube[LNW].upper.y += len.y;

			//subcube[LSW]
			subcube[LSW].lower = coorcenter - len;
			subcube[LSW].upper = coorcenter;

			//subcube[LSE]
			subcube[LSE].lower = coorcenter;
			subcube[LSE].lower.y -= len.y;
			subcube[LSE].lower.z -= len.z;
			subcube[LSE].upper = coorcenter;
			subcube[LSE].upper.x += len.x;


			ans = min(ans, chckIntersectWithLine(start, end, node->pChildNodes[UNE], subcube[UNE]));
			ans = min(ans, chckIntersectWithLine(start, end, node->pChildNodes[UNW], subcube[UNW]));
			ans = min(ans, chckIntersectWithLine(start, end, node->pChildNodes[USW], subcube[USW]));
			ans = min(ans, chckIntersectWithLine(start, end, node->pChildNodes[USE], subcube[USE]));
			ans = min(ans, chckIntersectWithLine(start, end, node->pChildNodes[LNE], subcube[LNE]));
			ans = min(ans, chckIntersectWithLine(start, end, node->pChildNodes[LNW], subcube[LNW]));
			ans = min(ans, chckIntersectWithLine(start, end, node->pChildNodes[LSW], subcube[LSW]));
			ans = min(ans, chckIntersectWithLine(start, end, node->pChildNodes[LSE], subcube[LSE]));
		}
		else//叶子节点
		{
			auto &cell_data = node->order_data;
			for (int i = 0; i < cell_data.size() ; i++) {//这里不取等号
				double p=OCT::tri_inter_with_a_line(pOctreeAgent->pEle[cell_data[i].id].indexs, pOctreeAgent->getNod(),start,end);
				if (ans < 1e-5)
					ans = 1;
				ans = min(ans, p);
			}
		}
	}

	return ans;
}

double Octree::chckIntersectWithLine(BLVector start, BLVector end)
{
	TreeNode* tmp = node_before;
	OCCUBE cube1 = tmp->getCube();
	if (!tmp)
	{
		tmp = getRootNode();
		cube1 = tmp->getCube();
	}
	while (tmp != pRoot && !pOctreeAgent->LineTotalInCube(start, end, cube1)) {

		tmp = tmp->getFather();
		cube1 = tmp->getCube();
	}
	//cout << "path=" << path << "layer=" << node_before->layer << endl;

	return chckIntersectWithLine(start, end, tmp, cube1);
}
void Octree::record_stone_position()
{
	//time_stone_ = time_machine;
}
inline int Octree::get_stone() {
	return time_stone_;
}
void Octree::sort_octants() {

	for (auto& target : leaf_nodes) {
		auto& cell_data = target->order_data;
		int count = 0;
		for (auto ele : cell_data) {
			cell_data[count] = ele;
			if (ele.valid) {
				count++;
			}
		}
		cell_data.resize(count);
		target->memory_queue.clear();

		sort(cell_data.begin(), cell_data.end());
		for (int i = 1; i < cell_data.size(); i++) {
			if (cell_data[i].max_value < cell_data[i - 1].max_value)
				cell_data[i].max_value = cell_data[i - 1].max_value;
		}

	}

}
void Octree::find_intersected_triangles() {
	if (sort_telerant) {
		sort_octants();
	}
	//currentlayer_cost = 0;
	possible_items.clear();
	for(auto& target: leaf_nodes)
	{
		if (!target->changed)
			continue;
		target->changed = false;
			const auto& cell_data = target->order_data;
			//currentlayer_cost += cell_data.size()*cell_data.size();
			for (int i = 1; i < cell_data.size(); i++) {
				const auto &triIdx2 = pOctreeAgent->pEle[cell_data[i].id].indexs;
				//for (int j = 0; j < cell_data.size(); j++) {
#if ONED_FILTER ==1
				for (int j = i - 1; j&&cell_data[j].max_value > cell_data[i].min_value; j--) {

#else
				for (int j = i - 1; j>=0; j--) {
#endif

					if (cell_data[i].box[1] >= cell_data[j].box[4] && cell_data[j].box[1] >= cell_data[i].box[4] && \
						cell_data[i].box[2] >= cell_data[j].box[5] && cell_data[j].box[2] >= cell_data[i].box[5] && \
						cell_data[i].box[0] >= cell_data[j].box[3] && cell_data[j].box[0] >= cell_data[i].box[3]) {
						const auto &pNods = pOctreeAgent->pNod;

						auto &triIdx1 = pOctreeAgent->pEle[cell_data[j].id].indexs;
						
						short share_node_num = 0, k,l;
						int share1, share2;

						for (k = 0; k < 3; k++) {
							for (l = 0; l < 3; l++) {
								if (triIdx1[k] == triIdx2[l])
								{
									share_node_num++;
									share1 = k;
									share2 = l;
								}
							}
						}


						if (share_node_num > 1) {
							continue;
						}



						if (!share_node_num)
						{
							if (tri_tri_inter(pNods[triIdx1[0]].coord, pNods[triIdx1[1]].coord, pNods[triIdx1[2]].coord, \
								pNods[triIdx2[0]].coord, pNods[triIdx2[1]].coord, pNods[triIdx2[2]].coord)) {
								possible_items.insert(triIdx1[0]);
								possible_items.insert(triIdx1[1]);
								possible_items.insert(triIdx1[2]);
								possible_items.insert(triIdx2[0]);
								possible_items.insert(triIdx2[1]);
								possible_items.insert(triIdx2[2]);

								//tri_tri_inter(pNods[triIdx1[0]].coord, pNods[triIdx1[1]].coord, pNods[triIdx1[2]].coord, \
								//	pNods[triIdx2[0]].coord, pNods[triIdx2[1]].coord, pNods[triIdx2[2]].coord);


								//OCT::printpls(pNods[triIdx1[0]].coord, pNods[triIdx1[1]].coord, pNods[triIdx1[2]].coord, \
								//	pNods[triIdx2[0]].coord, pNods[triIdx2[1]].coord, pNods[triIdx2[2]].coord,"111.pls");

								//exit(1);

	/*							for (int k = 0; k < 6; k++)
									cout << cell_data[i].box[k] << " ";
								cout << endl;
								for (int k = 0; k < 6; k++)
									cout << cell_data[j].box[k] << " ";
								cout << endl;
								cout << bool(cell_data[i].box[1] > cell_data[j].box[4]) << bool(cell_data[j].box[1] > cell_data[i].box[4]) << \
									    bool(cell_data[i].box[2] > cell_data[j].box[5]) << bool(cell_data[j].box[2] > cell_data[i].box[5]) << \
										bool(cell_data[i].box[0] > cell_data[j].box[3]) << bool(cell_data[j].box[0] > cell_data[i].box[3]) << endl;**/
								//break;
							}

						}

						else
						{
							if (OCT::share_one_node_tri_tri_inter(pNods[triIdx1[0]].coord, pNods[triIdx1[1]].coord, pNods[triIdx1[2]].coord, \
								pNods[triIdx2[0]].coord, pNods[triIdx2[1]].coord, pNods[triIdx2[2]].coord, share1, share2)) {
								possible_items.insert(triIdx1[0]);
								possible_items.insert(triIdx1[1]);
								possible_items.insert(triIdx1[2]);
								possible_items.insert(triIdx2[0]);
								possible_items.insert(triIdx2[1]);
								possible_items.insert(triIdx2[2]);
								//break;
							}
						}
						


					}
				}

			}
		}
	//cout << currentlayer_cost << " " << firstlayer_cost << endl;
	//static int rebuild_threadhold = 6;
	//if (currentlayer_cost > rebuild_threadhold * firstlayer_cost) {
	//	rebuild();
	//	rebuild_threadhold +=10;
	//}

}
bool Octree::check_intersection_in_set(int id)
{
	return possible_items.find(id) != possible_items.end();
	return false;
}
}
