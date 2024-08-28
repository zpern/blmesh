/* ----------------------------------------------------------------------------
 * ----------------------------------------------------------------------------
 *
 * 三维各向同性Delaunay网格生成器 (版本号：0.3)
 * 3D Isotropic Delaunay Mesh Generation (Version 0.3)
 *
 * 陈建军 中国 浙江大学工程与科学计算研究中心
 * 版权所有	  2005年9月15日
 * Chen Jianjun  Center for Engineering & Scientific Computation,
 * Zhejiang University, P. R. China
 * Copyright reserved, 2005, 10, 26
 * 
 * 联系方式
 *   电话：+86-571-87953165
 *   传真：+86-571-87953167
 *   邮箱：zdchenjj@yahoo.com.cn
 * For further information, please conctact
 *  Tel: +86-571-87953165
 *  Fax: +86-571-87953167
 * Mail: zdchenjj@yahoo.com.cn
 *
 * ------------------------------------------------------------------------------
 * ------------------------------------------------------------------------------*/
#include <spdlog/spdlog.h> 
 #include "earremoval.h"
#include "geom_func.h"
//#include "iso3d_define.h"
#include <assert.h>
#include <memory.h>
#ifndef __APPLE__
#include <malloc.h>
#endif
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

BOOL EarRemover::Left(double x1, double y1,
		  double x2, double y2,
		  double x3, double y3)
{
	double p1[2] = {x1, y1};
	double p2[2] = {x2, y2};
	double p3[2] = {x3, y3};
	return GEOM_FUNC::orient2d(p1, p2, p3) > 0.0;
#if 0
	double squad = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
//	assert(squad >= exp_abs);
	return ((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)) / squad >= exp_err;
#endif
}

BOOL EarRemover::LeftOn(double x1, double y1,
		    double x2, double y2,
		    double x3, double y3)
{
#if 0
	double squad = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
//	assert(squad >= exp_abs);
	return ((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)) / squad > -exp_err;
#endif
	double p1[2] = {x1, y1};
	double p2[2] = {x2, y2};
	double p3[2] = {x3, y3};
	return GEOM_FUNC::orient2d(p1, p2, p3) >= 0.0;
}

BOOL EarRemover::Collinear(double x1, double y1,
		       double x2, double y2,
		       double x3, double y3)
{
#if 0
	double squad = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
//	assert(squad >= exp_abs);
	return fabs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)) / squad < exp_err;
#endif
double p1[2] = {x1, y1};
	double p2[2] = {x2, y2};
	double p3[2] = {x3, y3};
	return GEOM_FUNC::orient2d(p1, p2, p3) == 0.0;
}


BOOL EarRemover::InCone(tVertex a, tVertex b)
{
	tVertex a0, a1; /* a0, a, a1 are consecutive vertices */

	a1 = a->next;
	a0 = a->prev;
	
	/* if a is a convex vertex ... */
	if (LeftOn(a->x, a->y, a1->x, a1->y, a0->x, a0->y))
		return Left(a->x, a->y, b->x, b->y, a0->x, a0->y) &&
		       Left(b->x, b->y, a->x, a->y, a1->x, a1->y);
	else
		return !(LeftOn(a->x, a->y, b->x, b->y, a1->x, a1->y) &&
		         LeftOn(b->x, b->y, a->x, a->y, a0->x, a0->y));

	return 0;
}

BOOL EarRemover::Xor(BOOL x, BOOL y)
{
	return !x ^ !y;
}

BOOL EarRemover::IntersectProp(double x1, double y1, double x2, double y2,
				   double x3, double y3, double x4, double y4)
{
	if (Collinear(x1, y1, x2, y2, x3, y3) ||
		Collinear(x1, y1, x2, y2, x4, y4) ||
		Collinear(x3, y3, x4, y4, x1, y1) ||
		Collinear(x3, y3, x4, y4, x2, y2))
		return FALSE;

	return Xor(Left(x1, y1, x2, y2, x3, y3), Left(x1, y1, x2, y2, x4, y4)) &&
		   Xor(Left(x3, y3, x4, y4, x1, y1), Left(x3, y3, x4, y4, x2, y2));
}

BOOL EarRemover::Between(double x1, double y1, double x2, double y2, double x3, double y3)
{
	if (!Collinear(x1, y1, x2, y2, x3, y3))
		return FALSE;

	/* if ab not vertical, check betweeness on x, else on y */
	if (fabs(x1 - x2) >= exp_abs)
		return ((x1 <= x3 && x3 <= x2) ||
		        (x1 >= x3 && x3 >= x2));
	else
		return ((y1 <= y3 && y3 <= y2) ||
		        (y1 >= y3 && y3 >= y2));
}

BOOL EarRemover::Intersect(double x1, double y1, double x2, double y2,
			   double x3, double y3, double x4, double y4)
{
	if (IntersectProp(x1, y1, x2, y2, x3, y3, x4, y4))
		return TRUE;
	else if (Between(x1, y1, x2, y2, x3, y3) ||
		     Between(x1, y1, x2, y2, x4, y4) ||
			 Between(x3, y3, x4, y4, x1, y1) ||
			 Between(x3, y3, x4, y4, x2, y2))
		return TRUE;
	else
		return FALSE;
}

BOOL EarRemover::Diagonalie(tVertex a, tVertex b)
{
	tVertex c, c1;

	/* For each edges (c, c1) of P */
	c = vertices;
	do
	{
		c1 = c->next;
		/* skip edges incident to a or b */
		if ((c != a) && (c1 != a) &&
			(c != b) && (c1 != b) &&
			Intersect(a->x, a->y, b->x, b->y, c->x, c->y, c1->x, c1->y))
			return FALSE;

		c = c->next;
	}
	while (c != vertices);
	
	return TRUE;
}

BOOL EarRemover::Diagonal(tVertex a, tVertex b)
{
	return InCone(a, b) && InCone(b, a) && Diagonalie(a, b);
}

void EarRemover::EarInit(void)
{
	tVertex v0, v1, v2; /* three consecutive vertices */
	
	/* Initialize v1->ear for all vertices */
	nears = 0;
	v1 = vertices;
	do
	{
		v2 = v1->next;
		v0 = v1->prev;
		v1->ear = nvertices > 3 ? Diagonal(v0, v2) : 0;
		if (v1->ear)
			nears++;
		v1->out_tri = -1;
		v1 = v1->next;

	}
	while (v1 != vertices);
}

bool EarRemover::TopuAllowed(int id1, int id2)
{
	int i, idx1, idx2;
	for (i = 0; i < g_colNDSize; i++)
	{
		idx1 = g_shellNDs_ColinIdx[i][id1];
		idx2 = g_shellNDs_ColinIdx[i][id2];
		if (idx1 >= 0 && idx2 >= 0 && abs(idx1 - idx2) > 1)
			return false;
	}

	return true;
}

int EarRemover::EarRemoval(void)
{
	tVertex v0, v1, v2, v3, v4; /* five consecutive vertices */
	int n = nvertices;

	EarInit();

	if (n > 3 && nears <= 0)
		return 0; /* cannot be triangulated */

	/* Each step of outer loop removes one ear */
	while ((n > 3 && nears > 0) || n == 3)
	{
		v2 = vertices;
		do
		{
			if (v2->ear || n == 3) 
			{/* Ear found. Fill variables. */
				v3 = v2->next;
				v4 = v3->next;

				v1 = v2->prev;
				v0 = v1->prev;

				/* 确定v3和v1不是拓扑不相容的点（共线且索引值相差超过1） */
				if (TopuAllowed(v1->vnum, v3->vnum))
				{
					/* (v1, v3) is a diagonal */
					(triangles[ntris].forms)[0] = v1->vnum;
					(triangles[ntris].forms)[1] = v2->vnum;
					(triangles[ntris].forms)[2] = v3->vnum;
					(triangles[ntris].neigs)[0] = 
						(triangles[ntris].neigs)[1] =
						(triangles[ntris].neigs)[2] = -1;
					if (v1->out_tri >= 0)
					{
						triangles[v1->out_tri].neigs[1] = ntris;
						triangles[ntris].neigs[2] = v1->out_tri;
					}
					if (v2->out_tri >= 0)
					{
						triangles[v2->out_tri].neigs[1] = ntris;
						triangles[ntris].neigs[0] = v2->out_tri;
					}
					if (v3->out_tri >= 0)
					{
						triangles[v3->out_tri].neigs[1] = ntris;
						triangles[ntris].neigs[1] = v3->out_tri;
					}
					v1->out_tri = ntris;
					if ((v2->vnum - v1->vnum + nvertices) % nvertices == 1)
						v1->prt = ntris;
					if ((v3->vnum - v2->vnum + nvertices) % nvertices == 1)
						v2->prt = ntris;
					if ((v1->vnum - v3->vnum + nvertices) % nvertices == 1)
						v3->prt = ntris;
					++ntris;


					/* update earity of diagonal end points */
					nears--; /* v2 does not exists */
					if (v1->ear)
						nears--;
					if (v3->ear)
						nears--;

					v1->ear = Diagonal(v0, v3);
					v3->ear = Diagonal(v1, v4);

					if (v1->ear)
						nears++;
					if (v3->ear)
						nears++;
				
					/* cut off the ear v2 */
					v1->next = v3;
					v3->prev = v1;
					vertices = v3; /* In case the head was v2 */
					n--;
					break;
				}
			} /* end if ear found */

			v2 = v2->next;
		}
		while (v2 != vertices);
	}

	return n < 3; /* triangulated or not */
}

void EarRemover::Init(double x[], double y[], int nn, int nf)
{
	int i;
	if (nn >= 3 && nf > 0)
	{
		vertices = (tVertex) malloc(sizeof(tsVertex) * nn);
		triangles = (Tri2D *) malloc(sizeof(Tri2D) * nf);
		ver_head = vertices;
		nvertices = nn;

		for (i = 0; i < nn; i++)
		{
			vertices[i].x = x[i];
			vertices[i].y = y[i];
			vertices[i].out_tri = -1;
			vertices[i].vnum = i;
			vertices[i].next = vertices[i].prev = nullptr;
			vertices[i].ear = FALSE;
			vertices[i].prt = -1;
		}
		vertices[0].prev = &vertices[nn-1];
		vertices[nn-1].next = &vertices[0];
		for (i = 1; i < nn; i++)
		{
			vertices[i].prev = &vertices[i - 1];
			vertices[i - 1].next = &vertices[i];
		}

		memset(triangles, -1, sizeof(Tri2D) * nf);
	}
}

void EarRemover::Unitize()
{
	double xmin, xmax, ymin, ymax, xcen, ycen;
	double scalex, scaley, scale;
	tVertex v;
	if (nvertices > 0)
	{
		v = vertices;
		assert(v);
		xmin = xmax = v->x;
		ymin = ymax = v->y;
		v = v->next;
		while (v != vertices)
		{
			if (v->x < xmin)
				xmin = v->x;
			if (v->x > xmax)
				xmax = v->x;
			if (v->y < ymin)
				ymin = v->y;
			if (v->y > ymax)
				ymax = v->y;

			v = v->next;
		}
		xcen = 0.5 * (xmin + xmax);
		ycen = 0.5 * (ymin + ymax);

		if (xmax - xmin > 0.0 && ymax - ymin > 0.0)
		{
			scalex = 1.0 / (xmax - xmin);
			scaley = 1.0 / (ymax - ymin);
			scale = scalex > scaley ? scaley : scalex;

			v = vertices;
			do
			{
				v->x = (v->x - xcen) * scale + 0.5;
				v->y = (v->y - ycen) * scale + 0.5;

				v = v->next;
			}
			while (v != vertices);
		}
	}
}

void EarRemover::Swap()
{
	int i, ii;
	int it1, it2;
	int m;
	int nSwap;

	/* void */
	for (i = 0; i < 100; i++)
	{
		nSwap = 0;
		for (ii = 0; ii < ntris; ii++)
		{
			it1 = ii;
			for (m = 0; m < 3; m++)
			{
				it2 = (triangles[it1].neigs)[m];
				if (it2 > it1)
				{
					if (SwapBase(it1, it2) != 0)
					{
						nSwap++;
						break;
					}
				}
			}
		}

		if (nSwap == 0)
			break;
	}
}

int EarRemover::SwapBase(int it1, int it2)
{
	int ia, ib1, ic1, ib2, ic2, id;
	int a, b, c, d;
	int nb1, nc1, nb2, nc2;
	int m;
	double doubArea, cenx, ceny, sqRad, cen_to_d;

	ia = -1;
	for (m = 0; m < 3; m++)
		if ((triangles[it1].neigs)[m] == it2)
		{
			ia = m;
			break;
		}
	if (ia < 0)
		return 0;
	ib1 = (ia + 1) % 3;
	ic1 = (ib1 + 1) % 3;

	id = -1;
	for (m = 0; m < 3; m++)
		if ((triangles[it2].neigs)[m] == it1)
		{
			id = m;
			break;
		}
	if (id < 0)
		return 0;
	ic2 = (id + 1) % 3;
	ib2 = (ic2 + 1) % 3;

	if ((triangles[it1].forms)[ib1] != (triangles[it2].forms)[ib2] || 
		(triangles[it1].forms)[ic1] != (triangles[it2].forms)[ic2])
		return 0;

	a = (triangles[it1].forms)[ia];
	b = (triangles[it1].forms)[ib1];
	c = (triangles[it1].forms)[ic1];
	d = (triangles[it2].forms)[id];

	/* ad是交换后形成的边，确定其是否拓扑允许 */
	if (!TopuAllowed(a, d))
		return 0;

	if (Intersect(ver_head[a].x, ver_head[a].y, ver_head[d].x, ver_head[d].y,
				  ver_head[b].x, ver_head[b].y, ver_head[c].x, ver_head[c].y) &&
		!Collinear(ver_head[a].x, ver_head[a].y, 
		           ver_head[d].x, ver_head[d].y, 
				   ver_head[b].x, ver_head[b].y) &&
		!Collinear(ver_head[a].x, ver_head[a].y, 
		           ver_head[d].x, ver_head[d].y, 
				   ver_head[c].x, ver_head[c].y)
		)
	{

		calcTriParams(ver_head[a].x, ver_head[a].y, 
			          ver_head[b].x, ver_head[b].y,
					  ver_head[c].x, ver_head[c].y,
					  &doubArea, &cenx, &ceny, &sqRad);

		cen_to_d = (cenx - ver_head[d].x) * (cenx - ver_head[d].x) + 
			       (ceny - ver_head[d].y) * (ceny - ver_head[d].y);
		if (cen_to_d < sqRad)
		{/* Swap it */
			nb1 = (triangles[it1].neigs)[ib1];
			nc1 = (triangles[it1].neigs)[ic1];
			nb2 = (triangles[it2].neigs)[ib2];
			nc2 = (triangles[it2].neigs)[ic2];
			
			(triangles[it1].forms)[0] = a;
			(triangles[it1].forms)[1] = d;
			(triangles[it1].forms)[2] = c;
			(triangles[it1].neigs)[0] = nb2;
			(triangles[it1].neigs)[1] = nb1;
			(triangles[it1].neigs)[2] = it2;
			if (nb2 >= 0)
			{
				for (m = 0; m < 3; m++)
				{
					if ((triangles[nb2].forms)[m] != c &&
						(triangles[nb2].forms)[m] != d)
					{
						(triangles[nb2].neigs)[m] = it1;
						break;
					}
				}
			}
			if (nb1 >= 0)
			{
				for (m = 0; m < 3; m++)
				{
					if ((triangles[nb1].forms)[m] != a &&
						(triangles[nb1].forms)[m] != c)
					{
						(triangles[nb1].neigs)[m] = it1;
						break;
					}
				}
			}
			
			(triangles[it2].forms)[0] = d;
			(triangles[it2].forms)[1] = a;
			(triangles[it2].forms)[2] = b;
			(triangles[it2].neigs)[0] = nc1;
			(triangles[it2].neigs)[1] = nc2;
			(triangles[it2].neigs)[2] = it1;
			if (nc1 >= 0)
			{
				for (m = 0; m < 3; m++)
				{
					if ((triangles[nc1].forms)[m] != a &&
						(triangles[nc1].forms)[m] != b)
					{
						(triangles[nc1].neigs)[m] = it2;
						break;
					}
				}
			}
			if (nc2 >= 0)
			{
				for (m = 0; m < 3; m++)
				{
					if ((triangles[nc2].forms)[m] != b &&
						(triangles[nc2].forms)[m] != d)
					{
						(triangles[nc2].neigs)[m] = it2;
						break;
					}
				}
			}

			return 1;

		}
	}

	return 0;
}

void EarRemover::calcTriParams(
				double x1, double y1, 
				double x2, double y2, 
				double x3, double y3,
				double *doubArea, 
				double *cenx, double *ceny,
				double *sqRad
				)
{
	double a, b, d, e, d1, d2, d3, c, f, te, t1, t2;
	a = 2.0 * (x3 - x2);
	b = 2.0 * (y3 - y2);
	d = 2.0 * (x1 - x2);
	e = 2.0 * (y1 - y2);
	d1 = x1 * x1 + y1 * y1;
	d2 = x2 * x2 + y2 * y2;
	d3 = x3 * x3 + y3 * y3;
	c = d3 - d2;
	f = d1 - d2;
	te = a * e - b * d;
	t1 = c * e - f * b;
	t2 = a * f - c * d;
	*doubArea = te;
	*cenx = t1 / te;
	*ceny = t2 / te;
	*sqRad = ((*cenx - x1) * (*cenx - x1) + (*ceny - y1) * (*ceny - y1));
}

void EarRemover::Output(int prt[], int forms[], int neigs[], int nf)
{
	int i, j;
	int m;

/*	for (i = 0; i < nvertices; i++)
	{
		prt[i] = ver_head[i].prt;
	}*/

	j = 0;
	for (i = 0; i < nf; i++, j++)
	{
		forms[3 * j] = (triangles[i].forms)[0];
		forms[3 * j + 1] = (triangles[i].forms)[1];
		forms[3 * j + 2] = (triangles[i].forms)[2];
		neigs[3 * j] = (triangles[i].neigs)[0];
		neigs[3 * j + 1] = (triangles[i].neigs)[1];
		neigs[3 * j + 2] = (triangles[i].neigs)[2];

		for (m = 0; m < 3; m++)
		{
			if ((forms[3*j + (m+1)%3] - forms[3*j + m] + nvertices) % nvertices == 1)
			{
				prt[forms[3*j + m]] = i;
			}
			else if ((forms[3*j + m] - forms[3*j + (m+1)%3] + nvertices) % nvertices == 1)
			{
				prt[forms[3*j + (m+1)%3]] = i;
			}
		}
	}
}

void EarRemover::FreeMemory()
{
	if (ver_head)
	{
		free(ver_head);
		ver_head = vertices = nullptr;
		nvertices = 0;
		nears = 0;
	}

	if (triangles)
	{
		free(triangles);
		triangles = nullptr;
		ntris = 0;
	}
}

//#ifdef __PolyTri_Test__
int EarRemover::readFr2(const char *fname, double x[], double y[], int *nn)
{
	int i, np, ne;
	int iTok, is, ie, ic, il;
	char chLine[512];
	FILE *fp = fopen(fname, "r");
	if (!fp)
	{
		spdlog::info("cannot open file %s\n", fname);
		return 0;
	}

	fgets(chLine, 512, fp);
	sscanf(chLine, "%d %d", &np, &ne);
	if (np != ne || np >= *nn)
	{
		fclose(fp);
		spdlog::info("error paramerters or file data: %s\n", fname);

		return 0;
	}
	
	*nn = np;
	for (i = 0; i < *nn; i++)
		fscanf(fp, "%d %lf %lf\n", &iTok, &x[i], &y[i]);

	for (i = 0; i < *nn; i++)
	{
		fscanf(fp, "%d %d %d %d %d\n", &iTok, &is, &ie, &ic, &il);
	}

	fclose(fp);

	return 1;
}

int EarRemover::writeFr2(const char *fname, double x[], double y[], int nn)
{
	int i;
	FILE *fp = fopen(fname, "w");
	if (!fp)
	{
		spdlog::info("cannot open file %s\n", fname);
		return 0;
	}

	fprintf(fp, "%d %d 0 0 0 0 0\n", nn, nn);
	for (i = 0; i < nn; i++)
		fprintf(fp, "%d %24.20g %24.20g\n", i+1, x[i], y[i]);
	for (i = 0; i < nn - 1; i++)
		fprintf(fp, "%d %d %d 1 1\n", i+1, i+1, i+2);
	fprintf(fp, "%d %d %d 1 1\n", nn, nn, 1);

	fclose(fp);

	return 1;
}

int EarRemover::writePL2(const char *fname, double x[], double y[], int prt[], int nn,
			 int forms[], int nf)
{
	int i;
	FILE *fp = fopen(fname, "w");
	if (!fp)
	{
		spdlog::info("cannot open file %s\n", fname);
		return 0;
	}

	fprintf(fp, "%d %d %d\n", nf, nn, nn);
	for (i = 0; i < nn; i++)
		fprintf(fp, "%d %f %f\n", i+1, x[i], y[i]);
	for (i = 0; i < nf; i++)
		fprintf(fp, "%d %d %d %d 0 0\n", i + 1, forms[3*i] + 1, forms[3*i+1] + 1, forms[3*i+2] + 1);
	for (i = 0; i < nn - 1; i++)
		fprintf(fp, "%d %d %d %d 0 0\n", i+1, i+1, i+2, prt[i] + 1);
	fprintf(fp, "%d %d %d %d 0 0\n", nn, nn, 1, prt[nn - 1] + 1);

	fclose(fp);

	return 1;
}

//#endif /* #ifdef __PolyTri_Test__ */

int EarRemover::triangulate(double x[], double y[], int nn,
	int prt[], int forms[], int neigs[], int nf, int shellNDs_ColinIdx[][MAX_SHELL_ED], int colNDSize)
{
	int nRet = 0;
	int i;

	g_colNDSize = colNDSize;
	for (i  = 0; i < colNDSize; i++)
		memcpy(g_shellNDs_ColinIdx[i], shellNDs_ColinIdx[i], sizeof(int)*MAX_SHELL_ED);

	Init(x, y, nn, nf);
	Unitize();
#ifdef _DEBUG
	writeFr2("shelltest.fr2", x, y, nn);
#endif
	if (EarRemoval())
	{
		Swap();
		Output(prt, forms, neigs, nf);
#ifdef _DEBUG
		writePL2("shelltest.pl2", x, y, prt, nn, forms, nf);
#endif
		nRet = 1;
	}
	FreeMemory();
	return nRet;
}

EarRemover::~EarRemover()
{
	FreeMemory();
}
