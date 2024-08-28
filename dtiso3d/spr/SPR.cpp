//#define _SPR_USING_STD_VECTOR

#ifndef WIN32
#include <sys/time.h>
#else
#include <windows.h>
#include <winbase.h>
#include <sys/timeb.h>
#endif

#include <spdlog/spdlog.h> 
 #include "spr.h"
#include "IntIntMap.h"
#include <vector>
#include <string.h>
#include <string>
#include <exception>

bool operator == (const divide_record & t1, const divide_record& t2)
{
	return t1.num == t2.num && memcmp(t1.div_rd, t2.div_rd, sizeof(int)*t1.num) == 0;
}

bool operator < (const divide_record & t1, const divide_record& t2)
{
	return (t1.num < t2.num) || ((t1.num == t2.num) && memcmp(t1.div_rd, t2.div_rd, sizeof(int)*t1.num) < 0);
}

bool operator > (const divide_record & t1, const divide_record& t2)
{
	return (t1.num > t2.num) || ((t1.num == t2.num) && memcmp(t1.div_rd, t2.div_rd, sizeof(int)*t1.num) > 0);
}

/* local functions */
inline void vec_2p(double *a, double *b, double *ab)
{
	int i;
	for (i = 0; i < 3; i++)
		ab[i] = b[i] - a[i];
}

inline void vec_crop(double *vector1, double *vector2, double *vector3)
{
	vector3[0] = vector1[1] * vector2[2] - vector2[1] * vector1[2];
	vector3[1] = vector2[0] * vector1[2] - vector1[0] * vector2[2];
	vector3[2] = vector1[0] * vector2[1] - vector2[0] * vector1[1];
}

inline double vec_val(double *vec)
{
	return(sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]));
}

inline void norm_3p(double *p1, double *p2, double *p3, double *normal)
{
	double v12[3], v13[3];
	vec_2p(p1, p2, v12);
	vec_2p(p1, p3, v13); // the direction of crop has changed ?!
	vec_crop(v12, v13, normal);
}

inline int vec_uni(double *vec)
{
	double len;
	int i;

	len = (double)sqrt(vec[0] * vec[0] + vec[1] * vec[1]
		+ vec[2] * vec[2]);
	if (len <= 0.0)
		return(0);        /* eps_2 zero vector bound :min_siz/2. */
	for (i = 0; i < 3; i++)
		vec[i] /= len;
	return(1);
}

inline double vec_dotp(double *vector1, double *vector2)
{
	return(vector1[0] * vector2[0] + vector1[1] * vector2[1] + vector1[2] * vector2[2]);
}

inline double comp2fAngle(double pa[3], double pb[3], double pc[3], double pd[3])
{
	double norm1[3], norm2[3];

	norm_3p(pa, pb, pc, norm1);
	norm_3p(pa, pb, pd, norm2);

#if 0
	if (vec_uni(norm1) == 0)
		fprintf(stderr, "err comp2fAngle, change mesh size and try again.\n ");
	if (vec_uni(norm2) == 0)
		fprintf(stderr, "err comp2fAngle, change mesh size and try again.\n ");
#endif

	return (vec_dotp(norm1, norm2));
}

inline double heightOfTriangle(double pa[3], double pb[3], double pc[3], double normal[3])
{
	double vab[3], vac[3], vbc[3], *vabc = normal;//vabc[3];
	double l1, l2, l3, lmax;
	vec_2p(pa, pb, vab);
	vec_2p(pa, pc, vac);
	vec_2p(pb, pc, vbc);

	l1 = vab[0] * vab[0] + vab[1] * vab[1] + vab[2] * vab[2];
	l2 = vac[0] * vac[0] + vac[1] * vac[1] + vac[2] * vac[2];
	l3 = vbc[0] * vbc[0] + vbc[1] * vbc[1] + vbc[2] * vbc[2];
	lmax = MAX_VALUE(MAX_VALUE(l1, l2), l3);

	vec_crop(vab, vac, vabc);

	return sqrt((vabc[0] * vabc[0] + vabc[1] * vabc[1] + vabc[2] * vabc[2]) / lmax);
}

inline double distance3D(double *v1, double *v2)
{
	return sqrt(
		(v1[0] - v2[0]) * (v1[0] - v2[0]) +
		(v1[1] - v2[1]) * (v1[1] - v2[1]) +
		(v1[2] - v2[2]) * (v1[2] - v2[2]));
}

inline void pointProjectToPlane(APOINT pa, double d, APOINT norm, APOINT pj)
{
	int i;
	double dp;

	dp = norm[0] * pa[0] + norm[1] * pa[1] + norm[2] * pa[2] + d; /* 点法式 */

	for (i = 0; i < 3; i++)
		pj[i] = pa[i] - dp * norm[i];
}

SPRImpl::SPRImpl()
	:dp_stkTreeSear(MAX_SPR_POLY_SIZE)
{
}

/* end local functions */

int SPRImpl::init_spr()
{
#ifdef _SPR_MAP_TETRA_VOLUME
	mapTetraVolume.clear();
#endif /* _SPR_MAP_TETRA_VOLUME */
	real = 0;
	setDivideRecord.clear();

	//spr_int_node_num = 0; /* 临时变量，记录新面片和约束相交的次数 */
	spr_no_more_recur = false; /* 临时变量，用来记录是否需要提前终止递归：
										* 情形1：超时 */

	spr_verbose_output = false;

	spr_optimal_type = MAX_MIN_CELL_QUALITY;
	spr_max_int_node = 0;

	return 1;
}

int SPRImpl::small_poly_reconn(polyhedron *poly, divide *divi, double iniQ,
	Constraints *csts, int *smallPolySize, int maxIntNum)
{
	int res = 0;
	divide_record divRd;
	int i;
	GridEdgeHash edgehash;
	int rtn;
#ifdef _SPR_USING_STD_VECTOR
	std::vector<struct polyhedron*> vec_simp_poly;
#else
	polyhedron *simpPolyArray[MAX_SIMP_POLY_NUM];
	int simpPolyNum = 0;
#endif
	bool freeSimpPoly = true;
	std::vector<int> vecETopus;
	std::vector<int> vecENodes;
	polyhedron *simpPoly = NULL;

	remove_overlap_faces(poly);

	//	spdlog::info("Begin to call small_poly_reconn(...).\n");
#if 0
//#ifdef _DEBUG
	savepls("sprtest.pls", poly, true);
	saveConstraints("sprtest.cst", csts, iniQ);
	if (csts)
		memcpy(&spr_csts, csts, sizeof(Constraints));
	//#endif
#endif

	if (poly->num_t < 4)
		return 0;

	/* 确保每个单元的第一个节点编号总是最小的 */
	rearrage(poly->t, poly->num_t);

	/* 大部分情况下，网格是满足二边流形准则的，所以可以更快确定相邻关系 */
	rtn = 0;
	if (!build_manifold_poly_neig(poly))
	{
		/* 这段代码挺耗时，以后需要改 */
		vecETopus.resize(poly->num_t);
		vecENodes.resize(3 * poly->num_t);

		count = 0;
		real = 0;

		for (i = 0; i < poly->num_t; i++)
		{
			vecETopus[i] = EEMAS::TRIANGLE;

			vecENodes[3 * i + 0] = poly->t[i][0];
			vecENodes[3 * i + 1] = poly->t[i][1];
			vecENodes[3 * i + 2] = poly->t[i][2];
		}

		edgehash.set_hash_type(GridEdgeHash::MIN_KEY);
		edgehash.enable_elem_record(true);
		edgehash.initialize(vecETopus, vecENodes);

		rtn = sort_comp_edges(poly, &edgehash);

		if (rtn >= 0)
		{
			build_poly_neig(poly, &edgehash);
			assert(assert_poly_neig(poly));

#ifdef _SPR_USING_STD_VECTOR
			if (rtn > 0)
			{
				divide_poly(poly, vec_simp_poly);
			}
			else
			{
				vec_simp_poly.push_back(poly);
				freeSimpPoly = false;
			}
#else
		}
	}
	if (rtn >= 0)
	{
		if (rtn > 0)
			divide_poly(poly, simpPolyArray, &simpPolyNum);
		else
		{
			simpPolyArray[0] = poly;
			simpPolyNum = 1;
			freeSimpPoly = false;
		}
#endif /* _SPR_USING_STD_VECTOR */

		if (smallPolySize)
			*smallPolySize = 0;

#ifdef _SPR_USING_STD_VECTOR
		for (i = 0; i < vec_simp_poly.size(); i++)
		{
			simpPoly = vec_simp_poly[i];
#else
		for (i = 0; i < simpPolyNum; i++)
		{
			simpPoly = simpPolyArray[i];
#endif
			if (smallPolySize && simpPoly->num_t > *smallPolySize)
				*smallPolySize = simpPoly->num_t;

			if (simpPoly->num_t >= MAX_ALLOWED_POLY_SIZE)
			{
				res = 0;
				goto END;
			}
		}

		spr_allowed_face_num = spr_shell_size - 1; /* 比shell size 要小才称为有效解 */

#ifdef _SPR_USING_STD_VECTOR
		res = optimal_twoStep_firstCall(iniQ, vec_simp_poly, divi, &divRd, csts, maxIntNum);
#else
		res = optimal_twoStep_firstCall(iniQ, simpPolyArray, simpPolyNum, divi, &divRd, csts, maxIntNum);
#endif /* _SPR_USING_STD_VECTOR */

		//spdlog::info("End to call small_poly_reconn(...).rtn = {}\n", rtn);
		if (res == 1)
		{
#if 0
			savepl3("sprtest.pl3", poly, divi);
#endif
			//			spdlog::info("An optimal tetrahedralization is found. minimal quality is {}.\n", divi->q);
		}
		else
		{
			//			spdlog::info("Fail to find an optimal tetrahedralization.\n");
		}
		}
	else
	{
#ifdef _VERBOSE
		spdlog::info("Unsupported complex polyhedron.\n");
#endif
		goto END;
	}

END:
#ifdef _SPR_USING_STD_VECTOR
	if (freeSimpPoly)
	{
		for (i = 0; i < vec_simp_poly.size(); i++)
			free(vec_simp_poly[i]);
	}
	vec_simp_poly.clear();
#else
	if (freeSimpPoly)
	{
		for (i = 0; i < simpPolyNum; i++)
			free(simpPolyArray[i]);
	}
	simpPolyNum = 0;
#endif /* _SPR_USING_STD_VECTOR */
	return res;
	}

bool SPRImpl::is_valid_poly(struct polyhedron *poly)
{
	int i, j, m, k;
	int facei1[3], facei2[3];
	double facep1[3][3], facep2[3][3];
	bool isInt = false;

	for (i = 0; i < poly->num_t; i++)
	{
		for (m = 0; m < 3; m++)
		{
			facei1[m] = poly->t[i][m];
			for (k = 0; k < 3; k++)
				facep1[m][k] = vertices[facei1[m]][k];
		}
		for (j = i + 1; j < poly->num_t; j++)
		{
			for (m = 0; m < 3; m++)
			{
				facei2[m] = poly->t[j][m];
				for (k = 0; k < 3; k++)
					facep2[m][k] = vertices[facei2[m]][k];
			}
			isInt = GEOM_FUNC::tri_tri_intersect3d(facei1, facei2, facep1, facep2);

			if (isInt)
				return false;
		}
	}

	return true;
}

int SPRImpl::remove_overlap_faces(struct polyhedron *poly)
{
	int i, j, m, idx;
	int facei1[3], facei2[3];
	int flag[MAX_SPR_POLY_SIZE];
	bool bChanged = false;

	memset(flag, 0, sizeof(int)*MAX_SPR_POLY_SIZE);
	for (i = 0; i < poly->num_t; i++)
	{
		for (m = 0; m < 3; m++)
		{
			facei1[m] = poly->t[i][m];
		}
		for (j = i + 1; j < poly->num_t; j++)
		{
			for (m = 0; m < 3; m++)
			{
				facei2[m] = poly->t[j][m];
				if (facei2[m] != facei1[0] &&
					facei2[m] != facei1[1] &&
					facei2[m] != facei1[2])
					break;
			}

			if (m >= 3)
			{
				flag[i] = 1;
				flag[j] = 1;
				poly->used_v[facei1[0]] -= 2;
				poly->used_v[facei1[1]] -= 2;
				poly->used_v[facei1[2]] -= 2;
				bChanged = true;
			}
		}
	}

	if (!bChanged)
		return 0;

	for (i = 0, j = 0; i < poly->num_t; i++)
	{
		if (flag[i] == 0 && i > j)
		{
			poly->t[j][0] = poly->t[i][0];
			poly->t[j][1] = poly->t[i][1];
			poly->t[j][2] = poly->t[i][2];
			poly->surface[j] = poly->surface[i];
			j++;
		}
	}
	poly->num_t = j;

	for (i = 0; i < poly->num_t; i++)
	{
		for (j = 0; j < 3; j++)
		{
			idx = poly->t[i][j];
			poly->vprt[idx] = (i << 2) | j;
		}
	}
	poly->num_vable = 0;
	for (i = 0; i < num_vertices; i++)
	{
		if (poly->used_v[i] != 0)
			poly->num_vable++;
	}

	return poly->num_t;
}

#if 0
void main(int argc, char **argv)
{
	double init_q;
	int res;
	struct polyhedron *poly;
	struct divide *div;
	divide_record divRd;
	int i;

	//     // 获取每秒多少CPU Performance Tick
	//     QueryPerformanceFrequency( &lv );
	//     // 转换为每个Tick多少秒
	//     secondsPerTick = 1.0 / lv.QuadPart;
	//     
	poly = (struct polyhedron*) malloc(sizeof(struct polyhedron));
	div = (struct divide*) malloc(sizeof(struct divide));

	count = 0;
	real = 0;

	init(&init_q, poly, div);
	memset(&divRd, 0, sizeof(divide_record));

#ifdef _OCTAHEDRON_EXAMPLE
	octahedron(P);
#endif

#ifdef _SPHERE_EXAMPLE
	openpls("sphere.pls", poly);
#endif 

#ifdef _FILE_EXAMPLE
	openpls("sphere.pls", poly);// true);
//	openpls("remesh_bnd.pls", poly);
#endif

	GridEdgeHash edgehash;
	int rtn;
	std::vector<struct polyhedron*> vec_simp_poly;
	bool freeSimpPoly = true;

	std::vector<int> vecETopus(poly->num_t);
	std::vector<int> vecENodes(3 * poly->num_t);

	for (i = 0; i < poly->num_t; i++)
	{
		vecETopus[i] = EEMAS::TRIANGLE;

		vecENodes[3 * i + 0] = poly->t[i][0];
		vecENodes[3 * i + 1] = poly->t[i][1];
		vecENodes[3 * i + 2] = poly->t[i][2];
	}

	edgehash.set_hash_type(GridEdgeHash::MIN_KEY);
	edgehash.enable_elem_record(true);
	edgehash.initialize(vecETopus, vecENodes);

	rtn = sort_comp_edges(poly, &edgehash);

	if (rtn >= 0)
	{
		build_poly_neig(poly, &edgehash);

		assert(assert_poly_neig(poly));

#ifdef _SORT_FACET
		divide_poly(poly, vec_simp_poly);
#else
		if (rtn > 0)
		{
			divide_poly(poly, vec_simp_poly);
		}
		else
		{
			vec_simp_poly.push_back(poly);
			freeSimpPoly = false;
		}
#endif

		start = SPRlogTime();
		//	res = optimal(init_q, vec_simp_poly, div, &divRd);
		res = optimal_twoStep(init_q, vec_simp_poly, div, &divRd);
		stop = SPRlogTime();

		duration = stop - start;
		if (res == 1)
		{
			//			spdlog::info("An optimal tetrahedralization is found. minimal quality is {}.\n", div->q);
			savepl3("B.pl3", poly, div);
		}
		else
		{
			//			spdlog::info("Fail to find an optimal tetrahedralization.\n");
		}
	}
	else
	{
		spdlog::info("Unspported complex polyhedron.\n");
		goto END;
	}

	spdlog::info("Duration: {}\n", duration);
	spdlog::info("\tInters. Calc.\t\t:%-20.6f(%-20.6f)\n", dura[3], dura[8]);
	spdlog::info("\tQuality Calc.\t\t:%-20.6f\n", dura[4]);
	spdlog::info("\tVisibi. Calc.\t\t:%-20.6f\n", dura[5]);
	spdlog::info("\tLineCrossFace\t\t:%-20.6f\n", dura[6]);
	spdlog::info("\tOrient3d\t\t:%-20.6f\n", dura[7]);
	spdlog::info("Orient3d_n\t:%ld\n", orient3d_n);
	spdlog::info("Call_V\t:%ld\n", Call_V);
	spdlog::info("Call_V_R\t:%ld\n", Call_V_R);
	spdlog::info("Call\t%ld\t%ld\t%ld\n", Call[0], Call[1], Call[2]);
	spdlog::info("Recursive Number.\t{}\n", real);
	spdlog::info("Divide Record.\t{}\n", setDivideRecord.size());

END:
	if (poly)
		free(poly);
	if (div)
		free(div);

	if (freeSimpPoly)
	{
		for (i = 0; i < vec_simp_poly.size(); i++)
			free(vec_simp_poly[i]);
	}
	vec_simp_poly.clear();

	return;
}
#endif

//现在整体的初始化的调试工作结束
//对数据进行初始化
void SPRImpl::init(double *q, struct polyhedron *P, struct divide *T)
{
	//	int i, j, k, m;
	*q = 0.0;
	init_poly(P);
	init_tetra(T);

	GEOM_FUNC::exactinit_threadSafe();

	// 	for(i=0;i<MAX_S;i++)
	// 		for(j=0;j<MAX_S;j++)
	// 			for(k=0;k<MAX_S;k++)
	// 				for(m=0;m<MAX_S;m++)
	/*					Visible[i][j][k][m] = MAX_VOLUME;*/
}

//对多面体结构进行初始化
void SPRImpl::init_poly(struct polyhedron *P)
{
	int i, j;
	for (i = 0; i < MAX_SPR_POLY_SIZE; i++)
	{
		// 		for(j=0;j<3;j++)
		// 			P->v[i][j]=0;	//把初始点初始化为0
		P->used_v[i] = 0;
	}
	/*	P->num_v = 0;*/
	for (i = 0; i < MAX_SPR_POLY_SIZE; i++)
	{
		for (j = 0; j < 3; j++)
			P->t[i][j] = -1;		//把边面三角形中的各点索引值初始化为-1
		P->surface[i] = 0;
	}
	P->num_t = 0;
	P->num_vable = 0;
}

//对四面体划分分法进行初始化
void SPRImpl::init_tetra(struct divide *T)
{
	int i, j;
	for (i = 0; i < MAX_SPR_POLY_SIZE; i++)
	{
		for (j = 0; j < 4; j++)
			T->te[i][j] = -1;		//把四面体划分的各点索引值初始化为-1
	}
	T->num = 0;
	T->q = 0;
}
//构造一个八面体
//5个点，8个面，这个还是相当的清爽的
void SPRImpl::octahedron(struct polyhedron *P)
{
	int i;
	vertices[0][0] = 0;
	vertices[0][1] = 0;
	vertices[0][2] = 1;
	for (i = 1; i <= 4; i++)
		vertices[i][2] = 0;
	vertices[1][0] = vertices[2][0] = sqrt(2.0) / 2;
	vertices[3][0] = vertices[4][0] = -1 * sqrt(2.0) / 2;
	vertices[1][1] = vertices[4][1] = -1 * sqrt(2.0) / 2;
	vertices[2][1] = vertices[3][1] = sqrt(2.0) / 2;
	vertices[5][0] = vertices[5][1] = 0;
	vertices[5][2] = -1;
	num_vertices = 6;
	P->num_vable = 6;

	for (i = 0; i <= 3; i++)
		P->t[i][0] = 0;
	P->t[0][1] = 2;
	P->t[0][2] = 1;
	P->t[1][1] = 3;
	P->t[1][2] = 2;
	P->t[2][1] = 4;
	P->t[2][2] = 3;
	P->t[3][1] = 1;
	P->t[3][2] = 4;
	P->t[4][0] = 1;
	P->t[4][1] = 2;
	P->t[4][2] = 5;
	P->t[5][0] = 2;
	P->t[5][1] = 3;
	P->t[5][2] = 5;
	P->t[6][0] = 3;
	P->t[6][1] = 4;
	P->t[6][2] = 5;
	P->t[7][0] = 1;
	P->t[7][1] = 5;
	P->t[7][2] = 4;

	for (i = 0; i <= 5; i++)
		P->used_v[i] = 4;
	P->num_t = 8;
}

//构造一个小一点的球体
//遵从右手法则，就意味着在设计程序的时候要保证按照右手法则这个方向始终是指向里的！！顺序很重要！
//为了后面编程容易，还有一个法则，对于任何一个三角形，最小的那个索引值必须放到第一位！
//现在的这个球体是14个点，24个面。

void SPRImpl::sphere(struct polyhedron *P)
{
	int i;
	double R, r;
	R = 1;
	vertices[0][2] = R;
	vertices[0][1] = 0;
	vertices[0][0] = 0;
	r = sqrt(pow(R, 2) - pow(0.5*R, 2));
	for (i = 1; i <= 4; i++)
	{
		vertices[i][2] = R / 2;
		vertices[i][0] = r * cos(PI / 2 * (i - 1));
		vertices[i][1] = r * sin(PI / 2 * (i - 1));
	}
	for (i = 5; i <= 8; i++)
	{
		vertices[i][2] = 0;
		vertices[i][0] = R * cos(PI / 2 * (i - 5) + PI / 4);		/*这样可以把点错开，方便一点 */
		vertices[i][1] = R * sin(PI / 2 * (i - 5) + PI / 4);
	}
	for (i = 9; i <= 12; i++)
	{
		vertices[i][2] = -1 * R / 2;
		vertices[i][0] = r * cos(PI / 2 * (i - 9));
		vertices[i][1] = r * sin(PI / 2 * (i - 9));
	}
	vertices[13][2] = -1 * R;
	vertices[13][1] = 0;
	vertices[13][0] = 0;
	num_vertices = 14;
	P->num_vable = 14;

	for (i = 0; i <= 2; i++)			/*现在的任务是构造出初始的三角形出来*/
	{
		P->t[i][0] = 0;
		P->t[i][1] = i + 2;
		P->t[i][2] = i + 1;
	}
	P->t[3][0] = 0;
	P->t[3][1] = 1;
	P->t[3][2] = 4;
	for (i = 4; i <= 6; i++)
	{
		P->t[i][0] = i - 2;
		P->t[i][1] = i + 2;
		P->t[i][2] = i + 1;
	}
	P->t[7][0] = 1;
	P->t[7][1] = 5;
	P->t[7][2] = 8;
	for (i = 8; i <= 10; i++)
	{
		P->t[i][0] = i - 7;
		P->t[i][1] = i - 6;
		P->t[i][2] = i - 3;
	}
	P->t[11][0] = 1;
	P->t[11][1] = 8;
	P->t[11][2] = 4;
	/*上半球已经构造完成，下面是构造下半球，基本上是类似，整个球体应该是48个面，一直到索引号47*/
	for (i = 12; i <= 14; i++)
	{
		P->t[i][0] = i - 7;
		P->t[i][1] = i - 2;
		P->t[i][2] = i - 3;
	}
	P->t[15][0] = 8;
	P->t[15][1] = 9;
	P->t[15][2] = 12;
	for (i = 16; i <= 18; i++)
	{
		P->t[i][0] = i - 11;
		P->t[i][1] = i - 10;
		P->t[i][2] = i - 6;
	}
	P->t[19][0] = 5;
	P->t[19][1] = 9;
	P->t[19][2] = 8;
	for (i = 20; i <= 22; i++)
	{
		P->t[i][0] = i - 11;
		P->t[i][1] = i - 10;
		P->t[i][2] = 13;
	}
	P->t[23][0] = 9;
	P->t[23][1] = 13;
	P->t[23][2] = 12;
	P->num_t = 24;

	P->used_v[0] = 4;
	for (i = 1; i <= 4; i++)
		P->used_v[i] = 5;
	for (i = 5; i <= 8; i++)
		P->used_v[i] = 6;
	for (i = 9; i <= 12; i++)
		P->used_v[i] = 5;
	P->used_v[13] = 4;
}

#ifdef _SPR_USING_STD_VECTOR
int SPRImpl::record_subproblems(double iniQ, struct polyhedron *poly, int f, int maxSize, divide_record *divR,
	std::vector<SubDivideProblem*>& vec_sub_prob,
	Constraints *csts, int maxIntNum)
#else
int SPRImpl::record_subproblems(double iniQ, struct polyhedron *poly, int f, int maxSize, divide_record *divR,
	SubDivideProblem *subProbArray[], int *subProbNum,
	Constraints *csts, int maxIntNum)
#endif

{
	int r = 1;
	struct polyhedron *comp_poly;
#ifdef _SPR_USING_STD_VECTOR
	std::vector<struct polyhedron*> vec_simp_poly;
#else
	static polyhedron *simpPolyArray[MAX_SIMP_POLY_NUM];
	int simpPolyNum = 0;
#endif
	struct divide *Tq = NULL;
	TETRA ELE;
	int i, j, k, j1, j2, j3, j4, fIter = 0, nIter = 0;
	int ipoly, ipoly2, validpoly_num = 0;
	double quality_ELE;
	int test, found;
	int key = 0;
#if 0
#ifdef _SPR_MAP_DIVIDE_VALID
	std::map<int, int>	  mapDivideValid;
	std::map<int, int>::iterator it;
#else
	IntIntMap mapDivideValid;
	int mapValue = -1;
#endif /* _SPR_MAP_DIVIDE_VALID */
#endif
	bool bFaceSucc;
	bool freeSimpPoly;
	divide_record cur_record;
	SubDivideProblem *subProb = NULL;
	int intNum = 0; /* 子问题和约束的交点个数 */

	assert(poly->num_t > 4);

#ifdef _DEBUG
	//	savepls("sphere-sample-output.pls", poly, true);
	//	savepls("sphere-sample-output-final.pls", poly, false);
#endif

	for (i = 0; i < 3; i++)
		ELE[i] = poly->t[f][i];

	bFaceSucc = false;

#ifdef _SPR_USING_STD_VECTOR
	for (j = 0; j < num_vertices && vec_sub_prob.size() < maxSize; j++)
	{
#else
	*subProbNum = 0; /* initialise */
	for (j = 0; j < num_vertices && *subProbNum < maxSize; j++)
	{
#endif
		if (poly->used_v[j] == 0 || ELE[0] == j || ELE[1] == j || ELE[2] == j)
			continue;

		ELE[3] = j;

		j1 = ELE[0];
		j2 = ELE[1];
		j3 = ELE[2];
		j4 = ELE[3];

		/* 如果这个单元不被欢迎，则直接略过 */
		if (csts && csts->badVolmCst_size > 0)
		{
			for (k = 0; k < csts->badVolmCst_size; k++)
			{
				if ((j1 == csts->badVolmCsts[k][0] || j1 == csts->badVolmCsts[k][1] ||
					j1 == csts->badVolmCsts[k][2] || j1 == csts->badVolmCsts[k][3]) &&
					(j2 == csts->badVolmCsts[k][0] || j2 == csts->badVolmCsts[k][1] ||
						j2 == csts->badVolmCsts[k][2] || j2 == csts->badVolmCsts[k][3]) &&
						(j3 == csts->badVolmCsts[k][0] || j3 == csts->badVolmCsts[k][1] ||
							j3 == csts->badVolmCsts[k][2] || j3 == csts->badVolmCsts[k][3]) &&
							(j4 == csts->badVolmCsts[k][0] || j4 == csts->badVolmCsts[k][1] ||
								j4 == csts->badVolmCsts[k][2] || j4 == csts->badVolmCsts[k][3]))
					break;
			}
			if (k < csts->badVolmCst_size)
				continue;
		}

#if 0
		sort_quad_pair(&j1, &j2, &j3, &j4);
		key = make_quad_pair(j1, j2, j3, j4);
#ifdef _SPR_MAP_DIVIDE_VALID
		it = mapDivideValid.find(key);
		if (it == mapDivideValid.end())
#else
		if (!mapDivideValid.get(key, &mapValue))
#endif
		{/* 当前层未被测试过 */
			/* 确定小空腔是否已经四面体化 */
			if (divR->num < 0)
			{/* 为防止记录过长，仅记录前面若干层子问题 */
				create_div_record(divR, key, &cur_record, &found);
				if (found == 0)
				{/* 新的子问题 */
					if (is_tested_div_record(&cur_record))
						continue;
					else
						add_div_record(&cur_record);
				}
				else	/* 已有的子问题 */
					continue;
			}
			else
				cur_record.num = divR->num + 1;
		}
		else /* 当前层已经测试过 */
		{
#ifdef _SPR_MAP_DIVIDE_VALID
			if (it->second)
				bFaceSucc = true;	/* 是一个有效的测试 */
#else
			if (mapValue)
				bFaceSucc = true;
#endif

			continue;
		}
#endif /* #if 0 */

		/* 判断这个四面体是不是一个有效的四面体 */
		freeSimpPoly = true;
		comp_poly = (struct polyhedron*)malloc(sizeof(struct polyhedron));
		intNum = 0; /* 初始化为 0 */
#ifdef _SPR_USING_STD_VECTOR
		test = istetra(ELE, poly, iniQ, &quality_ELE, comp_poly, vec_simp_poly, f, &freeSimpPoly, csts, maxIntNum, &intNum);
#else
		test = istetra(ELE, poly, iniQ, &quality_ELE, comp_poly, simpPolyArray, &simpPolyNum, f, &freeSimpPoly, csts, maxIntNum, &intNum);
#endif

#if 0
#ifdef _SPR_MAP_DIVIDE_VALID
		mapDivideValid[key] = test;
#else
		mapDivideValid.put(key, test ? 1 : 0);
#endif
#endif /* #if 0 */

		if (test)
		{/* 是一个有效的四面体, 获得一个新的子问题 */

#ifndef _SPR_USING_STD_VECTOR
			if (*subProbNum >= MAX_SUB_PROBLEM_NUM)
			{
				printf("SPR Error: Increase MAX_SUB_PROBLEM_NUM(%d), please.\n", MAX_SUB_PROBLEM_NUM);
				exit(1);
			}
#endif
			bFaceSucc = true; /* 面片有解 */

			subProb = new SubDivideProblem;
			assert(subProb);

			subProb->fIdx = f;
			subProb->nIdx = j;
			subProb->tetraID = key;
			subProb->initQ = iniQ;
			subProb->tetrQ = quality_ELE;
			subProb->poly = poly;
			subProb->intNum = intNum; /* 记录交点个数 */
			subProb->facNum = spr_subprob_face_num;	/* 记录包含待删除边的面片数目 */

#ifdef _SPR_USING_STD_VECTOR
			validpoly_num = vec_simp_poly.size();
			for (ipoly = 0; ipoly < vec_simp_poly.size(); ipoly++)
			{
				if (vec_simp_poly[ipoly]->num_t < 4)
					validpoly_num--;
			}
			assert(validpoly_num >= 1);
			subProb->vec_simp_poly.resize(validpoly_num);
			for (ipoly = 0, ipoly2 = 0; ipoly < vec_simp_poly.size(); ipoly++)
			{
				if (vec_simp_poly[ipoly]->num_t < 4)
				{
					spdlog::info("An incorrect subproblem model. Output its parent model now (%s).\n", "fail-parent.pls");
					spdlog::info("Labels of digging face & vertice: f:{};v:{}.\n", f, j);
					savepls("fail-parent.pls", poly);
					assert(vec_simp_poly[ipoly]->num_t == 2);
				}
				else
					subProb->vec_simp_poly[ipoly2++] = vec_simp_poly[ipoly];
			}
			vec_sub_prob.push_back(subProb);
#else
			validpoly_num = simpPolyNum;
			for (ipoly = 0; ipoly < simpPolyNum; ipoly++)
			{
				if (simpPolyArray[ipoly]->num_t < 4)
					validpoly_num--;
			}
			assert(validpoly_num >= 1);
			subProb->simpPolyNum = validpoly_num;
			for (ipoly = 0, ipoly2 = 0; ipoly < simpPolyNum; ipoly++)
			{
				if (simpPolyArray[ipoly]->num_t < 4)
				{
					spdlog::info("An incorrect subproblem model. Output its parent model now (%s).\n", "fail-parent.pls");
					spdlog::info("Labels of digging face & vertice: f:{};v:{}.\n", f, j);
					savepls("fail-parent.pls", poly);
					assert(simpPolyArray[ipoly]->num_t == 2);
				}
				else
					subProb->simpPolyArray[ipoly2++] = simpPolyArray[ipoly];
			}
			subProbArray[(*subProbNum)++] = subProb;
#endif
			memcpy(&subProb->divR, &cur_record, sizeof(divide_record));
		}
		if (freeSimpPoly)
		{/* 意味着comp_poly的内存未转储到vec_simp_poly中，应该删除掉 */
			free(comp_poly);
		}
#ifdef _SPR_USING_STD_VECTOR
		vec_simp_poly.clear();
#else
		simpPolyNum = 0;
#endif /* _SPR_USING_STD_VECTOR */
	}

	/* 如果某个面没有找到一个合适的四面体递归下去，终止整个过程
		注意，有可能前面已经成功过一次 所以返回值可以是1，也可以是0 */
	if (!bFaceSucc)
		r = 0;

#ifdef _SPR_USING_STD_VECTOR
	return vec_sub_prob.size() * r;
#else
	return (*subProbNum) * r;
#endif /* _SPR_USING_STD_VECTOR */
	}

/* ---------------------------------------------------------------------------------------
 * 分两步做递归
 * 1. 记录当前子问题的所有可能分解方案（增加1个四面体单元情形下）
 *
 * 老的策略/2012年12月30前
 * 1)计算所有面片的质量（基于AFT启发式准则，具体参见JF Liu在IJNME的文章
 * 2)选取质量最差的面片
 * 3)记录和这一面片相关的所有子问题
 * 4)如果子问题数目为0，则返回无解
 *
 * 新的策略/2012年12月30修改
 * 尝试结合LIU和Ziji Wu提出的Last_Resort算法（发表在Finite Element in Analysis and Design, 2010, 46:74-83
 * 1)计算所有面片的质量（基于AFT启发式准则，具体参见JF Liu在IJNME的文章
 * 2)按面片质量从小到大排序
 * 3)依此计算所有面片的子问题集合
 * 4)返回数量最小的子问题集合
 * 5)如果存在子问题集合大小为0的情形，则返回无解
 *
 * 参数说明
 * iniQ			输入			单元的最差质量限制，一个有效的子问题其单元质量最小值必须大于iniQ
 * poly			输入			多面体(可能是原始多面体的子多面体）
 * divR			输入			已经形成的四面体单元对应的记录形式(用来避免重复计算，实际这段逻辑没有起作用)
 * vec_sub_prob	输出			子问题集合
 * csts			输入			约束集合（要求必须存在的线和面）
 * maxIntNum    输入            子问题新增加的面允许的最大交点数（和约束）
 * --------------------------------------------------------------------------------------*/
#ifdef _SPR_USING_STD_VECTOR
int SPRImpl::record_subproblems(double iniQ, struct polyhedron *poly, divide_record *divR,
	std::vector<SubDivideProblem*>& vec_sub_prob,
	Constraints *csts, int maxIntNum)
#else
int SPRImpl::record_subproblems(double iniQ, struct polyhedron *poly, divide_record *divR,
	SubDivideProblem *subProbArray[], int *subProbNum,
	Constraints *csts, int maxIntNum)
#endif /* _SPR_USING_STD_VECTOR */
{
	int f;
#ifdef _SPR_USING_STD_VECTOR
	std::vector<SubDivideProblem*> vec_sub_prob_temp;
#else
	SubDivideProblem* subProbArray_temp[MAX_SUB_PROBLEM_NUM];
	int subProbNum_temp = 0;
#endif /* _SPR_USING_STD_VECTOR */
	struct divide dive_temp;
	divide_record divR_temp;
	int ret, size, min_size = MAX_SPR_POLY_SIZE + 1;
	SubDivideProblem *subProb = NULL;
	int i, j, k;
	const int SWITCH_POLY_SIZE_1 = 4;
	const int SWITCH_POLY_SIZE_2 = 120;//120;	/* 当多面体比较小时，可以选取比较简单的策略 default value = 12*/
	int polySize = 0;
	IDX_QUALITY fQuality[MAX_SPR_POLY_SIZE];
	const int MAX_COMPARISON_COUNT = 30000;
	int iterCount;

	polySize = poly->num_t;
	assert(polySize > 4);

#ifdef _DEBUG
	//	savepls("sphere-sample-output.pls", poly, true);
	//	savepls("sphere-sample-output-final.pls", poly, false);
#endif

	if (polySize <= SWITCH_POLY_SIZE_1)
	{
#ifdef _SPR_USING_STD_VECTOR
		return record_subproblems(iniQ, poly, 0, MAX_SPR_POLY_SIZE, divR, vec_sub_prob, csts, maxIntNum);
#else
		return record_subproblems(iniQ, poly, 0, MAX_SPR_POLY_SIZE, divR, subProbArray, subProbNum, csts, maxIntNum);
#endif
	}
	else if (polySize <= SWITCH_POLY_SIZE_2)
	{
		//		sort_facets(poly);
		//		f = get_sorted_facet(poly, 0);
		f = calc_worst_facet(poly);

#ifdef _SPR_USING_STD_VECTOR
		return record_subproblems(iniQ, poly, f, MAX_SPR_POLY_SIZE, divR, vec_sub_prob, csts, maxIntNum);
#else
		return record_subproblems(iniQ, poly, f, MAX_SPR_POLY_SIZE, divR, subProbArray, subProbNum, csts, maxIntNum);
#endif /* _SPR_USING_STD_VECTOR */
	}
	else
	{
		//	if (divR)
		//		memcpy(&divR_temp, divR, sizeof(divide_record));

#if 0
		for (f = poly->num_t - 1; f >= 0; f--)
		{
#else
		for (i = 0; i < poly->num_t; i++)
		{
			fQuality[i].idx = i;
			fQuality[i].quality = calc_facet_quality(poly, i);
		}
		qsort(fQuality, poly->num_t, sizeof(IDX_QUALITY), compare_idx_quality);

		iterCount = MIN_VALUE(poly->num_t, MAX_COMPARISON_COUNT);
		for (k = 0; k < poly->num_t; k++)
		{
			f = fQuality[k].idx;
#endif

#ifdef _SPR_USING_STD_VECTOR
			ret = record_subproblems(iniQ, poly, f, min_size, divR,
				//divR ? &divR_temp : NULL, 
				vec_sub_prob_temp, csts, maxIntNum);

			size = vec_sub_prob_temp.size();
#else
			ret = record_subproblems(iniQ, poly, f, min_size, divR,
				//divR ? &divR_temp : NULL, 
				subProbArray_temp, &subProbNum_temp, csts, maxIntNum);

			size = subProbNum_temp;
#endif
			if (size < min_size)
			{/* 当前子问题规模更小，则返回它 */
				/* 先删除vec_sub_poly对应的内存 */
#ifdef _SPR_USING_STD_VECTOR
				for (i = 0; i < vec_sub_prob.size(); i++)
				{
					subProb = vec_sub_prob[i];
					for (j = 0; j < subProb->vec_simp_poly.size(); j++)
						free(subProb->vec_simp_poly[j]);
					delete subProb;
				}
				vec_sub_prob.clear();

				for (i = 0; i < size; i++)
				{
					subProb = vec_sub_prob_temp[i];
					vec_sub_prob.push_back(subProb);
				}
				vec_sub_prob_temp.clear();
#else
				for (i = 0; i < *subProbNum; i++)
				{
					subProb = subProbArray[i];
					for (j = 0; j < subProb->simpPolyNum; j++)
						free(subProb->simpPolyArray[j]);
					delete subProb;
				}

				for (i = 0; i < size; i++)
					subProbArray[i] = subProbArray_temp[i];
				*subProbNum = size;
				subProbNum_temp = 0;
#endif
				min_size = size;
			}
			else
			{
#ifdef _SPR_USING_STD_VECTOR
				for (i = 0; i < size; i++)
				{
					subProb = vec_sub_prob_temp[i];
					for (j = 0; j < subProb->vec_simp_poly.size(); j++)
						free(subProb->vec_simp_poly[j]);
					delete subProb;
				}
				vec_sub_prob_temp.clear();
#else
				for (i = 0; i < size; i++)
				{
					subProb = subProbArray_temp[i];
					for (j = 0; j < subProb->simpPolyNum; j++)
						free(subProb->simpPolyArray[j]);
					delete subProb;
				}
				subProbNum_temp = 0;
#endif
			}

#ifdef _SPR_USING_STD_VECTOR
			if (vec_sub_prob.size() <= 0)
				return 0;
#else
			if (*subProbNum <= 0)
				return 0;
#endif
		}

#ifdef _SPR_USING_STD_VECTOR
		return vec_sub_prob.size();
#else
		return *subProbNum;
#endif
		}
	}

int SPRImpl::optimal_fourFace(double iniQ, polyhedron *poly, divide *dive, divide_record *divR, Constraints *csts, int maxIntNum)
{
	int i, k, ELE[4];
	double quality_cur, volume;
	int r = 0;
	int intNum;

	//spdlog::info("Begin Calling optimal_fourFace.\n");
	assert(poly->num_t == 4);

	for (i = 0; i < 3; i++)
		ELE[i] = poly->t[0][i];

	for (i = 0; i < 4; i++)
	{
		if (poly->t[1][i] != ELE[0] &&
			poly->t[1][i] != ELE[1] &&
			poly->t[1][i] != ELE[2])
		{
			ELE[3] = poly->t[1][i];
			break;
		}
	}
	assert(i < 4);

	/* 是不是不受欢迎的bad element */
	if (csts && csts->badVolmCst_size > 0)
	{
		for (k = 0; k < csts->badVolmCst_size; k++)
		{
			if ((ELE[0] == csts->badVolmCsts[k][0] || ELE[0] == csts->badVolmCsts[k][1] ||
				ELE[0] == csts->badVolmCsts[k][2] || ELE[0] == csts->badVolmCsts[k][3]) &&
				(ELE[1] == csts->badVolmCsts[k][0] || ELE[1] == csts->badVolmCsts[k][1] ||
					ELE[1] == csts->badVolmCsts[k][2] || ELE[1] == csts->badVolmCsts[k][3]) &&
					(ELE[2] == csts->badVolmCsts[k][0] || ELE[2] == csts->badVolmCsts[k][1] ||
						ELE[2] == csts->badVolmCsts[k][2] || ELE[2] == csts->badVolmCsts[k][3]) &&
						(ELE[3] == csts->badVolmCsts[k][0] || ELE[3] == csts->badVolmCsts[k][1] ||
							ELE[3] == csts->badVolmCsts[k][2] || ELE[3] == csts->badVolmCsts[k][3]))
				return 0;
		}
	}

	quality_cur = isVisible_Search(ELE[0], ELE[1], ELE[2], ELE[3], &volume);
	if (volume > SPR_UNALLOWED_QUALITY)
	{
		if (spr_optimal_type == MIN_INT_NODE_NUMBER)
		{

			if (quality_cur < SPR_ZERO_QUALITY)
				quality_cur = SPR_ZERO_QUALITY;

#if 0
			assert(quality_cur < 1.0);
			quality_cur = 1.0 / (1.0 - quality_cur);
#endif
		}
		if (maxIntNum > 0 || (maxIntNum == 0 && quality_cur - iniQ > 0.0))
		{/* 质量符合要求 */
			r = 1;

#if 0
			if (spr_optimal_type == MIN_INT_NODE_NUMBER)
			{
				intNum = dive->intNum;
				memset(dive, 0, sizeof(struct divide));
				dive->intNum = intNum;
			}
			else
				memset(dive, 0, sizeof(struct divide));
#endif

			memset(dive, 0, sizeof(struct divide));
			mergeT(dive, ELE);
			dive->q = quality_cur;
			dive->intNum = 0;
		}
	}

	//	spdlog::info("After Calling optimal_fourFace.\n");
	return r;
}


/* ---------------------------------------------------------------------------------------
 * 分两步做递归
 * 1. 记录当前子问题的所有可能分解方案（增加1个四面体单元情形下）
 * 2. 启用某种寻优策略来确定优先的分解方案
 * 由于现有横向的比较，而不是盲目地递归，有可能会获得更好的求解效率
 * --------------------------------------------------------------------------------------*/
#ifdef _SPR_USING_STD_VECTOR
int SPRImpl::optimal_twoStep_firstCall(double iniQ, std::vector<polyhedron *>& vec_simp_poly, divide *dive, divide_record *divR,
	Constraints *csts, int maxIntNum)
#else
int SPRImpl::optimal_twoStep_firstCall(double iniQ, polyhedron *simpPolyArray[], int simpPolyNum, divide *dive, divide_record *divR,
	Constraints *csts, int maxIntNum)
#endif /* _SPR_USING_STD_VECTOR */
{
	int i, res;
	divide simp_div;
	int allowedFaceNum = 0;
	polyhedron *simpPoly = NULL;
	memset(dive, 0, sizeof(divide));
	dive->q = 1.0e10;

	/* 初始值 */
	allowedFaceNum = spr_allowed_face_num;

#ifdef _SPR_USING_STD_VECTOR
	for (i = 0; i < vec_simp_poly.size(); i++)
	{
		simpPoly = vec_simp_poly[i];
#else
	for (i = 0; i < simpPolyNum; i++)
	{
		simpPoly = simpPolyArray[i];
#endif /* _SPR_USING_STD_VECTOR */

#ifdef _SPR_ERROR_CHK
		assert(assert_poly_neig(simpPoly));
		assert(assert_simp_poly_edge(simpPoly));
#endif /* _SPR_ERROR_CHK */

		memset(divR, 0, sizeof(divide_record));
		memset(&simp_div, 0, sizeof(divide));
		if (spr_time_limit > 0.0)
			spr_time_start = SPRlogTime(); /* 开始计时 */

		/* simp_div的当前facNum个数要叠加前面的facNum数目 */
		if (spr_optimal_type == EDGE_REMOVAL)
		{
			simpPoly->allowedFaceNum = allowedFaceNum;
		}

		res = optimal_twoStep(iniQ, simpPoly, &simp_div, divR, csts, maxIntNum);

		if (res == 0)
			break;

		mergeT(dive, &simp_div);
		/* dive的当前facNum个数要用simp_div相应数值更新 */
		if (spr_optimal_type == EDGE_REMOVAL)
		{
			dive->facNum += simp_div.facNum;
			allowedFaceNum -= simp_div.facNum;
		}

		/* 取最小值 */
		if (dive->q > simp_div.q)
			dive->q = simp_div.q;
		if (spr_optimal_type == MIN_INT_NODE_NUMBER)
		{
			dive->intNum += simp_div.intNum; /* 总共有多少个交点 */
			maxIntNum -= simp_div.intNum;	/* 减去已经存在的交点要减去 */
		}

	}

	return res;
	}

/* ---------------------------------------------------------------------------------------
 * 分两步做递归
 * 1. 记录当前子问题的所有可能分解方案（增加1个四面体单元情形下）
 * 2. 启用某种寻优策略来确定优先的分解方案
 * 由于现有横向的比较，而不是盲目地递归，有可能会获得更好的求解效率
 * --------------------------------------------------------------------------------------*/
#ifdef _SPR_USING_STD_VECTOR
int SPRImpl::optimal_twoStep(double iniQ, std::vector<polyhedron *>& vec_simp_poly, divide *dive, divide_record *divR,
	Constraints *csts, int maxIntNum)
#else
int SPRImpl::optimal_twoStep(double iniQ, polyhedron *simpPolyArray[], int simpPolyNum, divide *dive, divide_record *divR,
	Constraints *csts, int maxIntNum)
#endif
{
	int i, res;
	divide simp_div;
	int allowedFaceNum = 0;
	polyhedron *simpPoly = NULL;

	memset(dive, 0, sizeof(divide));
	dive->q = 1.0e10;

	/* 初始值 */
	allowedFaceNum = spr_allowed_face_num;

#ifdef _SPR_USING_STD_VECTOR
	for (i = 0; i < vec_simp_poly.size(); i++)
	{
		simpPoly = vec_simp_poly[i];
#else
	for (i = 0; i < simpPolyNum; i++)
	{
		simpPoly = simpPolyArray[i];
#endif /* _SPR_USING_STD_VECTOR */

#ifdef _SPR_ERROR_CHK
		assert(assert_poly_neig(simpPoly));
		assert(assert_simp_poly_edge(simpPoly));
#endif /* _SPR_ERROR_CHK

		memset(divR, 0, sizeof(divide_record));
		memset(&simp_div, 0, sizeof(struct divide));
		/* simp_div的当前facNum个数要叠加前面的facNum数目 */
		if (spr_optimal_type == EDGE_REMOVAL)
		{
			simpPoly->allowedFaceNum = allowedFaceNum;
		}

		res = optimal_twoStep(iniQ, simpPoly, &simp_div, divR, csts, maxIntNum);

		if (res == 0)
			break;

		mergeT(dive, &simp_div);
		if (dive->q > simp_div.q)
			dive->q = simp_div.q;
		/* dive的当前facNum个数要用simp_div相应数值更新 */
		if (spr_optimal_type == EDGE_REMOVAL)
		{
			dive->facNum += simp_div.facNum;
			allowedFaceNum -= simp_div.facNum;
		}

		if (spr_optimal_type == MIN_INT_NODE_NUMBER)
		{
			dive->intNum += simp_div.intNum; /* 总共有多少个交点 */
			maxIntNum -= simp_div.intNum;	/* 减去已经存在的交点要减去 */
		}
	}

	return res;
	}


inline int is_point_in_triangle(APOINT pa, APOINT p1, APOINT p2, APOINT p3, double dist)
{

	APOINT v1a, v12, v1a2, v2a, v23, v2a3, v3a, v31, v3a1, norm;

	norm_3p(p1, p2, p3, norm);

	vec_2p(p1, pa, v1a);
	vec_2p(p1, p2, v12);

	vec_crop(v1a, v12, v1a2);
	if (vec_dotp(norm, v1a2) > dist)
		return 0;

	vec_2p(p2, pa, v2a);
	vec_2p(p2, p3, v23);
	vec_crop(v2a, v23, v2a3);
	if (vec_dotp(norm, v2a3) > dist)
		return 0;

	vec_2p(p3, pa, v3a);
	vec_2p(p3, p1, v31);
	vec_crop(v3a, v31, v3a1);
	if (vec_dotp(norm, v3a1) > dist)
		return 0;

	return 1;
}

inline double distPointToSegm(APOINT p, APOINT ps, APOINT pe, APOINT pnear, int *ip)
{
	APOINT vsp, vse, pj;
	double dsp, dse, dep, dist, dsptose;
	int i;

	vec_2p(ps, p, vsp);
	vec_2p(ps, pe, vse);
	dsp = vec_val(vsp);
	dse = vec_val(vse);

	if (dse == 0.0)
	{
		if (pnear)
			memcpy(pnear, ps, 3 * sizeof(double));
		if (ip)
			*ip = 0;
		return dsp;
	}

	vec_uni(vse);
	dsptose = vec_dotp(vsp, vse);
	if (dsptose >= 0.0 && dsptose <= dse)
	{
		for (i = 0; i < 3; i++)
			pj[i] = ps[i] + dsptose * vse[i];

		if (pnear)
			memcpy(pnear, pj, 3 * sizeof(double));
		if (ip)
			*ip = 1;
		return distance3D(p, pj);
	}
	else
	{
		dep = distance3D(pe, p);
		dist = MIN_VALUE(dsp, dep);

		if (dist == dsp)
		{
			if (pnear)
				memcpy(pnear, ps, 3 * sizeof(double));
			if (ip)
				*ip = 0;
		}
		else
		{
			if (pnear)
				memcpy(pnear, pe, 3 * sizeof(double));
			if (ip)
				*ip = 2;
		}
	}

	return dist;
}

#if 0
double distPointToTriangle(APOINT p, APOINT p1, APOINT p2, APOINT p3, APOINT pnear, int *ip)
{
	int ip1, ip2, ip3;
	APOINT pn1, pn2, pn3;
	double norm[3], d, pj[3], d1, d2, d3, dist;

	norm_3p(p1, p2, p3, norm);
	vec_uni(norm);

	d = -vec_dotp(norm, p1);
	pointProjectToPlane(p, d, norm, pj); /* 获得投影点 */

	if (is_point_in_triangle(pj, p1, p2, p3, 0.0))
	{
		if (pnear)
			memcpy(pnear, pj, 3 * sizeof(double));
		if (ip)
			*ip = 6;	/* 代表投影点在三角形内部 */
		return distance3D(pj, p);
	}

	d1 = distPointToSegm(p, p1, p2, pn1, &ip1);
	d2 = distPointToSegm(p, p2, p3, pn2, &ip2);
	d3 = distPointToSegm(p, p3, p1, pn3, &ip3);
	dist = MIN_VALUE(MIN_VALUE(d1, d2), d3);

	if (pnear || ip)
	{
		if (dist == d1)
		{
			if (pnear != NULL)
				memcpy(pnear, pn1, 3 * sizeof(double));
			if (ip != NULL) *ip = ip1;
		}
		else if (dist == d2)
		{
			if (pnear != NULL)
				memcpy(pnear, pn2, 3 * sizeof(double));
			if (ip != NULL)
				*ip = ip2 + 2;
		}
		else
		{
			if (pnear != NULL)
				memcpy(pnear, pn3, 3 * sizeof(double));
			if (ip != NULL)
			{
				*ip = ip3 + 4;
				if (*ip == 6)
				{
					*ip = 0;
				}
			}
		}
	}

	return dist;
}
#else
typedef REAL PX_REAL;
/******************************************************************/
//  FUNCTION Definition: PXDistance2Triangle
int PXDistance2Triangle(PX_REAL const queryPoint[3], PX_REAL const triangle[9],
	PX_REAL*  resultPointRef, PX_REAL*  resultPoint,
	PX_REAL*  distance, int*  region0Flag)
{
	/* Determines the distance to an arbitrary triangle from a point in 3D by
	   analytically minimizing the square distance from a point to the triangle.
	   This distance is Q = |T(s,t) - queryPoint|^2
	   Q = a*s^2 + 2*b*s*t + c*t^2 + 2*d*s + 2*e*t + f
	   T(s,t) = B + s*E0 + t*E1;
	   B = v1; E0 = (v2-v1)/norm(v2-v1); E1 = (v3-v1)/norm(v3-v1) */
	   /* ASSUMPTION: input IS A TRIANGLE (else this will fail ungracefully) */

	   /* INPUTS */
	   /* queryPoint[3]: [x,y,z] coordinates of the query point
		  triangle[9]: vertices of the triangle in this order:
					   [v0x,v0y,v0z,v1x,v1y,v1z,v2x,v2y,v2z] */

					   /* OUTPUTS */
					   /* resultPointRef: [s,t]: the nearest point on the triangle in reference coords
						  resultPoint: [x,y,z]: the nearest point on the triangle in physical coords
						  distance: distance from input queryPoint to input triangle
						  region0Flag: PXE_True if queryPoint lies strictly inside region 0, false otherwise. */

						  /*  The algorithm is based on  */
						  /*  "David Eberly, 'Distance Between Point and Triangle in 3D', */
						  /*  Geometric Tools, LLC, (1999)" */

						  /*        ^t
							  \     |
							   \reg2|
								\   |
								 \  |
								  \ |
								   \|
									*v2
									|\
									| \
							  reg3  |  \ reg1
									|	\
									|reg0\
									|     \
									|      \ v1
							 -------*-------*------->s
									|v0      \
							  reg4  | reg5    \ reg6
						  */

	PX_REAL B[3], E0[3], E1[3], D[3];
	PX_REAL invnormE0, invnormE1;
	PX_REAL a, b, c, d, e, f, tmp0, tmp1, numer, denom, s, t, det, invdet, distSquared;

	/* Rewrite triangle in normal form */
	memcpy(B, triangle, 3 * sizeof(PX_REAL));
	memcpy(E0, triangle + 3, 3 * sizeof(PX_REAL));
	memcpy(E1, triangle + 6, 3 * sizeof(PX_REAL));
	E0[0] -= B[0]; E0[1] -= B[1]; E0[2] -= B[2];
	E1[0] -= B[0]; E1[1] -= B[1]; E1[2] -= B[2];
	/*   invnormE0 = 1/sqrt(E0[0]*E0[0] + E0[1]*E0[1] + E0[2]*E0[2]); */
	/*   invnormE1 = 1/sqrt(E1[0]*E1[0] + E1[1]*E1[1] + E1[2]*E1[2]); */
	/*   E0[0] *= invnormE0; E0[1] *= invnormE0; E0[2] *= invnormE0; */
	/*   E1[0] *= invnormE1; E1[1] *= invnormE1; E1[2] *= invnormE1; */
	D[0] = B[0] - queryPoint[0];
	D[1] = B[1] - queryPoint[1];
	D[2] = B[2] - queryPoint[2];

	a = E0[0] * E0[0] + E0[1] * E0[1] + E0[2] * E0[2];
	b = E0[0] * E1[0] + E0[1] * E1[1] + E0[2] * E1[2];
	c = E1[0] * E1[0] + E1[1] * E1[1] + E1[2] * E1[2];
	d = E0[0] * D[0] + E0[1] * D[1] + E0[2] * D[2];
	e = E1[0] * D[0] + E1[1] * D[1] + E1[2] * D[2];
	f = D[0] * D[0] + D[1] * D[1] + D[2] * D[2];

	det = a * c - b * b;
	s = b * e - c * d;
	t = b * d - a * e;

	if (region0Flag)
		(*region0Flag) = 0;

	if (s + t <= det)
	{
		if (s < 0)
		{
			if (t < 0)
			{ /* region 4*/
				if (d < 0)
				{
					t = 0.0;
					if (-d >= a)
					{
						s = 1.0;
						distSquared = a + 2 * d + f;
					}
					else
					{
						s = -d / a;
						distSquared = d * s + f;
					}
				}
				else
				{
					s = 0.0;
					if (e >= 0.0)
					{
						t = 0.0;
						distSquared = f;
					}
					else
					{
						if (-e >= c)
						{
							t = 1.0;
							distSquared = c + 2 * e + f;
						}
						else
						{
							t = -e / c;
							distSquared = e * t + f;
						}
					}
				}
			}
			else
			{ /*region 3 */
				/* F(t) = Q(0,t) = ct^2 + 2et + f */
				/* F'(t)/2 = ct+e */
				/* F'(T) = 0 when T = -e/c */
				/* Note that the logic for region 5 is similar */
				s = 0.0;
				if (e >= 0.0)
				{
					t = 0.0;
					distSquared = f;
				}
				else
				{
					if (-e >= c)
					{
						t = 1.0;
						distSquared = c + 2 * e + f;
					}
					else
					{
						t = -e / c;
						distSquared = e * t + f;
					}
				}
			}
		}
		else
		{
			if (t < 0)
			{ /* region 5*/
				t = 0.0;
				if (d >= 0)
				{
					s = 0.0;
					distSquared = f;
				}
				else
				{
					if (-d >= a)
					{
						s = 1.0;
						distSquared = a + d * s + f;
					}
					else
					{
						s = -d / a;
						distSquared = d * s + f;
					}
				}
			}
			else
			{ /*region 0*/
				if (region0Flag)
					(*region0Flag) = 0;
				invdet = 1 / det;
				s *= invdet;
				t *= invdet;
				distSquared = s * (a*s + b * t + 2 * d) + t * (b*s + c * t + 2 * e) + f;
			}
		}
	}
	else
	{
		if (s < 0)
		{ /*region 2*/
			tmp0 = b + d; tmp1 = c + e;
			if (tmp1 > tmp0)
			{ /* min is on edge s+t = 1 */
				numer = tmp1 - tmp0;
				denom = a - 2 * b + c;
				if (numer >= denom)
				{
					s = 1.0;
					t = 0.0;
					distSquared = a + 2 * b + f;
				}
				else
				{
					s = numer / denom;
					t = 1.0 - s;
					distSquared = s * (a*s + b * t + 2 * d) + t * (b*s + c * t + 2 * e) + f;
				}
			}
			else
			{ /* min is on s = 0 */
				s = 0.0;
				if (tmp1 <= 0.0)
				{
					t = 1.0;
					distSquared = c + 2 * e + f;
				}
				else
				{
					if (e >= 0.0)
					{
						t = 0.0;
						distSquared = f;
					}
					else
					{
						t = -e / c;
						distSquared = e * t + f;
					}
				}
			}
		}
		else
		{
			if (t < 0)
			{ /*region 6*/
				tmp0 = b + e; tmp1 = a + d;
				if (tmp1 > tmp0)
				{
					numer = tmp1 - tmp0;
					denom = a - 2 * b + c;
					if (numer >= denom)
					{
						t = 1.0;
						s = 0.0;
						distSquared = c + 2 * e + f;
					}
					else
					{
						t = numer / denom;
						s = 1.0 - t;
						distSquared = s * (a*s + b * t + 2 * d) + t * (b*s + c * t + 2 * e) + f;
					}
				}
				else
				{
					t = 0.0;
					if (tmp1 <= 0.0)
					{
						s = 1.0;
						distSquared = a + 2 * d + f;
					}
					else
					{
						if (d >= 0)
						{
							s = 0.0;
							distSquared = f;
						}
						else
						{
							s = -d / a;
							distSquared = d * s + f;
						}
					}
				}
			}
			else
			{ /*region 1*/
				/*  F(s) = Q(s,1-s) = (a-2b+c)s^2 + 2(b-c+d-e)s + (c+2e+f) */
				/*  F(s)/2 = (a-2b+c)s + (b-c+d-e) */
				/*  F(S) = 0 when S = (c+e-b-d)/(a-2b+c) */
				/*  a-2b+c = |E0-E1|^2 > 0, so only sign of c+e-b-d need be considered */
				numer = c + e - b - d;
				if (numer <= 0.0)
				{
					s = 0.0;
					t = 1.0;
					distSquared = c + 2 * e + f;
				}
				else
				{
					denom = a - 2 * b + c;
					if (numer >= denom)
					{
						s = 1.0;
						t = 0.0;
						distSquared = a + 2 * d + f;
					}
					else
					{
						s = numer / denom;
						t = 1.0 - s;
						distSquared = s * (a*s + b * t + 2 * d) + t * (b*s + c * t + 2 * e) + f;
					}
				}
			}
		}
	}
	if (distSquared < 0)
		distSquared = 0.0;

	if (distance != NULL)
		(*distance) = sqrt(distSquared);

	if (resultPointRef != NULL)
	{
		resultPointRef[0] = s;
		resultPointRef[1] = t;
	}

	if (resultPoint != NULL)
	{
		resultPoint[0] = B[0] + s * E0[0] + t * E1[0];
		resultPoint[1] = B[1] + s * E0[1] + t * E1[1];
		resultPoint[2] = B[2] + s * E0[1] + t * E1[2];
	}

	return 0;

}

double distPointToTriangle(APOINT p, APOINT p1, APOINT p2, APOINT p3, APOINT pnear, int *ip)
{
	double dist = 0.0, trianglePs[] = { p1[0], p1[1], p1[2], p2[0], p2[1], p2[2], p3[0], p3[1], p3[2] };
	PXDistance2Triangle(p, trianglePs, NULL, pnear, &dist, ip);
	return dist;
}
#endif



double distSegmToTriangle(double p0[3], double p1[3], double pa[3], double pb[3], double pc[3])
{
	double pm[3];
	const int n = 5;
	int i;
	double dlt = 1. / n, w1, w2;
	double dist, mindst;

	mindst = distPointToTriangle(p1, pa, pb, pc, 0, 0);

	for (i = 0; i < n; i++)
	{
		w2 = dlt * i;
		w1 = 1.0 - w2;

		pm[0] = w1 * p0[0] + w2 * p1[0];
		pm[1] = w1 * p0[1] + w2 * p1[1];
		pm[2] = w1 * p0[2] + w2 * p1[2];
		dist = distPointToTriangle(pm, pa, pb, pc, 0, 0);
		if (dist < mindst)
			mindst = dist;
	}
	return mindst;
}

// inline double compEdgeNarrowDist(polyhedron *simp_poly, int a, int b)
// {
// 	double dist, mindst = 1.0e20, lab, vab[3], norm1[3], norm2[3], norm3[3], vtad[3];
// 	int i, j, jIdx, i1, i2, i3, i4, numc;
// 	static double cos_th = 0.0;//sqrt(3)/2.0;
// 
// 	vec_2p(vertices[a], vertices[b], vab);
// 	lab = vec_val(vab);
// 	assert(lab > 0.0);
// 
// 	for (i = 0; i < simp_poly->num_t; i++)
// 	{
// 		numc = 0;
// 		for (j = 0; j < 3; j++)
// 		{
// 			i1 = simp_poly->t[i][j];
// 			if (i1 == a || i1 == b)
// 			{
// 				jIdx = j;
// 				numc++;
// 			}
// 		}
// 
// 		if (numc == 0)
// 		{/* 边和面没有公共节点 */
// 			dist=distSegmToTriangle(
// 				vertices[a], vertices[b],
// 				vertices[simp_poly->t[i][0]],
// 				vertices[simp_poly->t[i][1]],
// 				vertices[simp_poly->t[i][2]]) / lab;
// 
// 			if(dist < mindst) 
// 				mindst=dist;
// 		}
// 		else if (numc == 1)
// 		{/* 边和面有一个公共节点 */
// 			
// 			i1 = simp_poly->t[i][jIdx];
// 			i2 = simp_poly->t[i][(jIdx+1)%3];
// 			i3 = simp_poly->t[i][(jIdx+2)%3];
// 			i4 = a + b - i1;
// 			if (isVisible_Search(vertices[i1], vertices[i2], vertices[i3], vertices[i4]) <= 0.0)
// 			/* 如果边和面形成的是负体积，则不要计算 */
// 				continue;
// 			
// 			/* 确定两个二面角小于30度 */
// 			norm_3p(vertices[i1], vertices[i2], vertices[i3], norm1);
// 			norm_3p(vertices[i1], vertices[i4], vertices[i2], norm2);
// 			norm_3p(vertices[i1], vertices[i3], vertices[i4], norm3);
// 			vec_uni(norm1);
// 			vec_uni(norm2);
// 			vec_uni(norm3);
// 
// 			if (-vec_dotp(norm1, norm2) >= cos_th && 
// 				-vec_dotp(norm1, norm3) >= cos_th)
// 			{
// 				vec_2p(vertices[i1], vertices[i4], vtad);
// 				dist = vec_dotp(vtad, norm1);
// 				dist /= vec_val(vtad);
// 
// 				if (dist < mindst) 
// 					mindst=dist;
// 			}
// 		}
// 	}
// 
// 	return mindst;
// }

//计算新边a-b的距离   zhaodawei modify 2010-07-16
double SPRImpl::compEdgeNarrowDist(polyhedron *simp_poly, int a, int b, polyhedron *old_poly, double mindst)
{
	double dist, lab, lnorm1, lext, vab[3], norm1[3], norm2[3], norm3[3], vtad[3], minx, miny, minz, maxx, maxy, maxz;
	double *p1, *p2, *p3, *pa, *pb;
	int i, j, jIdx, i1, i2, i3, i4, numc;
	const double cos_th = 0.5;//sqrt(3)/2.0;

	vec_2p(vertices[a], vertices[b], vab);
	lab = vec_val(vab);  //计算边vab的长度
	assert(lab > 0.0);
	//这个循环是计算新边的距离
	for (i = 0; i < old_poly->num_t; i++)
	{
		numc = 0;
		for (j = 0; j < 3; j++)
		{
			i1 = old_poly->t[i][j];
			if (i1 == a || i1 == b)
			{
				jIdx = j;
				numc++;
			}
		}

		if (numc == 0)
		{/* 边和面没有公共节点 */
			//计算边a-b与面simp_poly->t[i][0]、simp_poly->t[i][1]、simp_poly->t[i][2]的距离
			p1 = vertices[old_poly->t[i][0]];
			p2 = vertices[old_poly->t[i][1]];
			p3 = vertices[old_poly->t[i][2]];
			pa = vertices[a];
			pb = vertices[b];
			/* 利用box做屏蔽 */
			lext = mindst * lab;
			minx = MIN_VALUE(p1[0], MIN_VALUE(p2[0], p3[0])) - lext;
			miny = MIN_VALUE(p1[1], MIN_VALUE(p2[1], p3[1])) - lext;
			minz = MIN_VALUE(p1[2], MIN_VALUE(p2[2], p3[2])) - lext;
			maxx = MAX_VALUE(p1[0], MAX_VALUE(p2[0], p3[0])) + lext;
			maxy = MAX_VALUE(p1[1], MAX_VALUE(p2[1], p3[1])) + lext;
			maxz = MAX_VALUE(p1[2], MAX_VALUE(p2[2], p3[2])) + lext;

			if (((pa[0] >= minx && pa[0] <= maxx) &&
				(pa[1] >= miny && pa[1] <= maxy) &&
				(pa[2] >= minz && pa[2] <= maxz)) ||
				((pb[0] >= minx && pb[0] <= maxx) &&
				(pb[1] >= miny && pb[1] <= maxy) &&
					(pb[2] >= minz && pb[2] <= maxz)))
			{/* a和b在盒子里 */
				dist = distSegmToTriangle(pa, pb, p1, p2, p3) / lab;   //这里为什么要除以lab

				if (dist < mindst)
					mindst = dist;
			}
		}
		else if (numc == 1)
		{/* 边和面有一个公共节点 */

			i1 = old_poly->t[i][jIdx];
			i2 = old_poly->t[i][(jIdx + 1) % 3];
			i3 = old_poly->t[i][(jIdx + 2) % 3];
			i4 = a + b - i1;

			norm_3p(vertices[i1], vertices[i2], vertices[i3], norm1);
			vec_2p(vertices[i1], vertices[i4], vtad);
			dist = vec_dotp(vtad, norm1);

			if (dist > 0.0)
			{/* 如果边和面形成的是负体积，则不要计算 */
				lnorm1 = vec_val(norm1);
				dist /= vec_val(vtad) * lnorm1;
				if (dist < mindst)
				{
					norm1[0] /= lnorm1;
					norm1[1] /= lnorm1;
					norm1[2] /= lnorm1;
					/* 确定两个二面角小于30度 */
					norm_3p(vertices[i1], vertices[i4], vertices[i2], norm2);
					vec_uni(norm2);
					if (-vec_dotp(norm1, norm2) >= cos_th)
					{
						norm_3p(vertices[i1], vertices[i3], vertices[i4], norm3);
						vec_uni(norm3);
						if (-vec_dotp(norm1, norm3) >= cos_th)
							mindst = dist;
					}
				}
			}
		}
	}

	return mindst;
}


// inline double compNarrowDist(polyhedron *simp_poly)
// {
// 	int i, j, fi, ng, a, b;
// 	double dist, mindst = 1.0e20;
// 
// 	/* 计算所有新边对应的最小距离 */
// 	for (i = 0; i < simp_poly->newf_ct; i++)
// 	{
// 		fi = simp_poly->newf_id[i];
// 		for (j = 0; j < 3; j++)
// 		{
// 			ng = simp_poly->n[fi][j];
// 
// 			if (simp_poly->newf_fg[ng] != 0 && ng > fi)
// 			{/* also a new face */
// 				a = simp_poly->t[fi][(j+1)%3];
// 				b = simp_poly->t[fi][(j+2)%3];
// 				dist = compEdgeNarrowDist(simp_poly, a, b);
// 
// 				if (dist < mindst)
// 					mindst = dist;
// 			}
// 		}
// 	}
// 
// 	/* 计算所有新面对应的最小距离 */
// 	for (i = 0; i < simp_poly->newf_ct; i++)
// 	{
// 		fi = simp_poly->newf_id[i];
// 		dist = calc_facet_quality_rel(simp_poly, fi);
// 
// 		if (dist < mindst)
// 			mindst = dist;
// 	}
// 
// 	return mindst;
// }

//zhaodawei modify  2010-07-16
double SPRImpl::compNarrowDist(polyhedron *simp_poly, polyhedron *old_poly)
{
#if 1
	double dist, mindst = 1.0e100;
	int a, b;
	int i, j, ng, fIdx;
#else
	double glob_mindst = 1.0e100, height;
	int c, d;
	int *poly_base_t = NULL, *poly_base_n = NULL, *poly_neig_t = NULL;
	double *pa, *pb, *pc, *pd, vad[3];
	const double cos_th = 0.5;//sqrt(3)/2.0;
	int worstf = 0;
	static double norms[MAX_SPR_POLY_SIZE][4]; /* record the normal vector and its magnitude**2 */
	static char faceNormalised[MAX_SPR_POLY_SIZE];
	double vab[3], vac[3], vbc[3], *vabc = NULL, *vadb = NULL;
	double l1, l2, l3, lmax;
#endif

	/* 计算所有新边对应的最小距离 */
	for (i = 0; i < simp_poly->newf_ct; i++)
	{
		fIdx = simp_poly->newf_id[i]; //得到新面的编号
		for (j = 0; j < 3; j++)
		{
			ng = simp_poly->n[fIdx][j]; //该新面的邻接单元编号

			if (simp_poly->newf_fg[ng] != 0 && ng > fIdx)  //等于1表示是新面 为什么必须有大于关系
			{/* also a new face */
				a = simp_poly->t[fIdx][(j + 1) % 3];  //a-b是新边
				b = simp_poly->t[fIdx][(j + 2) % 3];
				dist = compEdgeNarrowDist(simp_poly, a, b, old_poly, mindst);

				if (dist < mindst)
					mindst = dist;
			}
		}
	}


	/* 计算所有新面对应的最小距离 */
	for (i = 0; i < simp_poly->newf_ct; i++)
	{
		fIdx = simp_poly->newf_id[i];
#if 1
		dist = calc_facet_quality_rel(simp_poly, fIdx, old_poly, mindst);
#else
		/* 法向求出来了吗 */
		poly_base_t = &simp_poly->t[fIdx][0];
		a = poly_base_t[0];
		b = poly_base_t[1];
		c = poly_base_t[2];

		pa = vertices[a];
		pb = vertices[b];
		pc = vertices[c];

		vab[0] = pb[0] - pa[0];
		vab[1] = pb[1] - pa[1];
		vab[2] = pb[2] - pa[2];
		vac[0] = pc[0] - pa[0];
		vac[1] = pc[1] - pa[1];
		vac[2] = pc[2] - pa[2];
		//		vec_2p(pa,pb,vab);
		//		vec_2p(pa,pc,vac);
		vec_crop(vab, vac, vabc);
		vabc[3] = sqrt(vabc[0] * vabc[0] + vabc[1] * vabc[1] + vabc[2] * vabc[2]);
		if (vabc[3] != 0.0)
		{/* degenerate elements */
			vabc[0] /= vabc[3];
			vabc[1] /= vabc[3];
			vabc[2] /= vabc[3];
		}

		vbc[0] = pc[0] - pb[0];
		vbc[1] = pc[1] - pb[1];
		vbc[2] = pc[2] - pb[2];

		l1 = vab[0] * vab[0] + vab[1] * vab[1] + vab[2] * vab[2];
		l2 = vac[0] * vac[0] + vac[1] * vac[1] + vac[2] * vac[2];
		l3 = vbc[0] * vbc[0] + vbc[1] * vbc[1] + vbc[2] * vbc[2];
		lmax = MAX_VALUE(MAX_VALUE(l1, l2), l3);
		mindst = height = vabc[3] / sqrt(lmax);

		poly_base_n = &simp_poly->n[fIdx][0];
		for (i = 0; i < 3; i++)
		{
			ng = poly_base_n[i];
			poly_neig_t = &simp_poly->t[ng][0];

			/* 确保d点在上方 a, b等于poly_base_t[(i+1)%3]/poly_base_t[(i+2)%3]*/
			a = poly_base_t[(i + 1) % 3];
			b = poly_base_t[(i + 2) % 3];
			d = poly_neig_t[0] + poly_neig_t[1] + poly_neig_t[2] - a - b;

			/* 确定d在上方 */
			vad[0] = vertices[d][0] - vertices[a][0];
			vad[1] = vertices[d][1] - vertices[a][1];
			vad[2] = vertices[d][2] - vertices[a][2];
			height = vec_dotp(vad, vabc);

			if (height > 0.0 && height < MIN_VALUE(mindst, glob_mindst))
			{
				vadb = &norms[ng][0];
				if (faceNormalised[ng] == 0)
				{
					/* 求相邻面片的法向, 面片法向朝里 */
					vab[0] = vertices[b][0] - vertices[a][0];
					vab[1] = vertices[b][1] - vertices[a][1];
					vab[2] = vertices[b][2] - vertices[a][2];

					vec_crop(vad, vab, vadb);
					vadb[3] = sqrt(vadb[0] * vadb[0] + vadb[1] * vadb[1] + vadb[2] * vadb[2]);
					if (vadb[3] != 0.0)
					{/* degenerate elements */
						vadb[0] /= vadb[3];
						vadb[1] /= vadb[3];
						vadb[2] /= vadb[3];
					}
					faceNormalised[ng] = 1;
				}

				if (-vec_dotp(vabc, vadb) >= cos_th)
				{/* 只关心小于90度（cos_th = 0.0)或30度（cos_th = sqrt(3.0)/2的二面角 */
					mindst = height;
				}
			}
		}

		if (mindst < glob_mindst)
		{
			glob_mindst = mindst;
			worstf = fIdx;
		}
#endif
		if (dist < mindst)
			mindst = dist;
	}

	return mindst;
}

int SPRImpl::eval_one_subprob_LJF(SubDivideProblem *subProb)
{
	int i;
	polyhedron *simp_poly;
	double dist = 0.0, mindst = 1.0e20;

#ifdef _SPR_USING_STD_VECTOR
	for (i = 0; i < subProb->vec_simp_poly.size(); i++)
	{
		simp_poly = subProb->vec_simp_poly[i];
#else
	for (i = 0; i < subProb->simpPolyNum; i++)
	{
		simp_poly = subProb->simpPolyArray[i];
#endif
		assert(simp_poly);

#ifdef _DEBUG
		//		savepls("sphere-sample-sub-output.pls", simp_poly, true);
		//		savepls("sphere-sample-sub-output-final.pls", simp_poly, false);
#endif
		if (simp_poly->newf_ct > 0)
		{/* 有新的面 */
		/*	dist = compNarrowDist(simp_poly);*/
			dist = compNarrowDist(simp_poly, subProb->poly);   //zhaodawei modify 2010-07-16
			if (dist < mindst)
				mindst = dist;
		}
	}

	subProb->diviQ = mindst;

	return 1;
	}

int SPRImpl::eval_one_subprob(SubDivideProblem *subProb)
{
	int eval_method = 1;

	switch (eval_method)
	{
	case 0: /* 贪婪法 */
		subProb->diviQ = subProb->tetrQ;
		break;
	case 1:
		/* 刘剑飞的算法 */
		eval_one_subprob_LJF(subProb);
		break;
	}

	return 1;
}

#ifdef _SPR_USING_STD_VECTOR
int SPRImpl::eval_subproblems(std::vector<SubDivideProblem*>& vec_sub_prob)
#else
int SPRImpl::eval_subproblems(SubDivideProblem *subProbArray[], int subProbNum)
#endif
{
	int i;
	SubDivideProblem *subProb = NULL;

#ifdef _SPR_USING_STD_VECTOR
	for (i = 0; i < vec_sub_prob.size(); i++)
	{
		subProb = vec_sub_prob[i];
#else
	for (i = 0; i < subProbNum; i++)
	{
		subProb = subProbArray[i];
#endif /* _SPR_USING_STD_VECTOR */
		assert(subProb);
		eval_one_subprob(subProb);
	}

	return 1;
	}

bool SPRImpl::comp_subproblems(const SubDivideProblem *prob1, const SubDivideProblem *prob2)
{
	/* ----------------------------------------------
	 * 包含多个子多面体的子问题总是优先
	 * 在VS2010编译时，要采用所谓的strict weak ordering
	 * 准则，即相等时，返回false
	 * --------------------------------------------*/
#if 1
	int size1, size2;

#if 1
	if (spr_optimal_type == MIN_INT_NODE_NUMBER)
	{
		if (prob1->intNum != prob2->intNum)
			return prob1->intNum < prob2->intNum;
	}
#endif /* #if 1 */

#ifdef _SPR_USING_STD_VECTOR
	size1 = prob1->vec_simp_poly.size();
	size2 = prob2->vec_simp_poly.size();
#else
	size1 = prob1->simpPolyNum;
	size2 = prob2->simpPolyNum;
#endif /*  _SPR_USING_STD_VECTOR */
	if (size1 != size2)
	{
		return size1 > size2;
	}
	else
		return prob1->diviQ > prob2->diviQ;
#else
	return prob1->diviQ > prob2->diviQ;
#endif
}

bool SPRImpl::comp_subproblemsV2(const void *arg1, const void *arg2)
{
	/* ----------------------------------------------
	 * 包含多个子多面体的子问题总是优先
	 * 在VS2010编译时，要采用所谓的strict weak ordering
	 * 准则，即相等时，返回false
	 * --------------------------------------------*/
	SubDivideProblem *prob1 = *((SubDivideProblem **)arg1);
	SubDivideProblem *prob2 = *((SubDivideProblem **)arg2);
#if 1 
	int size1, size2;
	double diviQ1, diviQ2;

#if 1
	if (spr_optimal_type == MIN_INT_NODE_NUMBER)
	{
		if (prob1->intNum != prob2->intNum)
			return prob1->intNum <= prob2->intNum ? false : true;
	}
#endif

#ifdef _SPR_USING_STD_VECTOR
	size1 = prob1->vec_simp_poly.size();
	size2 = prob2->vec_simp_poly.size();
#else
	size1 = prob1->simpPolyNum;
	size2 = prob2->simpPolyNum;
#endif
	if (size1 != size2)
	{
		return size1 >= size2 ? false : true;
	}
	else
	{
		return prob1->diviQ >= prob2->diviQ ? false : true;
	}
#else
	return prob1->diviQ > prob2->diviQ ? -1 : 1;
#endif
}

//************************************************zhaodawei add  2010-07-16

void swap_nprall(int a, int b, std::vector<SubDivideProblem*>& vec_sub_prob)
{
	polyhedron *tmp1 = vec_sub_prob[a]->poly;
	int tmp2 = vec_sub_prob[a]->tetraID;
	int tmp3 = vec_sub_prob[a]->fIdx;
	int tmp4 = vec_sub_prob[a]->nIdx;
	double tmp5 = vec_sub_prob[a]->initQ;
	double tmp6 = vec_sub_prob[a]->tetrQ;
	double tmp7 = vec_sub_prob[a]->diviQ;
#ifdef _SPR_USING_STD_VECTOR
	std::vector<polyhedron*> tmp8 = vec_sub_prob[a]->vec_simp_poly;
#endif
	divide_record tmp9 = vec_sub_prob[a]->divR;

	vec_sub_prob[a]->poly = vec_sub_prob[b]->poly;
	vec_sub_prob[a]->tetraID = vec_sub_prob[b]->tetraID;
	vec_sub_prob[a]->fIdx = vec_sub_prob[b]->fIdx;
	vec_sub_prob[a]->nIdx = vec_sub_prob[b]->nIdx;
	vec_sub_prob[a]->initQ = vec_sub_prob[b]->initQ;
	vec_sub_prob[a]->tetrQ = vec_sub_prob[b]->tetrQ;
	vec_sub_prob[a]->diviQ = vec_sub_prob[b]->diviQ;
#ifdef _SPR_USING_STD_VECTOR
	vec_sub_prob[a]->vec_simp_poly = vec_sub_prob[b]->vec_simp_poly;
#endif
	vec_sub_prob[a]->divR = vec_sub_prob[b]->divR;


	vec_sub_prob[b]->poly = tmp1;
	vec_sub_prob[b]->tetraID = tmp2;
	vec_sub_prob[b]->fIdx = tmp3;
	vec_sub_prob[b]->nIdx = tmp4;
	vec_sub_prob[b]->initQ = tmp5;
	vec_sub_prob[b]->tetrQ = tmp6;
	vec_sub_prob[b]->diviQ = tmp7;
#ifdef _SPR_USING_STD_VECTOR
	vec_sub_prob[b]->vec_simp_poly = tmp8;
#endif
	vec_sub_prob[b]->divR = tmp9;
}

void sort_index(int start, int end, std::vector<SubDivideProblem*>& vec_sub_prob)
{
	int i = 0, j = 0;
	int init_n = vec_sub_prob[start]->nIdx;
	for (i = start; i <= end; i++)
	{
		init_n = vec_sub_prob[i]->nIdx;
		for (j = i + 1; j <= end; j++)
		{
			if (init_n > vec_sub_prob[j]->nIdx)
			{
				swap_nprall(j, i, vec_sub_prob);
				init_n = vec_sub_prob[j]->nIdx;
			}
		}

	}
}


//************************************************************
#ifdef _SPR_USING_STD_VECTOR
int SPRImpl::sort_subproblems(std::vector<SubDivideProblem*>& vec_sub_prob)
#else
int SPRImpl::sort_subproblems(SubDivideProblem* subProbArray[], int subProbNum)
#endif
{
	int i;
#ifdef _DEBUG
	static int layerID = 0;
#endif

	/* ----------------------------------------------------------------------------
	 * 最简单的，基于局部贪婪思想，哪些当前形成的四面体单元质量最佳，则其放在首位
	 * --------------------------------------------------------------------------*/
#ifdef _SPR_USING_STD_VECTOR
	eval_subproblems(vec_sub_prob);
	std::sort(vec_sub_prob.begin(), vec_sub_prob.end(), [this](const void *arg1, const void *arg2) {
		return this->comp_subproblemsV2(arg1, arg2);
	});
#else
	eval_subproblems(subProbArray, subProbNum);
	//#BUG? 从qsort改过来可能有问题
	// yhf 2022-02-15 注： 注意，自己写的程序两个相等时必须返回false，否则可能会出现invalid comparator
	std::sort(subProbArray, subProbArray + subProbNum, [this](const void *arg1, const void *arg2) {
		return this->comp_subproblemsV2(arg1, arg2);
	});
#endif

#if 0
	//*********************************zhaodawei add 2010-07-16
	std::vector<int> tmp_vec;
	tmp_vec.clear();
	int start, end;
	double qq = vec_sub_prob[0]->diviQ;
	start = 0;
	for (i = 0; i < vec_sub_prob.size(); i++)
	{
		qq = vec_sub_prob[i]->diviQ;
		tmp_vec.push_back(i);
		while ((i < vec_sub_prob.size() - 1) && (qq == vec_sub_prob[++i]->diviQ));
		if (i == vec_sub_prob.size() - 1)
		{
			if (vec_sub_prob[i]->diviQ == vec_sub_prob[i - 1]->diviQ)
				tmp_vec.push_back(i);
			else
				tmp_vec.push_back(i - 1);
			break;
		}
		tmp_vec.push_back(i - 1);
		i--;
	}
	int count_n = tmp_vec.size() / 2;
	if (count_n > 1)
	{
		for (i = 0; i < count_n; i++)
		{
			start = tmp_vec[2 * i];
			end = tmp_vec[2 * i + 1];
			if (start == end) continue;
			sort_index(start, end, vec_sub_prob);
		}
	}
#endif

	//***********************************************

#ifdef _DEBUG
#if 0
	if (vec_sub_prob.size() > 3)
	{
		spdlog::info("Layer {} Subproblems: size = {}.", ++layerID, vec_sub_prob.size());
		for (i = 0; i < vec_sub_prob.size(); i++)
		{
			if (i % 3 == 0)
			{
				spdlog::info("\n");
			}

			printf("(%.6f,%.6f); ", vec_sub_prob[i]->diviQ, vec_sub_prob[i]->tetrQ);
		}
		spdlog::info("\n");
	}
#endif
#endif

	return 1;
}

inline int decode_quality(double qual, int *nodeNum, double *tetraQ)
{
	assert(qual > 0.0);
	qual = 1.0 / qual;

	*nodeNum = floor(qual);
	*tetraQ = 1.0 - (qual - (double)(*nodeNum));

	assert(*nodeNum >= 0);
	assert(*tetraQ >= 0.0 && *tetraQ < 1.0);

	return 1;
}

/* ---------------------------------------------------------------------------------------
 * 分两步做递归
 * 1. 记录当前子问题的所有可能分解方案（增加1个四面体单元情形下）
 * 2. 启用某种寻优策略来确定优先的分解方案
 * 由于现有横向的比较，而不是盲目地递归，有可能会获得更好的求解效率
 * --------------------------------------------------------------------------------------*/
int SPRImpl::optimal_twoStep(double iniQ, polyhedron *poly, divide *dive, divide_record *divR,
	Constraints *csts, int maxIntNum)
{
	int i, j, r = 0, sr = 0;
#ifdef _SPR_USING_STD_VECTOR
	std::vector<SubDivideProblem*> vec_sub_prob;
#else
	SubDivideProblem *subProbArray[MAX_SUB_PROBLEM_NUM];
	int subProbNum = 0;
#endif
	SubDivideProblem *subProb = NULL;
	double quality_cur = iniQ;
	divide *sub_divide = NULL;
	int ELE[4];
	double time_end, time_elap;
	int intNodeNum1, intNodeNum2;
	double quality1, quality2, quality_sum;
	bool bValidSolution = false;
	int subP_maxIntNum = 0;
	int allowedFaceNum;
	int subProbPolySize = 0;

	/* ----------------------------------------
	* 在开始递归前，确定是否超时
	* --------------------------------------*/
	if (spr_no_more_recur)
		return 0;

	if (spr_time_limit > 0.0)
	{
		time_end = SPRlogTime();
		time_elap = time_end - spr_time_start;
		if (time_elap > spr_time_limit)
		{
			spr_no_more_recur = true;
			spdlog::info("Predefined time has been elapsed.({})\n", time_elap);
			//		exit(1);
			return 0;
		}
	}

	real++;
#ifdef _VERBOSE
	if (real % 1000 == 0)
		spdlog::info("recursive number == {}\n", real);
#endif

	if (poly->num_t == 4)
	{
		/* 不会增加任何交点，确保dive.numInt的值在调用前后不会有改变 */
		return optimal_fourFace(quality_cur, poly, dive, divR, csts, maxIntNum);
	}

	/* 每个子问题所引起的交点个数不能超过所能允许的最大值 */
#ifdef _SPR_USING_STD_VECTOR
	r = record_subproblems(quality_cur, poly, divR, vec_sub_prob, csts, maxIntNum);
#else
	r = record_subproblems(quality_cur, poly, divR, subProbArray, &subProbNum, csts, maxIntNum);
#endif

	if (r == 0)
	{
		goto FAIL;
	}
#ifdef _SPR_USING_STD_VECTOR
	if (vec_sub_prob.size() > 1)
		/* 如果希望尽量少的交点，交点数最少的要优先 */
		sort_subproblems(vec_sub_prob);
#else
	if (subProbNum > 1)
		/* 如果希望尽量少的交点，交点数最少的要优先 */
		sort_subproblems(subProbArray, subProbNum);
#endif /* _SPR_USING_STD_VECTOR */

	sub_divide = (divide*)malloc(sizeof(divide));
	assert(sub_divide);
	r = 0;
	allowedFaceNum = poly->allowedFaceNum;

#ifdef _SPR_USING_STD_VECTOR
	for (i = 0; i < vec_sub_prob.size(); i++)
	{
		subProb = vec_sub_prob[i];
		subProbPolySize = subProb->vec_simp_poly.size();
#else
	for (i = 0; i < subProbNum; i++)
	{
		subProb = subProbArray[i];
		subProbPolySize = subProb->simpPolyNum;
#endif /* _SPR_USING_STD_VECTOR */

		if (spr_optimal_type == EDGE_REMOVAL)
		{
			if (subProb->facNum > allowedFaceNum)
				continue;

			if (subProbPolySize == 1)
#ifdef _SPR_USING_STD_VECTOR
				subProb->vec_simp_poly[0]->allowedFaceNum = allowedFaceNum - subProb->facNum;
#else
				subProb->simpPolyArray[0]->allowedFaceNum = allowedFaceNum - subProb->facNum;
#endif
			else
				spr_allowed_face_num = allowedFaceNum - subProb->facNum;
		}
		else if (spr_optimal_type == MIN_INT_NODE_NUMBER)
		{
			subP_maxIntNum = maxIntNum - subProb->intNum;
			if (subP_maxIntNum < 0 || (subP_maxIntNum == 0 && subProb->tetrQ - quality_cur <= 0.0))
				continue;
		}
		else
		{
			if (subProb->tetrQ - quality_cur <= 0.0)//1.0e-6)
				continue;
		}

		memset(sub_divide, 0, sizeof(divide));
		if (subProbPolySize == 1)
		{
#ifdef _SPR_USING_STD_VECTOR
			sr = optimal_twoStep(quality_cur, subProb->vec_simp_poly[0], sub_divide, &subProb->divR,
				csts, subP_maxIntNum);		//递归调用
#else
			sr = optimal_twoStep(quality_cur, subProb->simpPolyArray[0], sub_divide, &subProb->divR,
				csts, subP_maxIntNum);		//递归调用
#endif /* _SPR_USING_STD_VECTOR */
		}
		else
		{
#ifdef _SPR_USING_STD_VECTOR
			sr = optimal_twoStep(quality_cur, subProb->vec_simp_poly, sub_divide, &subProb->divR,
				csts, subP_maxIntNum);		//递归调用
#else
			sr = optimal_twoStep(quality_cur, subProb->simpPolyArray, subProb->simpPolyNum, sub_divide, &subProb->divR,
				csts, subP_maxIntNum);		//递归调用
#endif /* _SPR_USING_STD_VECTOR */
		}

		if (sr)
		{
			bValidSolution = false;

			if (spr_optimal_type == EDGE_REMOVAL)
			{
				/* 得到一个最优解，需要确定这个最优解是符合我们要求的，并为后续解设定更高的标准 */
				quality_cur = sub_divide->q;
				if (subProb->tetrQ < quality_cur)
				{
					quality_cur = subProb->tetrQ;
				}
				bValidSolution = true;	/* 这是一个有效解 */
			}
			else if (spr_optimal_type != MIN_INT_NODE_NUMBER || !csts)
			{
				/* 更新quality_cur 和dive */
				quality_cur = sub_divide->q;
				if (subProb->tetrQ < quality_cur)
				{
					quality_cur = subProb->tetrQ;
				}
				bValidSolution = true;	/* 这是一个有效解 */
			}
			else
			{/* 此时，要叠加交点个数，并且根据新交点个数来确定子问题的质量 */
#if 0
				decode_quality(sub_divide->q, &intNodeNum1, &quality1);
				decode_quality(subProb->tetrQ, &intNodeNum2, &quality2);

				intNodeNum1 += intNodeNum2;
				if (quality1 > quality2)
					quality1 = quality2;

				assert(intNodeNum1 >= 0);
				assert(quality1 > 0.0 && quality1 < 1.0);

				quality_sum = 1.0 / (intNodeNum1 + (1.0 - quality1));
				if (quality_sum > quality_cur)
				{
					bValidSolution = true; /* 这是一个有效解 */
					quality_cur = quality_sum;
				}
#endif
				quality_cur = sub_divide->q;
				if (subProb->tetrQ < quality_cur)
				{
					quality_cur = subProb->tetrQ;
				}
				bValidSolution = true;	/* 这是一个有效解 */
			}

			if (bValidSolution)
			{
				r = 1;  /* 有一个子问题有解 */
				ELE[0] = poly->t[subProb->fIdx][0];
				ELE[1] = poly->t[subProb->fIdx][1];
				ELE[2] = poly->t[subProb->fIdx][2];
				ELE[3] = subProb->nIdx;

				memcpy(dive, sub_divide, sizeof(struct divide));
				mergeT(dive, ELE);
				dive->q = quality_cur;

				if (spr_optimal_type == MIN_INT_NODE_NUMBER)
				{/* 已经得到一个最优解，它所包含的最大交点数即为最优交点数 */
					maxIntNum = sub_divide->intNum + subProb->intNum;
					dive->intNum = sub_divide->intNum + subProb->intNum;
				}

				if (spr_optimal_type == EDGE_REMOVAL)
				{
					dive->facNum = sub_divide->facNum + subProb->facNum;
					allowedFaceNum = dive->facNum - 1; /* 必须比已经获得的解要小 */
				}

				if (!spr_opt_solution)
					break;
			}
		}
	}
	if (sub_divide)
	{
		free(sub_divide);
		sub_divide = NULL;
	}

	goto END;

FAIL:

END:
#ifdef _SPR_USING_STD_VECTOR
	for (i = 0; i < vec_sub_prob.size(); i++)
	{
		subProb = vec_sub_prob[i];
		for (j = 0; j < subProb->vec_simp_poly.size(); j++)
			free(subProb->vec_simp_poly[j]);
		delete subProb;
	}
#else
	for (i = 0; i < subProbNum; i++)
	{
		subProb = subProbArray[i];
		for (j = 0; j < subProb->simpPolyNum; j++)
			free(subProb->simpPolyArray[j]);
		delete subProb;
	}
#endif /* _SPR_USING_STD_VECTOR  */
	return r;
	}


/* ---------------------------------------------------------------------------------------
 * 深度优先递归
 * --------------------------------------------------------------------------------------*/
int SPRImpl::optimal(double init_q, struct polyhedron *P, struct divide *T, divide_record *divRd)
{
	int r = 0;
	struct polyhedron *Q;
	std::vector<struct polyhedron*> vec_simp_poly;
	struct divide *Tq = NULL;
	TETRA ELE;
	int i, j, j1, j2, j3, j4, fIter = 0, nIter = 0;
	int f;		//记录现在处理的面的位置
	int ipoly;
	double quality_ELE, quality_cur;
	int test;
	int key;
	std::map<int, int>	  mapDivideValid;
	std::map<int, int>::iterator it;
	bool bFaceSucc;
	bool freeSimpPoly;
	std::set<int> setMutexNodes;

	//#ifdef _TIMING_PERFORMANCE
	real++;
	//#endif
	divide_record cur_record;

#ifndef _DEBUG
	if (real % 5000 == 0)
		spdlog::info("recursive number == {}\n", real);
#endif

	/* 剩下最后四个单元时 */
#ifdef _DEBUG
//	savepls("sphere-sample-output.pls", P, true);
//	savepls("sphere-sample-output-final.pls", P, false);
#endif

	if (P->num_vable <= 4)
	{
		assert(P->num_vable == 4);
		for (i = 0; i < 3; i++)
			ELE[i] = P->t[0][i];

		for (i = 0; i < 4; i++)
		{
			if (P->t[1][i] != ELE[0] &&
				P->t[1][i] != ELE[1] &&
				P->t[1][i] != ELE[2])
			{
				ELE[3] = P->t[1][i];
				break;
			}
		}
		assert(i < 4);

		quality_cur = isVisible_Search(ELE[0], ELE[1], ELE[2], ELE[3]);

		//如果体积是小于零的话，就根本没有再考虑的必要，必然是不合理的四面体
		if (quality_cur > 0.0 && quality_cur - init_q > 1.0e-6)
		{
			r = 1;	/* 质量符合 */

			memset(T, 0, sizeof(struct divide));
			mergeT(T, ELE);
			T->q = quality_cur;
		}
	}
	else
	{
		quality_cur = init_q;
		Q = (struct polyhedron*)malloc(sizeof(struct polyhedron));
		Tq = (struct divide*)malloc(sizeof(struct divide));

#ifndef _SORT_FACET_AND_NODE
		sort_facets(P);
		for (fIter = 0; fIter < P->num_t; fIter++)
		{
			f = get_sorted_facet(P, fIter);

#else

		for (f = P->num_t - 1; f >= 0; f--)
			//	for (f = 0; f < P->num_t; f++)
		{
#endif
			// 			setMutexNodes.clear();
			// 			get_mutex_facet_node(P, f, setMutexNodes);

			for (i = 0; i < 3; i++)
				ELE[i] = P->t[f][i];

			bFaceSucc = false;

#ifdef _SORT_FACET_AND_NODE
			sort_nodes(P, f);
			for (nIter = 0; nIter < num_vertices; nIter++)
			{
				j = get_sorted_node(P, nIter);
#else
			for (j = 0; j < num_vertices; j++)
			{
#endif
				if (P->used_v[j] == 0)
					continue;
				if (ELE[0] == j || ELE[1] == j || ELE[2] == j)
					continue;
				ELE[3] = j;

				// 				if (setMutexNodes.find(j) != setMutexNodes.end())
				// 				{
				// 					continue;
				// 				}

#ifdef _TIMING_PERFORMANCE
				count++;		//记录函数运行次数
#endif

				//测试用代码
#ifdef _TIMING_PERFORMANCE
				starttmp[0] = SPRlogTime();
#endif
				j1 = ELE[0];
				j2 = ELE[1];
				j3 = ELE[2];
				j4 = ELE[3];

				sort_quad_pair(&j1, &j2, &j3, &j4);
				key = make_quad_pair(j1, j2, j3, j4);
				it = mapDivideValid.find(key);
				if (it == mapDivideValid.end())
				{

#ifndef _DIVIDE_RECORD
					if (divRd->num < 0)
					{
						/* 确定小空腔是否已经四面体化 */
// 						if (!is_nonoverlap_tet(divRd, key))
// 							continue;

						int found = 0;
						create_div_record(divRd, key, &cur_record, &found);

						if (found == 0)
						{
							if (is_tested_div_record(&cur_record))
								continue;
							else
								add_div_record(&cur_record);
						}
						else
							continue;
					}
					else
						cur_record.num = divRd->num;
#endif

#ifdef _DIVIDE_POLY
					freeSimpPoly = false;
					test = istetra(ELE, P, quality_cur, &quality_ELE, Q, vec_simp_poly, f, &freeSimpPoly);
#else
					test = istetra(ELE, P, quality_cur, &quality_ELE, Q, f);
					vec_simp_poly.push_back(Q);
#endif
					mapDivideValid[key] = test;
				}
				else
					continue;

#ifdef _TIMING_PERFORMANCE
				stoptmp[0] = SPRlogTime();
				dura[0] += (stoptmp[0] - starttmp[0]);
#endif

				if (test)	//判断这个四面体是不是一个有效地四面体
				{
					bFaceSucc = true;

					if (vec_simp_poly.size() == 1)
					{
#ifdef _DIVIDE_POLY
						assert(assert_poly_neig(Q));
						assert(assert_simp_poly_edge(Q));
#endif
						if (optimal(quality_cur, Q, Tq, &cur_record))		//递归调用
						{
							r = 1;

							/* 更新quality_cur 和 T */
							quality_cur = Tq->q;
							if (quality_ELE < quality_cur)
							{
								quality_cur = quality_ELE;
							}

							/* 如果成功了，记录T */
							memcpy(T, Tq, sizeof(struct divide));
							//	mergeT(T, Tq);
							mergeT(T, ELE);
							T->q = quality_cur;
						}
					}
					else
					{
						if (optimal(quality_cur, vec_simp_poly, Tq, &cur_record))		//递归调用
						{
							r = 1;

							/* 更新quality_cur 和 T */
							quality_cur = Tq->q;
							if (quality_ELE < quality_cur)
							{
								quality_cur = quality_ELE;
							}

							/* 如果成功了，记录T */
							memcpy(T, Tq, sizeof(struct divide));
							//	mergeT(T, Tq);
							mergeT(T, ELE);
							T->q = quality_cur;
						}
					}
				}

				if (freeSimpPoly)
				{
					for (ipoly = 0; ipoly < vec_simp_poly.size(); ipoly++)
						free(vec_simp_poly[ipoly]);
				}
				vec_simp_poly.clear();
			}

			/* 如果某个面没有找到一个合适的四面体递归下去，终止整个过程
			 * 注意，有可能前面已经成功过一次 所以返回值可以是1，也可以是0 */
			if (!bFaceSucc)
			{
				break;
			}
			}
		free(Q);//及时释放内存
		free(Tq);
		}

#ifdef _DEBUG
	static bool first = true;
	if (r == 0 && P->num_vable > 0 && first)
	{
		savepls("fail.pls", P);
		first = false;
	}
#endif 

	return r;
		}

int SPRImpl::optimal(double init_q, std::vector<struct polyhedron*> &vec_simp_poly, struct divide *T, divide_record *divRd)
{
	int i, res;
	struct divide simp_div;

	memset(T, 0, sizeof(divide));
	T->q = 1.0e10;
	for (i = 0; i < vec_simp_poly.size(); i++)
	{
		assert(assert_poly_neig(vec_simp_poly[i]));
		assert(assert_simp_poly_edge(vec_simp_poly[i]));

		res = optimal(init_q, vec_simp_poly[i], &simp_div, divRd);

		if (res == 0)
			break;

		mergeT(T, &simp_div);
		if (T->q > simp_div.q)
			T->q = simp_div.q;
	}

	return res;
}

int SPRImpl::istetra(int ELE[], struct polyhedron *P, double q, double *quality_ELE,
	struct polyhedron *Q, int f, Constraints *csts)
{
	int newCnt = 0;

#ifdef _TIMING_PERFORMANCE
	starttmp[5] = SPRlogTime();
	Call[0]++;
#endif
	*quality_ELE = isVisible_Search(ELE[0], ELE[1], ELE[2], ELE[3]);

#ifdef _TIMING_PERFORMANCE
	stoptmp[5] = SPRlogTime();
	dura[5] += (stoptmp[5] - starttmp[5]);
#endif

	if (*quality_ELE <= 0.0 || *quality_ELE - q <= 0.00001)
		return 0;	/* 质量不符合 */

#ifdef _TIMING_PERFORMANCE
	starttmp[1] = SPRlogTime();
#endif

	memset(Q, 0, sizeof(struct polyhedron));
	if (!delELE(ELE, P, Q, f, &newCnt))		//将ELE从P中分离出来，得到Q
		return 0;

#ifdef _TIMING_PERFORMANCE
	stoptmp[1] = SPRlogTime();
	dura[1] += (stoptmp[1] - starttmp[1]);
#endif

	return isNewFacetsValid_BruteForce(Q, newCnt);
}

int SPRImpl::isNewFacetsValid_BFS(polyhedron *poly, int newCnt)
{
	int i, j, m, iElem, iSrch;
	std::queue<int> queTreeSear;
	int test_flag[MAX_SPR_POLY_SIZE];

	int *face1 = NULL, *face2 = NULL;
	APOINT facep1[3], facep2[3];
	double test;

	//程序此时需要判断可能增加的三个面是否都是可行的
	for (i = poly->num_t - newCnt; i < poly->num_t; i++)
	{
		memset(test_flag, 0, sizeof(int)*MAX_SPR_POLY_SIZE);

		face1 = poly->t[i];		//face1是遍历现在已有的所有的面
		for (j = 0; j < 3; j++)
		{
			for (m = 0; m < 3; m++)
				facep1[j][m] = vertices[face1[j]][m];
		}
		test_flag[i] = 1;

		for (m = 0; m <= 2; m++)
		{
			iSrch = poly->n[i][m];
			queTreeSear.push(iSrch);
			test_flag[iSrch] = 1;
		}

		do
		{
			/* visit the elements first */
			iElem = queTreeSear.front();
			queTreeSear.pop();
			face2 = poly->t[iElem];		//face1是遍历现在已有的所有的面
			for (j = 0; j < 3; j++)
			{
				for (m = 0; m < 3; m++)
					facep2[j][m] = vertices[face2[j]][m];
			}

			//以下为测试用代码
#ifdef _TIMING_PERFORMANCE
			starttmp[3] = SPRlogTime();
#endif

#ifdef _OPT_INTERSEC_ALGO
			test = isintersect(face1, face2, facep1, facep2, poly);		//将face1中的每一个面和face2比较，看看是否相交
#else
			test = isintersect_old(face1, face2, facep1, facep2, poly);
#endif 

#ifdef _TIMING_PERFORMANCE
			stoptmp[3] = SPRlogTime();
			dura[3] += (stoptmp[3] - starttmp[3]);
#endif

			if (test)
			{
				return 0;
			}

			for (m = 0; m <= 2; m++)
			{
				iSrch = poly->n[iElem][m];
				if (iSrch >= 0 && test_flag[iSrch] == 0)
				{
					queTreeSear.push(iSrch);
					test_flag[iSrch] = 1;
				}
			}
		} while (!queTreeSear.empty());
	}

	return 1;
}

int SPRImpl::isNewFacetsValid_BruteForce(polyhedron *poly, int newCnt, int newIdx[],
	Constraints *csts)
{
	int i, j, k, m;
	int *face1 = NULL, *face2 = NULL;
	APOINT facep1[3], facep2[3];
	double test;
	int test_flag[MAX_SPR_POLY_SIZE];

	memset(test_flag, 0, sizeof(int)*MAX_SPR_POLY_SIZE);
	for (i = 0; i < newCnt; i++)
		test_flag[newIdx[i]] = 1;

	//程序此时需要判断可能增加的三个面是否都是可行的
	for (i = 0; i < newCnt; i++)
	{
		face1 = poly->t[i];		//face1是遍历现在已有的所有的面
		for (j = 0; j < 3; j++)
		{
			for (m = 0; m < 3; m++)
				facep1[j][m] = vertices[face1[j]][m];
		}

		/* 判断新增加的面是否和这些直线相交 */

		for (k = 0; k < poly->num_t; k++)
		{
			if (test_flag[k] == 1)
				continue;

			face2 = poly->t[k];		//face1是遍历现在已有的所有的面
			for (j = 0; j < 3; j++)
			{
				for (m = 0; m < 3; m++)
					facep2[j][m] = vertices[face2[j]][m];
			}

			//以下为测试用代码
#ifdef _TIMING_PERFORMANCE
			starttmp[3] = SPRlogTime();
#endif

#ifdef _OPT_INTERSEC_ALGO
			test = isintersect(face1, face2, facep1, facep2, poly);		//将face1中的每一个面和face2比较，看看是否相交
#else
			test = isintersect_old(face1, face2, facep1, facep2, poly);
#endif 

#ifdef _TIMING_PERFORMANCE
			stoptmp[3] = SPRlogTime();
			dura[3] += (stoptmp[3] - starttmp[3]);
#endif

			if (test)
			{
				return 0;
			}
		}
	}

	return 1;
}

int SPRImpl::calcInstPointNumber_LineCst(polyhedron *poly, int newCnt,
	Constraints *csts, int maxIntNum, int *intNum)
{
	typedef struct NewEdge
	{
		int ipnt[2];/* 边的2个端点 */
		int ifac[2];/* 包含边的个面 */
	} NewEdge;
	int i, j, k, m, i1, i2, if1, if2;
	int *face1 = NULL, *face2 = NULL;
	APOINT facep1[3], facep2[3], linep[2];
	int iface1[3], iline[2];
	int test;
	int lineSameCnt = 0;
	NewEdge newEdges[3]; /* 最多只有3条新边 */
	int newENum = 0;
	double proj1[2], proj2[2], proj3[2], proj4[2];
	bool isInt = 0;
	double ort1, ort2, ort3, ort4;
	int iNewFaces[3] = { 0, 0, 0 }; /* 3个新面是否需要判断相交 */

	assert(spr_optimal_type == MIN_INT_NODE_NUMBER);
	assert(csts->lineCst_size == 1);
	assert(csts->faceCst_size == 0);

	iline[0] = csts->lineCsts[0][0];
	iline[1] = csts->lineCsts[0][1];
	for (j = 0; j < 2; j++)
	{
		for (m = 0; m < 3; m++)
			linep[j][m] = vertices[iline[j]][m];
	}

	//spr_int_node_num = 0; /* 记录和约束相交的面片的个数 */
	*intNum = 0;
	/* 计算线和单元新面的交点时，有以下一些准则
	 * 1. 忽略和线共面的新面
	 * 2. 如果和面的交点在边上（不可能在点上），且边为新边，则交点数加1
	 * 3. 有包含某条和线约束相交的边面片只需与约束线求交1次
	 */

	 /* 构建3条新边 */
	for (i = poly->num_t - newCnt, newENum = 0; i < poly->num_t - 1; i++)
	{
		face1 = poly->t[i];		//face1是遍历现在已有的所有的面
		for (j = i + 1; j < poly->num_t; j++)
		{
			face2 = poly->t[j];		//face1是遍历现在已有的所有的面

			for (m = 0; m < 3; m++)
			{
				i1 = face1[(m + 1) % 3];
				i2 = face1[(m + 2) % 3];

				if ((i1 == face2[0] || i1 == face2[1] || i1 == face2[2]) &&
					(i2 == face2[0] || i2 == face2[1] || i2 == face2[2]))
				{
					newEdges[newENum].ipnt[0] = i1;
					newEdges[newENum].ipnt[1] = i2;
					newEdges[newENum].ifac[0] = i + newCnt - poly->num_t;
					newEdges[newENum].ifac[1] = j + newCnt - poly->num_t;
					newENum++;
					assert(newENum <= 3);
				}
			}
		}
	}

	/* 看看新边是否和约束线相交 */
	for (i = 0; i < newENum; i++)
	{
		i1 = newEdges[i].ipnt[0];
		i2 = newEdges[i].ipnt[1];
		if1 = newEdges[i].ifac[0];
		if2 = newEdges[i].ifac[1];

		if (i1 == iline[0] || i2 == iline[0] ||
			i1 == iline[1] || i2 == iline[1])
			continue; /* 有一个点相同时，则不相交 */

		/* 不共面时，必不想交 */
		if (GEOM_FUNC::orient3d(vertices[i1], vertices[i2], vertices[iline[0]], vertices[iline[1]]) != 0.0)
			continue;

		/* 共面了，先转化为二维问题 */
		GEOM_FUNC::proj_four_coplanr_points(vertices[iline[0]], vertices[iline[1]], vertices[i1], vertices[i2],
			proj1, proj2, proj3, proj4);

		isInt = false;
		/* proj3和proj4 不可能是线段上的点,所以省略了相应的共面判断 */
		ort1 = GEOM_FUNC::orient2d(proj3, proj4, proj1);
		if (ort1 == 0.0)
		{
			isInt = (proj1[0] > proj3[0] && proj1[0] < proj4[0]) ||
				(proj1[0] < proj3[0] && proj1[0] > proj4[0]);
		}
		else
		{
			ort2 = GEOM_FUNC::orient2d(proj3, proj4, proj2);
			if (ort2 == 0.0)
				isInt = (proj2[0] > proj3[0] && proj2[0] < proj4[0]) ||
				(proj2[0] < proj3[0] && proj2[0] > proj4[0]);
			else if (ort1*ort2 < 0.0)
			{
				ort3 = GEOM_FUNC::orient2d(proj1, proj2, proj3);
				ort4 = GEOM_FUNC::orient2d(proj1, proj2, proj4);
				isInt = ort3 * ort4 < 0.0;
			}
		}

		if (isInt)
		{
			//spr_int_node_num++;
			(*intNum)++;

			/* 以下新面不需要再判断了 */
			iNewFaces[if1] = 1;
			iNewFaces[if2] = 1;

			if (*intNum > maxIntNum)
			{
				return 0;	/* 交点已超过允许的最大值，则提前结束 */
			}
		}
	}

	/* 最多只可能有2个交点 */
	if (*intNum > 2)
		return 0;

	/* 确定新面和约束线的关系 */
	for (i = poly->num_t - newCnt; i < poly->num_t; i++)
	{
		if1 = i + newCnt - poly->num_t;
		if (iNewFaces[if1] == 1)
			continue;

		face1 = poly->t[i];		//face1是遍历现在已有的所有的面
		for (j = 0; j < 3; j++)
		{
			for (m = 0; m < 3; m++)
				facep1[j][m] = vertices[face1[j]][m];
		}

		/* 忽略和约束线共面的面 */
		if (GEOM_FUNC::orient3d(facep1[0], facep1[1], facep1[2], linep[0]) == 0.0 &&
			GEOM_FUNC::orient3d(facep1[0], facep1[1], facep1[2], linep[1]) == 0.0)
			continue;

		/* 如果面的边和线相交，也应该忽略 */
		if (GEOM_FUNC::orient3d(facep1[0], facep1[1], linep[0], linep[1]) == 0.0 ||
			GEOM_FUNC::orient3d(facep1[1], facep1[2], linep[0], linep[1]) == 0.0 ||
			GEOM_FUNC::orient3d(facep1[2], facep1[0], linep[0], linep[1]) == 0.0)
			continue;

		/* 允许和约束的相交，且统计新面和约束相交的个数 */
		test = isintersect(face1, facep1, csts);
		if (test)
		{
			//spr_int_node_num++;
			(*intNum)++;
			if (*intNum > maxIntNum)
			{
				return 0;	/* 交点已超过允许的最大值，则提前结束 */
			}
		}
	}

	return 1;
}

int SPRImpl::isNewFacetsValid_BruteForce(polyhedron *poly, int newCnt,
	Constraints *csts, int maxIntNum, int *intNum, InsertedTetra *instTetra, polyhedron *prtpoly)
{
	int i, j, k, m;
	int *face1 = NULL, *face2 = NULL;
	APOINT facep1[3], facep2[3], linep[2];
	double *facep1_fast[3], *facep2_fast[3];
	int iface1[3], iline[2];
	int test;
	int lineSameCnt = 0;
	int nRet = 0;

	if (spr_optimal_type == MIN_INT_NODE_NUMBER)
	{
		nRet = calcInstPointNumber_LineCst(poly, newCnt, csts, maxIntNum, intNum);
		if (nRet == 0)
			return 0;
	}
	else
	{
		/* 确定新面包含待删除边的个数小于允许的最大值 */
		if (spr_optimal_type == EDGE_REMOVAL)
		{
			assert(edgeConsArraySize > 0);
			spr_subprob_face_num = 0;	/* 记录子问题实际包含的面片数，临时变量 */
			for (i = poly->num_t - newCnt; i < poly->num_t; i++)
			{
				face1 = poly->t[i];		//face1是遍历现在已有的所有的面
#if 1
				/* 新面不应包含任何type = 0的dependent edge */
				for (j = 0; j < spr_dependent_edges_cnt; j++)
				{
					if (true && /*spr_dependent_edges[j][2] >= 0 && /* 没有顶点为壳边端点 */
						(face1[0] == spr_dependent_edges[j][0] || face1[1] == spr_dependent_edges[j][0] || face1[2] == spr_dependent_edges[j][0]) &&
						(face1[0] == spr_dependent_edges[j][1] || face1[1] == spr_dependent_edges[j][1] || face1[2] == spr_dependent_edges[j][1]))
						return 0;
				}
#endif
				if ((face1[0] == spr_shell_edge_nodeS || face1[1] == spr_shell_edge_nodeS || face1[2] == spr_shell_edge_nodeS) &&
					(face1[0] == spr_shell_edge_nodeE || face1[1] == spr_shell_edge_nodeE || face1[2] == spr_shell_edge_nodeE))
				{/* 面包含壳边 */
					spr_subprob_face_num++;
					if (spr_subprob_face_num > prtpoly->allowedFaceNum)
						return 0;
				}
#if 0
				else
				{/* 面不包含壳边 */
					/* 不包括shellEdge的新面不应包含type = 1的dependent edge */
					for (j = 0; j < spr_dependent_edges_cnt; j++)
					{
						if (spr_dependent_edges[j][2] == 1 && /* 没有顶点为壳边端点 */
							(face1[0] == spr_dependent_edges[j][0] || face1[1] == spr_dependent_edges[j][0] || face1[2] == spr_dependent_edges[j][0]) &&
							(face1[0] == spr_dependent_edges[j][1] || face1[1] == spr_dependent_edges[j][1] || face1[2] == spr_dependent_edges[j][1]))
							return 0;
					}
				}
#endif
			}
		}

		if (csts || edgeConsArraySize > 0)
		{
			/* 确定新面和约束线的关系 */
			for (i = poly->num_t - newCnt; i < poly->num_t; i++)
			{
				face1 = poly->t[i];		//face1是遍历现在已有的所有的面
				for (j = 0; j < 3; j++)
				{
					for (m = 0; m < 3; m++)
						facep1[j][m] = vertices[face1[j]][m];
				}

				/* 如果约束不为空，需要确定新面和约束不相交 */
				if (csts)
				{
					test = isintersect(face1, facep1, csts);
					if (test)
						return 0;
				}

				if (edgeConsArraySize > 0)
				{
					test = isIntersectEdgeConstraints(face1, facep1);
					if (test)
						return 0;
				}
			}
		}
	}

	//程序此时需要判断可能增加的三个面是否都是可行的
	for (i = poly->num_t - newCnt; i < poly->num_t; i++)
	{
		face1 = poly->t[i];		//face1是遍历现在已有的所有的面
#if 0
		for (j = 0; j < 3; j++)
		{
			for (m = 0; m < 3; m++)
				facep1[j][m] = vertices[face1[j]][m];
		}
#else
		facep1_fast[0] = vertices[face1[0]];
		facep1_fast[1] = vertices[face1[1]];
		facep1_fast[2] = vertices[face1[2]];
#endif

		for (k = 0; k < poly->num_t - newCnt; k++)
		{
			face2 = poly->t[k];		//face1是遍历现在已有的所有的面
#if 0
			for (j = 0; j < 3; j++)
			{
				for (m = 0; m < 3; m++)
					facep2[j][m] = vertices[face2[j]][m];
			}
#else
			facep2_fast[0] = vertices[face2[0]];
			facep2_fast[1] = vertices[face2[1]];
			facep2_fast[2] = vertices[face2[2]];
#endif

			//以下为测试用代码
#ifdef _TIMING_PERFORMANCE
			starttmp[3] = SPRlogTime();
#endif

#ifdef _OPT_INTERSEC_ALGO
			//			test = isintersect(face1, face2, facep1, facep2, poly);		//将face1中的每一个面和face2比较，看看是否相交
			//			test = GEOM_FUNC::tri_tri_intersect3d(face1, face2, facep1, facep2);
			test = GEOM_FUNC::tri_tri_intersect3d_fast(face1, face2, facep1_fast, facep2_fast);
#else
			test = isintersect_old(face1, face2, facep1, facep2, poly);
#endif 

#ifdef _TIMING_PERFORMANCE
			stoptmp[3] = SPRlogTime();
			dura[3] += (stoptmp[3] - starttmp[3]);
#endif

			if (test)
			{
				return 0;
			}
		}
	}

#if 0
	if (spr_optimal_type == EDGE_REMOVAL && spr_dependent_edges_cnt > 0)
	{
		for (i = poly->num_t - newCnt; i < poly->num_t; i++)
		{
			face1 = poly->t[i];		//face1是遍历现在已有的所有的面

			/* 新面不应包含dependent edge */
			for (j = 0; j < spr_dependent_edges_cnt; j++)
			{
				if ((face1[0] == spr_dependent_edges[j][0] || face1[1] == spr_dependent_edges[j][0] || face1[2] == spr_dependent_edges[j][0]) &&
					(face1[0] == spr_dependent_edges[j][1] || face1[1] == spr_dependent_edges[j][1] || face1[2] == spr_dependent_edges[j][1]))
				{
					return 0;
				}
			}
		}
	}
#endif


	return 1;
}

#ifdef _SPR_USING_STD_VECTOR
int SPRImpl::istetra(int ELE[], struct polyhedron *P, double q, double *quality_ELE,
	struct polyhedron *Q,
	std::vector<struct polyhedron *>& vec_simp_poly, int f, bool *freeSimpPoly,
	Constraints *csts, int maxIntNum, int *intNum)
#else
int SPRImpl::istetra(int ELE[], struct polyhedron *P, double q, double *quality_ELE,
	struct polyhedron *Q,
	polyhedron *simpPolyArray[], int *simpPolyNum, int f, bool *freeSimpPoly,
	Constraints *csts, int maxIntNum, int *intNum)
#endif
{
	double tetraQ = 0.0, volume;
	int nRet = 0;

#ifdef _TIMING_PERFORMANCE
	starttmp[5] = SPRlogTime();
	Call[0]++;
#endif
	tetraQ = isVisible_Search(ELE[0], ELE[1], ELE[2], ELE[3], &volume);

#ifdef _TIMING_PERFORMANCE
	stoptmp[5] = SPRlogTime();
	dura[5] += (stoptmp[5] - starttmp[5]);
#endif
	if (spr_optimal_type == EDGE_REMOVAL)
	{
		if (volume <= SPR_UNALLOWED_QUALITY)
			return 0;	/* 质量不符合 */

#if 1
		if (tetraQ < SPR_ZERO_QUALITY)
			tetraQ = SPR_ZERO_QUALITY;
#endif

		/* -----------------------------------------------------------
		 * 确认相邻3个面和待插入点的关系
		 * ---------------------------------------------------------*/
		if (!isNeigborAllowed(f, ELE[3], P))
			return 0;

#ifdef _SPR_USING_STD_VECTOR
		nRet = delELE(ELE, P, Q, vec_simp_poly, f, freeSimpPoly, csts, maxIntNum, intNum);//将ELE从P中分离出来，得到Q
#else
		nRet = delELE(ELE, P, Q, simpPolyArray, simpPolyNum, f, freeSimpPoly, csts, maxIntNum, intNum);//将ELE从P中分离出来，得到Q
#endif

		if (nRet)
		{/* 成功了，继续检查单元质量 */

			/* 计算质量，确定它是否合适 */
			assert(tetraQ > 0.0 && tetraQ < 1.0);
			*quality_ELE = tetraQ;
		}

		return nRet;
	}
	else if (spr_optimal_type != MIN_INT_NODE_NUMBER || !csts)
	{
		*quality_ELE = tetraQ;
		if (volume <= SPR_UNALLOWED_QUALITY || *quality_ELE - q <= 0.0)
			return 0;	/* 质量不符合 */

		/* -----------------------------------------------------------
		 * 确认相邻3个面和待插入点的关系
		 * ---------------------------------------------------------*/
		if (!isNeigborAllowed(f, ELE[3], P))
			return 0;

#ifdef _SPR_USING_STD_VECTOR
		return delELE(ELE, P, Q, vec_simp_poly, f, freeSimpPoly, csts, maxIntNum, intNum);//将ELE从P中分离出来，得到Q
#else
		return delELE(ELE, P, Q, simpPolyArray, simpPolyNum, f, freeSimpPoly, csts, maxIntNum, intNum);//将ELE从P中分离出来，得到Q
#endif
	}
	else
	{
		if (volume <= SPR_UNALLOWED_QUALITY)
			return 0;	/* 质量不符合 */

#if 1
		if (tetraQ < SPR_ZERO_QUALITY)
			tetraQ = SPR_ZERO_QUALITY;
#endif

		/* -----------------------------------------------------------
		 * 确认相邻3个面和待插入点的关系
		 * ---------------------------------------------------------*/
		if (!isNeigborAllowed(f, ELE[3], P))
			return 0;

#if 1
		if (maxIntNum == 0 && tetraQ - q <= 0.0)
			return 0;
#endif

#ifdef _SPR_USING_STD_VECTOR
		nRet = delELE(ELE, P, Q, vec_simp_poly, f, freeSimpPoly, csts, maxIntNum, intNum);//将ELE从P中分离出来，得到Q
#else
		nRet = delELE(ELE, P, Q, simpPolyArray, simpPolyNum, f, freeSimpPoly, csts, maxIntNum, intNum);//将ELE从P中分离出来，得到Q
#endif

		if (nRet)
		{/* 成功了，继续检查单元质量 */

			/* 计算质量，确定它是否合适 */
			assert(tetraQ > 0.0 && tetraQ < 1.0);
			//assert(spr_int_node_num >= 0);
#if 0
			* quality_ELE = 1.0 / ((1.0 - tetraQ) + spr_int_node_num);
#endif
			*quality_ELE = tetraQ;

			assert(*intNum <= maxIntNum);
			if (*intNum >= maxIntNum && *quality_ELE - q <= 0.0)
				return 0;
		}

		return nRet;
	}
}

//几个原则：四面体中新增加的线可以和别的三角形共享
int SPRImpl::isintersect_old(TRI face1, TRI face2, APOINT facep1[], APOINT facep2[], struct polyhedron *P)
{
	//	double arear1, arear2;
	int res[6];
	int i;

	face_n++;

	// 	if (face1[0] == 3 && face1[1] == 5 && face1[2] == 7 &&
	// 		face2[0] == 0 && face2[1] == 4 && face2[2] == 5)
		/*	face_n = face_n;*/

	for (i = 0; i < 6; i++)
		res[i] = 0;

	//分别考虑每一个边与面的相交问题，就有六个相交判断了，其中要先做整数判断，后做浮点判断

	if (isintercross(face1, face2, facep1, facep2, 0, 1, &res[0], P))
		return 1;
	else if (isintercross(face1, face2, facep1, facep2, 1, 2, &res[1], P))
		return 1;
	else if (isintercross(face1, face2, facep1, facep2, 0, 2, &res[2], P))
		return 1;
	else if (isintercross(face2, face1, facep2, facep1, 0, 1, &res[3], P))
		return 1;
	else if (isintercross(face2, face1, facep2, facep1, 1, 2, &res[4], P))
		return 1;
	else if (isintercross(face2, face1, facep2, facep1, 0, 2, &res[5], P))
		return 1;

	if (!(isboxinter(facep1, facep2)))		//如果在box的层面上发现这个是不相交的话，那就肯定不会相交了
	{
		box_n++;
		return 0;
	}
	/*
		arear2=areas(facep2[0][0], facep2[1][0], facep2[2][0],
					 facep2[0][1], facep2[1][1], facep2[2][1],
					 facep2[0][2], facep2[1][2], facep2[2][2]);
	*/
	float_count++;
	if (islinecrossface(face1, face2, facep1, facep2, 0, 1, res[0], P))	//再做六个线面的浮点计算判断
		return 1;
	if (islinecrossface(face1, face2, facep1, facep2, 1, 2, res[1], P))	//这里有一个res来作为信息的传递
		return 1;
	if (islinecrossface(face1, face2, facep1, facep2, 0, 2, res[2], P))
		return 1;
	/*
		arear1=areas(facep1[0][0], facep1[1][0], facep1[2][0],
					facep1[0][1], facep1[1][1], facep1[2][1],
					facep1[0][2], facep1[1][2], facep1[2][2]);
	*/
	if (islinecrossface(face2, face1, facep2, facep1, 0, 1, res[3], P))
		return 1;
	if (islinecrossface(face2, face1, facep2, facep1, 1, 2, res[4], P))
		return 1;
	if (islinecrossface(face2, face1, facep2, facep1, 0, 2, res[5], P))
		return 1;
	return 0;
}


// 几个原则：四面体中新增加的线可以和别的三角形共享
// inline int isintersect(TRI face1, TRI face2, APOINT facep1[], APOINT facep2[], struct polyhedron *P)
// {
// 	int sameVtxCnt = 0, i;
// 	int nRet = 0;
// 
// 	face_n++;
// 
// 	for (i = 0; i < 3; i++)
// 	{
// 		if (face1[i] == face2[0] ||
// 			face1[i] == face2[1] ||
// 			face1[i] == face2[2])
// 			sameVtxCnt++;
// 	}
// 
// 	if (sameVtxCnt == 0)
// 	{
// 		nRet = tri_tri_overlap_test_3d(facep1[0], facep1[1], facep1[2],
// 			facep2[0], facep2[1], facep2[2]);
// 	}
// 	else if (sameVtxCnt == 1)
// 	{
// #ifdef _TIMING_PERFORMANCE
// 		starttmp[8] = SPRlogTime();
// #endif
// 		nRet = isintersect_old(face1, face2, facep1, facep2, P);
// 
// #ifdef _TIMING_PERFORMANCE
// 		stoptmp[8] = SPRlogTime();
// 		dura[8] += (stoptmp[8] - starttmp[8]);
// #endif
// 	}
// 
// 	return nRet;
// }

//*************************************zhaodawei add 2010-07-16

void swap_apoint(APOINT facep1, APOINT facep2)
{
	APOINT tmp;
	tmp[0] = facep1[0];
	tmp[1] = facep1[1];
	tmp[2] = facep1[2];

	facep1[0] = facep2[0];
	facep1[1] = facep2[1];
	facep1[2] = facep2[2];

	facep2[0] = tmp[0];
	facep2[1] = tmp[1];
	facep2[2] = tmp[2];
}

/* iline[0]为共享点 */
inline int isintersect_oneSharePoint(int iface[], int iline[], APOINT facep[], APOINT linep[])
{
	int isInt = 0;
	double s1, s2;
	double facepj[3][2], linepj[3][2] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };	/* 共面情形的投影点 */
	double normalf[3];
	double nx, ny, nz;
	int i, m, i1, i2;
	double vect0[2], vect1[2], vect2[2], scale1, scale2;

	/* 不共面，必不相交 */
	s1 = GEOM_FUNC::orient3d(facep[0], facep[1], facep[2], linep[1]);
	if (s1 == 0.0)
	{/* 共面? */
		for (m = 0; m < 3; m++)
		{
			if (iface[m] == iline[0])
				break;
		}
		assert(m < 3);

		/* 求得平面的法向，沿使得平面投影面积最大的方向投影平面和直线 */
		norm_3p(facep[0], facep[1], facep[2], normalf);
		nx = ((normalf[0] < 0.0) ? -normalf[0] : normalf[0]);
		ny = ((normalf[1] < 0.0) ? -normalf[1] : normalf[1]);
		nz = ((normalf[2] < 0.0) ? -normalf[2] : normalf[2]);

		if (nx > 0.0 || ny > 0.0 || nz > 0.0)
		{
			if ((nx > nz) && (nx >= ny))
			{ //往YOZ平面投影
				i1 = 1;
				i2 = 2;
			}
			else if ((ny > nz) && (ny >= nx))
			{//往ZOX平面投影
				i1 = 2;
				i2 = 0;
			}
			else
			{//往XOY平面投影
				i1 = 0;
				i2 = 1;
			}

			facepj[0][0] = facep[0][i1];
			facepj[0][1] = facep[0][i2];
			facepj[1][0] = facep[1][i1];
			facepj[1][1] = facep[1][i2];
			facepj[2][0] = facep[2][i1];
			facepj[2][1] = facep[2][i2];
			linepj[0][0] = linep[0][i1];
			linepj[0][1] = linep[0][i2];
			linepj[1][0] = linep[1][i1];
			linepj[1][1] = linep[1][i2];

			/* 判断linepj[1]是否为三角形顶点facepj[m]对应的内角所包含 */
			if (GEOM_FUNC::orient2d(facepj[0], facepj[1], facepj[2]) >= 0.0)
			{
				isInt = GEOM_FUNC::orient2d(facepj[m], facepj[(m + 1) % 3], linepj[1]) >= 0.0 &&
					GEOM_FUNC::orient2d(facepj[(m + 2) % 3], facepj[m], linepj[1]) >= 0.0;
			}
			else
			{
				isInt = GEOM_FUNC::orient2d(facepj[m], facepj[(m + 2) % 3], linepj[1]) >= 0.0 &&
					GEOM_FUNC::orient2d(facepj[(m + 1) % 3], facepj[m], linepj[1]) >= 0.0;
			}
		}
		else
		{/* 平面退化成一条直线 */
			/* 确定直线和2条以facepj[m]为顶点的直线不重叠 */
			isInt = 0;
			vect0[0] = linepj[1][0] - linepj[0][0];
			vect0[1] = linepj[1][1] - linepj[0][1];
			vect1[0] = facepj[(m + 1) % 3][0] - facepj[m][0];
			vect1[1] = facepj[(m + 1) % 3][1] - facepj[m][1];
			vect2[0] = facepj[(m + 2) % 3][0] - facepj[m][0];
			vect2[1] = facepj[(m + 2) % 3][1] - facepj[m][1];

			scale1 = vect0[0] * vect1[0] + vect0[1] * vect1[1];
			scale2 = vect0[0] * vect2[0] + vect0[1] * vect2[1];

			if (scale1 > 0.0)
			{
				isInt = GEOM_FUNC::orient2d(facepj[m], facepj[(m + 1) % 3], linepj[1]) == 0.0;
			}
			if (isInt == 0 && scale2 > 0.0)
			{
				isInt = GEOM_FUNC::orient2d(facepj[m], facepj[(m + 2) % 3], linepj[1]) == 0.0;
			}
		}
		if (isInt == 1)
			isInt = isInt;
	}

	return isInt;
}

int SPRImpl::isintersect(int iface[], int iline[], APOINT facep[], APOINT linep[], int colineNodes[], int colineNSize)
{
	int isInt = 0;
	int retValue, val, retValueCpy;
	APOINT pnt;
	double s1, s2;
	bool ln1Similar, ln2Similar;
	double facepj[3][2], linepj[3][2] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };	/* 共面情形的投影点 */
	double normalf[3];
	double nx, ny, nz;
	int i, m, i1, i2;
	double vect0[2], vect1[2], vect2[2], scale1, scale2;
	int fcSimilarCnt = 0, iSimilarIdx;
	int ilineCpy[2];
	APOINT linepCpy[2];
	int iTemp;
	double dTemp;

	/* 面上有几个节点是共面的 */
	for (m = 0, fcSimilarCnt = 0; m < 3; m++)
	{
		if (iface[m] == iline[0] || iface[m] == iline[1])
		{
			iSimilarIdx = m;
			fcSimilarCnt++;
		}
		else
		{
			for (i = 0; i < colineNSize; i++)
				if (colineNodes[i] == iface[m])
				{
					iSimilarIdx = m;
					fcSimilarCnt++;
					break;
				}
		}
	}

	switch (fcSimilarCnt)
	{
	case 0:
		/* 暂时借用成熟的函数，后面再统一处理 */
#ifdef _DEBUG
		retValue = GEOM_FUNC::lnFacInt(linep[0], linep[1], facep[0], facep[1], facep[2], pnt, &val);
		retValueCpy = GEOM_FUNC::lnFacInt2(linep[0], linep[1], facep[0], facep[1], facep[2],
			pnt, &val, -1, -1, -1, NULL, NULL, NULL, false);
		if (retValue != retValueCpy)
			retValue = retValueCpy;
#else
		retValue = GEOM_FUNC::lnFacInt2(linep[0], linep[1], facep[0], facep[1], facep[2],
			pnt, &val, -1, -1, -1, NULL, NULL, NULL, false);
#endif
		isInt = (retValue == -1) ? 0 : 1;
		break;
	case 1:
		ln1Similar = (iline[0] == iface[0] ||
			iline[0] == iface[1] ||
			iline[0] == iface[2]);
		ln2Similar = (iline[1] == iface[0] ||
			iline[1] == iface[1] ||
			iline[1] == iface[2]);

		if (ln1Similar || ln2Similar)
		{
			if (ln2Similar)
			{
				assert(!ln1Similar);
				iTemp = iline[0];
				iline[0] = iline[1];
				iline[1] = iTemp;
				for (m = 0; m < 3; m++)
				{
					dTemp = linep[0][m];
					linep[0][m] = linep[1][m];
					linep[1][m] = dTemp;
				}
			}
			else
				assert(!ln2Similar);
			isInt = isintersect_oneSharePoint(iface, iline, facep, linep);
		}
		else
		{
			/* 线段分成2段，调用isIntersect */
			ilineCpy[0] = iface[iSimilarIdx];
			ilineCpy[1] = iline[1];
			for (m = 0; m < 3; m++)
			{
				linepCpy[0][m] = vertices[ilineCpy[0]][m];
				linepCpy[1][m] = vertices[ilineCpy[1]][m];
			}
			isInt = isintersect(iface, ilineCpy, facep, linepCpy, NULL, 0);

			if (!isInt)
			{
				ilineCpy[1] = iface[iSimilarIdx];
				ilineCpy[0] = iline[0];
				for (m = 0; m < 3; m++)
				{
					linepCpy[0][m] = vertices[ilineCpy[0]][m];
					linepCpy[1][m] = vertices[ilineCpy[1]][m];
				}
				isInt = isintersect(iface, ilineCpy, facep, linepCpy, NULL, 0);
			}
		}
		break;
	case 2:
		isInt = 0; /* 边为面片的子线段 */
		break;
	case 3:
		isInt = 1; /* 3点共线的面片，自然不允许 */
		break;
	}

	return isInt;
}

int SPRImpl::isintersect(int iface[], APOINT facep[], Constraints *cst)
{
	int ipnt[3];
	APOINT pnt[3];
	int iline[2];
	APOINT linep[2];
	int isInt = 0;
	int i, m;
	bool noFaceTest = false;
	bool fc;
	int diffVtxCnt = 0;
	int diffVtx[3];

	/* 线约束 */
	for (i = 0; i < cst->lineCst_size; i++)
	{
		ipnt[0] = cst->lineCsts[i][0];
		ipnt[1] = cst->lineCsts[i][1];
		for (m = 0; m < 3; m++)
		{
			pnt[0][m] = vertices[ipnt[0]][m];
			pnt[1][m] = vertices[ipnt[1]][m];
		}
		isInt = isintersect(iface, ipnt, facep, pnt, cst->colineNodes, cst->colineN_size);
		if (isInt)
			return isInt;
	}

	/* 面约束 */
	if (cst->faceCst_size > 0 && cst->cofaceN_size > 0)
	{
		/* 如果3个顶点均在共面顶点序列，则直接返回不相交 */
		for (m = 0; m < 3; m++)
		{
			for (i = 0; i < cst->cofaceN_size; i++)
			{
				if (iface[m] == cst->cofaceNodes[i])
					break;
			}
			if (i >= cst->cofaceN_size)
			{
				diffVtx[diffVtxCnt++] = m;//这个顶点不在共面顶点列表
			}
		}
	}

	if (cst->cofaceN_size > 0 && diffVtxCnt > 0)
	{
		for (i = 0; i < cst->faceCst_size; i++)
		{
			ipnt[0] = cst->faceCsts[i][0];
			ipnt[1] = cst->faceCsts[i][1];
			ipnt[2] = cst->faceCsts[i][2];

			for (m = 0; m < 3; m++)
			{
				pnt[0][m] = vertices[ipnt[0]][m];
				pnt[1][m] = vertices[ipnt[1]][m];
				pnt[2][m] = vertices[ipnt[2]][m];
			}

			if (diffVtxCnt == 3)
			{
				//			isInt = isintersect(iface, ipnt, facep, pnt, NULL);
				isInt = GEOM_FUNC::tri_tri_intersect3d(iface, ipnt, facep, pnt);
			}
			else if (diffVtxCnt == 1)
			{/* 这个点不在平面上即不相交 */
				isInt = GEOM_FUNC::orient3d(pnt[0], pnt[1], pnt[2], facep[diffVtx[0]]) == 0.0;
			}
			else if (diffVtxCnt == 2)
			{/* 判断这条线和面是否相交 */
				iline[0] = iface[diffVtx[0]];
				iline[1] = iface[diffVtx[1]];
				for (m = 0; m < 3; m++)
				{
					linep[0][m] = vertices[iline[0]][m];
					linep[1][m] = vertices[iline[1]][m];
				}
				isInt = isintersect(ipnt, iline, pnt, linep, NULL, 0);
			}
			if (isInt)
				return isInt;
		}
	}

	return isInt;
}

int SPRImpl::isIntersectEdgeConstraints(int iface[], APOINT facep[])
{
	int i, m, k, simCnt = 0;
	InstEntList *pList = NULL;
	APOINT p1, p2, pnt[2];
	int isInt = 0;
	int ipnt[2];

	if (edgeConsArraySize <= 0)
		return 0;

	/* 如果面本身包含必须删除的边或面，返回1 */
	if (spr_optimal_type != EDGE_REMOVAL)
	{
		if (edgeConsArray[0].j3 < 0)
		{/* 一条边 */
			for (m = 0, simCnt = 0; m < 3; m++)
			{
				if (iface[m] == edgeConsArray[0].j1 || iface[m] == edgeConsArray[0].j2)
					simCnt++;
			}
			if (simCnt == 2)
				return 1;
		}
		else
		{/* 一个面 */
			for (m = 0, simCnt = 0; m < 3; m++)
			{
				if (iface[m] == edgeConsArray[0].j1 || iface[m] == edgeConsArray[0].j2 || iface[m] == edgeConsArray[0].j3)
					simCnt++;
			}
			if (simCnt == 3)
				return 1;
		}
	}
	else
	{
		assert(spr_shell_edge_nodeS >= 0);
		assert(spr_shell_edge_nodeE >= 0);
	}

	/*
	 * 如果面是允许的边或面, 则返回0；一个面称之为允许的面，需满足以下条件：
	 * 1. 该面在允许列表中；或
	 * 2. 该面有2条边在允许列表中；或
	 * 3. 该面有1条边在允许列表中：
	 *    a.另1个点为约束线的起点或终点；或
	 *    b.且该面不和约束线共面
	 */

	ipnt[0] = edgeConsArray[0].i1;
	ipnt[1] = edgeConsArray[0].i2;
	for (k = 0; k < 3; k++)
	{
		pnt[0][k] = vertices[edgeConsArray[0].i1][k];
		pnt[1][k] = vertices[edgeConsArray[0].i2][k];
	}
	pList = edgeConsArray[0].pListAllowed;
	if (pList)
	{
		if (pList->findFace(iface[0], iface[1], iface[2]) > 0)
			return 0;
		else
		{
			for (m = 0, simCnt = 0; m < 3; m++)
			{
				if (pList->findEdge(iface[m], iface[(m + 1) % 3]) > 0)
				{
					if (iface[(m + 2) % 3] == edgeConsArray[0].i1 ||
						iface[(m + 2) % 3] == edgeConsArray[0].i2)
						return 0;
					simCnt++;
				}
			}
			if (simCnt >= 2)
				return 0;
			else if (simCnt == 1)
			{
				if (GEOM_FUNC::orient3d(facep[0], facep[1], facep[2], pnt[0]) != 0.0 ||
					GEOM_FUNC::orient3d(facep[0], facep[1], facep[2], pnt[1]) != 0.0)
					return 0;
			}
		}
	}


	/* 调用线面相交算法，确定2者是否相交 */
	isInt = isintersect(iface, ipnt, facep, pnt, NULL, 0);
	if (isInt)
		return isInt;

	/* 看看是否和已经恢复的边相交 */
	for (i = 1; i < edgeConsArraySize; i++)
	{
		assert(edgeConsArray[i].type == MUST_EXIST);
		ipnt[0] = edgeConsArray[i].i1;
		ipnt[1] = edgeConsArray[i].i2;
		for (k = 0; k < 3; k++)
		{
			pnt[0][k] = vertices[edgeConsArray[i].i1][k];
			pnt[1][k] = vertices[edgeConsArray[i].i2][k];
		}
		isInt = isintersect(iface, ipnt, facep, pnt, NULL, 0);
		if (isInt)
			return isInt;
	}

	return 0;
}

//*****************************************************

//几个原则：四面体中新增加的线可以和别的三角形共享
//返回0表示有效
//zhaodawei modify 2010-07-16
int SPRImpl::isintersect(TRI face1, TRI face2, APOINT facep1[], APOINT facep2[], struct polyhedron *P)
{
	int sameVtxCnt, i;
	int nRet = 0;
	int first_index, second_index;

	face_n++;

	for (i = 0, sameVtxCnt = 0; i < 3; i++)
	{
		if ((face1[i] == face2[0]) ||
			(face1[i] == face2[1]) ||
			(face1[i] == face2[2]))
		{
			sameVtxCnt++;
		}
	}
	if (sameVtxCnt == 1)
	{
		for (i = 0; i < 3; i++)
		{
			if (facep1[i][0] == facep2[0][0]
				&& facep1[i][1] == facep2[0][1]
				&& facep1[i][2] == facep2[0][2])  //zhaodawei
			{
				first_index = i;
				second_index = 0;
			}
			else if (facep1[i][0] == facep2[1][0]
				&& facep1[i][1] == facep2[1][1]
				&& facep1[i][2] == facep2[1][2])
			{
				first_index = i;
				second_index = 1;
			}
			else if (facep1[i][0] == facep2[2][0]
				&& facep1[i][1] == facep2[2][1]
				&& facep1[i][2] == facep2[2][2])
			{
				first_index = i;
				second_index = 2;
			}
		}
	}


	if (sameVtxCnt == 0)
	{
		nRet = GEOM_FUNC::tri_tri_overlap_test_3d(facep1[0], facep1[1], facep1[2],
			facep2[0], facep2[1], facep2[2]);
	}
	else if (sameVtxCnt == 1)
	{
#ifdef _TIMING_PERFORMANCE
		starttmp[8] = SPRlogTime();
#endif
		nRet = isintersect_old(face1, face2, facep1, facep2, P);
		if (first_index == 1)
		{
			swap_apoint(facep1[0], facep1[1]);
		}
		else if (first_index == 2)
		{
			swap_apoint(facep1[0], facep1[2]);
		}

		if (second_index == 1)
		{
			swap_apoint(facep2[0], facep2[1]);
		}
		else if (second_index == 2)
		{
			swap_apoint(facep2[0], facep2[2]);
		}
		nRet = GEOM_FUNC::one_node_same_tri_tri_overlap_3d(facep1[0], facep1[1], facep1[2],
			facep2[0], facep2[1], facep2[2]);

#ifdef _TIMING_PERFORMANCE
		stoptmp[8] = SPRlogTime();
		dura[8] += (stoptmp[8] - starttmp[8]);
#endif
	}

	return nRet;
}

int SPRImpl::isintercross(TRI face1, TRI face2, APOINT facep1[], APOINT facep2[], int m, int n, int *res, struct polyhedron *P)
{
	int tmp;
	int l, k;
	//	double V[2];

	tmp = isborder(face1, face2, m, n, &l, &k);

	if (tmp == 2)	//这条线是和三角形共线
	{
		*res = 1;
		return 0;
	}
	else if (tmp == 1)	//有一个点是相同的，这样就比较麻烦，如果是异面那可以，如果是共面就避免相交
	{
		/*
				V[0]=fabs(isVisible(facep2[0], facep2[1], facep2[2], facep1[m]));
				V[1]=fabs(isVisible(facep2[0], facep2[1], facep2[2], facep1[n]));
				if(V[0] < SPR_SPR_EPS_ZERO && V[1] < SPR_SPR_EPS_ZERO)		//这时认为这条线是和这个三角形共面
		*/
		Call[1]++;
		if (fabs(isVisible_Search(face2[0], face2[1], face2[2], face1[m + n - k])) < SPR_EPS_ZERO)
		{

			if (islinecross(face1[m], face1[n], face2[(l + 1) % 3], face2[(l + 2) % 3], facep1[m], facep1[n], facep2[(l + 1) % 3], facep2[(l + 2) % 3], P))		//如果这两条线相交的话
				return 1;
			/*
						if(l == 0)
						{
							if(islinecross(facep1[m], facep1[n], facep2[1], facep2[2]))		//如果这两条线相交的话
								return 1;
						}
						else if(l == 1)
						{
							if(islinecross(facep1[m], facep1[n], facep2[0], facep2[2]))
								return 1;
						}
						else if(l == 2)
						{
							if(islinecross(facep1[m], facep1[n], facep2[0], facep2[1]))
								return 1;
						}
			*/
			*res = 1;
			return 0;
		}
		else	//如果是有一个点相同，然后又是异面，那么一定是可以的
		{
			*res = 1;
			return 0;
		}
	}
	else
	{
		*res = 0;
		return 0;
	}
}

int SPRImpl::islinecrossface(TRI face1, TRI face2, APOINT facep1[], APOINT facep2[], int m, int n, int res, struct polyhedron *P)
{
	double V[2];
	int result;

	/*
	//以下代码测试使用
	int tmp;
	double arear;
	APOINT pnt;
	int val;
	double w1, w2, w3;

	arear=areas(facep2[0][0],facep2[0][1],facep2[0][2],
				facep2[1][0],facep2[1][1],facep2[1][2],
				facep2[2][0],facep2[2][1],facep2[2][2]);
	*/

	if (res == 1)
		return 0;

	line_n++;
	if (!(isboxinter_l(facep1[m], facep1[n], facep2)))
	{
		box_l_n++;
		return 0;
	}

	Call[2] += 2;
	V[0] = fabs(isVisible_Search(face2[0], face2[1], face2[2], face1[m]));
	V[1] = fabs(isVisible_Search(face2[0], face2[1], face2[2], face1[n]));

	if (V[0] < SPR_EPS_ZERO && V[1] < SPR_EPS_ZERO)		//这时认为这条线是和这个三角形共面
	{
		if (islinecross(face1[m], face1[n], face2[0], face2[1], facep1[m], facep1[n], facep2[0], facep2[1], P) ||
			islinecross(face1[m], face1[n], face2[1], face2[2], facep1[m], facep1[n], facep2[1], facep2[2], P) ||
			islinecross(face1[m], face1[n], face2[0], face2[2], facep1[m], facep1[n], facep2[0], facep2[2], P))
			return 1;
		else
			return 0;
	}
	else
	{
#ifdef _TIMING_PERFORMANCE
		starttmp[6] = SPRlogTime();
#endif

		result = lnFacInt_own(face1[m], face1[n], face2, facep1[m], facep1[n], facep2[0], facep2[1], facep2[2], P);

#ifdef _TIMING_PERFORMANCE
		stoptmp[6] = SPRlogTime();
		dura[6] += (stoptmp[6] - starttmp[6]);
#endif


		/*
		tmp=lnFacInt(facep1[m], facep1[n], facep2[0], facep2[1], facep2[2],
						distance(facep1,m,n), arear,
						pnt, &val, &w1, &w2, &w3);
		if((tmp != -1 && result == 0) || (tmp == -1 && result == 1))
			spdlog::info("Wrong!");
		*/

		// 		if(res == 1)
		// 			return 1;
		// 		else
		// 			return 0;
		/*
				res=lnFacInt(facep1[m], facep1[n], facep2[0], facep2[1], facep2[2],
								distance(facep1,m,n), arear,
								pnt, &val, &w1, &w2, &w3);
				if(res == -1)
					return 0;
				else
					return 1;
		*/
		return result;
	}
}



int SPRImpl::isboxinter(APOINT facep1[], APOINT facep2[])
{
	int i;
	for (i = 0; i < 3; i++)
	{
		if (MIN_VALUE(MIN_VALUE(facep1[0][i], facep1[1][i]), facep1[2][i]) >
			MAX_VALUE(MAX_VALUE(facep2[0][i], facep2[1][i]), facep2[2][i]) ||
			MAX_VALUE(MAX_VALUE(facep1[0][i], facep1[1][i]), facep1[2][i]) <
			MIN_VALUE(MIN_VALUE(facep2[0][i], facep2[1][i]), facep2[2][i]))
			return 0;
	}
	return 1;
}



int SPRImpl::isboxinter_l(APOINT A, APOINT B, APOINT facep2[])
{
	int i;

	for (i = 0; i < 3; i++)
	{
		if (MIN_VALUE(A[i], B[i]) > MAX_VALUE(MAX_VALUE(facep2[0][i], facep2[1][i]), facep2[2][i]) ||
			MAX_VALUE(A[i], B[i]) < MIN_VALUE(MIN_VALUE(facep2[0][i], facep2[1][i]), facep2[2][i]))
			return 0;
	}
	return 1;
}

//分析face1中的m和n所代表的边是否就是face2中的边，如果是的话，返回2，否则返回有多少个点是相同的
//在face1和face2内的数据没有错误的前提下（很重要），返回值只可能是2，1，0
//l和k是用来记录是那个点相同
int SPRImpl::isborder(TRI face1, TRI face2, int m, int n, int *l, int *k)
{
	int i;
	int a, b;
	int flag;
	a = face1[m];
	b = face1[n];
	flag = 0;
	for (i = 0; i < 3; i++)
	{
		if (a == face2[i])
		{
			*l = i;
			*k = m;
			flag++;
		}
		else if (b == face2[i])
		{
			*l = i;
			*k = n;
			flag++;
		}
	}

	return flag;
}

int SPRImpl::islinecross(int index1, int index2, int index3, int index4, APOINT A, APOINT B, APOINT C, APOINT D, struct polyhedron *P)
{
	int res[2];
	//	res[0] = islinecross_new(A, B, C, D);
	//	res[1] = islinecross_old(A, B, C, D);
	res[0] = islinecross_new(A, B, C, D);//islinecross_newer(index1, index2, index3, index4, P);

	return res[0];
}


int SPRImpl::islinecross_new(APOINT A, APOINT B, APOINT C, APOINT D)
{
	double res[2];

	res[0] = GEOM_FUNC::orient2d(A, B, C);
	res[1] = GEOM_FUNC::orient2d(A, B, D);

	if (res[0] * res[1] < 0)
	{
		res[0] = GEOM_FUNC::orient2d(C, D, A);
		res[1] = GEOM_FUNC::orient2d(C, D, B);
		if (res[0] * res[1] < 0)
			return 1;
	}
	return 0;
}

//先找到一个点不在这个面上面，这样，可以用这个点作为一个基准来思考，现在关键是找到一个合理的方法来找到这个点
//可以保证速度，同时，确保鲁棒性
int SPRImpl::islinecross_newer(int index1, int index2, int index3, int index4, struct polyhedron *P)
{
	double res[2];
	int i, j;
	double max, ab_m, ac_m, ad_m, bc_m, bd_m;
	APOINT AB, AC, AD;

	res[0] = 0;
	res[1] = 0;

	for (i = 0; i < num_vertices && fabs(res[0]) < SPR_EPS_ZERO; i++)
		res[0] = isVisible_Search(index1, index2, index3, i);

	if (fabs(res[0]) < SPR_EPS_ZERO)
	{
		for (i = 0; i < num_vertices && fabs(res[1]) < SPR_EPS_ZERO; i++)
			res[1] = isVisible_Search(index1, index2, index4, i);
		if (fabs(res[1]) < SPR_EPS_ZERO)		//说明AB、CD共线，通过向量的方法来判别
		{
			for (i = 0; i < 3; i++)
			{
				AB[i] = vertices[index2][i] - vertices[index1][i];
				AC[i] = vertices[index3][i] - vertices[index1][i];
				AD[i] = vertices[index4][i] - vertices[index1][i];
			}

			max = SPR_EPS_ZERO;

			for (i = 0; i < 3; i++)
			{
				if (AB[i] > max)
				{
					max = AB[i];
					j = i;
				}
			}

			ab_m = AB[j];
			ac_m = AC[j];
			ad_m = AD[j];

			if (ac_m * ad_m < -1 * SPR_EPS_ZERO)		//如果两者异号，则自然相交
				return 1;
			else if (ac_m * ab_m < -1 * SPR_EPS_ZERO)	//如果ac，ad是同号，而ac与ab是异号，则表明不相交
				return 0;
			else
			{
				bc_m = vertices[index3][j] - vertices[index2][j];
				bd_m = vertices[index4][j] - vertices[index2][j];
				if (bc_m * bd_m < -1 * SPR_EPS_ZERO)		//如果bc，bd是异号，则表明是相交
					return 1;
				else if (ac_m * bc_m < -1 * SPR_EPS_ZERO)	//如果ac，bc是异号，那么表明是相交
					return 1;
				else
					return 0;	//此时，cd在ab的同侧，自然是不相交
			}
			//至此，向量判断结束，只是有一个符号比较，无需比较大小，不会有误差
		}
	}

	i--;	//把i值恢复正常

	res[1] = isVisible_Search(index1, index2, index4, i);

	if (res[0] * res[1] < SPR_EPS_ZERO)		//这里很重要，因为res[0]*res[1]可能为零，如果是零的话说明有共线发生，自然是相交关系
	{
		res[0] = isVisible_Search(index3, index4, index1, i);
		res[1] = isVisible_Search(index3, index4, index2, i);
		if (res[0] * res[1] < SPR_EPS_ZERO)
			return 1;
	}
	return 0;
}

// 以下为原来的做法，不严密
int SPRImpl::islinecross_old(APOINT A, APOINT B, APOINT C, APOINT D)
{
	int i;
	APOINT AB, AC, AD, CD, CA, CB;
	APOINT res1, res2;
	double a, b, c;

	//该函数为自己编写，主要的思想是向量叉乘，再比较方向

	for (i = 0; i < 3; i++)	//获得AB、AC、AD向量
	{
		AB[i] = B[i] - A[i];
		AC[i] = C[i] - A[i];
		AD[i] = D[i] - A[i];
	}
	res1[0] = AB[1] * AC[2] - AB[2] * AC[1];	//得到叉乘的结果
	res1[1] = AB[2] * AC[0] - AB[0] * AC[2];
	res1[2] = AB[0] * AC[1] - AB[1] * AC[0];

	res2[0] = AB[1] * AD[2] - AB[2] * AD[1];
	res2[1] = AB[2] * AD[0] - AB[0] * AD[2];
	res2[2] = AB[0] * AD[1] - AB[1] * AD[0];

	a = res1[0] * res2[0];	//分析符号是否相反
	b = res1[1] * res2[1];
	c = res1[2] * res2[2];

	if (a < SPR_EPS_ZERO && b < SPR_EPS_ZERO && c < SPR_EPS_ZERO)		//如果符号相反的话，那么就说明CD两点位于AB的两侧
	{
		for (i = 0; i < 3; i++)		//然后再判断AB在CD的两侧，原理相同
		{
			CD[i] = D[i] - C[i];
			CA[i] = A[i] - C[i];
			CB[i] = B[i] - C[i];
		}
		res1[0] = CD[1] * CA[2] - CD[2] * CA[1];
		res1[1] = CD[2] * CA[0] - CD[0] * CA[2];
		res1[2] = CD[0] * CA[1] - CD[1] * CA[0];

		res2[0] = CD[1] * CB[2] - CD[2] * CB[1];
		res2[1] = CD[2] * CB[0] - CD[0] * CB[2];
		res2[2] = CD[0] * CB[1] - CD[1] * CB[0];

		a = res1[0] * res2[0];
		b = res1[1] * res2[1];
		c = res1[2] * res2[2];

		if (a < SPR_EPS_ZERO && b < SPR_EPS_ZERO && c < SPR_EPS_ZERO)		//此时说明两个线段之间相交
			return 1;
	}
	return 0;
}

inline double distance_square(APOINT p1, APOINT p2)
{
	double dx = p2[0] - p1[0];
	double dy = p2[1] - p1[1];
	double dz = p2[2] - p1[2];
	return dx * dx + dy * dy + dz * dz;
}

/*
几个原则，当012本身是指向四面体内部时，那么，321、213、132；230，302，023；310，103，031也都是指向内部
肯定要一个默认的规则，不然会很乱，约定第一数字要是最小的，这样判断起来容易些
*/
//该函数已经整体调试通过，又解决了一个大的问题

int SPRImpl::delELE(int ELE[4], struct polyhedron *P, struct polyhedron *Q, int f, int *newCnt)
{
	int i, j, k, m;
	int add[4][3], del[4] = { -1, -1, -1, -1 };	//这个del是记录是否删除这个面
	int src_del[MAX_SPR_POLY_SIZE];

	//下面的代码主要是在确定各个面的法线方向的指向
	if (ELE[3] < ELE[2] && ELE[3] < ELE[1])
	{/* ELE[3] 最小 */
		add[0][0] = ELE[3];
		add[0][1] = ELE[2];
		add[0][2] = ELE[1];
	}
	else if (ELE[2] < ELE[1] && ELE[2] < ELE[3])
	{/* ELE[2] 最小 */
		add[0][0] = ELE[2];
		add[0][1] = ELE[1];
		add[0][2] = ELE[3];
	}
	else
	{/* ELE[1] 最小 */
		add[0][0] = ELE[1];
		add[0][1] = ELE[3];
		add[0][2] = ELE[2];
	}

	if (ELE[3] < ELE[2] && ELE[3] < ELE[0])
	{/* ELE[3] 最小 */
		add[1][0] = ELE[3];
		add[1][1] = ELE[0];
		add[1][2] = ELE[2];
	}
	else if (ELE[2] < ELE[0] && ELE[2] < ELE[3])
	{/* ELE[2] 最小 */
		add[1][0] = ELE[2];
		add[1][1] = ELE[3];
		add[1][2] = ELE[0];
	}
	else
	{/* ELE[0] 最小 */
		add[1][0] = ELE[0];
		add[1][1] = ELE[2];
		add[1][2] = ELE[3];
	}

	if (ELE[3] < ELE[1] && ELE[3] < ELE[0])
	{/* ELE[3] 最小 */
		add[2][0] = ELE[3];
		add[2][1] = ELE[1];
		add[2][2] = ELE[0];
	}
	else if (ELE[1] < ELE[0] && ELE[1] < ELE[3])
	{/* ELE[1] 最小 */
		add[2][0] = ELE[1];
		add[2][1] = ELE[0];
		add[2][2] = ELE[3];
	}
	else
	{/* ELE[0] 最小 */
		add[2][0] = ELE[0];
		add[2][1] = ELE[3];
		add[2][2] = ELE[1];
	}

	memset(src_del, -1, sizeof(int)*MAX_SPR_POLY_SIZE);
	//上面的判断添加的代码已经调试通过，运行效果符合预期
	for (i = 0; i < P->num_t; i++)
	{
		for (j = 0; j < 3; j++)
		{
			if (del[j] == -1 &&
				P->t[i][0] == add[j][0] &&
				P->t[i][1] == add[j][1] &&
				P->t[i][2] == add[j][2])	//如果本来就有这个面，那么就有必要将其删掉
			{
				src_del[i] = j;
				del[j] = i;	//记录要删除的面的索引，和对应的点的坐标
			}
		}

		if (del[0] != -1 && del[1] != -1 && del[2] != -1)
			break;
	}
	del[3] = f;
	src_del[f] = 3;

	for (i = 0; i < 3; i++)
		add[3][i] = ELE[i];	//把基面的坐标点暂存到add数组中，便于以后编程

	//将点的使用情况的数据转存到used_v数组中，防止修改P中的数据
	for (i = 0; i < MAX_SPR_POLY_SIZE; i++)
		Q->used_v[i] = P->used_v[i];

	/* 将P中的面拷贝到Q中 */
	for (m = 0, k = 0; m < P->num_t; m++)
	{
		if (src_del[m] >= 0)
		{
			for (j = 0; j < 3; j++)
				Q->used_v[add[src_del[m]][j]]--;		//对于要删去的面，不仅仅不将其转入Q中，反倒是降低其引用次数
		}
		else
		{
			for (j = 0; j < 3; j++)
				Q->t[k][j] = P->t[m][j];		//把这个面加入到Q中，同时，移动写入的位置
			k++;
		}
	}

	/* 将新增加的三个面拷贝到Q中 */
	*newCnt = 0;
	for (i = 0; i < 3; i++)
	{
		if (del[i] == -1)		//说明这个面需要加进去，但是，要注意旋转方向
		{
			Q->t[k][0] = add[i][0];
			Q->t[k][1] = add[i][2];	//旋转方向，这个相当的重要
			Q->t[k][2] = add[i][1];

			Q->used_v[add[i][0]]++;	//增加对应面的点的引用次数
			Q->used_v[add[i][1]]++;
			Q->used_v[add[i][2]]++;
			k++;
			(*newCnt)++;
		}
	}
	Q->num_t = k;

	//处理点的问题, 记录可用点的数目
	//直接把所有的点都复制到Q里面去，然后再通过点的引用次数来表征这个点是否还存在
	for (i = 0, k = 0; i < num_vertices; i++)
	{
		if (Q->used_v[i] > 0)
			k++;
	}
	Q->num_vable = k;		//这里存的就是现在还可以使用的点的数目

	return 1;
}

#ifdef _SPR_USING_STD_VECTOR
int SPRImpl::delELE(int ELE[4], struct polyhedron *P, struct polyhedron *Q,
	std::vector<struct polyhedron *>& vec_simp_poly, int f, bool *freeSimpPoly,
	Constraints *csts, int maxIntNum, int *intNum)
#else
int SPRImpl::delELE(int ELE[4], struct polyhedron *P, struct polyhedron *Q,
	polyhedron * simpPolyArray[], int *simpPolyNum, int f, bool *freeSimpPoly,
	Constraints *csts, int maxIntNum, int *intNum)
#endif /* _SPR_USING_STD_VECTOR */
{

	InsertedTetra instTetra;
	int i;

	if (fill_inserted_tetra(ELE, f, &instTetra, P, Q))
	{
		if (isNewFacetsValid_BruteForce(Q, 4 - instTetra.num_del, csts, maxIntNum, intNum, &instTetra, P))
		{
			Q->newf_ct = 4 - instTetra.num_del;	/* 新增加的面 */
			memset(Q->newf_fg, 0, sizeof(int)*Q->num_t);
			for (i = 0; i < Q->newf_ct; i++)
			{
				Q->newf_fg[Q->num_t - 1 - i] = 1;
				Q->newf_id[i] = Q->num_t - 1 - i;
			}

			//		instTetra.complex = false;
			if (instTetra.complex)
			{
				*freeSimpPoly = true;
#ifdef _SPR_USING_STD_VECTOR
				divide_poly(Q, vec_simp_poly);
#else
				divide_poly(Q, simpPolyArray, simpPolyNum);
#endif
			}
			else
			{
				*freeSimpPoly = false;
#ifdef _SPR_USING_STD_VECTOR
				vec_simp_poly.push_back(Q);
#else
				simpPolyArray[0] = Q;
				*simpPolyNum = 1;
#endif
			}
			return 1;
		}
	}
	return 0;
}

//函数是把ELE这个四面体给整合到Tc里面，函数会返回一个指针，即指向Tc
inline void mergeT(struct divide *Tq, int ELE[])
{
	int j;

	assert(Tq->num >= 0);
	for (j = 0; j < 4; j++)
		Tq->te[Tq->num][j] = ELE[j];		//再把ELE加入到Tc中
	Tq->num++;							//增加Tc的数目
}

inline void mergeT(struct divide *sT, struct divide *tT)
{
	int j;

	for (j = 0; j < tT->num; j++)
		memcpy(sT->te[sT->num++], tT->te[j], sizeof(TETRA));
}

//函数把Tc中的内容整体复制到T中去，内容比较简单
inline void copyT(struct divide *Tc, struct divide *T)
{
	int i, j;
	for (i = 0; i < MAX_SPR_POLY_SIZE; i++)
	{
		for (j = 0; j < 4; j++)
			T->te[i][j] = Tc->te[i][j];
	}
	T->num = Tc->num;
	T->q = Tc->q;
}

/*
 * 三点面积绝对值
 * absolute value of area of a triangle
 这里是只要给定x，y，z的三点坐标，就可以求出这个三角形的面积
 */
inline double areas(double xa, double ya, double za,
	double xb, double yb, double zb,
	double xc, double yc, double zc)
{
	double aa = 0.0;
	double pq1, pq2, pq3, pr1, pr2, pr3, aq1, aq2, aq3;

	pq1 = xb - xa;
	pq2 = yb - ya;
	pq3 = zb - za;
	pr1 = xc - xa;
	pr2 = yc - ya;
	pr3 = zc - za;
	aq1 = pq2 * pr3 - pq3 * pr2;
	aq2 = -(pq1*pr3 - pq3 * pr1);
	aq3 = pq1 * pr2 - pq2 * pr1;
	aa = 0.5*sqrt(aq1*aq1 + aq2 * aq2 + aq3 * aq3);
	return aa;
}


/*
 * p4点对三角面p1p2p3是否可视，即p1p2p3p4的体积是否为正
 * check if p4 is visible to p1p2p3, i.e. the volume of p1p2p3p4 is positive
 遵从右手法则，就意味着在设计程序的时候要保证按照右手法则这个方向始终是指向里的！！
 这里对原来的代码有修改，使得函数直接返回这个四面体的体积
 */

inline double isVisible_SPR(APOINT p1, APOINT p2, APOINT p3, APOINT p4)
{
	double aa, bb, cc, id1, id2, d3, d4, dd, hh, al, ee, ff, gg, ai, aj, ak;
	double s1, s2, s3, s4, s5, s6;
	double te;

	aa = 2 * (p2[0] - p1[0]);
	bb = 2 * (p2[1] - p1[1]);
	cc = 2 * (p2[2] - p1[2]);
	id1 = p1[0] * p1[0] + p1[1] * p1[1] + p1[2] * p1[2];
	id2 = p2[0] * p2[0] + p2[1] * p2[1] + p2[2] * p2[2];
	d3 = p3[0] * p3[0] + p3[1] * p3[1] + p3[2] * p3[2];
	d4 = p4[0] * p4[0] + p4[1] * p4[1] + p4[2] * p4[2];
	dd = id2 - id1;
	hh = d3 - id1;
	al = d4 - id1;
	ee = 2 * (p3[0] - p1[0]);
	ff = 2 * (p3[1] - p1[1]);
	gg = 2 * (p3[2] - p1[2]);
	ai = 2 * (p4[0] - p1[0]);
	aj = 2 * (p4[1] - p1[1]);
	ak = 2 * (p4[2] - p1[2]);
	s1 = ff * ak - aj * gg;
	s2 = ee * ak - ai * gg;
	s3 = ee * aj - ai * ff;
	s4 = hh * ak - al * gg;
	s5 = hh * aj - al * ff;
	s6 = ee * al - hh * ai;
	te = aa * s1 - bb * s2 + cc * s3;

	return te;
}

#if 0
inline double isVisible_Search(int index1, int index2, int index3, int index4)
{
	double res[2];

#ifdef _TIMING_PERFORMANCE
	Call_V++;
#endif

	if (Visible[index1][index2][index3][index4] == MAX_SPR_POLY_SIZE)	//表明没有进行过计算 A bug
	{
#ifdef _TIMING_PERFORMANCE	
		orient3d_n++;
		starttmp[7] = SPRlogTime();
#endif

		res[0] = orient3d(vertices[index1], vertices[index3], vertices[index2], vertices[index4]);		//orient3d的定义与isVisible是相反的，这点很重要

#ifdef _TIMING_PERFORMANCE
		stoptmp[7] = SPRlogTime();
		dura[7] += (stoptmp[7] - starttmp[7]);
#endif

#if 1
		Visible[index1][index2][index3][index4] = res[0];

		res[1] = -res[0];	//res[1]是res[0]的相反数

		//后面是用穷举法弄得，先并未考虑更好的表达形式
		//校对过没有笔误，调试下面的代码的工作也完全正常

		Visible[index1][index2][index4][index3] = res[1];
		Visible[index1][index3][index2][index4] = res[1];
		Visible[index1][index3][index4][index2] = res[0];
		Visible[index1][index4][index3][index2] = res[1];
		Visible[index1][index4][index2][index3] = res[0];

		Visible[index2][index1][index3][index4] = res[1];
		Visible[index2][index1][index4][index3] = res[0];
		Visible[index2][index3][index1][index4] = res[0];
		Visible[index2][index3][index4][index1] = res[1];
		Visible[index2][index4][index1][index3] = res[1];
		Visible[index2][index4][index3][index1] = res[0];

		Visible[index3][index1][index2][index4] = res[0];
		Visible[index3][index1][index4][index2] = res[1];
		Visible[index3][index2][index1][index4] = res[1];
		Visible[index3][index2][index4][index1] = res[0];
		Visible[index3][index4][index2][index1] = res[1];
		Visible[index3][index4][index1][index2] = res[0];

		Visible[index4][index1][index3][index2] = res[0];
		Visible[index4][index1][index2][index3] = res[1];
		Visible[index4][index2][index1][index3] = res[0];
		Visible[index4][index2][index3][index1] = res[1];
		Visible[index4][index3][index2][index1] = res[0];
		Visible[index4][index3][index1][index2] = res[1];
#endif

		return res[0];
	}
	else
	{
#ifdef _TIMING_PERFORMANCE
		Call_V_R++;
#endif

		return Visible[index1][index2][index3][index4];
	}
}
#else
double SPRImpl::isVisible_Search(int j1, int j2, int j3, int j4, double *volume)
{
	int key;
	double v = 0.0;
	double AB, BC, AC, AD, BD, CD;

#ifdef _SPR_MAP_TETRA_VOLUME
	std::map<int, double>::iterator it;
	int ort = 1;
#endif

#ifdef _TIMING_PERFORMANCE
	Call_V++;
#endif

#ifdef _SPR_MAP_TETRA_VOLUME
	ort = sort_quad_pair(&j1, &j2, &j3, &j4);
	key = make_quad_pair(j1, j2, j3, j4);

	it = mapTetraVolume.find(key);

	if (it == mapTetraVolume.end())
#else
	if (true)
#endif
	{
#ifdef _TIMING_PERFORMANCE	
		orient3d_n++;
		starttmp[7] = SPRlogTime();
#endif
		v = GEOM_FUNC::orient3d(vertices[j4], vertices[j1], vertices[j2], vertices[j3]);

#ifdef _SPR_MAP_TETRA_VOLUME
		if (volume)
			*volume = v * ort;
#else
		if (volume)
			*volume = v;
#endif /* _SPR_MAP_TETRA_VOLUME */

		/* 计算质量，确定是否可行 */
		AB = distance_square(vertices[j1], vertices[j2]);
		BC = distance_square(vertices[j2], vertices[j3]);
		AC = distance_square(vertices[j1], vertices[j3]);
		AD = distance_square(vertices[j1], vertices[j4]);
		BD = distance_square(vertices[j2], vertices[j4]);
		CD = distance_square(vertices[j3], vertices[j4]);

		v /= pow(AB + BC + AC + AD + BD + CD, 1.5);	//72*sqrt(3)*算出这个四面体的质量


#ifdef _TIMING_PERFORMANCE
		stoptmp[7] = SPRlogTime();
		dura[7] += (stoptmp[7] - starttmp[7]);
#endif

#ifdef _SPR_MAP_TETRA_VOLUME
		mapTetraVolume[key] = v;
#endif
	}
#ifdef _SPR_MAP_TETRA_VOLUME
	else
	{
#ifdef _TIMING_PERFORMANCE
		Call_V_R++;
#endif
		v = it->second;
		if (volume)
			*volume = ort * GEOM_FUNC::orient3d(vertices[j4], vertices[j1], vertices[j2], vertices[j3]);
	}
#endif /* _SPR_MAP_TETRA_VOLUME */

#ifdef _SPR_MAP_TETRA_VOLUME
	return v * ort;
#else
	return v;
#endif /* _SPR_MAP_TETRA_VOLUME */
}
#endif

double SPRImpl::isVisible_Search(APOINT p1, APOINT p2, APOINT p3, APOINT p4)
{
	return  GEOM_FUNC::orient3d(p4, p1, p2, p3);
}
/* 我们要求可能新形成的3个面片位于被删除面片和其相邻面片形成的内角中间 */
int SPRImpl::isNeigborAllowed(int fidx, int nidx, polyhedron *poly)
{
	int neig, m, i1, i2, i3, i4;
	double quality;

	assert(poly->t[fidx][0] != nidx && poly->t[fidx][1] != nidx && poly->t[fidx][2] != nidx);
	for (m = 0; m < 3; m++)
	{
		neig = poly->n[fidx][m];
		i1 = poly->t[fidx][(m + 1) % 3];
		i2 = poly->t[fidx][(m + 2) % 3];
		i3 = poly->t[fidx][m];
		i4 = poly->t[neig][0] + poly->t[neig][1] + poly->t[neig][2] - i1 - i2;

		assert(i4 != i1 && i4 != i2 && i4 != i3);

		if (nidx == i4)
			continue; /* 这也是一个旧面片，所以略过 */

		if (isVisible_Search(vertices[i1], vertices[i2], vertices[i3], vertices[i4]) >= 0.0)
		{/* 内角张开的是一个小于180度的角 */
			if (isVisible_Search(vertices[i1], vertices[i2], vertices[i3], vertices[nidx]) <= 0.0 ||
				isVisible_Search(vertices[i1], vertices[i2], vertices[nidx], vertices[i4]) <= 0.0)
				return 0;
		}
		else
		{/* 内角张开的是一个大于180度的角 */
			if (isVisible_Search(vertices[i1], vertices[i2], vertices[i4], vertices[nidx]) >= 0.0 &&
				isVisible_Search(vertices[i1], vertices[i2], vertices[nidx], vertices[i3]) >= 0.0)
				return 0;
		}
	}

	return 1;
}

//dont be used
#if 0
namespace SPR
{
	REAL tetquality(APOINT vtx1, APOINT vtx2, APOINT vtx3, APOINT vtx4, int qualmeasure)
	{
		REAL quality = 0.0; /* the quality of this tetrahedron */

		switch (qualmeasure)
		{
			//case QUALMINSINE:
			//	quality = SPR::minsine(vtx1, vtx2, vtx3, vtx4);
			//	break;
			//case QUALMEANSINE:
			//	quality = SPR::meansine(vtx1, vtx2, vtx3, vtx4);
			//	break;
			//case QUALMINSINEANDEDGERATIO:
			//	quality = SPR::minsineandedgeratio(vtx1, vtx2, vtx3, vtx4);
			//	break;
			//case QUALRADIUSRATIO:
			//	quality = SPR::radiusratio(vtx1, vtx2, vtx3, vtx4);
			//	break;
			//case QUALVLRMS3RATIO:
			//	quality = SPR::vlrms3ratio(vtx1, vtx2, vtx3, vtx4);
			//	break;
		case QUALWARPEDMINSINE:
			quality = warpedminsine(vtx1, vtx2, vtx3, vtx4);
			break;
			//case QUALMINANGLE:
			//	quality = SPR::minmaxangle(vtx1, vtx2, vtx3, vtx4, false);
			//	break;
			//case QUALMAXANGLE:
			//	quality = SPR::minmaxangle(vtx1, vtx2, vtx3, vtx4, true);
			//	break;
		default:
			spdlog::info("I don't know what quality measure {} is. Dying...\n", qualmeasure);
			exit(1);
		}

		return quality;
	}

	/* compute the (square) of the minimum sine
	of all the dihedral angles in the tet defined
	by the four vertices (vtx1, vtx2, vtx3, vtx4)
	*/
	REAL warpedminsine(APOINT vtx1, APOINT vtx2, APOINT vtx3, APOINT vtx4)
	{
		REAL *point[4] = { NULL, NULL, NULL, NULL };      /* tet vertices */
		REAL edgelength[3][4]; /* the lengths of each of the edges of the tet */
		REAL facenormal[4][3]; /* the normals of each face of the tet */
		REAL dx, dy, dz;       /* intermediate values of edge lengths */
		REAL facearea2[4];     /* areas of the faces of the tet */
		REAL pyrvolume;        /* volume of tetrahedron */
		REAL sine2, minsine2;  /* the sine (squared) of the dihedral angle */
		int i, j, k, l;          /* loop indices */
		REAL E[3][3] = { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0} };
		const REAL HUGEFLOAT = 1.0e100; /* 这个值要尽量大，否则不保险 */

		/* get tet vertices */
		point[0] = vtx1;
		point[1] = vtx2;
		point[2] = vtx3;
		point[3] = vtx4;

		/* calculate the volume*6 of the tetrahedron */
		pyrvolume = (REAL)GEOM_FUNC::orient3d(point[0], point[1], point[2], point[3]);

		/* if the volume is zero, the quality is zero, no reason to continue */
		if (pyrvolume <= 0.0)
		{
			return 0.0;
		}

		/* for each vertex/face of the tetrahedron */
		for (i = 0; i < 4; i++) {
			j = (i + 1) & 3;
			if ((i & 1) == 0) {
				k = (i + 3) & 3;
				l = (i + 2) & 3;
			}
			else {
				k = (i + 2) & 3;
				l = (i + 3) & 3;
			}

			/* compute the normal for each face */
			facenormal[i][0] =
				(point[k][1] - point[j][1]) * (point[l][2] - point[j][2]) -
				(point[k][2] - point[j][2]) * (point[l][1] - point[j][1]);
			facenormal[i][1] =
				(point[k][2] - point[j][2]) * (point[l][0] - point[j][0]) -
				(point[k][0] - point[j][0]) * (point[l][2] - point[j][2]);
			facenormal[i][2] =
				(point[k][0] - point[j][0]) * (point[l][1] - point[j][1]) -
				(point[k][1] - point[j][1]) * (point[l][0] - point[j][0]);

			/* compute (2 *area)^2 for this face */
			facearea2[i] = facenormal[i][0] * facenormal[i][0] +
				facenormal[i][1] * facenormal[i][1] +
				facenormal[i][2] * facenormal[i][2];

			/* compute edge lengths (squared) */
			for (j = i + 1; j < 4; j++) {
				dx = point[i][0] - point[j][0];
				dy = point[i][1] - point[j][1];
				dz = point[i][2] - point[j][2];
				edgelength[i][j] = dx * dx + dy * dy + dz * dz;
			}
		}

		minsine2 = HUGEFLOAT;     /* start with absurdly big value for sine */

		/* for each edge in the tetrahedron */
		for (i = 0; i < 3; i++) {
			for (j = i + 1; j < 4; j++) {
				k = (i > 0) ? 0 : (j > 1) ? 1 : 2;
				l = 6 - i - j - k;

				/* compute the expression for minimum sine, squared, over 4
				   The reason it's over 4 is because the area values we have
				   are actually twice the area squared */
				   /* if either face area is zero, the sine is zero */
				if (facearea2[k] > 0 && facearea2[l] > 0)
				{
					sine2 = edgelength[i][j] / (facearea2[k] * facearea2[l]);
				}
				else
				{
					/* Encountered zero-area face */
					sine2 = 0.0;
				}

				/* check whether this angle is obtuse */
				if (facenormal[k][0] * facenormal[l][0] + facenormal[k][1] * facenormal[l][1] + facenormal[k][2] * facenormal[l][2] > 0)
				{
					/* if so, warp it down */
					//sine2 = warpsine(sqrt(sine2));
					sine2 = sqrt(sine2)*WRAPPED_SINE_RATIO;
					sine2 *= sine2;
				}

				/* update minimum sine */
				if (sine2 < minsine2)
				{
					minsine2 = sine2;
				}
			}
		}

		return sqrt(minsine2) * pyrvolume;
	}
}
#endif

inline int lnFacInt_own(int line1, int line2, TRI face, APOINT ln1, APOINT ln2, APOINT fac1, APOINT fac2, APOINT fac3, struct polyhedron *Poly)
{
#if 1
	APOINT pnt;
	int val = 0;
	return GEOM_FUNC::lnFacInt(ln1, ln2, fac1, fac2, fac3, pnt, &val) > 1;
#else
	int i;
	double O[5];
	double T[4];
	APOINT P, Q;
	int P_i, Q_i;

	O[0] = isVisible_Search(face[0], face[1], face[2], line1);
	O[1] = isVisible_Search(face[0], face[1], face[2], line2);

	T[0] = O[0] * O[1];

	if (T[0] > SPR_EPS_ZERO)
		return 0;
	else if (T[0] < -1 * SPR_EPS_ZERO)
	{
		O[2] = isVisible_Search(face[0], face[1], line1, line2);
		T[1] = O[0] * O[2];
		if (T[1] > SPR_EPS_ZERO)
			return 0;
		O[3] = isVisible_Search(face[1], face[2], line1, line2);
		T[2] = O[0] * O[3];
		if (T[2] > SPR_EPS_ZERO)
			return 0;
		O[4] = isVisible_Search(face[2], face[0], line1, line2);
		T[3] = O[0] * O[4];
		if (T[3] > SPR_EPS_ZERO)
			return 0;
		//此时可以认定线和面绝对是异面相交

		return 1;
	}
	else	//此时有一个点与ABC是共面，所以需要判断具体是哪个点
	{
		if (O[0] < SPR_EPS_ZERO && O[0] > -1 * SPR_EPS_ZERO)
		{
			for (i = 0; i < 3; i++)
			{
				P[i] = ln1[i];
				Q[i] = ln2[i];
			}
			P_i = line1;
			Q_i = line2;
		}
		else
		{
			for (i = 0; i < 3; i++)
			{
				P[i] = ln2[i];
				Q[i] = ln1[i];
			}
			P_i = line2;
			Q_i = line1;
		}
		O[2] = isVisible_Search(face[0], face[1], P_i, Q_i);
		if (O[2] < -1 * SPR_EPS_ZERO)
			return 0;
		O[3] = isVisible_Search(face[1], face[2], P_i, Q_i);
		if (O[3] < -1 * SPR_EPS_ZERO)
			return 0;
		O[4] = isVisible_Search(face[2], face[0], P_i, Q_i);
		if (O[4] < -1 * SPR_EPS_ZERO)
			return 0;

		return 1;
	}
#endif
}


/*这里面的ln1和ln2是线的两个端点，fac1，fac2，fac3是那个面的三个坐标点，d是这条线的长度，areaa是这个面的面积，
通过参量可以返回的值是这个pnt，返回的是如果线和面相交，那么交点的情况。val是相交的具体情况
w1，w2，w3具体是指什么？
在函数中，当返回一个NOD（1）的时候，表明这条线和这个面相交的情况是这条线穿过了面的端点
当返回一个EDG（2）的时候，表明这条线和这个面相交的情况是这条线穿过了面的边界线
当返回一个FAC（3）的时候，表明这条线和这个面直接相交
当返回-1的时候，表明这条线和这个面不相交
当返回-2的时候，表明程序出现严重错误，但是具体还是不大懂*/

inline int lnFacInt(APOINT ln1, APOINT ln2,
	APOINT fac1, APOINT fac2, APOINT fac3,
	double d, double areaa,
	APOINT pnt, int *val,
	double *w1, double *w2, double *w3)
{
	double xa, ya, za, xb, yb, zb;
	double x1, y1, z1, x2, y2, z2, x3, y3, z3;
	double xi, yi, zi;
	double a1, b1, c1, d1;
	double xconst, yconst, zconst;
	double e1, f1, g1, h1;
	double det, det1, det2, det3;
	double area1, area2, area3, suma, diff;
	double dist1, dist2;
	double area1n, area2n, area3n;
	int ntype = -1;

	xa = ln1[0];
	ya = ln1[1];
	za = ln1[2];
	xb = ln2[0];
	yb = ln2[1];
	zb = ln2[2];

	x1 = fac1[0];
	y1 = fac1[1];
	z1 = fac1[2];
	x2 = fac2[0];
	y2 = fac2[1];
	z2 = fac2[2];
	x3 = fac3[0];
	y3 = fac3[1];
	z3 = fac3[2];
	a1 = (y2 - y1)*(z3 - z1) - (z2 - z1)*(y3 - y1);
	b1 = -((x2 - x1)*(z3 - z1) - (z2 - z1)*(x3 - x1));
	c1 = (x2 - x1)*(y3 - y1) - (y2 - y1)*(x3 - x1);
	d1 = a1 * x1 + b1 * y1 + c1 * z1;

	xconst = fabs(xb - xa) / d;
	yconst = fabs(yb - ya) / d;
	zconst = fabs(zb - za) / d;

	if (xconst < SPR_EPS_ZERO && yconst < SPR_EPS_ZERO &&
		zconst < SPR_EPS_ZERO)
	{
		spdlog::info("A line is degenerated into a APOINT\n");
		return -2;
	}

	if (xconst < SPR_EPS_ZERO && yconst < SPR_EPS_ZERO)
	{
		if (fabs(c1) < SPR_EPS_ZERO_SQ)
			return -1;

		xi = xa;
		yi = ya;
		zi = (d1 - a1 * xa - b1 * ya) / c1;
	}
	else if (yconst < SPR_EPS_ZERO && zconst < SPR_EPS_ZERO)
	{
		if (fabs(a1) < SPR_EPS_ZERO_SQ)
			return -1;

		xi = (d1 - b1 * ya - c1 * za) / a1;
		yi = ya;
		zi = za;
	}
	else if (zconst < SPR_EPS_ZERO && xconst < SPR_EPS_ZERO)
	{
		if (fabs(b1) < SPR_EPS_ZERO_SQ)
			return -1;

		xi = xa;
		yi = (d1 - a1 * xa - c1 * za) / b1;
		zi = za;
	}
	else if (xconst < SPR_EPS_ZERO)
	{
		e1 = -(yb - ya) / (zb - za);
		f1 = ya + e1 * za;
		det = b1 * e1 - c1;

		if (fabs(det / areaa) < SPR_EPS_ZERO)
			return -1;

		det1 = (d1 - a1 * xa)*e1 - c1 * f1;
		det2 = b1 * f1 - (d1 - a1 * xa);

		xi = xa;
		yi = det1 / det;
		zi = det2 / det;
	}
	else if (yconst < SPR_EPS_ZERO)
	{
		g1 = -(xb - xa) / (zb - za);
		h1 = xa + g1 * za;
		det = a1 * g1 - c1;

		if (fabs(det / areaa) < SPR_EPS_ZERO)
			return -1;

		det1 = (d1 - b1 * ya)*g1 - h1 * c1;
		det2 = a1 * h1 - (d1 - b1 * ya);
		xi = det1 / det;
		yi = ya;
		zi = det2 / det;
	}
	else if (zconst < SPR_EPS_ZERO)
	{
		e1 = -(xb - xa) / (yb - ya);
		f1 = xa + e1 * ya;
		det = a1 * e1 - b1;

		if (fabs(det / areaa) < SPR_EPS_ZERO)
			return -1;

		det1 = (d1 - c1 * za)*e1 - b1 * f1;
		det2 = a1 * f1 - (d1 - c1 * za);

		xi = det1 / det;
		yi = det2 / det;
		zi = za;
	}
	else
	{
		e1 = -(xb - xa) / (yb - ya);
		f1 = xa + e1 * ya;
		g1 = -(xb - xa) / (zb - za);
		h1 = xa + g1 * za;

		det = a1 * e1*g1 - b1 * g1 - c1 * e1;
		det1 = d1 * e1*g1 - b1 * f1*g1 - c1 * e1*h1;
		det2 = a1 * f1*g1 - d1 * g1 + c1 * (h1 - f1);
		det3 = a1 * e1*h1 - b1 * (h1 - f1) - d1 * e1;

		if (fabs(det / areaa) < SPR_EPS_ZERO)
			return -1;

		xi = det1 / det;
		yi = det2 / det;
		zi = det3 / det;
	}

	area1 = areas(xi, yi, zi, x1, y1, z1, x2, y2, z2);
	area2 = areas(xi, yi, zi, x2, y2, z2, x3, y3, z3);
	area3 = areas(xi, yi, zi, x3, y3, z3, x1, y1, z1);
	suma = area1 + area2 + area3;
	diff = fabs(suma - areaa);
	if (fabs(diff / areaa) > SPR_EPS_ZERO)
		return -1;

	dist1 = sqrt((xi - xa) * (xi - xa) +
		(yi - ya) * (yi - ya) +
		(zi - za) * (zi - za));
	dist2 = sqrt((xi - xb) * (xi - xb) +
		(yi - yb) * (yi - yb) +
		(zi - zb) * (zi - zb));
	if (fabs((dist1 + dist2 - d) / d) >= SPR_EPS_ZERO)
	{
		return -1;
	}


	pnt[0] = xi;
	pnt[1] = yi;
	pnt[2] = zi;

	area1n = area1 / areaa;
	area2n = area2 / areaa;
	area3n = area3 / areaa;

	if (fabs(area1n) < SPR_EPS_ZERO &&
		fabs(area2n) < SPR_EPS_ZERO &&
		fabs(area3n) < SPR_EPS_ZERO)
	{
		spdlog::info("Face area is near zero!\n");
		return -2;
	}

	if (fabs(area1n) < SPR_EPS_ZERO &&
		fabs(area2n) < SPR_EPS_ZERO)
	{
		ntype = NOD;
		*val = 1;
	}
	else if (fabs(area2n) < SPR_EPS_ZERO &&
		fabs(area3n) < SPR_EPS_ZERO)
	{
		ntype = NOD;
		*val = 2;
	}
	else if (fabs(area3n) < SPR_EPS_ZERO &&
		fabs(area1n) < SPR_EPS_ZERO)
	{
		ntype = NOD;
		*val = 0;
	}
	else if (fabs(area1n) < SPR_EPS_ZERO)
	{
		ntype = EDG;
		*val = 0;
	}
	else if (fabs(area2n) < SPR_EPS_ZERO)
	{
		ntype = EDG;
		*val = 1;
	}
	else if (fabs(area3n) < SPR_EPS_ZERO)
	{
		ntype = EDG;
		*val = 2;
	}
	else
	{
		ntype = FAC;
	}

	if (w1)
		*w1 = area1n;
	if (w2)
		*w2 = area2n;
	if (w3)
		*w3 = area3n;

	return ntype;
}

int SPRImpl::saveConstraints(const char *fname, Constraints *csts, double iniQ)
{
	FILE *fp = NULL;
	int i, size;

	if (!(fp = fopen(fname, "w")))
	{
		spdlog::info("Cannot open file: %s.\n", fname);
		return -1;
	}

	/* 单元的初始质量值 */
	fprintf(fp, "%24.20g\n", iniQ);

	/* 输出线约束 */
	size = csts ? csts->lineCst_size : 0;
	fprintf(fp, "%d\n", size);
	for (i = 0; i < size; i++)
		fprintf(fp, "%d %d\n", csts->lineCsts[i][0], csts->lineCsts[i][1]);
	size = csts ? csts->colineN_size : 0;
	fprintf(fp, "%d\n", size);
	for (i = 0; i < size; i++)
		fprintf(fp, "%d\n", csts->colineNodes[i]);

	/* 输出面约束 */
	size = csts ? csts->faceCst_size : 0;
	fprintf(fp, "%d\n", size);
	for (i = 0; i < size; i++)
		fprintf(fp, "%d %d %d\n", csts->faceCsts[i][0], csts->faceCsts[i][1], csts->faceCsts[i][2]);
	size = csts ? csts->cofaceN_size : 0;
	fprintf(fp, "%d\n", size);
	for (i = 0; i < size; i++)
		fprintf(fp, "%d\n", csts->cofaceNodes[i]);

	/* 输出体约束 */
	size = csts ? csts->badVolmCst_size : 0;
	fprintf(fp, "%d\n", size);
	for (i = 0; i < size; i++)
		fprintf(fp, "%d %d %d %d\n", csts->badVolmCsts[i][0], csts->badVolmCsts[i][1],
			csts->badVolmCsts[i][2], csts->badVolmCsts[i][3]);

	fclose(fp);
}

int SPRImpl::readConstraints(const char *fname, Constraints *csts, double *iniQ)
{
	FILE *fp = NULL;
	int i, size;

	memset(csts, 0, sizeof(Constraints));
	if (!(fp = fopen(fname, "r")))
	{
		spdlog::info("Cannot open file: %s.\n", fname);
		return -1;
	}

	/* 单元的初始质量值 */
	fscanf(fp, "%lf", iniQ);

	/* 输出线约束 */
	fscanf(fp, "%d", &size);
	csts->lineCst_size = size;
	for (i = 0; i < size; i++)
		fscanf(fp, "%d%d", &csts->lineCsts[i][0], &csts->lineCsts[i][1]);
	fscanf(fp, "%d", &size);
	csts->colineN_size = size;
	for (i = 0; i < size; i++)
		fscanf(fp, "%d", &csts->colineNodes[i]);

	/* 输出面约束 */
	fscanf(fp, "%d\n", &size);
	csts->faceCst_size = size;
	for (i = 0; i < size; i++)
		fscanf(fp, "%d %d %d\n", &csts->faceCsts[i][0], &csts->faceCsts[i][1], &csts->faceCsts[i][2]);
	fscanf(fp, "%d\n", &size);
	csts->cofaceN_size = size;
	for (i = 0; i < size; i++)
		fscanf(fp, "%d\n", &csts->cofaceNodes[i]);

	fclose(fp);
}

//这个代码是把输入的数据翻转了一下再来计算的
int SPRImpl::openpls(char *address, struct polyhedron *P, bool reverse)
{
	FILE *pf;
	int ne, np;
	int i, j, index;
	APOINT v[MAX_SPR_POLY_SIZE];
	TRI t[MAX_SPR_POLY_SIZE];
	int surface[MAX_SPR_POLY_SIZE];
	if (!(pf = fopen(address, "r")))
	{
		spdlog::info("Can't open file: %s!\n", address);
		return 0;
	}
	fscanf(pf, "\t%d\t%d\t0\t0\t0\t0\n", &ne, &np);
	P->num_t = ne;
	num_vertices = np;

	for (i = 1; i <= np; i++)
		fscanf(pf, "\t%d\t%lf\t%lf\t%lf\n", &index, &v[i - 1][0], &v[i - 1][1], &v[i - 1][2]);
	for (i = 0; i < num_vertices; i++)
	{
		for (j = 0; j < 3; j++)
			vertices[i][j] = v[i][j];
	}

	if (reverse)
	{
		for (i = 1; i <= ne; i++)
			fscanf(pf, "\t%d\t%d\t%d\t%d\t%d\n", &index, &t[i - 1][0], &t[i - 1][2], &t[i - 1][1], &surface[i - 1]);
	}
	else
	{
		for (i = 1; i <= ne; i++)
			fscanf(pf, "\t%d\t%d\t%d\t%d\t%d\n", &index, &t[i - 1][0], &t[i - 1][1], &t[i - 1][2], &surface[i - 1]);
	}

	rearrage(t, P->num_t);

	for (i = 0; i < P->num_t; i++)
	{
		for (j = 0; j < 3; j++)
		{
			P->t[i][j] = t[i][j] - 1;
			P->used_v[P->t[i][j]]++;
		}
		P->surface[i] = surface[i];

		for (j = 0; j < 3; j++)
		{
			index = P->t[i][j];
			P->vprt[index] = (i << 2) | j;
		}
	}

	for (i = 0; i < num_vertices; i++)
	{
		if (P->used_v[i] != 0)
			P->num_vable++;
	}
	fclose(pf);
	return 1;
}

int SPRImpl::savepls(char *address, struct polyhedron *P, bool bAllNode)
{
	FILE *pf;
	int i, j;
	int g2l[MAX_SPR_POLY_SIZE];

	if (!(pf = fopen(address, "w")))
	{
		spdlog::info("Can't open file: %s!\n", address);
		return 0;
	}

	fprintf(pf, "\t%d\t%d\t0\t0\t0\t0\n", P->num_t, bAllNode ? num_vertices : P->num_vable);

	if (bAllNode)
	{
		for (i = 0; i < num_vertices; i++)
		{
			//fprintf(pf,"\t%d\t%lf\t%lf\t%lf\n", i+1, vertices[i][0], vertices[i][1], vertices[i][2]);
			fprintf(pf, "\t%d\t%14.14lf\t%14.14lf\t%14.14lf\n", spr_l2g[i] - 7, vertices[i][0], vertices[i][1], vertices[i][2]);
		}

		for (i = 0; i < P->num_t; i++)
			fprintf(pf, "\t%d\t%d\t%d\t%d\t%d\n", i + 1,
				P->t[i][0] + 1,
				P->t[i][2] + 1,
				P->t[i][1] + 1,
				P->surface[i]);
	}
	else
	{
		for (i = 0, j = 0; i < num_vertices; i++)
		{
			if (P->used_v[i] > 0)
			{
				g2l[i] = j++;
				fprintf(pf, "\t%d\t%14.14lf\t%14.14lf\t%14.14lf\n", j, vertices[i][0], vertices[i][1], vertices[i][2]);
			}
		}

		for (i = 0; i < P->num_t; i++)
			fprintf(pf, "\t%d\t%d\t%d\t%d\t%d\n", i + 1,
				g2l[P->t[i][0]] + 1,
				g2l[P->t[i][2]] + 1,
				g2l[P->t[i][1]] + 1,
				P->surface[i]);
	}
	fclose(pf);
	return 1;
}

int SPRImpl::saveneg(char *address, struct polyhedron *P)
{
	FILE *pf;
	int i;

	if (!(pf = fopen(address, "w")))
	{
		spdlog::info("Can't open file: %s!\n", address);
		return 0;
	}

	fprintf(pf, "\t%d\n", P->num_t);

	for (i = 0; i < P->num_t; i++)
		fprintf(pf, "\t%d\t%d\t%d\t%d\n", i,
			P->n[i][0],
			P->n[i][1],
			P->n[i][2]);

	fclose(pf);
	return 1;
}

void SPRImpl::rearrage(TRI t[], int n)
{
	int i, j;
	int tmp[3];
	for (i = 0; i < n; i++)
	{
		if (t[i][0] < t[i][1] && t[i][0] < t[i][2])
			;
		else if (t[i][1] < t[i][0] && t[i][1] < t[i][2])
		{
			for (j = 0; j < 3; j++)
				tmp[j] = t[i][j];
			t[i][0] = tmp[1];
			t[i][1] = tmp[2];
			t[i][2] = tmp[0];
		}
		else if (t[i][2] < t[i][0] && t[i][2] < t[i][1])
		{
			for (j = 0; j < 3; j++)
				tmp[j] = t[i][j];
			t[i][0] = tmp[2];
			t[i][1] = tmp[0];
			t[i][2] = tmp[1];
		}
	}
}


int SPRImpl::savepl3(char *address, struct polyhedron *P, struct divide *T)
{
	int i;
	//	int parent[MAX_SPR_POLY_SIZE];
	FILE *pf;

	if (!(pf = fopen(address, "wb")))
	{
		spdlog::info("Can't open file: %s!", address);
		return 0;
	}

	//	findparent(parent, P, T);
	fprintf(pf, "\t%d\t%d\t%d\n", T->num, num_vertices, P->num_t);
	for (i = 0; i < num_vertices; i++)
		fprintf(pf, "%d\t%lf\t%lf\t%lf\n", i + 1, vertices[i][0], vertices[i][1], vertices[i][2]);
	for (i = 0; i < T->num; i++)
		fprintf(pf, "%d\t%d\t%d\t%d\t%d\n", i + 1, T->te[i][0] + 1, T->te[i][1] + 1, T->te[i][2] + 1, T->te[i][3] + 1);
	for (i = 0; i < P->num_t; i++)
		fprintf(pf, "%d\t%d\t%d\t%d\t%d\t%d\n", i + 1, P->t[i][0] + 1, P->t[i][1] + 1, P->t[i][2] + 1,
			0, 0);//parent[i], P->surface[i]);
	fclose(pf);
	return 1;
}

void SPRImpl::findparent(int parent[], struct polyhedron *P, struct divide *T)
{
	int i, j, k, m;
	int flag;
	TRI face[6];
	TRI tetraface[MAX_SPR_POLY_SIZE * 4];
	k = 0;
	flag = 0;
	for (i = 0; i < T->num; i++)
	{
		tetraface[k][0] = T->te[i][0];
		tetraface[k][1] = T->te[i][1];
		tetraface[k][2] = T->te[i][2];
		k++;

		tetraface[k][0] = T->te[i][0];
		tetraface[k][1] = T->te[i][1];
		tetraface[k][2] = T->te[i][3];
		k++;

		tetraface[k][0] = T->te[i][0];
		tetraface[k][1] = T->te[i][2];
		tetraface[k][2] = T->te[i][3];
		k++;

		tetraface[k][0] = T->te[i][1];
		tetraface[k][1] = T->te[i][2];
		tetraface[k][2] = T->te[i][3];
		k++;
	}
	//这样k便可以代表tetraface的数目
	for (i = 0; i < P->num_t; i++)
	{
		face[0][0] = P->t[i][0];
		face[0][1] = P->t[i][1];
		face[0][2] = P->t[i][2];
		face[1][0] = P->t[i][0];
		face[1][1] = P->t[i][2];
		face[1][2] = P->t[i][1];
		face[2][0] = P->t[i][1];
		face[2][1] = P->t[i][0];
		face[2][2] = P->t[i][2];
		face[3][0] = P->t[i][1];
		face[3][1] = P->t[i][2];
		face[3][2] = P->t[i][0];
		face[4][0] = P->t[i][2];
		face[4][1] = P->t[i][0];
		face[4][2] = P->t[i][1];
		face[5][0] = P->t[i][2];
		face[5][1] = P->t[i][1];
		face[5][2] = P->t[i][0];
		for (j = 0; j < k; j++)
		{
			for (m = 0; m < 6; m++)
			{
				if (face[m][0] == tetraface[j][0] && face[m][1] == tetraface[j][1] && face[m][2] == tetraface[j][2])
				{
					parent[i] = (int)(j / 4) + 1;
					flag = 1;
				}
			}
			if (flag == 1)
			{
				flag = 0;
				break;
			}
		}
	}
}

#if 0
inline void sort_quad_pair(int i1, int i2, int i3, int i4,
	int *j1, int *j2, int *j3, int *j4)
{
	int max1, max2, min1, min2;

	max1 = MAX_VALUE(i1, i2);
	max2 = MAX_VALUE(i3, i4);
	min1 = MIN_VALUE(i1, i2);
	min2 = MIN_VALUE(i3, i4);

	*j1 = MIN_VALUE(min1, min2);
	*j4 = MAX_VALUE(max1, max2);
	*j3 = MIN_VALUE(max1, max2);
	*j2 = MAX_VALUE(min1, min2);
}
#endif

#if 1
inline int sort_quad_pair(int *i1, int *i2, int *i3, int *i4)
{
	int t;
	int cnt = 0;
	if (*i1 > *i2)
	{
		SWAP_VALUE(*i1, *i2, t)
			cnt++;
	}
	if (*i3 > *i4)
	{
		SWAP_VALUE(*i3, *i4, t)
			cnt++;
	}
	if (*i2 > *i4)
	{
		SWAP_VALUE(*i2, *i4, t)
			cnt++;
	}
	if (*i1 > *i3)
	{
		SWAP_VALUE(*i1, *i3, t)
			cnt++;
	}
	if (*i2 > *i3)
	{
		SWAP_VALUE(*i2, *i3, t)
			cnt++;
	}

	return cnt % 2 ? -1 : 1;
}
#endif

inline int make_quad_pair(int j1, int j2, int j3, int j4)
{
	return j4 << 24 | j3 << 16 | j2 << 8 | j1;
}

/* ----------------------------------------------------------------------
 * 建立复杂多面体的相邻关系，如果GridEdgeHash存在复边，假定复边对应的
 * 单元是有序的
 * -------------------------------------------------------------------*/
int build_poly_neig(struct polyhedron *poly, GridEdgeHash *eHash)
{
	int i, j, k, i1, i2, t, neig, idx;
	std::vector<int> vecNeigElems;
	const int tri_en[3][2] = {
		{1, 2}, {2, 0}, {0, 1}
	};

	memset(poly->n, -1, sizeof(NEG)*MAX_SPR_POLY_SIZE);

	for (i = 0; i < poly->num_t; i++)
	{
		for (j = 0; j < 3; j++)
		{
			if (poly->n[i][j] >= 0)
				continue;

			i1 = poly->t[i][tri_en[j][0]];
			i2 = poly->t[i][tri_en[j][1]];
			if (i1 > i2)
			{
				t = i1;
				i1 = i2;
				i2 = t;
			}

			vecNeigElems.clear();
			eHash->find_edge_elems(i1, i2, vecNeigElems);

			assert(vecNeigElems.size() > 0 && vecNeigElems.size() % 2 == 0);
			for (k = 0; k < vecNeigElems.size(); k++)
			{
				if (vecNeigElems[k] == i)
				{
					idx = k % 2 ? k - 1 : k + 1;
					neig = vecNeigElems[idx];
					break;
				}
			}
			assert(k < vecNeigElems.size());

			poly->n[i][j] = neig;
			for (k = 0; k < 3; k++)
			{
				if (poly->t[neig][k] != i1 &&
					poly->t[neig][k] != i2)
				{
					poly->n[neig][k] = i;
					break;
				}
			}
		}
	}

	return 1;
}

int build_manifold_poly_neig(struct polyhedron *poly)
{
	int i, j, k, i1, i2, t, neig, idx, key, value, count;
	std::vector<int> vecNeigElems;
	const int tri_en[3][2] = {
		{1, 2}, {2, 0}, {0, 1}
	};
	IntIntMap edgeMap(poly->num_t * 3 / 2); /* 边的数目不是很多 */

	memset(poly->n, -1, sizeof(NEG)*MAX_SPR_POLY_SIZE);

	for (i = 0; i < poly->num_t; i++)
	{
		for (j = 0; j < 3; j++)
		{
			if (poly->n[i][j] >= 0)
				continue;

			i1 = poly->t[i][tri_en[j][0]];
			i2 = poly->t[i][tri_en[j][1]];
			if (i1 > i2)
			{
				t = i1;
				i1 = i2;
				i2 = t;
			}
			/* i1/i2不是很大，可将小的压到后16位，大的压到前16位 */
			key = (i1 << 16) | i2;
			if (edgeMap.get(key, &value))
			{/* 前面已经碰到这条边 */
				count = value & 0x3;
				if (count > 1)
				{/* 非二边流形，返回0 */
					return 0;
				}
				assert(count == 1);
				neig = value >> 4; /* 右移4位 */
				idx = (value & 0xc) >> 2;
				poly->n[i][j] = neig;
				poly->n[neig][idx] = i;

				value++; /* 多一次引用 */
				edgeMap.put(key, value);
			}
			else
			{
				value = (i << 4) | (j << 2) | 0x1; /* count = 1 */
				edgeMap.put(key, value);
			}
		}
	}

	return 1;
}
/* ----------------------------------------------------------------------
 * 分解复杂多面体为多个简单多面体
 * 基于相邻关系的着色算法
 * -------------------------------------------------------------------*/
#ifdef _SPR_USING_STD_VECTOR
int SPRImpl::divide_poly(struct polyhedron *cmpx_poly,
	std::vector<struct polyhedron*>& vec_simp_poly)
#else
int SPRImpl::divide_poly(struct polyhedron *cmpx_poly,
	polyhedron *simpPolyArray[], int *simpPolyNum)
#endif /* _SPR_USING_STD_VECTOR */
{
	int i, j, k, iElem, iSrch, neig, nidx;
	int m, newf_l, newf_g;
	//	std::stack<int> stkTreeSear;

	char test_flag[MAX_SPR_POLY_SIZE], add_fg[MAX_SPR_POLY_SIZE];
	struct polyhedron *simp_poly = NULL;

#ifndef _SPR_USING_STD_VECTOR 
	*simpPolyNum = 0; /* initialise */
#endif

#ifdef _SPR_ERROR_CHK
//	savepls("fail.pls", cmpx_poly);
//	saveneg("fail.neg", cmpx_poly);
	assert(assert_poly_neig(cmpx_poly));
#endif /* _SPR_ERROR_CHK */

	assert(dp_stkTreeSear.empty());

	memset(test_flag, 0, sizeof(char)*cmpx_poly->num_t);
	memset(dp_g2l, -1, sizeof(int)*cmpx_poly->num_t);
	memset(add_fg, 0, sizeof(char)*cmpx_poly->num_t);

	for (i = 0; i < cmpx_poly->num_t; i++)
	{
		if (test_flag[i] != 0)
			continue;

		simp_poly = (struct polyhedron*)malloc(sizeof(struct polyhedron));
		memset(simp_poly, 0, sizeof(struct polyhedron));

		memcpy(simp_poly->t[simp_poly->num_t], cmpx_poly->t[i], sizeof(TRI));
		memcpy(simp_poly->n[simp_poly->num_t], cmpx_poly->n[i], sizeof(NEG));
		dp_g2l[i] = simp_poly->num_t;
		simp_poly->num_t++;
		test_flag[i] = 1;
		dp_stkTreeSear.push(i);

		/* find a seed, perform tree search */
		while (!dp_stkTreeSear.empty())
		{
			iElem = dp_stkTreeSear.top();
			dp_stkTreeSear.pop();
			for (m = 0; m <= 2; m++)
			{
				iSrch = cmpx_poly->n[iElem][m];
				if (iSrch >= 0 && test_flag[iSrch] == 0)
				{
					memcpy(simp_poly->t[simp_poly->num_t], cmpx_poly->t[iSrch], sizeof(TRI));
					memcpy(simp_poly->n[simp_poly->num_t], cmpx_poly->n[iSrch], sizeof(NEG));
					dp_g2l[iSrch] = simp_poly->num_t;
					simp_poly->num_t++;
					test_flag[iSrch] = 1;
					dp_stkTreeSear.push(iSrch);
				}
			}
		}

		/* ----------------------------------
		 * 记录新增加的面
		 * --------------------------------*/
		for (j = 0; j < cmpx_poly->newf_ct; j++)
		{
			newf_g = cmpx_poly->newf_id[j];
			if (add_fg[newf_g] == 0)
			{
				newf_l = dp_g2l[newf_g];
				if (newf_l >= 0)
				{
					simp_poly->newf_id[simp_poly->newf_ct++] = newf_l;
					simp_poly->newf_fg[newf_l] = 1;

					add_fg[newf_g] = 1;
				}
			}
		}

#ifdef _SPR_USING_STD_VECTOR
		vec_simp_poly.push_back(simp_poly);
#else
		if (*simpPolyNum >= MAX_SIMP_POLY_NUM)
		{
			spdlog::info("SPR Error: Increase MAX_SIMP_POLY_NUM({}).\n", MAX_SIMP_POLY_NUM);
			exit(1);
		}
		simpPolyArray[(*simpPolyNum)++] = simp_poly;
#endif
	}

#ifdef _SPR_USING_STD_VECTOR
	for (i = 0; i < vec_simp_poly.size(); i++)
	{
		simp_poly = vec_simp_poly[i];
#else
	for (i = 0; i < *simpPolyNum; i++)
	{
		simp_poly = simpPolyArray[i];
#endif
		//		memset(simp_poly->used_v, 0, sizeof(int)*MAX_SPR_POLY_SIZE);
		//		simp_poly->num_vable = 0;
		for (j = 0; j < simp_poly->num_t; j++)
		{
			for (k = 0; k < 3; k++)
			{
				neig = simp_poly->n[j][k];
				if (neig >= 0)
					simp_poly->n[j][k] = dp_g2l[neig];

				nidx = simp_poly->t[j][k];
				if (simp_poly->used_v[nidx] == 0)
					simp_poly->num_vable++;
				simp_poly->used_v[nidx]++;
				simp_poly->vprt[nidx] = (j << 2) | k;
			}
		}
	}

	return 1;
	}

/* ----------------------------------------------------------------------
 * 获得首个包含某节点的父亲单元
 * -------------------------------------------------------------------*/
inline void get_node_first_parent(struct polyhedron *poly,
	int nidx, int *tidx, int *tcod)
{
	*tidx = poly->vprt[nidx] >> 2;
	*tcod = poly->vprt[nidx] & 0x3;
}

inline int node_face_code(struct polyhedron *poly, int nidx, int fidx)
{
	int i;
	for (i = 0; i < 3; i++)
	{
		if (poly->t[fidx][i] == nidx)
			return i;
	}

	return -1;
}

inline int edge_face_code(struct polyhedron *poly,
	int idx1, int idx2, int fidx,
	int *ecod, int *eort)
{
	int i, i1, i2;
	const int tri_en[3][2] = {
		{1, 2}, {2, 0}, {0, 1}
	};

	*ecod = -1;
	for (i = 0; i < 3; i++)
	{
		i1 = poly->t[fidx][tri_en[i][0]];
		i2 = poly->t[fidx][tri_en[i][1]];
		if (i1 == idx1 && i2 == idx2)
		{
			*ecod = i;
			*eort = 1;
			break;
		}
		else if (i2 == idx1 && i1 == idx2)
		{
			*ecod = i;
			*eort = 0;
			break;
		}
	}

	return *ecod;
}

inline int add_face_edge_pair(struct polyhedron *poly,
	int idx1, int idx2, int fidx,
	std::vector<int>& vecp)
{
	int eprt, m;
	int ecod, eort;

	if (edge_face_code(poly, idx1, idx2, fidx, &ecod, &eort) >= 0)
	{
		eprt = (fidx << 3) & (ecod << 1) & eort;

		vecp.push_back(eprt);

		fidx = poly->n[fidx][ecod];
		for (m = 0; m <= 2; m++)
		{
			if (poly->t[fidx][m] != idx1 &&
				poly->t[fidx][m] != idx2)
				break;
		}
		if (m <= 2)
		{
			ecod = m;
			eort = !eort;

			eprt = (fidx << 3) & (ecod << 1) & eort;

			vecp.push_back(eprt);
		}

		return vecp.size();
	}

	return 0;
}
/* ----------------------------------------------------------------------
 * 获得所有包含某节点的父亲单元
 * -------------------------------------------------------------------*/
 // inline void get_node_parents(struct polyhedron *poly, int nidx, 
 // 				std::vector<int>& vecp)
 // {
 // 	int vprt, fidx, fcod, sidx, scod;
 // 	int m;
 // 	std::stack<int> stkTreeSear;
 // 	int test_flag[MAX_SPR_POLY_SIZE];
 // 
 // 	memset(test_flag, 0, sizeof(int)*MAX_SPR_POLY_SIZE);
 // 
 // 	vprt = poly->vprt[nidx];
 // 	fidx = vprt >> 2;
 // 	fcod = vprt & 0x3;
 // 	vecp.push_back(fidx);
 // 	stkTreeSear.push(vprt);
 // 	test_flag[fidx] = 1;
 // 
 // 	while (!stkTreeSear.empty())
 // 	{
 // 		vprt = stkTreeSear.top();
 // 		stkTreeSear.pop();
 // 		fidx = vprt >> 2;
 // 		fcod = vprt & 0x3;
 // 
 // 		for (m = 0; m <= 2; m++)
 // 		{
 // 			if (m == fcod)
 // 				continue;
 // 
 // 			sidx = poly->n[fidx][m];
 // 			if (sidx >= 0 && test_flag[sidx] == 0)
 // 			{
 // 				test_flag[sidx] = 1;
 // 				if ((scod = node_face_code(poly, nidx, sidx)) >= 0)
 // 				{
 // 					vprt = (sidx << 2) | scod;
 // 					stkTreeSear.push(vprt);
 // 					vecp.push_back(sidx);
 // 				}
 // 			}
 // 		}
 // 	}
 // 
 // #ifdef _DEBUG
 // 	if (vecp.size() != poly->used_v[nidx])
 // 	{
 // 		m = vecp.size();
 // 	}
 // 	assert(vecp.size() == poly->used_v[nidx]);
 // #endif
 // 
 //}
inline void get_node_parents(struct polyhedron *poly, int nidx, int *vecp, int *num)
{
	int i, j, lmt = poly->used_v[nidx];

	for (i = 0, j = 0; i < poly->num_t && j < lmt; i++)
	{
		if (poly->t[i][0] == nidx ||
			poly->t[i][1] == nidx ||
			poly->t[i][2] == nidx)
		{
			if (j >= *num)
			{
				spdlog::info("error in get_node_parents(..): Input array is too small.\n");
				exit(1);
			}
			vecp[j++] = i;
		}
	}
	*num = j;
#ifdef _DEBUG
	assert(*num == poly->used_v[nidx]);
#endif
}

#if 1
inline void get_node_parents(struct polyhedron *poly, int nidx,
	std::vector<int>& vecp)
{
	int i, j, lmt = poly->used_v[nidx];

	for (i = 0, j = 0; i < poly->num_t && j < lmt; i++)
	{
		if (poly->t[i][0] == nidx ||
			poly->t[i][1] == nidx ||
			poly->t[i][2] == nidx)
		{
			++j;
			vecp.push_back(i);
		}
	}

#ifdef _DEBUG
	if (vecp.size() != poly->used_v[nidx])
	{
		i = vecp.size();
	}
	assert(vecp.size() == poly->used_v[nidx]);
#endif

}
#else
inline void get_node_parents(struct polyhedron *poly, int nidx,
	std::vector<int>& vecp)
{
	int vprt, fidx, fcod, sidx, scod;
	int i, m;
	std::stack<int> stkTreeSear;
	int test_flag[MAX_SPR_POLY_SIZE];

	memset(test_flag, 0, sizeof(int)*MAX_SPR_POLY_SIZE);

	vprt = poly->vprt[nidx];
	fidx = vprt >> 2;
	fcod = vprt & 0x3;
	vecp.push_back(fidx);
	//	stkTreeSear.push(vprt);
	test_flag[fidx] = 1;

	do
	{
		for (m = 0; m <= 2; m++)
		{
			if (m == fcod)
				continue;

			sidx = poly->n[fidx][m];
			if (sidx >= 0 && test_flag[sidx] == 0)
			{
				test_flag[sidx] = 1;
				if ((scod = node_face_code(poly, nidx, sidx)) >= 0)
				{
					vprt = (sidx << 2) | scod;
					stkTreeSear.push(vprt);
					vecp.push_back(sidx);
					if (vecp.size() == poly->used_v[nidx])
						return;
				}
			}
		}

		if (stkTreeSear.empty())
			break;
		else
		{
			vprt = stkTreeSear.top();
			stkTreeSear.pop();
			fidx = vprt >> 2;
			fcod = vprt & 0x3;
		}

	} while (true);

	assert(vecp.size() <= poly->used_v[nidx]);

	if (vecp.size() != poly->used_v[nidx])
	{
		for (i = 0; i < poly->num_t; i++)
		{
			if (test_flag[i] == 0 &&
				(poly->t[i][0] == nidx ||
					poly->t[i][1] == nidx ||
					poly->t[i][2] == nidx))
				vecp.push_back(i);
		}
	}

#ifdef _DEBUG
	if (vecp.size() != poly->used_v[nidx])
	{
		m = vecp.size();
	}
	assert(vecp.size() == poly->used_v[nidx]);
#endif

}
#endif

inline int get_edge_parents(int i1, int i2, std::vector<int>& vecp,
	polyhedron *poly, int *face1, int *face2,
	int *code1, int *code2)
{
	int i, m, idx, num = 0, j1, j2;
	const int tri_en[3][2] = {
		{1, 2}, {2, 0}, {0, 1}
	};

	for (i = 0; i < vecp.size() && num < 2; i++)
	{
		idx = vecp[i];

		for (m = 0; m < 3 && num < 2; m++)
		{
			j1 = poly->t[idx][tri_en[m][0]];
			j2 = poly->t[idx][tri_en[m][1]];

			if (j1 == i1 && j2 == i2)
			{
				*face1 = idx;
				*code1 = m;
				num++;
			}
			else if (j1 == i2 && j2 == i1)
			{
				*face2 = idx;
				*code2 = m;
				num++;
			}
		}
	}

	assert(num == 0 || num == 2);
	return num;
}


inline int get_edge_parents(int i1, int i2, std::vector<int>& vecp,
	polyhedron *poly, int *face1, int *face2,
	int *code1, int *code2, int *nume)
{
	int i, m, k, idx, neg, num = 0, j1, j2;
	const int tri_en[3][2] = {
		{1, 2}, {2, 0}, {0, 1}
	};
	int test[MAX_SPR_POLY_SIZE];

	memset(test, 0, sizeof(int)*MAX_SPR_POLY_SIZE);
	for (i = 0; i < vecp.size(); i++)
	{
		idx = vecp[i];

		if (test[idx] == 1)
			continue;

		for (m = 0; m < 3 && num < 2; m++)
		{
			j1 = poly->t[idx][tri_en[m][0]];
			j2 = poly->t[idx][tri_en[m][1]];

			if (j1 == i1 && j2 == i2)
			{
				face1[num] = idx;
				code1[num] = m;

				/* for its neighbor */
				neg = face2[num] = poly->n[idx][m];
				for (k = 0; k < 3; k++)
					if (poly->n[neg][k] == idx)
					{
						code2[num] = k;
						break;
					}
				if (!(k < 3 && poly->t[neg][tri_en[k][0]] == i2 &&
					poly->t[neg][tri_en[k][1]] == i1))
				{
					k = k;
					spdlog::info("Error find in get_edge_parents.\n");
					exit(1);
				}
				//			assert(k < 3 && poly->t[neg][tri_en[k][0]] == i2 &&
				//				poly->t[neg][tri_en[k][1]] == i1);

				test[idx] = 1;
				test[neg] = 1;

				num++;
			}
			else if (j1 == i2 && j2 == i1)
			{
				face2[num] = idx;
				code2[num] = m;

				/* for its neighbor */
				neg = face1[num] = poly->n[idx][m];
				for (k = 0; k < 3; k++)
					if (poly->n[neg][k] == idx)
					{
						code1[num] = k;
						break;
					}
				assert(k < 3 && poly->t[neg][tri_en[k][0]] == i1 &&
					poly->t[neg][tri_en[k][1]] == i2);

				test[idx] = 1;
				test[neg] = 1;

				num++;
			}
		}
	}

	*nume = num;

	return num;
}

inline int get_edge_parents(int i1, int i2, int nodeParents[], int numParents,
	polyhedron *poly, int *face1, int *face2,
	int *code1, int *code2, int *nume)
{
	int i, m, k, idx, neg, num = 0, j1, j2;
	const int tri_en[3][2] = {
		{1, 2}, {2, 0}, {0, 1}
	};
	int test[MAX_SPR_POLY_SIZE];

	memset(test, 0, sizeof(int)*MAX_SPR_POLY_SIZE);
	for (i = 0; i < numParents; i++)
	{
		idx = nodeParents[i];

		if (test[idx] == 1)
			continue;

		for (m = 0; m < 3 && num < 2; m++)
		{
			j1 = poly->t[idx][tri_en[m][0]];
			j2 = poly->t[idx][tri_en[m][1]];

			if (j1 == i1 && j2 == i2)
			{
				face1[num] = idx;
				code1[num] = m;

				/* for its neighbor */
				neg = face2[num] = poly->n[idx][m];
				for (k = 0; k < 3; k++)
					if (poly->n[neg][k] == idx)
					{
						code2[num] = k;
						break;
					}
				if (!(k < 3 && poly->t[neg][tri_en[k][0]] == i2 &&
					poly->t[neg][tri_en[k][1]] == i1))
				{
					k = k;
					spdlog::info("Error find in get_edge_parents.\n");
					exit(1);
				}
				//			assert(k < 3 && poly->t[neg][tri_en[k][0]] == i2 &&
				//				poly->t[neg][tri_en[k][1]] == i1);

				test[idx] = 1;
				test[neg] = 1;

				num++;
			}
			else if (j1 == i2 && j2 == i1)
			{
				face2[num] = idx;
				code2[num] = m;

				/* for its neighbor */
				neg = face1[num] = poly->n[idx][m];
				for (k = 0; k < 3; k++)
					if (poly->n[neg][k] == idx)
					{
						code1[num] = k;
						break;
					}
				assert(k < 3 && poly->t[neg][tri_en[k][0]] == i1 &&
					poly->t[neg][tri_en[k][1]] == i2);

				test[idx] = 1;
				test[neg] = 1;

				num++;
			}
		}
	}

	*nume = num;

	return num;
}
/* ----------------------------------------------------------------------
 * 获得所有包含某条边的所有父亲单元
 * -------------------------------------------------------------------*/
inline void get_edge_parents(struct polyhedron *poly, int idx1, int idx2,
	std::vector<int>& vecp)
{
	int vprt, fidx, fcod, sidx, scod;
	int m;
	std::stack<int> stkTreeSear;
	int test_flag[MAX_SPR_POLY_SIZE];
	const int tri_en[3][2] = {
		{1, 2}, {2, 0}, {0, 1}
	};

	memset(test_flag, 0, sizeof(int)*MAX_SPR_POLY_SIZE);

	vprt = poly->vprt[idx1];
	fidx = vprt >> 2;
	fcod = vprt & 0x3;
	stkTreeSear.push(vprt);
	test_flag[fidx] = 1;

	if (add_face_edge_pair(poly, idx1, idx2, fidx, vecp) > 0)
		return;

	while (!stkTreeSear.empty())
	{
		vprt = stkTreeSear.top();
		stkTreeSear.pop();
		fidx = vprt >> 2;
		fcod = vprt & 0x3;

		for (m = 0; m <= 2; m++)
		{
			if (m == fcod)
				continue;

			sidx = poly->n[fidx][m];
			if (sidx >= 0 && test_flag[sidx] == 0)
			{
				test_flag[sidx] = 1;
				if ((scod = node_face_code(poly, idx1, sidx)) >= 0)
				{
					if (add_face_edge_pair(poly, idx1, idx2, fidx, vecp) > 0)
						return;

					vprt = (sidx << 2) | scod;
					stkTreeSear.push(vprt);
				}
			}
		}
	}
}

int SPRImpl::judge_pair_from_four_facet(
	int edgeS, int edgeE,	 /* 公共边的首末点 */
	int othe[], bool orit[], /* othe: 四个面片中除首末点外的另外一个点，
							  * orit: 边在面中的绕向 */
	int retn[]
)
{
	int g[2] = { -1, -1 }, h, i, j;
	double v1, v2, v4;
	//	double v3;

	retn[0] = 0;
	for (i = 1, j = 0; i < 4; i++)
	{
		if (orit[0] ^ orit[i])
		{
			g[j] = i;
			j++;
			assert(j <= 2);
		}
		else
		{
			h = i;
		}
	}
	assert(j == 2);

	v1 = isVisible_Search(vertices[edgeS], vertices[edgeE], vertices[othe[0]], vertices[othe[g[0]]]);
	v2 = isVisible_Search(vertices[edgeS], vertices[edgeE], vertices[othe[0]], vertices[othe[g[1]]]);

	if (!orit[0])
	{
		v1 = -v1;
		v2 = -v2;
	}

#if 0
	if (v1 >= 0.0 && v2 >= 0.0)
	{
		//	v3 = v1;
			//v3 = isVisible_Search(edgeS, edgeE, othe[g[0]], othe[0]);
		v4 = isVisible_Search(vertices[edgeS], vertices[edgeE], vertices[othe[g[0]]], vertices[othe[g[1]]]);
		v4 = orit[0] ? -v4 : v4;
		//	if ((v3 <= 0.0 && v4 >= 0.0) || (v3 >= 0.0 && v4 <= 0.0))
		if (v4 <= 0.0)

		{
			retn[1] = g[0];
			retn[2] = g[1];
		}
		else
		{
			retn[1] = g[1];
			retn[2] = g[0];
		}
		retn[3] = h;
	}
	else if (v1 <= 0.0 && v2 <= 0.0)
	{
		//	v3 = isVisible_Search(edgeS, edgeE, othe[g[0]], othe[0]);
	//		v3 = orit[0] ? -v1 : v1;
		v4 = isVisible_Search(vertices[edgeS], vertices[edgeE], vertices[othe[g[0]]], vertices[othe[g[1]]]);
		v4 = orit[0] ? -v4 : v4;

		//		if ((v3 <= 0.0 && v4 >= 0.0) || (v3 >= 0.0 && v4 <= 0.0))
		if (v4 >= 0.0)
		{
			retn[1] = g[1];
			retn[2] = g[0];
		}
		else
		{
			retn[1] = g[0];
			retn[2] = g[1];
		}
	}
	else
	{
		if (v1 > 0.0)
		{
			retn[1] = g[0];
			retn[2] = g[1];
		}
		else
		{
			retn[1] = g[1];
			retn[2] = g[0];
		}
	}
#endif

	if (v1*v2 >= 0.0)
	{
		v4 = isVisible_Search(vertices[edgeS], vertices[edgeE], vertices[othe[g[0]]], vertices[othe[g[1]]]);
		v4 = orit[0] ? -v4 : v4;

		if (v4 >= 0.0)
		{
			retn[1] = g[1];
			retn[2] = g[0];
		}
		else
		{
			retn[1] = g[0];
			retn[2] = g[1];
		}
	}
	else
	{
		if (v1 > 0.0)
		{
			retn[1] = g[0];
			retn[2] = g[1];
		}
		else
		{
			retn[1] = g[1];
			retn[2] = g[0];
		}
	}
	retn[3] = h;

	return 1;
}

/* ----------------------------------------------------------------------
 * 根据单元信息填充一个InsertedTetra对象
 * -------------------------------------------------------------------*/
int SPRImpl::fill_inserted_tetra(
	int ELE[], int f,			/* 节点和基面 */
	InsertedTetra *instTetra,
	polyhedron *poly,
	polyhedron *sub_poly)
{
	int i, j, k, m, idx, num_keep;
	int i1, i2, i3, i4, i5, i6, f1, f2, f3, f4, if1, if2, numf = 0;
	int cod1, cod2, cod3, cod4, loc, locNeig;
	int neig, newNeig;
	//	std::vector<int> vecp;
	int numParents = MAX_SPR_POLY_SIZE;
	const int edg_fac[3][2] = { {2, 1}, {0, 2}, {1, 0} };	/* 和新增3条边相应的面对，正前负后 */

	/* 确定4个面的信息, 保证索引最小的点是第1个点，且
	 * 保证每个面的法向指向区域内部 */
	fill_inst_tetra_facets(ELE, instTetra);

	numParents = MAX_SPR_POLY_SIZE;
	get_node_parents(poly, ELE[3], fit_nodeParents, &numParents);

	memset(instTetra->idx, -1, sizeof(int) * 4);
	memset(instTetra->flg, -1, sizeof(int)*poly->num_t);
	instTetra->idx[3] = f;
	instTetra->flg[f] = 3;
	instTetra->num_del = 1;
	for (i = 0; i < numParents && instTetra->num_del < 4; i++)
	{
		idx = fit_nodeParents[i];
		for (j = 0; j < 3; j++)
		{
			if (poly->t[idx][0] == instTetra->t[j][0] &&
				poly->t[idx][1] == instTetra->t[j][1] &&
				poly->t[idx][2] == instTetra->t[j][2])	//如果本来就有这个面，那么就有必要将其删掉
			{
				instTetra->num_del++;
				instTetra->idx[j] = idx;
				instTetra->flg[idx] = j;

				break;
			}
		}
	}

	if (instTetra->num_del < 4) /* 需要增加某个面 */
	{
		for (i = 0, j = 0; i < 4; i++)
		{
			if (instTetra->idx[i] < 0)
			{
				instTetra->idx[i] = (instTetra->num_del - poly->num_t) - j;	/* 用负值以区分 */
				j++;
			}
		}
	}

	/* 形成新的子多面体. 将poly中的面拷贝到sub_poly中 */
	memcpy(sub_poly->used_v, poly->used_v, sizeof(int)*num_vertices);
	for (i = 0, j = 0; i < poly->num_t; i++)
	{
		loc = instTetra->flg[i];

		if (loc >= 0)
		{
			for (m = 0; m < 3; m++)
			{
				sub_poly->used_v[instTetra->t[loc][m]]--;		//对于要删去的面，不仅仅不将其转入Q中，反倒是降低其引用次数
			}
		}
		else
		{
			memcpy(sub_poly->t[j], poly->t[i], sizeof(TRI));
			memcpy(sub_poly->n[j], poly->n[i], sizeof(NEG));

			fit_g2l[i] = j++;
		}
	}
	num_keep = j;

	/* 将新增加的面拷贝到sub_poly中 */
	for (i = 0; i < 3; i++)
	{
		if (instTetra->idx[i] < 0)
		{/* 说明这个面需要加进去，但是，要注意旋转方向 */
			sub_poly->t[j][0] = instTetra->t[i][0];
			sub_poly->t[j][1] = instTetra->t[i][2];
			sub_poly->t[j][2] = instTetra->t[i][1];

			sub_poly->used_v[sub_poly->t[j][0]]++;	//增加对应面的点的引用次数
			sub_poly->used_v[sub_poly->t[j][1]]++;
			sub_poly->used_v[sub_poly->t[j][2]]++;

			j++;
		}
	}
	sub_poly->num_t = j;

	/* --------------------------------------------
	 * 理顺相邻关系
	 * 1. 进行索引映射
	 * 2. 将被删除的面的相邻关系导引到新增加的面上
	 * -------------------------------------------*/
	for (i = 0; i < num_keep; i++)
	{
		for (j = 0; j < 3; j++)
		{
			neig = sub_poly->n[i][j];
			loc = instTetra->flg[neig];
			if (loc >= 0)
			{/* 指向一个被删除的面 */
				i1 = sub_poly->t[i][edg_fac[j][0]];
				i2 = sub_poly->t[i][edg_fac[j][1]];

				for (k = 0; k < 3; k++)
				{
					locNeig = instTetra->n[loc][k];
					newNeig = -instTetra->idx[locNeig];
					if (newNeig > 0)
					{
						assert(newNeig < sub_poly->num_t);

						for (m = 0; m < 3; m++)
						{
							if (sub_poly->t[newNeig][edg_fac[m][0]] == i2 &&
								sub_poly->t[newNeig][edg_fac[m][1]] == i1)
							{
								sub_poly->n[i][j] = newNeig;
								sub_poly->n[newNeig][m] = i;
								break;
							}
						}

						if (m < 3)
							break;
					}
				}
				if (m >= 3)
				{
					/* 此时提前返回，是因为已检测出新增加的四面体单元的2个侧面同为旧面片，但
					 * 它们并不是面片对，因此，它们构成的内角必然和区域外部存在重叠，为无效情形 */
					return 0;
				}
				assert(k < 3 && m < 3);
			}
			else
			{
				/* 指向一个已有面 */
				sub_poly->n[i][j] = fit_g2l[neig];
			}
		}
	}

	/* ----------------------------------------------------
	 * 理顺相邻关系
	 * 确认可能需要新增加的3条边的情况 3-->0, 3-->1, 3-->2
	 * ----------------------------------------------------*/
	i1 = ELE[3];
	instTetra->complex = false;
	for (i = 0; i < 3; i++)
	{
		i2 = ELE[i];

		if1 = edg_fac[i][0];
		if2 = edg_fac[i][1];
		f1 = -instTetra->idx[if1];
		f2 = -instTetra->idx[if2];

		if (f1 > 0 && f2 > 0)
		{/* 新增加一条边，如果这条边存在于原来的多面体中，这条边将成为复边 */
			if (get_edge_parents(i1, i2, fit_nodeParents, numParents, poly, fit_fac3_arr, fit_fac4_arr, fit_cod3_arr, fit_cod4_arr, &numf) > 0)
			{/* f3和f2互为邻居，f4和f1互为邻居 */
				instTetra->complex = true;

				for (j = 0; j < numf; j++)
				{
					f3 = fit_g2l[fit_fac3_arr[j]];
					cod3 = fit_cod3_arr[j];
					assert(f3 >= 0 && f3 <= sub_poly->num_t);

					for (m = 0; m < 3; m++)
					{
						if (sub_poly->t[f2][m] != i1 && sub_poly->t[f2][m] != i2)
						{
							cod2 = m;
							break;
						}
					}
					assert(m < 3);

					/* 这段代码似乎没有道理 */
//					if (isVisible_Search(vertices[i1], vertices[i2], 
//					vertices[sub_poly->t[f3][cod3]], vertices[sub_poly->t[f2][cod2]]) <= 0.0)
//						continue;

					f4 = fit_g2l[fit_fac4_arr[j]];
					cod4 = fit_cod4_arr[j];
					assert(f4 >= 0 && f4 <= sub_poly->num_t);

					for (m = 0; m < 3; m++)
					{
						if (sub_poly->t[f1][m] != i1 && sub_poly->t[f1][m] != i2)
						{
							cod1 = m;
							break;
						}
					}
					assert(m < 3);

					/* 这段代码似乎没有道理 */
	//				if (isVisible_Search(vertices[i2], vertices[i1], vertices[sub_poly->t[f4][cod4]], vertices[sub_poly->t[f1][cod1]]) <= 0.0)
	//					continue;

					/* 需确保新生成的两个面f1和f2在f3和f4的张角内 */
					i3 = sub_poly->t[f3][cod3];
					i4 = sub_poly->t[f4][cod4];
					i5 = sub_poly->t[f1][cod1];
					i6 = sub_poly->t[f2][cod2];

#ifdef _DEBUG
					for (m = 0; m < 3; m++)
					{
						if (sub_poly->t[f1][m] == i1)
						{
							assert(sub_poly->t[f1][(m + 1) % 3] == i2);
							assert(sub_poly->t[f1][(m + 2) % 3] == i5);
						}
						if (sub_poly->t[f3][m] == i1)
						{
							assert(sub_poly->t[f3][(m + 1) % 3] == i2);
							assert(sub_poly->t[f3][(m + 2) % 3] == i3);
						}
						if (sub_poly->t[f2][m] == i2)
						{
							assert(sub_poly->t[f2][(m + 1) % 3] == i1);
							assert(sub_poly->t[f2][(m + 2) % 3] == i6);
						}
						if (sub_poly->t[f4][m] == i2)
						{
							assert(sub_poly->t[f4][(m + 1) % 3] == i1);
							assert(sub_poly->t[f4][(m + 2) % 3] == i4);
						}
					}
#endif
					if (isVisible_Search(vertices[i1], vertices[i2], vertices[i3], vertices[i4]) >= 0.0)
					{/* 内部张角为小于180度的角，则i5和i6和（i1,i2,i3）可视，
						(i1, i2, i5(或i6))和i4可视 */
						if (isVisible_Search(vertices[i1], vertices[i2], vertices[i3], vertices[i5]) <= 0.0 ||
							isVisible_Search(vertices[i1], vertices[i2], vertices[i3], vertices[i6]) <= 0.0 ||
							isVisible_Search(vertices[i1], vertices[i2], vertices[i5], vertices[i4]) <= 0.0 ||
							isVisible_Search(vertices[i1], vertices[i2], vertices[i6], vertices[i4]) <= 0.0)
							continue;
					}
					else
					{/* 内部张角为大于180度的角，则i5和i6必和(i1,i2,i3)不可视，
						或(i1,i2,i5(或i6))和i4不可视 */
						if (isVisible_Search(vertices[i1], vertices[i2], vertices[i4], vertices[i5]) >= 0.0 &&
							isVisible_Search(vertices[i1], vertices[i2], vertices[i5], vertices[i3]) >= 0.0 ||
							isVisible_Search(vertices[i1], vertices[i2], vertices[i4], vertices[i6]) >= 0.0 &&
							isVisible_Search(vertices[i1], vertices[i2], vertices[i6], vertices[i3]) >= 0.0)
							continue;
					}

					sub_poly->n[f3][cod3] = f2;
					sub_poly->n[f4][cod4] = f1;
					sub_poly->n[f2][cod2] = f3;
					sub_poly->n[f1][cod1] = f4;
					break;
				}

				if (j >= numf)
				{/* unqualified decomposition */
					/* 此时提前返回，是因为已检测出新增加的面片不在区域内部 */
					return 0;
				}
			}
			else
			{/* f1和f2互为邻居 */

				for (m = 0; m < 3; m++)
				{
					if (sub_poly->t[f2][m] != i1 && sub_poly->t[f2][m] != i2)
					{
						sub_poly->n[f2][m] = f1;
						break;
					}
				}

				for (m = 0; m < 3; m++)
				{
					if (sub_poly->t[f1][m] != i1 && sub_poly->t[f1][m] != i2)
					{
						sub_poly->n[f1][m] = f2;
						break;
					}
				}
			}
		}
	}

	//处理点的问题, 记录可用点的数目
	//直接把所有的点都复制到Q里面去，然后再通过点的引用次数来表征这个点是否还存在
	for (i = 0, j = 0; i < num_vertices; i++)
	{
		if (sub_poly->used_v[i] > 0)
			j++;
	}
	sub_poly->num_vable = j;		//这里存的就是现在还可以使用的点的数目

	for (i = 0; i < sub_poly->num_t; i++)
	{
		for (m = 0; m < 3; m++)
		{
			i1 = sub_poly->t[i][m];
			sub_poly->vprt[i1] = (i << 2) | m;
		}
	}

	return 1;
}

/* 确定4个面的信息, 保证索引最小的点是第1个点，且保证每个面的法向指向区域内部 */
void SPRImpl::fill_inst_tetra_facets(int ELE[], InsertedTetra *instTetra)
{
	if (ELE[3] < ELE[2] && ELE[3] < ELE[1])
	{/* ELE[3] 最小 */
		instTetra->t[0][0] = ELE[3];
		instTetra->t[0][1] = ELE[2];
		instTetra->t[0][2] = ELE[1];
		instTetra->n[0][0] = 3;
		instTetra->n[0][1] = 1;
		instTetra->n[0][2] = 2;
	}
	else if (ELE[2] < ELE[1] && ELE[2] < ELE[3])
	{/* ELE[2] 最小 */
		instTetra->t[0][0] = ELE[2];
		instTetra->t[0][1] = ELE[1];
		instTetra->t[0][2] = ELE[3];
		instTetra->n[0][0] = 2;
		instTetra->n[0][1] = 3;
		instTetra->n[0][2] = 1;
	}
	else
	{/* ELE[1] 最小 */
		instTetra->t[0][0] = ELE[1];
		instTetra->t[0][1] = ELE[3];
		instTetra->t[0][2] = ELE[2];
		instTetra->n[0][0] = 1;
		instTetra->n[0][1] = 2;
		instTetra->n[0][2] = 3;
	}

	if (ELE[3] < ELE[2] && ELE[3] < ELE[0])
	{/* ELE[3] 最小 */
		instTetra->t[1][0] = ELE[3];
		instTetra->t[1][1] = ELE[0];
		instTetra->t[1][2] = ELE[2];
		instTetra->n[1][0] = 3;
		instTetra->n[1][1] = 2;
		instTetra->n[1][2] = 0;
	}
	else if (ELE[2] < ELE[0] && ELE[2] < ELE[3])
	{/* ELE[2] 最小 */
		instTetra->t[1][0] = ELE[2];
		instTetra->t[1][1] = ELE[3];
		instTetra->t[1][2] = ELE[0];
		instTetra->n[1][0] = 2;
		instTetra->n[1][1] = 0;
		instTetra->n[1][2] = 3;
	}
	else
	{/* ELE[0] 最小 */
		instTetra->t[1][0] = ELE[0];
		instTetra->t[1][1] = ELE[2];
		instTetra->t[1][2] = ELE[3];
		instTetra->n[1][0] = 0;
		instTetra->n[1][1] = 3;
		instTetra->n[1][2] = 2;
	}

	if (ELE[3] < ELE[1] && ELE[3] < ELE[0])
	{/* ELE[3] 最小 */
		instTetra->t[2][0] = ELE[3];
		instTetra->t[2][1] = ELE[1];
		instTetra->t[2][2] = ELE[0];
		instTetra->n[2][0] = 3;
		instTetra->n[2][1] = 0;
		instTetra->n[2][2] = 1;
	}
	else if (ELE[1] < ELE[0] && ELE[1] < ELE[3])
	{/* ELE[1] 最小 */
		instTetra->t[2][0] = ELE[1];
		instTetra->t[2][1] = ELE[0];
		instTetra->t[2][2] = ELE[3];
		instTetra->n[2][0] = 1;
		instTetra->n[2][1] = 3;
		instTetra->n[2][2] = 0;
	}
	else
	{/* ELE[0] 最小 */
		instTetra->t[2][0] = ELE[0];
		instTetra->t[2][1] = ELE[3];
		instTetra->t[2][2] = ELE[1];
		instTetra->n[2][0] = 0;
		instTetra->n[2][1] = 1;
		instTetra->n[2][2] = 3;
	}

	instTetra->t[3][0] = ELE[0];
	instTetra->t[3][1] = ELE[1];
	instTetra->t[3][2] = ELE[2];
	instTetra->n[3][0] = 0;
	instTetra->n[3][1] = 1;
	instTetra->n[3][2] = 2;
}


/* ----------------------------------------------------------------------
 * 对复边进行排序
 * -------------------------------------------------------------------*/
int SPRImpl::sort_comp_edges(struct polyhedron *poly, GridEdgeHash *eHash)
{
	int i, fidx, ecod, eort;
	int edgeS, edgeE;
	int othe[4], retn[4];
	bool orit[4];
	GEdgeHashNode *hNode;
	std::multimap<int, GEdgeHashNode*>::iterator it = eHash->begin();
	std::vector<int> vecElemCpy;
	int rtn = 0;
	char filename[512];

	while (it != eHash->end())
	{
		hNode = it->second;

		if (hNode->vecElems.size() == 4)
		{
			edgeS = it->first;
			edgeE = hNode->node2;

			for (i = 0; i < 4; i++)
			{
				fidx = hNode->vecElems[i];
				edge_face_code(poly, edgeS, edgeE, fidx, &ecod, &eort);

				othe[i] = poly->t[fidx][ecod];
				orit[i] = eort != 0 ? true : false;
			}

			judge_pair_from_four_facet(edgeS, edgeE, othe, orit, retn);

			vecElemCpy.clear();
			for (i = 0; i < 4; i++)
				vecElemCpy.push_back(hNode->vecElems[retn[i]]);

			std::copy(vecElemCpy.begin(), vecElemCpy.end(), hNode->vecElems.begin());

			rtn = 1;
		}
		else if (hNode->vecElems.size() > 4)
		{
#ifdef _DEBUG
			sprintf(filename, "complex-poly-%d.pls", sce_compPolyID++);
			savepls(filename, poly, true);
#endif
			return -1;
		}
		else
		{
			if (hNode->vecElems.size() != 2)
			{
#ifdef _DEBUG
				sprintf(filename, "odd-poly-%d.pls", sce_oddPolyID++);
				savepls(filename, poly, true);
#endif
				return -1;
			}
			//		assert(hNode->vecElems.size() == 2);
		}

		++it;
	}

	return rtn;
}

#ifdef _WIN32
double secondsPerTick = 0.0;
#endif
#ifndef _WIN32
double SPRlogTime(void)
{
	double my_time;

#ifndef WIN32
	struct timeval tp;
	gettimeofday(&tp, NULL);
	my_time = ((double)tp.tv_sec) + (1.0e-6) * (double)(tp.tv_usec);
#else
	struct _timeb _tp;
	_ftime(&_tp);
	my_time = ((double)_tp.time) + (1.0e-3) * (_tp.millitm);
#endif
	return(my_time);
}
#else
double SPRlogTime(void)
{
	LARGE_INTEGER lv;

	// 获取CPU运行到现在的Tick数
	QueryPerformanceCounter(&lv);

	// 计算CPU运行到现在的时间
	// 比GetTickCount和timeGetTime更加精确
	//return lv.QuadPart;
	return secondsPerTick * lv.QuadPart;
}
#endif

void SPRlogTime_Init()
{
#ifdef WIN32
	LARGE_INTEGER ts;
	QueryPerformanceFrequency(&ts);
	secondsPerTick = 1.0 / ts.QuadPart;
#endif // WIN32
}

int SPRImpl::sort_quad_pair_test()
{
	int i, j, k, m, j1, j2, j3, j4, ort;

	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
		{
			if (j == i)
				continue;

			for (k = 0; k < 4; k++)
			{
				if (k == i || k == j)
					continue;

				m = 6 - i - j - k;
				j1 = i;
				j2 = j;
				j3 = k;
				j4 = m;
				ort = sort_quad_pair(&j1, &j2, &j3, &j4);
				printf("(%d,%d,%d,%d) %d (%d,%d,%d,%d)\n",
					i, j, k, m, ort, j1, j2, j3, j4);
			}
		}

	return 1;
}

int SPRImpl::four_facet_pair_test(struct polyhedron* poly)
{
	int othe[4], retn[4];
	bool orit[4];
	int i, j, i1, i2;

	othe[0] = 2;
	othe[1] = 3;
	othe[2] = 4;
	othe[3] = 5;

	assert(poly->num_t == 4);
	for (i = 0; i < 4; i++)
	{
		for (j = 0; j < 3; j++)
		{
			if (poly->t[i][j] == 0)
				i1 = j;
			if (poly->t[i][j] == 1)
				i2 = j;
		}
		if ((i2 - i1) == 1 || (i2 + 3 - i1) == 1)
		{
			orit[i] = true;
		}
		else
			orit[i] = false;
	}

	judge_pair_from_four_facet(0, 1, othe, orit, retn);

	printf("Pair Facet: (%d, %d) and (%d, %d).\n", retn[0], retn[1], retn[2], retn[3]);

	return 0;
}


int SPRImpl::assert_poly_neig(struct polyhedron *poly)
{
	int iElem, iNeig;
	int ret = 0;
	int m, l, i1, i2, i3, i4;
	const int tri_en[3][2] = {
		{1, 2}, {2, 0}, {0, 1}
	};

	for (iElem = 0; iElem < poly->num_t; iElem++)
	{
		for (m = 0; m <= 2; m++)
		{
			iNeig = poly->n[iElem][m];
			if (iNeig < 0 || iNeig >= poly->num_t)
				return 0;

			//ensure there are ia neighbor for pNeig POINTing to pElem
			for (l = 0; l <= 2; l++)
			{
				if (poly->n[iNeig][l] == iElem)
					break;
			}

			if (l > 2)
				return 0;

			i1 = poly->t[iElem][tri_en[m][0]];
			i2 = poly->t[iElem][tri_en[m][1]];
			i3 = poly->t[iNeig][tri_en[l][0]];
			i4 = poly->t[iNeig][tri_en[l][1]];

			if (!((i1 == i3 && i2 == i4) || (i2 == i3 && i1 == i4)))
				return 0;
		}
	}

	return 1;
}

int SPRImpl::assert_simp_poly_edge(struct polyhedron *simp_poly)
{
	int i;
	GridEdgeHash edgehash;
	std::vector<struct polyhedron*> vec_simp_poly;
	GEdgeHashNode *hNode;
	std::multimap<int, GEdgeHashNode*>::iterator it;
	int rtn = 1;

	std::vector<int> vecETopus(simp_poly->num_t);
	std::vector<int> vecENodes(3 * simp_poly->num_t);
	for (i = 0; i < simp_poly->num_t; i++)
	{
		vecETopus[i] = EEMAS::TRIANGLE;

		vecENodes[3 * i + 0] = simp_poly->t[i][0];
		vecENodes[3 * i + 1] = simp_poly->t[i][1];
		vecENodes[3 * i + 2] = simp_poly->t[i][2];
	}

	edgehash.set_hash_type(GridEdgeHash::MIN_KEY);
	edgehash.enable_elem_record(true);
	edgehash.initialize(vecETopus, vecENodes);

	it = edgehash.begin();
	while (it != edgehash.end())
	{
		hNode = it->second;

		if (hNode->vecElems.size() != 2)
		{
			rtn = 1;
			//		savepls("fail.pls", simp_poly);
			break;
		}
		++it;
	}

	return rtn;
}

/* ----------------------------------------------------------------------
 * 对中间结果进行记录
 * -------------------------------------------------------------------*/
int SPRImpl::is_nonoverlap_tet(divide_record *oldRd, int key, int *pos)
{
	int i;

	for (i = 0; i < oldRd->num; i++)
		if (oldRd->div_rd[i] == key)
			return 0;

	return 1;
}

#if 0
inline void create_div_record(divide_record *oldRd, int key, divide_record *newRd, int pos)
{
	int lft = 0, rgt = oldRd->num - 1, mid, pos = 0;
	int *opt = oldRd->div_rd, *npt = newRd->div_rd;

	newRd->num = oldRd->num + 1;
	/* find the position */
	while (lft <= rgt)
	{
		if (key < opt[i])
		{
			pos = i;
			break;
		}
		else if (key > opt[j])
		{
			pos = j + 1;
			break;
		}
		else
		{
			if (j == i + 1)
			{
				pos = j;
				break;
			}
			else
			{
				mid = (i + j) / 2;

				if (key > opt[mid])
					i = mid;
				else
					j = mid;
			}
		}
	}

	if (pos >= 1)
		memcpy(npt, opt, sizeof(int)*pos);
	npt[pos] = key;
	if (pos < oldRd->num)
		memcpy(&npt[pos + 1], &opt[pos], sizeof(int)*(oldRd->num - pos));
}
#else
int SPRImpl::create_div_record(divide_record *oldRd, int key, divide_record *newRd, int *found)
{
	int pos, *opt = oldRd->div_rd, *npt = newRd->div_rd;

	binarySearch(&pos, found, key, opt, oldRd->num - 1, 0);

	if (*found == 0)
	{
		if (pos >= 1)
			memcpy(npt, opt, sizeof(int)*pos);
		npt[pos] = key;
		if (pos < oldRd->num)
			memcpy(&npt[pos + 1], &opt[pos], sizeof(int)*(oldRd->num - pos));
		newRd->num = oldRd->num + 1;
	}

	return pos;
}
#endif

void SPRImpl::binarySearch(int * loc, int * found, int Item, int Array[], int right, int left)
{
	int middle = (right + left) / 2;
	(*loc) = right;
	if (*loc < 0)
		(*loc) = 0;
	(*found) = 0;
	if ((middle >= 0) && (Array[middle] == Item))
	{
		(*loc) = middle;
		(*found) = 1;
	}
	else if (right > left)
	{
		if (Array[middle] < Item)
			binarySearch(loc, found, Item, Array, right, middle + 1);
		else
			binarySearch(loc, found, Item, Array, middle - 1, left);
	}
	else if (Array[middle] < Item)
	{
		(*loc) = right + 1;
	}
}

int SPRImpl::is_tested_div_record(divide_record *cur_record)
{
	return setDivideRecord.find(*cur_record) != setDivideRecord.end();
}

void SPRImpl::add_div_record(divide_record *cur_record)
{
	setDivideRecord.insert(*cur_record);

#ifdef _DEBUG
	if (setDivideRecord.size() % 1000 == 0)
		spdlog::info("setDivideRecord.size() == {}\n", setDivideRecord.size());
#endif
}

/* ----------------------------------------------------------------------
 * 对面片和节点进行排序
 * -------------------------------------------------------------------*/
int compare_idx_quality(const void *arg1, const void *arg2)
{
	double val1 = ((IDX_QUALITY*)arg1)->quality;
	double val2 = ((IDX_QUALITY*)arg2)->quality;

	if (val1 < val2)		 return -1;
	else if (val1 > val2) return  1;
	else				  	 return  0;
}

double SPRImpl::calc_facet_quality(polyhedron *poly, int fIdx)
{
	double height, mindst;
	int a, b, c, d;
	int i, j, ng;
	int *poly_base_t = NULL, *poly_base_n = NULL, *poly_neig_t = NULL;
	double *pa, *pb, *pc, *pd, vad[3], norm1[3], norm2[3];
	const double cos_th = 0.0;//sqrt(3)/2.0;

	poly_base_t = &poly->t[fIdx][0];
	a = poly_base_t[0];
	b = poly_base_t[1];
	c = poly_base_t[2];

	pa = vertices[a];
	pb = vertices[b];
	pc = vertices[c];

	mindst = height = heightOfTriangle(pa, pb, pc, norm1);

	//	norm_3p(pa, pb, pc, norm1);
	vec_uni(norm1);

	poly_base_n = &poly->n[fIdx][0];
	for (i = 0; i < 3; i++)
	{
		ng = poly_base_n[i];
		poly_neig_t = &poly->t[ng][0];

		/* 确保d点在上方 a, b等于poly_base_t[(i+1)%3]/poly_base_t[(i+2)%3]*/
		a = poly_base_t[(i + 1) % 3];
		b = poly_base_t[(i + 2) % 3];
		d = poly_neig_t[0] + poly_neig_t[1] + poly_neig_t[2] - a - b;

		/* 确定d在上方 */
		vad[0] = vertices[d][0] - vertices[a][0];
		vad[1] = vertices[d][1] - vertices[a][1];
		vad[2] = vertices[d][2] - vertices[a][2];
		height = vec_dotp(vad, norm1);

		if (height > 0.0 && height < mindst)
		{
			/* 求相邻面片的法向, 面片法向朝里 */
			norm_3p(vertices[a], vertices[d], vertices[b], norm2);
			vec_uni(norm2);

			if (-vec_dotp(norm1, norm2) >= cos_th)
			{/* 只关心小于90度（cos_th = 0.0)或30度（cos_th = sqrt(3.0)/2的二面角 */
				mindst = height;
			}
		}
	}

	return mindst;
}

/* 计算相对值 */
double SPRImpl::calc_facet_quality_rel(polyhedron *poly, int fIdx, double mindst)
{
	double height;
	int a, b, c, d, e;
	int i, j, ng;
	double *pa, *pb, *pc, *pd, *pe, ved[3], norm1[3], norm2[3];
	const double cos_th = 0.0;//sqrt(3)/2.0;

	a = poly->t[fIdx][0];
	b = poly->t[fIdx][1];
	c = poly->t[fIdx][2];

	pa = vertices[a];
	pb = vertices[b];
	pc = vertices[c];

	norm_3p(pa, pb, pc, norm1);
	vec_uni(norm1);

	for (i = 0; i < 3; i++)
	{
		ng = poly->n[fIdx][i];
		e = poly->t[fIdx][i];
		pe = vertices[e];

		norm_3p(vertices[poly->t[ng][0]],
			vertices[poly->t[ng][1]],
			vertices[poly->t[ng][2]], norm2);
		vec_uni(norm2);

		for (j = 0; j < 3; j++)
		{
			if (poly->n[ng][j] == fIdx)
			{
				d = poly->t[ng][j];
				pd = vertices[d];

				if (isVisible_Search(vertices[a], vertices[b], vertices[c], vertices[d]) > 0.0 &&
					-vec_dotp(norm1, norm2) >= cos_th)
				{/* 小于30度的二面角 */
					vec_2p(pe, pd, ved);
					height = vec_dotp(ved, norm1) / vec_val(ved);

					if (height < mindst)
						mindst = height;
				}

				break;
			}
		}
		assert(j < 3);
	}

	return mindst;
}

#if 0
//zhaodawei add 2010-07-16
/* 计算相对值 */     //zhaodawei
inline double calc_facet_quality_rel(polyhedron *poly, int fIdx, polyhedron *pld_poly, double mindst)
{
	double height;
	int a, b, c, d, e;
	int i, j, ng;
	double *pa, *pb, *pc, *pd, *pe, ved[3], norm1[3], norm2[3], *t_1, *t_2, *t_3, *tmp_t;
	int *pld_poly_base = NULL;
	const double cos_th = 0.5;//sqrt(3)/2.0;
	int i_same = 0;
	double d_value1, d_value2;
	char b_f[3] = { 0, 0, 0 };
	int index = 0;
	double *pp_t1, *pp_t2, *pp_t3;
	int i_tmp1, i_tmp2;

	a = poly->t[fIdx][0];
	b = poly->t[fIdx][1];
	c = poly->t[fIdx][2];

	pa = vertices[a];
	pb = vertices[b];
	pc = vertices[c];

	norm_3p(pa, pb, pc, norm1);
	vec_uni(norm1);

	for (i = 0; i < pld_poly->num_t; i++)
	{
		b_f[0] = b_f[1] = b_f[2] = 0;
		pld_poly_base = pld_poly->t[i];
		for (j = 0; j < 3; j++)
		{
			e = pld_poly_base[j];
			if (e == a)
				b_f[0] = 1;
			else if (e == b)
				b_f[1] = 1;
			else if (e == c)
				b_f[2] = 1;
			else
				d = e;
		}
		if (b_f[0] + b_f[1] + b_f[2] != 2)
			continue;

		if (b_f[0] == false)
		{
			e = a;
			index = 0;
			pp_t1 = pb;
			pp_t2 = pc;
			pp_t3 = pa;
		}
		else if (b_f[1] == false)
		{
			e = b;
			index = 1;
			pp_t1 = pc;
			pp_t2 = pa;
			pp_t3 = pb;
		}
		else if (b_f[2] == false)
		{
			e = c;
			index = 2;
			pp_t1 = pa;
			pp_t2 = pb;
			pp_t3 = pc;
		}
		pe = vertices[e];
		pd = vertices[d];
		t_1 = vertices[pld_poly_base[0]];
		t_2 = vertices[pld_poly_base[1]];
		t_3 = vertices[pld_poly_base[2]];
		norm_3p(t_1, t_2, t_3, norm2);
		vec_uni(norm2);

		if (comp2fAngle(pp_t1, pp_t2, pp_t3, pd) <= cos_th) continue;

		vec_2p(vertices[i_tmp1], vertices[i_tmp2], ved);

		if (vec_dotp(ved, norm2) <= 0)
			continue;
		height = fabs(vec_dotp(ved, norm2));
		height /= vec_val(ved);
		if (height < mindst)
			mindst = height;
	}

	return mindst;
}
#else
double SPRImpl::calc_facet_quality_rel(polyhedron *poly, int fIdx, polyhedron *pld_poly, double mindst)
{
	double height;
	int a, b, c, d, e;
	int i0, i1, i2, i, j, ng;
	double *pa, *pb, *pc, *pd, *pe, ved[3], norm1[3], norm2[3];
	int *pld_poly_base = NULL;
	const double cos_th = 0.0;//0.5;//sqrt(3)/2.0;
	int i_same = 0;
	char b_f[3] = { 0, 0, 0 };

	a = poly->t[fIdx][0];
	b = poly->t[fIdx][1];
	c = poly->t[fIdx][2];

	pa = vertices[a];
	pb = vertices[b];
	pc = vertices[c];

	norm_3p(pa, pb, pc, norm1);
	vec_uni(norm1);

	for (i = 0; i < pld_poly->num_t; i++)
	{
		b_f[0] = b_f[1] = b_f[2] = 0;
		pld_poly_base = pld_poly->t[i];
		i0 = pld_poly_base[0];
		i1 = pld_poly_base[1];
		i2 = pld_poly_base[2];
		if (i0 == a || i0 == b || i0 == c)
			b_f[0] = 1;
		if (i1 == a || i1 == b || i1 == c)
			b_f[1] = 1;
		if (i2 == a || i2 == b || i2 == c)
			b_f[2] = 1;

		if (b_f[0] + b_f[1] + b_f[2] != 2)
			continue;
		d = b_f[0] ? (b_f[1] ? i2 : i1) : i0;
		e = a + b + c + d - i0 - i1 - i2;
		pe = vertices[e];
		pd = vertices[d];
		vec_2p(pe, pd, ved);
		height = vec_dotp(ved, norm1);
		if (height > 0.0)
		{
			height /= vec_val(ved);
			if (height < mindst)
			{
				norm_3p(vertices[i0], vertices[i1], vertices[i2], norm2);
				vec_uni(norm2);

				if (-vec_dotp(norm1, norm2) >= cos_th)
				{/* 只关心小于90度（cos_th = 0.0)或30度（cos_th = sqrt(3.0)/2的二面角 */
					mindst = height;
				}
			}
		}
	}

	return mindst;
}
#endif

//*****************************zhaodawei add 2010-07-16

float change_3int_float(int a, int b, int c)
{
	float aa, bb, cc;
	aa = (float)a;
	bb = ((float)b) / 10;
	cc = ((float)c) / 100;
	return aa + bb + cc;
}

void sort3int(int &a, int &b, int &c)
{
	int t;
	if (a > b)
	{
		t = a; a = b; b = t;
	}
	if (b > c)
	{
		t = b; b = c; c = t;
	}
	if (a > b)
	{
		t = a; a = b; b = t;
	}
}

void sort_tri(polyhedron *poly, int *pp)
{
	int mm[MAX_SPR_POLY_SIZE][3];
	float tt[MAX_SPR_POLY_SIZE];
	int i = 0;
	int j = 0;
	float min;
	int min_re;
	int aaa = -1;

	for (i = 0; i < poly->num_t; i++)
	{
		mm[i][0] = poly->t[i][0];
		mm[i][1] = poly->t[i][1];
		mm[i][2] = poly->t[i][2];
		sort3int(mm[i][0], mm[i][1], mm[i][2]);
		tt[i] = change_3int_float(mm[i][0], mm[i][1], mm[i][2]);
	}

	/* 按tt[i]增序排列 */
	for (i = 0; i < poly->num_t; i++)
	{
		min = tt[i];
		min_re = i;
		for (j = 0; j < poly->num_t; j++)
		{
			if (min > tt[j])
			{
				min = tt[j];
				min_re = j;
			}
		}

		pp[i] = min_re;
		tt[min_re] = (float)INT_MAX;
	}
}


//******************************************************

//zhaodawei modify 2010-07-16
void SPRImpl::sort_facets(polyhedron *poly)
{
	int i;
	IDX_QUALITY fQuality[MAX_SPR_POLY_SIZE];
	int min_num = 0;
	int pp[MAX_SPR_POLY_SIZE];
	double init_t;

	for (i = 0; i < poly->num_t; i++)
	{
		fQuality[i].idx = i;
		fQuality[i].quality = calc_facet_quality(poly, i);
	}

	// 	qsort(fQuality, poly->num_t, sizeof(IDX_QUALITY), compare_idx_quality);
	// 
	// 	for (i = 0; i < poly->num_t; i++)
	// 		poly->sort_f[i] = fQuality[i].idx;

#if 0
	sort_tri(poly, pp);

	init_t = fQuality[pp[0]].quality;
	min_num = pp[0];
	for (i = 1; i < poly->num_t; i++)
	{
		if (fQuality[pp[i]].quality < init_t)
		{
			init_t = fQuality[pp[i]].quality;
			min_num = pp[i];
		}
	}
#else
	init_t = fQuality[0].quality;
	min_num = 0;
	for (i = 1; i < poly->num_t; i++)
	{
		if (fQuality[i].quality < init_t)
		{
			init_t = fQuality[i].quality;
			min_num = i;
		}
	}
#endif

	poly->sort_f[0] = min_num;
}

/* compute the (square) of the minimum sine
   of all the dihedral angles in the tet defined
   by the four vertices (vtx1, vtx2, vtx3, vtx4)
*/
REAL spr_minsine(APOINT vtx1, APOINT vtx2, APOINT vtx3, APOINT vtx4)
{
	REAL *point[4] = { NULL, NULL, NULL, NULL };      /* tet vertices */
	REAL edgelength[3][4]; /* the lengths of each of the edges of the tet */
	REAL facenormal[4][3]; /* the normals of each face of the tet */
	REAL dx, dy, dz;       /* intermediate values of edge lengths */
	REAL facearea2[4];     /* areas of the faces of the tet */
	REAL pyrvolume;        /* volume of tetrahedron */
	REAL sine2, minsine2;  /* the sine (squared) of the dihedral angle */
	int i, j, k, l;          /* loop indices */
	REAL E[3][3] = { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0} };
	const REAL HUGEFLOAT = 1.0e100;
	double sinEdge01 = -1.0;

	/* get tet vertices */
	point[0] = vtx1;
	point[1] = vtx2;
	point[2] = vtx3;
	point[3] = vtx4;

	/* calculate the volume*6 of the tetrahedron */
	pyrvolume = (REAL)GEOM_FUNC::orient3d(point[0], point[1], point[2], point[3]);

	/* if the volume is zero, the quality is zero, no reason to continue */
	if (pyrvolume == 0.0)
	{
		return 0.0;
	}

	/* for each vertex/face of the tetrahedron */
	for (i = 0; i < 4; i++) {
		j = (i + 1) & 3;
		if ((i & 1) == 0) {
			k = (i + 3) & 3;
			l = (i + 2) & 3;
		}
		else {
			k = (i + 2) & 3;
			l = (i + 3) & 3;
		}

		/* compute the normal for each face */
		facenormal[i][0] =
			(point[k][1] - point[j][1]) * (point[l][2] - point[j][2]) -
			(point[k][2] - point[j][2]) * (point[l][1] - point[j][1]);
		facenormal[i][1] =
			(point[k][2] - point[j][2]) * (point[l][0] - point[j][0]) -
			(point[k][0] - point[j][0]) * (point[l][2] - point[j][2]);
		facenormal[i][2] =
			(point[k][0] - point[j][0]) * (point[l][1] - point[j][1]) -
			(point[k][1] - point[j][1]) * (point[l][0] - point[j][0]);

		/* compute (2 *area)^2 for this face */
		facearea2[i] = facenormal[i][0] * facenormal[i][0] +
			facenormal[i][1] * facenormal[i][1] +
			facenormal[i][2] * facenormal[i][2];

		/* compute edge lengths (squared) */
		for (j = i + 1; j < 4; j++) {
			dx = point[i][0] - point[j][0];
			dy = point[i][1] - point[j][1];
			dz = point[i][2] - point[j][2];
			edgelength[i][j] = dx * dx + dy * dy + dz * dz;
		}
	}

	minsine2 = HUGEFLOAT;     /* start with absurdly big value for sine */

	/* for each edge in the tetrahedron */
	for (i = 0; i < 3; i++) {
		for (j = i + 1; j < 4; j++) {
			k = (i > 0) ? 0 : (j > 1) ? 1 : 2;
			l = 6 - i - j - k;

			/* compute the expression for minimum sine, squared, over 4
			   The reason it's over 4 is because the area values we have
			   are actually twice the area squared */
			   /* if either face area is zero, the sine is zero */
			if (facearea2[k] > 0 && facearea2[l] > 0)
			{
				sine2 = edgelength[i][j] / (facearea2[k] * facearea2[l]);
			}
			else
			{
				/* Encountered zero-area face */
				sine2 = 0.0;
			}

			/* update minimum sine */
			if (sine2 < minsine2)
			{
				minsine2 = sine2;
			}

			if (k + l == 5)
			{
				sinEdge01 = sine2;
			}
		}
	}

	//  return sqrt(minsine2) * pyrvolume;
	return sqrt(sinEdge01) * pyrvolume;
}

#if 1
int SPRImpl::calc_worst_facet(polyhedron *poly)
{
	double height, mindst, glob_mindst = 1.0e100;
	int a, b, c, d;
	int i, j, ng, fIdx;
	int *poly_base_t = NULL, *poly_base_n = NULL, *poly_neig_t = NULL;
	double *pa, *pb, *pc, *pd, vad[3];
	const double cos_th = 0.8660254037844386;//0.0;//sqrt(3)/2.0;
	int worstf = 0;
	double vab[3], vac[3], vbc[3], *vabc = NULL, *vadb = NULL;
	double l1_sq, l2_sq, l3_sq, lmax_sq;

	memset(cwf_faceNormalised, 0, sizeof(char)*poly->num_t);
	for (fIdx = 0; fIdx < poly->num_t; fIdx++)
	{
		/* 法向求出来了吗 */
		poly_base_t = &poly->t[fIdx][0];
		a = poly_base_t[0];
		b = poly_base_t[1];
		c = poly_base_t[2];

		pa = vertices[a];
		pb = vertices[b];
		pc = vertices[c];

		vab[0] = pb[0] - pa[0];
		vab[1] = pb[1] - pa[1];
		vab[2] = pb[2] - pa[2];
		vac[0] = pc[0] - pa[0];
		vac[1] = pc[1] - pa[1];
		vac[2] = pc[2] - pa[2];
		//		vec_2p(pa,pb,vab);
		//		vec_2p(pa,pc,vac);
		vabc = &cwf_norms[fIdx][0];
		if (cwf_faceNormalised[fIdx] == 0)
		{/* 求出法向量 */
			vec_crop(vab, vac, vabc);
			vabc[3] = sqrt(vabc[0] * vabc[0] + vabc[1] * vabc[1] + vabc[2] * vabc[2]);
			if (vabc[3] != 0.0)
			{/* degenerate elements */
				vabc[0] /= vabc[3];
				vabc[1] /= vabc[3];
				vabc[2] /= vabc[3];
			}
			cwf_faceNormalised[fIdx] = 1;
		}
		vbc[0] = pc[0] - pb[0];
		vbc[1] = pc[1] - pb[1];
		vbc[2] = pc[2] - pb[2];
		//	vec_2p(pb,pc,vbc);

		l1_sq = vab[0] * vab[0] + vab[1] * vab[1] + vab[2] * vab[2];
		l2_sq = vac[0] * vac[0] + vac[1] * vac[1] + vac[2] * vac[2];
		l3_sq = vbc[0] * vbc[0] + vbc[1] * vbc[1] + vbc[2] * vbc[2];
		lmax_sq = MAX_VALUE(MAX_VALUE(l1_sq, l2_sq), l3_sq);
		mindst = height = vabc[3] / sqrt(lmax_sq);

		poly_base_n = &poly->n[fIdx][0];
		for (i = 0; i < 3; i++)
		{
			ng = poly_base_n[i];
			poly_neig_t = &poly->t[ng][0];

			/* 确保d点在上方 a, b等于poly_base_t[(i+1)%3]/poly_base_t[(i+2)%3]*/
			a = poly_base_t[(i + 1) % 3];
			b = poly_base_t[(i + 2) % 3];
			d = poly_neig_t[0] + poly_neig_t[1] + poly_neig_t[2] - a - b;

			/* 确定d在上方 */
			vad[0] = vertices[d][0] - vertices[a][0];
			vad[1] = vertices[d][1] - vertices[a][1];
			vad[2] = vertices[d][2] - vertices[a][2];
			height = vec_dotp(vad, vabc);

			if (height > 0.0 && height < MIN_VALUE(mindst, glob_mindst))
			{/* height > 0.0: 形成一个小于180度的角 */
				vadb = &cwf_norms[ng][0];
				if (cwf_faceNormalised[ng] == 0)
				{
					/* 求相邻面片的法向, 面片法向朝里 */
					vab[0] = vertices[b][0] - vertices[a][0];
					vab[1] = vertices[b][1] - vertices[a][1];
					vab[2] = vertices[b][2] - vertices[a][2];

					vec_crop(vad, vab, vadb);
					vadb[3] = sqrt(vadb[0] * vadb[0] + vadb[1] * vadb[1] + vadb[2] * vadb[2]);
					if (vadb[3] != 0.0)
					{/* degenerate elements */
						vadb[0] /= vadb[3];
						vadb[1] /= vadb[3];
						vadb[2] /= vadb[3];
					}
					cwf_faceNormalised[ng] = 1;
				}

#ifdef _DEBUG
				//算得对不对，需要核对
				pa = vertices[a];
				pb = vertices[b];
				pd = vertices[d];
				c = poly_base_t[0] + poly_base_t[1] + poly_base_t[2] - a - b;
				pc = vertices[c];
				double sinVal = spr_minsine(pa, pb, pd, pc);
				double cosVal = -vec_dotp(vabc, vadb);
				if (fabs(1.0 - (sinVal*sinVal + cosVal * cosVal)) > 1.0e-4)
				{
					spdlog::info("Incosistent dihedral angle computation. sinV = {}; cosV = {}; \n", sinVal, cosVal);
					//		exit(1);
				}
#endif
				if (-vec_dotp(vabc, vadb) >= cos_th)
				{/* 只关心小于90度（cos_th = 0.0)或30度（cos_th = sqrt(3.0)/2的二面角 */
					mindst = height;
				}
			}
		}

		if (mindst < glob_mindst)
		{
			glob_mindst = mindst;
			worstf = fIdx;
		}
	}

	return worstf;
}
#else
inline int calc_worst_facet(polyhedron *poly)
{
	double height, maxcos, glob_maxcos = -2.0, cosVal;
	int a, b, c, d;
	int i, j, ng, fIdx;
	int *poly_base_t = NULL, *poly_base_n = NULL, *poly_neig_t = NULL;
	double *pa, *pb, *pc, *pd, vad[3];
	const double cos_th = 0.8660254037844386;//0.0;//sqrt(3)/2.0;
	int worstf = 0;
#error 需要去除静态变量
	static double norms[MAX_SPR_POLY_SIZE][4]; /* record the normal vector and its magnitude**2 */
	static char faceNormalised[MAX_SPR_POLY_SIZE];
	double vab[3], vac[3], vbc[3], *vabc = NULL, *vadb = NULL;
	double l1_sq, l2_sq, l3_sq, lmax_sq, lmin_sq, lmid_sq;

	memset(faceNormalised, 0, sizeof(char)*poly->num_t);
	for (fIdx = 0; fIdx < poly->num_t; fIdx++)
	{
		/* 法向求出来了吗 */
		poly_base_t = &poly->t[fIdx][0];
		a = poly_base_t[0];
		b = poly_base_t[1];
		c = poly_base_t[2];

		pa = vertices[a];
		pb = vertices[b];
		pc = vertices[c];

		vab[0] = pb[0] - pa[0];
		vab[1] = pb[1] - pa[1];
		vab[2] = pb[2] - pa[2];
		vac[0] = pc[0] - pa[0];
		vac[1] = pc[1] - pa[1];
		vac[2] = pc[2] - pa[2];
		//		vec_2p(pa,pb,vab);
		//		vec_2p(pa,pc,vac);
		vabc = &norms[fIdx][0];
		if (faceNormalised[fIdx] == 0)
		{/* 求出法向量 */
			vec_crop(vab, vac, vabc);
			vabc[3] = sqrt(vabc[0] * vabc[0] + vabc[1] * vabc[1] + vabc[2] * vabc[2]);
			if (vabc[3] != 0.0)
			{/* degenerate elements */
				vabc[0] /= vabc[3];
				vabc[1] /= vabc[3];
				vabc[2] /= vabc[3];
			}
			faceNormalised[fIdx] = 1;
		}
		vbc[0] = pc[0] - pb[0];
		vbc[1] = pc[1] - pb[1];
		vbc[2] = pc[2] - pb[2];
		//	vec_2p(pb,pc,vbc);

		l1_sq = vab[0] * vab[0] + vab[1] * vab[1] + vab[2] * vab[2];
		l2_sq = vac[0] * vac[0] + vac[1] * vac[1] + vac[2] * vac[2];
		l3_sq = vbc[0] * vbc[0] + vbc[1] * vbc[1] + vbc[2] * vbc[2];
		lmax_sq = MAX_VALUE(MAX_VALUE(l1_sq, l2_sq), l3_sq);
		height = vabc[3] / sqrt(lmax_sq);

		/* 尝试用角度代替长度 */
		lmin_sq = MIN_VALUE(MIN_VALUE(l1_sq, l2_sq), l3_sq);
		lmid_sq = l1_sq + l2_sq + l3_sq - lmin_sq - lmax_sq;
		maxcos = sqrt(MAX_VALUE(0.0, 1.0 - height * height / lmid_sq));
		//	maxcos = -1.0;

		poly_base_n = &poly->n[fIdx][0];
		for (i = 0; i < 3; i++)
		{
			ng = poly_base_n[i];
			poly_neig_t = &poly->t[ng][0];

			/* 确保d点在上方 a, b等于poly_base_t[(i+1)%3]/poly_base_t[(i+2)%3]*/
			a = poly_base_t[(i + 1) % 3];
			b = poly_base_t[(i + 2) % 3];
			d = poly_neig_t[0] + poly_neig_t[1] + poly_neig_t[2] - a - b;

			/* 确定d在上方 */
			vad[0] = vertices[d][0] - vertices[a][0];
			vad[1] = vertices[d][1] - vertices[a][1];
			vad[2] = vertices[d][2] - vertices[a][2];
			height = vec_dotp(vad, vabc);

			if (height > 0.0)
			{/* height > 0.0: 形成一个小于180度的角 */
				vadb = &norms[ng][0];
				if (faceNormalised[ng] == 0)
				{
					/* 求相邻面片的法向, 面片法向朝里 */
					vab[0] = vertices[b][0] - vertices[a][0];
					vab[1] = vertices[b][1] - vertices[a][1];
					vab[2] = vertices[b][2] - vertices[a][2];

					vec_crop(vad, vab, vadb);
					vadb[3] = sqrt(vadb[0] * vadb[0] + vadb[1] * vadb[1] + vadb[2] * vadb[2]);
					if (vadb[3] != 0.0)
					{/* degenerate elements */
						vadb[0] /= vadb[3];
						vadb[1] /= vadb[3];
						vadb[2] /= vadb[3];
					}
					faceNormalised[ng] = 1;
				}

				cosVal = -vec_dotp(vabc, vadb);
#ifdef _DEBUG
				//算得对不对，需要核对
				c = poly_base_t[0] + poly_base_t[1] + poly_base_t[2] - a - b;
				double sinVal = spr_minsine(vertices[a], vertices[b], vertices[d], vertices[c]);
				if (fabs(1.0 - (sinVal*sinVal + cosVal * cosVal)) > 1.0e-4)
				{
					spdlog::info("Incosistent dihedral angle computation. sinV = {}; cosV = {}; \n", sinVal, cosVal);
					//			exit(1);
				}
#endif
				if (cosVal > maxcos)
					maxcos = cosVal;
			}
		}

		if (maxcos > glob_maxcos)
		{
			glob_maxcos = maxcos;
			worstf = fIdx;
		}
	}

	return worstf;
}
#endif

int SPRImpl::get_sorted_facet(polyhedron *poly, int fIter)
{
	return poly->sort_f[fIter];
}
void SPRImpl::sort_nodes(polyhedron *poly, int fDig)
{
	return;
}
int SPRImpl::get_sorted_node(polyhedron *poly, int nIter)
{
	return nIter;
}

/* ----------------------------------------------------------------------
 * 为防止面片和复边但非对面上的点形成对应关系，排除这些点
 * -------------------------------------------------------------------*/
inline int is_mutex_facet_node(polyhedron *poly, int fIdx, int nIdx)
{
	/* There is a compiling error when using VS2010, so we disable this code */
#if 0
	std::set<int> setAroundNodes, setNeighbNodes;
	std::set<int>::iterator it;

	get_around_nodes(poly, fIdx, setAroundNodes);
	get_neighb_nodes(poly, fIdx, setNeighbNodes);

	if (setAroundNodes.size() != setNeighbNodes.size())
	{
		it = setNeighbNodes.begin();
		while (it != setNeighbNodes.end())
		{
			std::remove(setAroundNodes.begin(), setAroundNodes.end(), *it);
			++it;
		}

		return setAroundNodes.find(nIdx) != setAroundNodes.end();
	}
#endif
	assert(false);
	return 0;
}

inline int get_mutex_facet_node(polyhedron *poly, int fIdx, std::set<int>& setMutexNodes)
{
	std::set<int> setNeighbNodes;
	std::set<int>::iterator it;

	get_around_nodes(poly, fIdx, setMutexNodes);
	get_neighb_nodes(poly, fIdx, setNeighbNodes);

	if (setMutexNodes.size() != setNeighbNodes.size())
	{
		it = setNeighbNodes.begin();
		while (it != setNeighbNodes.end())
		{
			setMutexNodes.erase(*it);
			++it;
		}
	}
	else
		setMutexNodes.clear();

	return setMutexNodes.size();
}

inline int get_around_nodes(polyhedron *poly, int fIdx, std::set<int>& setAroundNodes)
{
	int i, a1, b1, c1, a2, b2, c2;

	a1 = poly->t[fIdx][0];
	b1 = poly->t[fIdx][1];
	c1 = poly->t[fIdx][2];

	for (i = 0; i < poly->num_t; i++)
	{
		if (i == fIdx)
			continue;

		a2 = poly->t[i][0];
		b2 = poly->t[i][1];
		c2 = poly->t[i][2];

		if ((b2 == a1 || b2 == b1 || b2 == c1) &&
			(c2 == a1 || c2 == b1 || c2 == c1) &&
			(a2 != a1 && a2 != b1 && a2 != c1))
			setAroundNodes.insert(a2);
		else if (
			(c2 == a1 || c2 == b1 || c2 == c1) &&
			(a2 == a1 || a2 == b1 || a2 == c1) &&
			(b2 != a1 && b2 != b1 && b2 != c1))
			setAroundNodes.insert(b2);
		else if (
			(a2 == a1 || a2 == b1 || a2 == c1) &&
			(b2 == a1 || b2 == b1 || b2 == c1) &&
			(c2 != a1 && c2 != b1 && c2 != c1))
			setAroundNodes.insert(c2);
	}

	return setAroundNodes.size();
}

inline int get_neighb_nodes(polyhedron *poly, int fIdx, std::set<int>& setNeighbNodes)
{
	int i, j, ng;

	for (i = 0; i < 3; i++)
	{
		ng = poly->n[fIdx][i];

		for (j = 0; j < 3; j++)
		{
			if (poly->n[ng][j] == fIdx)
			{
				setNeighbNodes.insert(poly->t[ng][j]);

				break;
			}
		}
		assert(j < 3);
	}

	return setNeighbNodes.size();
}

/* ------------------------------------------------------------------------------
 * 需要考虑如下特殊情形
 * 如果最终选定的iAdjNode是边界点，但重连后的网格质量低于预期，则该方案不可接受
 * 应该考虑次优方案
 * poly				in   多面体
 * divi				out  四面体化结果
 * minQ_forAll		in	 所有四面体化结果必须满足的最小质量指标
 * minQ_forBnd		in   所有以边界点为iAdjNode的四面体化结果必须满足的最小质量指标
 * nodeFlags		in	 节点标志，nodeFlags[i] = 0: 节点i为非边界点；
 *                                 nodeFlags[i] = 1: 节点i为边界点
 * -----------------------------------------------------------------------------*/
int SPRImpl::small_ball_reconn(polyhedron *poly, divide *divi, double minQ_forAll, double minQ_forBnd, int *nodeFlags, int *iAdjNode)
{
	int i, iNode, iFace, i1, i2, i3, i4;
	double v, AB, BC, AC, AD, BD, CD, qThreshold = minQ_forAll, qCurrMin = 0.0;
	divide currDivide;

	divi->num = 0;
	*iAdjNode = -1;
	/* loop over all of the boundary nodes */
	for (iNode = 0; iNode < num_vertices; iNode++)
	{
		if (poly->used_v[iNode] == 0)
			continue;

		currDivide.q = 1.0e10;
		currDivide.num = 0;
		/* loop over all of the boundary faces */
		for (iFace = 0; iFace < poly->num_t; iFace++)
		{
			i4 = iNode;
			i1 = poly->t[iFace][0];
			i2 = poly->t[iFace][1];
			i3 = poly->t[iFace][2];
			if (i1 == iNode || i2 == iNode || i3 == iNode)
				continue;

			v = GEOM_FUNC::orient3d(vertices[i4], vertices[i1], vertices[i2], vertices[i3]);
			if (v <= 0.0)
				break;

			/* 计算质量，确定是否可行 */
			AB = distance_square(vertices[i1], vertices[i2]);
			BC = distance_square(vertices[i2], vertices[i3]);
			AC = distance_square(vertices[i1], vertices[i3]);
			AD = distance_square(vertices[i1], vertices[i4]);
			BD = distance_square(vertices[i2], vertices[i4]);
			CD = distance_square(vertices[i3], vertices[i4]);

			v /= pow(AB + BC + AC + AD + BD + CD, 1.5);	//72*sqrt(3)*算出这个四面体的质量

			if (v <= qThreshold || (nodeFlags[iNode] == 1 && v <= minQ_forBnd))
				break;

			if (v < currDivide.q)
				currDivide.q = v;

			currDivide.te[currDivide.num][0] = i1;
			currDivide.te[currDivide.num][1] = i2;
			currDivide.te[currDivide.num][2] = i3;
			currDivide.te[currDivide.num][3] = i4;
			currDivide.num++;
		}

		if (iFace >= poly->num_t)
		{
			for (i = 0; i < currDivide.num; i++)
			{
				divi->te[i][0] = currDivide.te[i][0];
				divi->te[i][1] = currDivide.te[i][1];
				divi->te[i][2] = currDivide.te[i][2];
				divi->te[i][3] = currDivide.te[i][3];
			}
			divi->num = currDivide.num;
			qThreshold = divi->q = currDivide.q;
			*iAdjNode = iNode;
		}
	}

	return divi->num > 0;
}

int SPRImpl::small_ball_reconn(polyhedron *poly, divide *divi, double iniQ, int *iAdjNode)
{
	int i, iNode, iFace, i1, i2, i3, i4;
	double v, AB, BC, AC, AD, BD, CD, qThreshold = iniQ, qCurrMin = 0.0;
	divide currDivide;

	divi->num = 0;
	*iAdjNode = -1;
	/* loop over all of the boundary nodes */
	for (iNode = 0; iNode < num_vertices; iNode++)
	{
		if (poly->used_v[iNode] == 0)
			continue;

		currDivide.q = 1.0e10;
		currDivide.num = 0;
		/* loop over all of the boundary faces */
		for (iFace = 0; iFace < poly->num_t; iFace++)
		{
			i4 = iNode;
			i1 = poly->t[iFace][0];
			i2 = poly->t[iFace][1];
			i3 = poly->t[iFace][2];
			if (i1 == iNode || i2 == iNode || i3 == iNode)
				continue;

			v = GEOM_FUNC::orient3d(vertices[i4], vertices[i1], vertices[i2], vertices[i3]);
			if (v <= 0.0)
				break;

			/* 计算质量，确定是否可行 */
			AB = distance_square(vertices[i1], vertices[i2]);
			BC = distance_square(vertices[i2], vertices[i3]);
			AC = distance_square(vertices[i1], vertices[i3]);
			AD = distance_square(vertices[i1], vertices[i4]);
			BD = distance_square(vertices[i2], vertices[i4]);
			CD = distance_square(vertices[i3], vertices[i4]);

			v /= pow(AB + BC + AC + AD + BD + CD, 1.5);	//72*sqrt(3)*算出这个四面体的质量

			if (v < qThreshold)
				break;

			if (v < currDivide.q)
				currDivide.q = v;

			currDivide.te[currDivide.num][0] = i1;
			currDivide.te[currDivide.num][1] = i2;
			currDivide.te[currDivide.num][2] = i3;
			currDivide.te[currDivide.num][3] = i4;
			currDivide.num++;
		}

		if (iFace >= poly->num_t)
		{
			for (i = 0; i < currDivide.num; i++)
			{
				divi->te[i][0] = currDivide.te[i][0];
				divi->te[i][1] = currDivide.te[i][1];
				divi->te[i][2] = currDivide.te[i][2];
				divi->te[i][3] = currDivide.te[i][3];
			}
			divi->num = currDivide.num;
			qThreshold = divi->q = currDivide.q;
			*iAdjNode = iNode;
		}
	}

	return divi->num > 0;
}
