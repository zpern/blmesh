/** ----------------------------------------------------------------------------
	* TIGER, a Trustable Intelligent GridER
	*
	* TIGER的四面体网格生成API集成接口定义 (版本号：1.0)
	* Definitions of APIs for Surface Mesh Generation (Version 1.0)
	*
	* 陈建军 中国 浙江大学航空航天学院
	* 版权所有	  2020年2月6日
	* Chen Jianjun  School of Aeronautics & Astronautics,
	* Zhejiang University, P. R. China
	* Copyright reserved, Feb 6th, 2020
	*
	* 联系方式
	*   电话：+86-571-87951883
	*   传真：+86-571-87953167
	*   邮箱：chenjj@zju.edu.cn
	* For further information, please conctact
	*  Tel: +86-571-87951883
	*  Fax: +86-571-87953167
	*  Mail: chenjj@zju.edu.cn
	*
	* ------------------------------------------------------------------------------*/
#ifndef __tiger_volmtetrahedralizer_h__
#define __tiger_volmtetrahedralizer_h__

#ifdef WIN32
#ifdef libvolmt_EXPORTS 
#define DECL_VOLTET __declspec(dllexport)
#else
#define DECL_VOLTET //__declspec(dllimport)
#endif libvolmt_EXPORTS
#else
#define DECL_VOLTET //__declspec(dllimport)
#endif WIN32

/** -------------------------------------------------------------------------------
		* Given a surface triangulation defining a single-connected domain and a user
		* defined sizing function for element spacing, try to discretize the domain into a
		* combination of tetrahedra
		*
		* The inputs include three parts:
		* The first part is the group of boundary facets with all their normals pointing to
		* the interior of the domain
		*
		*	The second part is the interior constraints, i.e., points, edges and/or facets (if any).
		*
		* The third part of the input is a user defined sizing function, which is expected to
		* output a sizing value for an arbitrary point inside the problem domain
		*
		* The outputs are a sequence of tetrahedra filling in the domain, with the input boundary
		* & interior constraints respected
		*
		* Input Paramerters:
		*	bndPtNum		# of boundary points defining the exterior boundary & interior constraints.
		*	bndPts			the coordinates of boundary points. dim = 3*bndPtNum
		*	bndFctNum		# of exterior boundary facets
		* bndFcts			the connections of exterior boundary facets. dim = 3*bndFctNum
		* intCstNum		# of interior constraints. dim=3.
		*							intCstNum[0] -- # of interior points
		*							intCstNum[1] -- # of interior edges
		*							intCstNum[2] -- # of interior facets
		*	intCsts			the connections of interior constraints
		*							dim = intCstNum[0] + 2 * intCstNum[1] + 3*intCstNum[2]
		* sizingFunc	the user specified sizing function (if any)
		* Output Paramerters:
		*	objVolmTet		The tetrahedra object
		* Returned value:
		*	0				Completely succeed
		*	1				Partially succeed
		*	otherwise		Fail. Error info. will be defined later
		* ------------------------------------------------------------------------------*/
DECL_VOLTET int API_Tetrahedralize_Single_Domain(
	int bndPtNum, double bndPts[],
	int bndFctNum, int bndFcts[],
	int intCstNum[], int intCsts[],
	double(*sizingFunc)(double x, double y, double z),
	int *objVolmTet
);

DECL_VOLTET int getTestObj(int* obj1, int* obj2);

//add by zyj
DECL_VOLTET int API_Tetrahedralize_Single_Domain(
	int bndPtNum, double bndPts[],
	int bndFctNum, int bndFcts[],
	int intCstNum[], int intCsts[],
	int optimizeCount,
	double(*sizingFunc)(double x, double y, double z),
	int *objVolmTet
);

//add by zyj
#ifdef USE_TETGEN
DECL_VOLTET int API_Tetrahedralize_Multi_Domain(
	int bndPtNum, double bndPts[],
	int bndFctNum, int bndFcts[],
	int *objVolmTet
);
#endif
//add by zyj 
DECL_VOLTET int API_Merge_Multi_Domain(
	int domainNum, int domainVolmTets[],
	int * bndPtNum,		// DIM = domainNum
	int * ltowBndPtMap,	//ltowMap[domainIndex][localPtIndex] = worldPtIndex
	int *objVolmTet
);

/** -------------------------------------------------------------------------------
	* input Parameter:
	*	objVolmTet     The tetrahedra object
	* Returned value: Number of mesh points
* ------------------------------------------------------------------------------*/
DECL_VOLTET int API_GetTetrahedraPntNum(int objVolmTet);

/** -------------------------------------------------------------------------------
	* input Parameter:
	*	objVolmTet     The tetrahedra object
	* Returned value: Number of volume elements
* ------------------------------------------------------------------------------*/
DECL_VOLTET int API_GetTetrahedraElemNum(int objVolmTet);

/** -------------------------------------------------------------------------------
	* input Parameter:
	*	objVolmTet     The tetrahedra object
	*   index          The index of point (start from 1 !)
	* Output Paramerters:
	*   x,y,z          The coordinate of mesh point
	* Returned value:
	*	0				Completely succeed
	*	1				Partially succeed
	*	otherwise		Fail. Error info. will be defined later
* ------------------------------------------------------------------------------*/
DECL_VOLTET int API_GetTetrahedraPointCoord(int objVolmTet, int index, double &x, double &y, double &z);

/** -------------------------------------------------------------------------------
	* input Parameter:
	*	objVolmTet     The tetrahedra object
	*   index          The index of element (start from 1 !)
	* Output Paramerters:
	*   a,b,c,d        The index of mesh point
	* Returned value:
	*	0				Completely succeed
	*	1				Partially succeed
	*	otherwise		Fail. Error info. will be defined later
* ------------------------------------------------------------------------------*/
DECL_VOLTET int API_GetTetrahedraElemPntIdx(int objVolmTet, int index, int &a, int &b, int &c, int &d);

/** -------------------------------------------------------------------------------
	* input Parameter:
	*	objVolmTet     The tetrahedra object
	* Returned value:
	*	0				Completely succeed
	*	otherwise		Fail. Error info. will be defined later
* ------------------------------------------------------------------------------*/
DECL_VOLTET int API_DelTetrahedraObj(int objVolmTet);

/** -------------------------------------------------------------------------------
	* Returned value:
	*	n				Number of tetrahedra objects
* ------------------------------------------------------------------------------*/
DECL_VOLTET int API_TetrahedraObjNum();

/* --------------------------------------------------------------------------------
 * Compute the volume of the polyhedron defined by a set of facets
 * -------------------------------------------------------------------------------*/
DECL_VOLTET int API_Compute_Poly_Volume(
	int bndPtNum, double bndPts[],
	int bndFctNum, int bndFcts[],
	double *vol
);
#endif /* __tiger_volmtetrahedralizer_h__ */
