#ifndef __BLMesh_BLMesh_h__
#define __BLMesh_BLMesh_h__

#pragma once

#include "blmeshapi.h"
#include "BLMesh_define.h"
#include "BLFrontList.h"
#include "PotentialBEM.h"
#include "BLVector.h"
#include <map>
#include <cstring>
#include <cstdlib>
#include <cfloat>
#include "coord.h"
#include "Octree.h"
#include "OctreeAgent.h"
#include "symmetry.h"

#include "common.h"
#include "BLDomain.h"
#include <iostream>
#include "geom_func.h"
#include <memory>
#include <iomanip>
#include <string>
#include <array>
#include <unordered_set>
#include "NormalSmoothStrategy.h"
#include "IterationNormalSmoothStrategy.h"
#include "SimpleNormalSmoothStrategy.h"
#include "PrismStorehouse.h"
#include "../include/singlesymmetry.h"
//#include "DynamicArray.h"

#ifdef  MEMORY_DEBUG
#include "vld.h"
#endif
using std::cout;
using namespace OCT;
typedef struct
{
	int nelm;
	int npt;
	int* belm;
	double* bpt;
	double* bc;
	double* norm;
	double* u;
}BCArg;
using INPUTFORMAT = std::tuple<std::string, double*, int*, int*, std::vector<double>>;
class BLMesh
{
public:
	BLMesh()=delete;
	void DeleteOctree();
	BLMesh(int totalayer, bool bmulti = false, double bisotropic = 0.0)
		:m_blNxtFList(nullptr), m_blFrontList(nullptr),m_TriElm(), transfer_time(0)
	{
		generate_pyramid = true;
		expan_ratio = 1.1;
		process_id = 0;
		m_pBElm = nullptr;
		m_pU = nullptr;
		first_layer_FrtNods = 0;
		m_nTotalLayer = totalayer;
		m_nCurrLayer = 0;
		generate_pyramid_time = 0;
		check_prism_quality_time = 0;

		max_depth_=9;
		max_obj_=120;


		m_nNodes = 0;
		m_nElems = 0;
		m_pNodes = nullptr;
		m_pElems = nullptr;
		
		m_nAllocNodes = 0;
		m_nAllocElems = 0;
		m_nNewTri = 0;
		m_nOutTri = 0;
		m_bMultiNormal = bmulti;
		m_bIsotropicStop = bisotropic;

		m_nPyramid = 0;
		m_nTet = 0;
		m_nPrism = 0;
		m_nTri = 0;
		m_nQuad = 0;

		m_dStpAgl = 79;
		m_dFstLyHig = 0.036;
		m_sysValue = 0.0;
		m_symVals = nullptr;
		m_symFidx = nullptr;

		m_nBdryPnt = 0;
		m_pBdryPnt = nullptr;

		m_nSymBdrys = 0;
		m_pSymBdrys = nullptr;
		m_nAllocSymBdrys = 0;

		m_pPntIdx = nullptr;
		m_pPntElm = nullptr;
		m_pNorm = nullptr;
		m_pAddBdry = nullptr;
		m_nAddBdry = 0;


		m_nAllocTriElems = 0;
		m_ocAgent = nullptr;
		m_ocTree = nullptr;
		check_prism_time = 0;
		smooth_time=0;
		insert_time = 0;
		m_pDomain = nullptr;
	}

	~BLMesh(void);

public:
	//
	void InitBLMesh();

	double  GetAvergeLayer();
	//


	int SetBoundary(INPUTFORMAT filee);
	int ReadBoundary(const INPUTFORMAT file,bool clear=false);

	void destroyNode();
	void destroyFront();
	void CalZeroNorm();
	//



	// VM SaveBLMesh();

	void SaveBLMeshAndFreeMemory(VM & v);

	int SaveTetMesh(char * filename);
#ifdef _DEBUG
	void SaveCurrentLayerBLMesh(std::string filename, int layer);
#endif
	//
	int SaveBLMesh(char* filename);

	int SaveBdry(VM &mesh);
	int SaveBdry(char* filename);



	int SaveSurfGrid(char* filename);

	int SaveSurfGrid2(char * filename);

	void MeshTopoStatics();

	//
	int GetOuterBoundary(int* npt, int* nlem, double** pt, int** elm, int** l_to_g,bool add_symm=true);
	int OutputOuterBoundary(string fn, int npt, int nlem, double* pt, int* elm);
	void OutputFr2(string fn, int npt, int nlem, double* pt, int* elm);
	void OutputPls(string fn, int npt, int nlem, double* pt, int* elm);
	void resetZeroHeight(int NumNodes,BLNode ** nodes);
	//
	int VTK2Pls(char* filename);

	//
	void PreparePotential();

	//
	void GenerateBLMesh();
	
	//
	double GetDistance(BLNode* nod);

	//
	void PropagateNode(BLNode* blNod, BLVector normal, int iLayer);
	void StopPropagateNode(BLNode* blNod);


	bool CheckPrismIntersection(BLFront * blFront);

	bool CheckPrismValid(BLFront * blFront);

	void SmoothHeightRatio(BLNode * blNod, MBLNode * pNodes);
	//
	// void CollapseNode(BLNode* blNode, double* pBC, double* pU);

	void CheckInsertSideSuface(BLFront * blFront);

	void CheckInsertSuface(BLFront * blFront);
	void insertAndRmTriInOctree(BLFront* blFront);

	void Propagate();

	void getConnection(BLFront * blFront);

	void PreCheckPrismValid(BLFront * blFront);

	//
	void prePropagate(BLFront* blFront);

	//

	//
	void CalFrontSize(BLFront* blFront);

	BLVector CalNorm(double coord0[3], double coord1[3], double coord2[3]);
	BLVector TransNorm(BLVector norm, double* thd = NULL);
	double NormAngle(BLVector norm);

	double AngleNorm(BLVector norm1, BLVector norm2);
	void OutputAngle(char* filename, double* angle, int cnt, int idx);


	//
	//void CalVolField();
     int AddNode(BLVector pnt, double space);

	inline int AddElem(int nconn, int* conn, BLEntityTopology topu, int igoem = -1);
	


	int AddSymBdry(int nconn, int* conn, int topu);
	inline INDEX_TYPE AddTriElem(std::array<int, 3>& conn);
	inline INDEX_TYPE AddTriElem(int nconn, int* conn);
	bool CheckTri(int nconn, int* conn, double coord[3][3]);

	int RmvSymBdry(int p1, int p2);
	void SetSymBdryDelete(int i);
	bool IsSymBdryDelete(int i);

	vector<int> dSet;
	//check whether the node is a convex node.
	bool IsConvexNode(double coord0[3], double coord1[3], double coord2[3]);
	
	void test();


	int ElmBdryPtCnt(int eidx);
	void UpdateBdryNorm();


	void smoothPostProcess(BLNode* blNod, MBLNode* pNodes, BLVector originalNormal);
	void SmoothNodeNormAndRatio(BLNode* blNod, MBLNode* pNodes);
	void SmoothVirtualFrontHeight(BLNode* blNod, MBLNode* pNodes);
	void SmoothHorNodeNorm(BLNode* blNod, MBLNode* pNodes);
	//check whether is it need to create pyramid elements in the neighbor of blFront


	//void PostChckIntersect(BLFront* blFront);
	double CheckPyramidVolumn(double coordinates[][3]);
	bool CheckPyramidSkewness(double coordinates[][3]);
    bool CheckPyramidOrth(double coordinates[][3], double neighcenter[3]);


	bool ChckIntersectforTransit(BLFront* blFront);
	void CreatePyramid(BLFront* blFront);
	/**
	* @brife create
	*/
	void CreateTransitionElements();
	void UpdateSymmetry();


	bool CheckIsotroStop(BLNode * blNod);

	bool CheckStop(BLNode* blNod, BLNode* blNodNew, int iLayer);
	bool CheckPrismEveryVolumn(int nconn, int * conn);
	bool CheckPrismVolumn(int nconn, int * conn);

	bool CheckPrismValidity(int nconn, int* conn, int* idx = NULL);
	inline bool CheckPrismSkewness(int nconn, int* conn);
    inline bool BLMesh::CheckPrismOrth(BLFront *baseFront,bool include_top_face);
	/*
	* @brife 求交算法
	*/
	bool CheckIntersection(BLFront* blFront, bool* flag, int* tridx);

	void CreateTriangles(BLFront* blFront, int &itrix);

	void FixHightRatio(BLNode * blNod, MBLNode * pNodes);

	void PstprecsMergedElm();

	void CheckElemNum(int itmck);

	std::vector<int> RemvNodElm();
	void FreeMemory();

	void setpntelm(int npt, int nElm, int* pElm, int** pntidx, int** pntelm);

	void setpntelm3d(int npt, int nElm, Elem* pElem, int** pntidx, int** pntelm);

	void updateneig(int nElm, Elem* pElem, int* pntidx, int* pntelm, vector<std::array<int, 3>>& neigh);

	bool iselmedge(int eidx, int idx1, int idx2, int* pidx, Elem* pElm);

	void UpdateFrontNeig();

	void OrientLoop(int npt, int nbdry, double* pt, int*  elm, int** elmN, int** ortEl,vector<std::array<int, 2>>& subholes);
	void RemoveBox(int npt, int &nbdry, double * pt, int * elm);
	void OrientFace(int npt, int nbdry, double* pt, int* elm, int** elmN, int** ortEl,vector<double>& hole);

	//void GenOuterMesh(int smooth_attempt);

	void ReadAddBdry(char* filename, bool bfirst = false);
	int GetAddBoundary(int* npt, int* nlem, double** pt, int** elm, int** l_to_g);

	void GenTopMesh(VM & v,bool add_symm=true);

	void RemoveOverlapNodeAndElement();
	void FreeMemoryInFrontAndNode();
	
	void NarrowPointperturbation();

	std::set<int> GetNarrowConstrainedPointidx();

	void GenOuterMesh(int smooth_attempt);


	/*
	* @brife ： 创建初始尺寸场 
	*/
	void CalOriginSize(int npt, int nlem, double *pt, int *elm, double* &pt_size,int *l_to_g);


	void UpdateDomainGrid(int ngp, int nbp, int nel, double* g_coordx,double *g_coordy,double* g_coordz, int* gtopu, int** l_to_g);
	//multi domain
	void UpdateDomainGrid(int ngp, int nbp, int nel, double* g_coord, int* gtopu, int** l_to_g, int index);
	void UpdateSymplnGrid(int ngp, int nbp, int nel, double* g_coord, int* gtopu, int** l_to_g, int ifc = 0);

	void RemoveNonManifoldEdgeInSymFace();
	
	/*
	* @brife：暴力不鲁棒的算法，如果还是有非二边流型可以尝试使用这个算法
	*/
	void RemoveNonManifoldFrontByForce();
	void RemoveNonManifoldFront();
	void RemoveOverlapFace(int* bndry, int nbdry);
	void RemoveOutbdry(int nRmBdry, int* pRmBdry);
	int GetNumTriangleByEdge(int p1, int p2);
	bool IsOutbdry(int conn[3]);
	inline void GetFacIdx(int* bndry, int i, int* i1, int* i2, int* i3);
	void RmvUperNeigFrontsAndFreeNode(BLNode* blNod);
	void RmvUperNeigFronts(BLNode* blNod, bool bnext = true);

	void SetElmDelete(int i) { m_pElems[i].nconn=-1; }
	inline bool IsElmDelete(int i) { return m_pElems[i].nconn<0; }

	//Only symmetry plane paralleled to axis supported currently
	void UpdateSymnode();

	//void UpdateSymNormal(SymmetryPlane symPln, double symVal, BLVector& norm);

	bool IsSymLine(int p1, int p2);
	void CalSymplnBdry(int* nBdry, int** pBdry, int ifc = 0, bool add_symm=true);
	void RemvDeletedSymBdry();
	void OutputSymplnBdry(string fn, int nBdry, int* pBdry, int ifc = 0);
	void CalSymplnMsh(int nBdry, int* pBdry, int* nBdryElm, int** pBdryElm, double symval = 0.0, int ifc = 0);
	void removeNonManifoldPoint(std::set<std::array<int, 2>>& symlines);
	void SetSizingFunction(std::function<double(std::array<double, 3>)> func);
	static std::function<double(std::array<double, 3>)>& GetSizingFunction();
	void GlobalOptmize();
	std::vector<double> recommand_length_calculation();
	/*
	*  @note  重写 原暴力算法效率太低，
	*  @author yhf
	*  @brife 该函数返回front1中的第三点的index
	*/
	void NeighIdx(BLFront*& ft1, BLFront*& ft2, int* idx)
	{
		for (*idx = 0; *idx < DIM3; (*idx)++)
		{
			if (!ft2->IncludeNode(ft1->m_pBLNods[*idx])) {
				return;
			}
		}
	}
 
	void SetFrontSymm();
	void AddSymmSegment(BLFront* blFront);
	void SetBdryFront();
#ifdef _DEBUG
	void RemoveSymLineByManifoldCretiria();
	bool closeToPoint(int node_index,BLVector coord);
#endif
public:
    vector<int> boundary_to_delete_;

    double transfer_time;
    ConfigArgc cf;
    double check_prism_time;
    double check_prism_quality_time;
    double generate_pyramid_time;
    double smooth_time;
    double insert_time;

    double expan_ratio;

private:

	static std::function<double(std::array<double, 3>)> sizefuntion;

	int first_layer_FrtNods;
	int m_nTotalLayer;
	int m_nCurrLayer;


	//symmetry
	SymmetryJudger sj;
	std::map<int,TiGER::SymmetryPlane> faceid2sp;

	std::shared_ptr<BLFrontList> m_blFrontList;
	std::shared_ptr<BLFrontList> m_blNxtFList;

	std::vector<BLFront*> front_to_delete_;
	std::vector<BLNode*> node_to_delete_;
	vector<std::shared_ptr<BLFrontList>> m_blFrontListAll;
	queue<BLFront*> createPyramid_queue;
	int m_nNodes;
	MBLNode* m_pNodes;

	int m_nElems;
	Elem* m_pElems;

	int m_nAllocNodes;
	int m_nAllocElems;

	//for symmetry boundaries
	int m_nSymBdrys;
	int m_nTtlInitSymBdrys;
	Elem* m_pSymBdrys;
	int m_nAllocSymBdrys;
	PrismStorehouse prismhouse;///结构化存储附面层元件
public:
	int max_depth_;
	int max_obj_;
	bool generate_pyramid;
	std::unordered_set<int> insert_container;
	std::unordered_set<int> rm_container;
	std::unordered_set<BLFront*> scheck;/// ???
	std::unordered_set<int> new_inserted;
	deque<BLFront*> inser_queue;
	void insertSymFace(pair<int, int> p);
private:
	//std::vector<int> rm_flag;///存储在生成金字塔时的临时存储对象
	std::map<int, std::set<int>> m_mapSymline;///存储可能的边界信息，这里原来是用multimap的，这玩意儿有bug
	std::set<std::array<int, 2>> m_pyramid_symline;///金字塔单独存储
	map<int, double> distance_ratio;
	int m_nSurfNodes;
	int m_nSurfElems;

	//
	int m_nNewTri;

	//
	int m_nOutTri;

	//
	bool m_bMultiNormal;

	//
	double m_bIsotropicStop;


	double ave_front_size;
	int* m_pBElm;
	double* m_pU;
	double* m_pNorm;
	int process_id;
	// outerboundary
	int* outbdry;
    int *outbdryfarfield;
	int noutbdry;
	std::vector<BLFront*> m_vBdyFront;//理论上存储着在表面的单元，理论上应该是大小和面网格大小一样

	//for multidomain
	int* m_pAddBdry;
	int m_nAddBdry;

	double m_averU;
	double m_averUb;

	double averFontSize;
	public:
	int m_nPyramid;
	int m_nTet;
	int m_nPrism;
	int m_nTri;
	int m_nQuad;
	private:

	//for surface mesh
	int m_nBdryPnt;
	int* m_pBdryPnt;

	//for surface mesh
	int* m_pPntIdx;
	int* m_pPntElm;

	//user parameters

	// angle = 0 means ideal, angle = 90 means the opposite extreme
	double m_dStpAgl;
	double m_dFstLyHig;
	SymmetryPlane m_sysPlane;
	double m_sysValue;
	double* m_symVals;
	int m_nSymLoop;
	int* m_symFidx;	//face id

	//intersection parameters
	std::vector<int> m_vecData;
	std::vector<int> m_vecData_symm;
	OCCUBE m_cbCube;
	//std::vector<std::array<double, 6>> boxs_;//快速求包围盒
	int box_offset_;//记录boxs_数组和三角形数组之间的index差距
	//int* m_pTriElm;
	//int m_nTriElm;
	
	DynamicArray< SearchTriangle > m_TriElm;
	int m_nAllocTriElems;
	OctreeAgent* m_ocAgent;
	OCT::Octree* m_ocTree;

    DynamicArray<SearchTriangle> m_TriElm_symm;

	OctreeAgent* m_ocAgent_symm;
    OCT::Octree* m_ocTree_symm;

	std::vector<int> temp_to_del;
	//fmm related
	BLDomain* m_pDomain;
	
	BLVector model_centor;
};

#endif
