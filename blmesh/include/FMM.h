#ifndef __FMM_H__
#define __FMM_H__
#include <cstdlib>
#include <cstring>
#include "algebra.h"
#include "BLDomain.h"
#include "BLMultipole_define.h"
#include "PotentialBEM3D.h"
#include "blmeshapi.h"
class BLDomain;
class CMultipole
{
public:
	TREECUBE rootCube;
	FMMTREENODE* fmmTreeRoot;
	BLDomain* pDomain;

	double* x;    //迭代向量/结果向量
//	double* x0;   //初始向量
	double* b;    //右端向量

	complex<double>** RSmn;
	int p;                //多极展开的阶数(number of terms in multipole expansion)
	int max_nodes;        //叶子节点实际含有最大节点数
	int min_nodes;		  //叶子节点实际含有最小节点数
	int max_nearNodes;    //叶子节点实际含有最大邻居节点数
	int min_nearNodes;	  //叶子节点实际含有最小邻居节点数
	long tlLeaf;          //叶子节点数
	long tlLeafguess;	  //
	FMMTREENODE** arryLeafPoint;//指向叶子节点的指针数组(叶子节点数)

	PotentialBEM3D* m_pPotentialBEM;

	CMultipole()
	{
		pDomain=NULL;
		fmmTreeRoot=NULL;
		x=NULL;
//		x0=NULL;
		b=NULL;
		arryLeafPoint=NULL;
		RSmn=NULL;

		m_pPotentialBEM = NULL;
		g2l = NULL;
		l2g = NULL;
		m_pPts = NULL;
		m_pElms = NULL;
	}
	~CMultipole(){};

	int createMultipole(BLDomain* pDom);
	void initMultipole(BLDomain* pDom);
	void deleteMultipole();

public:
	long tlNode;
	long tlUnknown;

	int MAX_nodes;		  //叶子节点限定最大节点数given
	long tlFmmTreeNode;
//	FMMTREENODE** arryNodeToLeaf;//各单元所属叶子节点 (tlNode)

	//GMRES
//	int flag;
//    flag:
//    0 GMRES converged to the desired tolerance TOL within MAXIT iterations.
//    1 GMRES iterated MAXIT times but did not converge.
//    2 preconditioner M was ill-conditioned.
//    3 GMRES stagnated (two consecutive iterates were the same).
	double tolerance;     //设定迭代误差

	void ErrDisplay( char* message );
	int IT_inner;         //设定内迭代次数given
	int IT_outer;         //设定外迭代次数given
public:
	int it_inner;         //实际内迭代次数
	int it_outer;         //实际外迭代次数
	int max_level;
	int min_level;

public:
	int initFmmTreeRoot();
	void initTreeNode(FMMTREENODE* fmmTreeNode);
	int createFmmTree(FMMTREENODE* fmmTreeNode);
	void createLeafPoint(FMMTREENODE* fmmTreeNode);
	void countLeafPoint(FMMTREENODE* fmmTreeNode);
	void checkLeafNods();
	void getTreeNodeElms(int nNeig, FMMTREENODE** nodeNeigs, int *npt, double **pts, int *nelm, int **elms);
	int createNeighbor(FMMTREENODE* fmmTreeNode);
	int createInteraction(FMMTREENODE* fmmTreeNode);
	void allocateMoments(FMMTREENODE* fmmTreeNode);
	void freeMoments(FMMTREENODE* fmmTreeNode);
	int initALeaf_PotentialAndFlux(int k, complex<double> *uu,complex<double> *qq, complex<double>*** INT_RmnU, complex<double>*** INT_RmnQ);
	void allocatePreCOND();
	void freePreCOND();
	void getExpansionR(D3POINT pDt0, D3POINT pDt1);
	int build_AX_constantQ(double* bt);
	int build_AX_constantU(double* x0,double* bt);
	int build_preCOND_constant();
	int fmm_GMRES_constant();
	
	void ComputeTransform();
	int GetPtLeaf(Node *pNod, int *k);
	FMMTREENODE* GetPtTreeNode(Node* pNod, FMMTREENODE* fmmTreeNode);
	int GetResultPt(double *pt, double *val, int ilyer = -1);
	void getExpansionR_Local(D3POINT pDt0, D3POINT pDt, complex<double> RSmn0[21][21]);
	void getExpansionS_Local(D3POINT pDt0, D3POINT pDt1, complex<double> RSmn0[21][21]);
	void PrintRSmn(complex<double> RSmn0[21][21], char *filename);

	//for real problem
	void build_B();     //建立右端向量B
	int build_AX();     //建立矩阵与向量相乘A*X
	int build_AX_2();   //called by CDomain::solveUnKnown() 
//	int build_preCOND();//求解快速算法的预处理矩阵--multidomain
//	void preconditioning();//vtQ=M*vtU. vtU, vtQ used here temporarily

private:
	int isInCube(Node* pNode,TREECUBE* pCube);
	void resetFmmTree(FMMTREENODE* fmmTreeNode);
	void deleteFmmTree(FMMTREENODE* fmmTreeNode);
	void deleteList(FMMTREENODE* fmmTreeNode);
	void addNeighbor(FMMTREENODE* fmmTreeNode1,FMMTREENODE* fmmTreeNode2);
	void addInteraction(FMMTREENODE* fmmTreeNode1,FMMTREENODE* fmmTreeNode2);
	int countInteraction(FMMTREENODE* fmmTreeNode);
	int isNeighbor(FMMTREENODE* fmmTreeNode1,FMMTREENODE* fmmTreeNode2);
	int setLeafSequence();
	void getExpansionS(D3POINT pDt0, D3POINT pDt1);
	void addSequence(FMMTREENODE* fmmTreeNode);
	void initLeafs(double* xx);
	void SonToFatherTrans(FMMTREENODE* fmmTreeNode);
	void SonToFatherTrans_para(FMMTREENODE* fmmTreeNode);
	void UpWardPass(FMMTREENODE* fmmTreeNode);
	void UpWardPass_para(FMMTREENODE* fmmTreeNode);
	void FarToNearField(FMMTREENODE* fmmTreeNode1, FMMTREENODE* fmmTreeNode2);
	void FarToNearField_para(FMMTREENODE* fmmTreeNode1,FMMTREENODE* fmmTreeNode2);
	void FarToNear(FMMTREENODE* fmmTreeNode1,FMMTREENODE* fmmTreeNode2);
	void FarToNear1(FMMTREENODE* fmmTreeNode1,FMMTREENODE* fmmTreeNode2);
	int isFarField(D3POINT* pDt, FMMTREENODE* fmmTreeNode);
	void FatherToSonTrans(FMMTREENODE* fmmTreeNode);
	void FatherToSonTrans_para(FMMTREENODE* fmmTreeNode);
	void GetOneByPIRSmn(D3POINT pNode1,D3POINT pNode2);
	void DownWardPass(FMMTREENODE* fmmTreeNode);
	void DownWardPass_para(FMMTREENODE* fmmTreeNode);
	void deleteTemps(double* r, double* vh, double* u, double* u2, double* q, 
					double* h, double* f, double** V, double** QT, double** R, double** W);
	void preconditioning_constant(double* v, double* v0);

	void getNodePoint(FMMTREENODE* fmmTreeNode);
	void getNodePoint(FMMTREENODE* fmmTreeNode, int *tlTreeNods, int *pTreeNods);
	void getNeigValue(FMMTREENODE* fmmTreeNode, double *pt, double *val);
	void getChildValue(FMMTREENODE* fmmTreeNode, double *pt, double *val);
	void getNodNeigs(FMMTREENODE* fmmTreeNode, int *nNeig, FMMTREENODE** NodNeigs);
	bool isLeafNode(FMMTREENODE* fmmTreeNode);
	bool isNeighbor(FMMTREENODE* fmmTreeNode1,FMMTREENODE* fmmTreeNode2, bool caled);
	bool isChild(FMMTREENODE* fmmTreeNode1,FMMTREENODE* fmmTreeNode2);

	void PrintLmn(FMMTREENODE* fmmTreeRoot);
	void PrintMmn(FMMTREENODE* fmmTreeRoot);

private://temporary varibles
//	long tempSequence[2600];
	long *tempSequence;
	long tempNum;


	//
	//PotentialBEM3D* m_pPotentialBEM;
	int tlTreeNodePts;
	int *m_pTreeNodePts;
	int *g2l;
	int *l2g;
	int *eg2l;
	int *el2g;
	double *m_pPts;
	int *m_pElms;
};


#endif
