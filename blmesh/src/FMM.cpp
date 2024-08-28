////////////////////////////////////////////////////////////////
//
// Written by Jianming Zhang (zhangjm@homer.shinshu-u.ac.jp)
// Department of Mechanical Systems engineering, Shinshu University
// Copyright (c) September 20, 2003
//
// Edited by Zhoufang Xiao(xiaozf@zju.edu.cn), Jan 5, 2017
///////////////////////////////////////////////////////////////////
#include <spdlog/spdlog.h> 
 #include "FMM.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <set>
#ifndef __APPLE__
#include <omp.h>
#endif

// Output error messeage
void CMultipole::ErrDisplay( char *message )
{
	spdlog::info("Error info:%\n%s\n", message);
}

void CMultipole::initMultipole(BLDomain* pDom)
{
	int i;
	pDomain=pDom;
	rootCube.center.x=pDomain->m_dmCen.coord[0];
	rootCube.center.y=pDomain->m_dmCen.coord[1];
	rootCube.center.z=pDomain->m_dmCen.coord[2];

	rootCube.halfD=pDomain->m_dLength/2.0;
	tlNode=pDomain->tlNode;
	tlUnknown=pDomain->tlUnknown;

// 	m_pPotentialBEM = new PotentialBEM3D;
// 	m_pPotentialBEM->AllocSurfElmInfo(tlNode);
	m_pElms = new int[tlNode*DIM3];
	m_pPts = new double[pDomain->m_nNods*DIM3];

	//local for each thread
	m_pTreeNodePts = new int[tlNode];
	g2l = new int[pDomain->m_nNods];
	l2g = new int[pDomain->m_nNods];
	el2g = new int[tlNode];
	eg2l = new int[tlNode];

	tlFmmTreeNode=0;
	MAX_nodes=45;       //叶子节点含有最大单元数
//	MAX_nodes=100;       //叶子节点含有最大单元数
	tlLeaf=0;            //叶子节点数
	tlLeafguess = 0;

	//arryLeafPoint=new FMMTREENODE*[tlNode/(MAX_nodes/8)]; //100/8:=12.5
	//tlLeafguess = tlNode*2/MAX_nodes;
	//arryLeafPoint=new FMMTREENODE*[tlNode/MAX_nodes];
//	arryNodeToLeaf=new FMMTREENODE*[tlNode];//各单元所属叶子节点 (tlNode)

	//p=10;                //多极展开的阶数
	p=3;
//	IT_inner=100;         //设定内迭代次数
	IT_inner=25;         //设定内迭代次数
	IT_outer=50;         //设定外迭代次数
	tolerance=1.0e-5;    //设定迭代误差1.0e-6

	fmmTreeRoot=new FMMTREENODE;

	//local for each thread
	RSmn=new  complex<double>*[2*p+1];
	for(i=0;i<=2*p;i++) RSmn[i]=new  complex<double>[i+1];

	x=new double[tlUnknown];    //迭代向量/结果向量
//	x0=new double[tlUnknown];   //初始向量
}

////////////////////////////////////
int CMultipole::createMultipole(BLDomain *pDom)
{
	initMultipole(pDom);
	if(initFmmTreeRoot())
	{
		ErrDisplay("InitFmmTreeRoot wrong");
		return 1;
	}
	if(createFmmTree(fmmTreeRoot))
	{
		ErrDisplay("Create FMM tree wrong");
		return 1;
	}
	countLeafPoint(fmmTreeRoot);
	createLeafPoint(fmmTreeRoot);
	//checkLeafNods();
	if(createNeighbor(fmmTreeRoot))
	{
		ErrDisplay("Create Neighbour wrong");
		return 1;
	}
	if(setLeafSequence())
	{
		ErrDisplay("AllocateCoeffients wrong");
		return 1;
	}
	if(createInteraction(fmmTreeRoot))
	{
		ErrDisplay("Create interaction list wrong");
		return 1;
	}

	allocateMoments(fmmTreeRoot);
	resetFmmTree(fmmTreeRoot);

	//deleteList(fmmTreeRoot);

	return 0;
}
////////////////////////////////////
void CMultipole::deleteMultipole()
{
	int i;
	if(fmmTreeRoot) deleteFmmTree(fmmTreeRoot);
	if(RSmn)
	{
		for(i=0;i<=2*p;i++) delete[] RSmn[i];
		delete RSmn;
	}
	delete[] arryLeafPoint;
	delete[] x;
//	delete[] x0;
}

//判断CNode是否在树节点里，是返回1，否返回0。
inline int CMultipole::isInCube(Node* pNode,TREECUBE* pCube)
{
//	double x=fabs(pNode->D3Center.x-pCube->center.x);
//	if(fabs(pNode->D3Center.x-pCube->center.x) <= pCube->halfD && 
//	   fabs(pNode->D3Center.y-pCube->center.y) <= pCube->halfD &&
//     fabs(pNode->D3Center.z-pCube->center.z) <= pCube->halfD)
// 	if(fabs(pNode->coord[0]-pCube->center.x) < pCube->halfD+TOLERANCE && 
// 	   fabs(pNode->coord[1]-pCube->center.y) < pCube->halfD+TOLERANCE &&
//       fabs(pNode->coord[2]-pCube->center.z) < pCube->halfD+TOLERANCE)
	if((pNode->coord[0]-pCube->center.x <= pCube->halfD)&&(pCube->center.x-pNode->coord[0]< pCube->halfD) && 
		(pNode->coord[1]-pCube->center.y <= pCube->halfD)&&(pCube->center.y-pNode->coord[1]< pCube->halfD) &&
		(pNode->coord[2]-pCube->center.z <= pCube->halfD)&&(pCube->center.z-pNode->coord[2]< pCube->halfD))
		return 1;
	else
		return 0;
}

int CMultipole::initFmmTreeRoot()
{
	if(fmmTreeRoot==0)
	{
		ErrDisplay("初始化根节点失败");
		return 1;
	}

	fmmTreeRoot->level=1;
    fmmTreeRoot->father=0;
    fmmTreeRoot->cube.center.x=rootCube.center.x;
    fmmTreeRoot->cube.center.y=rootCube.center.y;
    fmmTreeRoot->cube.center.z=rootCube.center.z;
    fmmTreeRoot->cube.halfD=rootCube.halfD;
    fmmTreeRoot->tlNode=tlNode;

    NODELIST* tempt;
    fmmTreeRoot->listNodes=new NODELIST;
    fmmTreeRoot->listNodes->nodeNo_=0;//从0开始
    fmmTreeRoot->listNodes->next=0;

	if(!isInCube(&(pDomain->m_pElmCen[0]), &(fmmTreeRoot->cube)))
	{
		ErrDisplay("包络方块太小");
		return 1;
	}

	tempt=fmmTreeRoot->listNodes;
    for(long i=1; i<tlNode; i++)
	{
	    tempt->next=new NODELIST;
	    tempt->next->nodeNo_=i;
	    tempt=tempt->next;
	    tempt->next=0;
		if(!isInCube(&(pDomain->m_pElmCen[i]),&(fmmTreeRoot->cube)))
		{
			ErrDisplay("包络方块太小");
			return 1;
		}
	}

	return 0;
}

void CMultipole::initTreeNode(FMMTREENODE* fmmTreeNode)
{
	fmmTreeNode->tlNode=0; //number of nodes included within this treenode


	fmmTreeNode->listNeighbors=0;
	fmmTreeNode->listInteract=0;

	fmmTreeNode->Mmn=0;//远端系数
	fmmTreeNode->Lmn=0;//近端系数

	fmmTreeNode->tlLocalNode=0;//for leaf treenode only
	fmmTreeNode->sequence=0;
}

int CMultipole::createFmmTree(FMMTREENODE* fmmTreeNode)
{
	//fmmTreeNode:树的节点
	//tlFmmTreeNode:总的节点数，初值赋0
	//MAX_nodes:叶子节点含有的最大单元数

	int i;
	if(fmmTreeNode==0)
	{
		ErrDisplay("初始化根节点失败");
		return 1;
	}

	tlFmmTreeNode++;

	if(fmmTreeNode->tlNode<=MAX_nodes /*&& fmmTreeNode->level>=4*/)//是叶子节点
	//if(fmmTreeNode->tlNode<=MAX_nodes && fmmTreeNode->level>=5)
	{
		for(i=0; i<8; i++) fmmTreeNode->child[i]=0;
		fmmTreeNode->listNeighbors=0;
		fmmTreeNode->listInteract=0;
		return 0;
	}

	for(i=0; i<8; i++) fmmTreeNode->child[i]=new FMMTREENODE;

	double halfD=fmmTreeNode->cube.halfD/2.0;
	int level;
	level = fmmTreeNode->level;
	for(i=0; i<8; i++) fmmTreeNode->child[i]->cube.halfD=halfD;

	fmmTreeNode->child[0]->cube.center.x=fmmTreeNode->cube.center.x+halfD;
	fmmTreeNode->child[3]->cube.center.x=fmmTreeNode->cube.center.x+halfD;
	fmmTreeNode->child[4]->cube.center.x=fmmTreeNode->cube.center.x+halfD;
	fmmTreeNode->child[7]->cube.center.x=fmmTreeNode->cube.center.x+halfD;
	fmmTreeNode->child[1]->cube.center.x=fmmTreeNode->cube.center.x-halfD;
	fmmTreeNode->child[2]->cube.center.x=fmmTreeNode->cube.center.x-halfD;
	fmmTreeNode->child[5]->cube.center.x=fmmTreeNode->cube.center.x-halfD;
	fmmTreeNode->child[6]->cube.center.x=fmmTreeNode->cube.center.x-halfD;

	fmmTreeNode->child[2]->cube.center.y=fmmTreeNode->cube.center.y-halfD;
	fmmTreeNode->child[3]->cube.center.y=fmmTreeNode->cube.center.y-halfD;
	fmmTreeNode->child[6]->cube.center.y=fmmTreeNode->cube.center.y-halfD;
	fmmTreeNode->child[7]->cube.center.y=fmmTreeNode->cube.center.y-halfD;
	fmmTreeNode->child[0]->cube.center.y=fmmTreeNode->cube.center.y+halfD;
	fmmTreeNode->child[1]->cube.center.y=fmmTreeNode->cube.center.y+halfD;
	fmmTreeNode->child[4]->cube.center.y=fmmTreeNode->cube.center.y+halfD;
	fmmTreeNode->child[5]->cube.center.y=fmmTreeNode->cube.center.y+halfD;

	fmmTreeNode->child[0]->cube.center.z=fmmTreeNode->cube.center.z+halfD;
	fmmTreeNode->child[1]->cube.center.z=fmmTreeNode->cube.center.z+halfD;
	fmmTreeNode->child[2]->cube.center.z=fmmTreeNode->cube.center.z+halfD;
	fmmTreeNode->child[3]->cube.center.z=fmmTreeNode->cube.center.z+halfD;
	fmmTreeNode->child[4]->cube.center.z=fmmTreeNode->cube.center.z-halfD;
	fmmTreeNode->child[5]->cube.center.z=fmmTreeNode->cube.center.z-halfD;
	fmmTreeNode->child[6]->cube.center.z=fmmTreeNode->cube.center.z-halfD;
	fmmTreeNode->child[7]->cube.center.z=fmmTreeNode->cube.center.z-halfD;
	/*
	for(i=0; i<8; i++) fmmTreeNode->child[i]->level=fmmTreeNode->level+1;

	for(i=0; i<8; i++) fmmTreeNode->child[i]->listNodes=0;

	for(i=0; i<8; i++) fmmTreeNode->child[i]->father=fmmTreeNode;

	for(i=0; i<8; i++) fmmTreeNode->child[i]->tlNode=0;*/
	for(i=0; i<8; i++)
	{
		fmmTreeNode->child[i]->level=fmmTreeNode->level+1;
		fmmTreeNode->child[i]->listNodes=0;
		fmmTreeNode->child[i]->father=fmmTreeNode;
		fmmTreeNode->child[i]->tlNode=0;

		fmmTreeNode->child[i]->Mmn = NULL;
		fmmTreeNode->child[i]->Lmn = NULL;
	}

	NODELIST* tempt;
	while(fmmTreeNode->listNodes)
	{
		bool bIn;
		bIn=false;
		for(i=0; i<8; i++)
		{
			if(isInCube(&(pDomain->m_pElmCen[fmmTreeNode->listNodes->nodeNo_]),&(fmmTreeNode->child[i]->cube)))
			{
				if(fmmTreeNode->child[i]->listNodes)
				{
					tempt=fmmTreeNode->child[i]->listNodes;
					fmmTreeNode->child[i]->listNodes=fmmTreeNode->listNodes;
					fmmTreeNode->listNodes=fmmTreeNode->listNodes->next;
					fmmTreeNode->child[i]->listNodes->next=tempt;
					fmmTreeNode->child[i]->tlNode++;
				}
				else
				{
					fmmTreeNode->child[i]->listNodes=fmmTreeNode->listNodes;
					fmmTreeNode->listNodes=fmmTreeNode->listNodes->next;
					fmmTreeNode->child[i]->listNodes->next=0;
					fmmTreeNode->child[i]->tlNode++;
				}
				bIn=true;
				break;
			}
		}
		if(bIn==false)
		{
			ErrDisplay("包络方盒尺寸小");
			return 1;
		}
	}
	tempt=0;

	//annotated by xzf
	
	for(i=0; i<8; i++)
	{
		if(fmmTreeNode->child[i]->listNodes==0 /*&& fmmTreeNode->level>=6*/)
		{
			delete fmmTreeNode->child[i];
			fmmTreeNode->child[i]=0;
		}
	}

    fmmTreeNode->listNeighbors=0;
    fmmTreeNode->listInteract=0;

	for(i=0; i<8; i++)
		if(fmmTreeNode->child[i])
			if(createFmmTree(fmmTreeNode->child[i])) return 1;

	return 0;
}

void CMultipole::createLeafPoint(FMMTREENODE* fmmTreeNode)
{
	if(!arryLeafPoint)
		arryLeafPoint=new FMMTREENODE*[tlLeafguess];

	//arryLeafPoint:指向所有叶子节点的指针
	if(fmmTreeNode->child[0]!=0||fmmTreeNode->child[1]!=0
		||fmmTreeNode->child[2]!=0||fmmTreeNode->child[3]!=0
		||fmmTreeNode->child[4]!=0||fmmTreeNode->child[5]!=0
		||fmmTreeNode->child[6]!=0||fmmTreeNode->child[7]!=0)
	{
		for(int i=0; i<8; i++)
			if(fmmTreeNode->child[i]!=0) createLeafPoint(fmmTreeNode->child[i]);
	}
	else
	{
		arryLeafPoint[tlLeaf]=fmmTreeNode;
		tlLeaf++;

		if(tlLeaf > tlLeafguess)
		{
			spdlog::info("Error: the allocated size for leaves is not enough!\n");
			exit(0);
		}
	}
}

void CMultipole::countLeafPoint(FMMTREENODE* fmmTreeNode)
{
	//arryLeafPoint:指向所有叶子节点的指针
	if(fmmTreeNode->child[0]!=0||fmmTreeNode->child[1]!=0
		||fmmTreeNode->child[2]!=0||fmmTreeNode->child[3]!=0
		||fmmTreeNode->child[4]!=0||fmmTreeNode->child[5]!=0
		||fmmTreeNode->child[6]!=0||fmmTreeNode->child[7]!=0)
	{
		for(int i=0; i<8; i++)
			if(fmmTreeNode->child[i]!=0) countLeafPoint(fmmTreeNode->child[i]);
	}
	else
	{
		//arryLeafPoint[tlLeaf]=fmmTreeNode;
		tlLeafguess++;

// 		if(tlLeaf > tlLeafguess)
// 		{
// 			spdlog::info("Error: the allocated size for leaves is not enough!\n");
// 			exit(0);
// 		}
	}
}

void CMultipole::checkLeafNods()
{
	int i, cnt = 0;
	for (i=0; i<tlLeaf; i++)
	{
		cnt += arryLeafPoint[i]->tlNode;
	}
	spdlog::info("Total elements: {}\n", cnt);
}

void CMultipole::getTreeNodeElms(int nNeig, FMMTREENODE** nodeNeigs, int *npt, double **pts, int *nelm, int **elms)
{
	int i, j, eidx, pidx, tlTreeNods, *pTreeNods = NULL;
	std::set<int> stp;
	std::set<int>::iterator siter;
	pTreeNods = new int[tlNode];

#if 0
	tlTreeNodePts = 0;

	for(i=0; i<nNeig; i++)
		getNodePoint(nodeNeigs[i]);

	for (i=0; i<tlTreeNodePts; i++)
	{
		eidx = m_pTreeNodePts[i];
		for (j=0; j<DIM3; j++)
			stp.insert(pDomain->m_pFEle[eidx*DIM3 + j]);
	}
#else
	tlTreeNods = 0;
	for(i=0; i<nNeig; i++)
		getNodePoint(nodeNeigs[i], &tlTreeNods, pTreeNods);

	for (i=0; i<tlTreeNods; i++)
	{
		eidx = pTreeNods[i];
		(*elms)[i] = eidx;
	}
	*nelm = tlTreeNods;
#endif
#if 0
	*npt = stp.size();
	*nelm = tlTreeNods;
	//*nelm = tlTreeNodePts;
	//*pts = new double[(*npt)*DIM3];
	//*elms = new int [(*nelm)*DIM3];

	pidx = 0;	//begin form 1
	siter = stp.begin();
	while(siter != stp.end())
	{
		g2l[*siter] = pidx;
		l2g[pidx] = *siter;

		(*pts)[pidx*DIM3 + 0] = pDomain->m_pFNod[*siter].coord[0];
		(*pts)[pidx*DIM3 + 1] = pDomain->m_pFNod[*siter].coord[1];
		(*pts)[pidx*DIM3 + 2] = pDomain->m_pFNod[*siter].coord[2];

		++pidx;
		++siter;
	}

	for (i=0; i<tlTreeNods/*tlTreeNodePts*/; i++)
	{
		//eidx = m_pTreeNodePts[i];
		eidx = pTreeNods[i];
		el2g[i] = eidx;
		eg2l[eidx] = i;

		(*elms)[i*DIM3 + 0] = g2l[pDomain->m_pFEle[eidx*DIM3 + 0]]+1;
		(*elms)[i*DIM3 + 1] = g2l[pDomain->m_pFEle[eidx*DIM3 + 1]]+1;
		(*elms)[i*DIM3 + 2] = g2l[pDomain->m_pFEle[eidx*DIM3 + 2]]+1;
	}
#endif

	if(pTreeNods)
	{
		delete []pTreeNods;
		pTreeNods = NULL;
	}
}

void CMultipole::getNodePoint(FMMTREENODE* fmmTreeNode)
{
	if(fmmTreeNode->child[0]!=0||fmmTreeNode->child[1]!=0
		||fmmTreeNode->child[2]!=0||fmmTreeNode->child[3]!=0
		||fmmTreeNode->child[4]!=0||fmmTreeNode->child[5]!=0
		||fmmTreeNode->child[6]!=0||fmmTreeNode->child[7]!=0)
	{
		for(int i=0; i<8; i++)
			if(fmmTreeNode->child[i]!=0) getNodePoint(fmmTreeNode->child[i]);
	}
	else
	{
		NODELIST* tempt = fmmTreeNode->listNodes;
		while(tempt)
		{
			m_pTreeNodePts[tlTreeNodePts] = tempt->nodeNo_;
			tlTreeNodePts++;
			tempt = tempt->next;
		}
	}
}

void CMultipole::getNodePoint(FMMTREENODE* fmmTreeNode, int *tlTreeNods, int *pTreeNods)
{
	if(fmmTreeNode->child[0]!=0||fmmTreeNode->child[1]!=0
		||fmmTreeNode->child[2]!=0||fmmTreeNode->child[3]!=0
		||fmmTreeNode->child[4]!=0||fmmTreeNode->child[5]!=0
		||fmmTreeNode->child[6]!=0||fmmTreeNode->child[7]!=0)
	{
		for(int i=0; i<8; i++)
			if(fmmTreeNode->child[i]!=0) getNodePoint(fmmTreeNode->child[i], tlTreeNods, pTreeNods);
	}
	else
	{
		NODELIST* tempt = fmmTreeNode->listNodes;
		while(tempt)
		{
			pTreeNods[*tlTreeNods] = tempt->nodeNo_;
			(*tlTreeNods)++;
			tempt = tempt->next;
		}
	}
}

void CMultipole::resetFmmTree(FMMTREENODE* fmmTreeNode)
{
	int i, j, k;

	if(fmmTreeNode->level>=2)
	{
		for (k=0; k<DIM3; k++)
		{
			for(i=0; i<=p; i++)
				for(j=0; j<=i; j++)
				{
					fmmTreeNode->Mmn[k][i][j]=0.0;
				}

			if(fmmTreeNode->level>=3)
			{			
				for(i=0; i<=p; i++)
					for(j=0; j<=i; j++)
					{
						fmmTreeNode->Lmn[k][i][j]=0.0;
					}
			}
		}
	}
	for(i=0; i<8; i++)
		if(fmmTreeNode->child[i]!=0) resetFmmTree(fmmTreeNode->child[i]);
}

void CMultipole::deleteFmmTree(FMMTREENODE* fmmTreeNode)
{
	for(int i=0; i<8; i++)
		if(fmmTreeNode->child[i]!=0) deleteFmmTree(fmmTreeNode->child[i]);

	if(fmmTreeNode->tlNode <= max_nodes) //is a leaf
		delete[] fmmTreeNode->sequence;

	delete fmmTreeNode;
}

void CMultipole::deleteList(FMMTREENODE* fmmTreeNode)
{
	for(int i=0; i<8; i++)
		if(fmmTreeNode->child[i]!=0) deleteList(fmmTreeNode->child[i]);

	NODELIST* tempt;
	while(fmmTreeNode->listNodes)
	{
		tempt=fmmTreeNode->listNodes;
		fmmTreeNode->listNodes=fmmTreeNode->listNodes->next;
		delete tempt;
	}

	CUBENODELIST* neighbor;
	while(fmmTreeNode->listNeighbors)
	{
		neighbor=fmmTreeNode->listNeighbors;
		fmmTreeNode->listNeighbors=fmmTreeNode->listNeighbors->next;
		delete neighbor;
	}
}

void CMultipole::allocateMoments(FMMTREENODE* fmmTreeNode)
{
	int i, j;
	if(fmmTreeNode->level>=2)
	{
		fmmTreeNode->Mmn=new  complex<double>**[DIM3];
		for (j=0; j<DIM3; j++)
		{
			fmmTreeNode->Mmn[j]=new  complex<double>*[p+1];
			for(i=0;i<=p;i++) fmmTreeNode->Mmn[j][i]=new  complex<double>[2*i+1];
		}

		if(fmmTreeNode->level>=3)
		{
			fmmTreeNode->Lmn=new  complex<double>**[DIM3];
			for (j=0; j<DIM3; j++)
			{
				fmmTreeNode->Lmn[j]=new  complex<double>*[p+1];
				for(i=0;i<=p;i++) fmmTreeNode->Lmn[j][i]=new  complex<double>[2*i+1];
			}
		}
		
	}

	for(i=0; i<8; i++)
		if(fmmTreeNode->child[i])
			allocateMoments(fmmTreeNode->child[i]);
}

void CMultipole::freeMoments(FMMTREENODE* fmmTreeNode)
{
	int i, j;
	if(fmmTreeNode->level>=2)
	{
		for (j=0; j<DIM3; j++)
		{
			for(i=0;i<=p;i++) 
				delete[] fmmTreeNode->Mmn[j][i];

			delete [] fmmTreeNode->Mmn[j];
		}
		delete[] fmmTreeNode->Mmn;

		if(fmmTreeNode->level>=3)
		{
			for (j=0; j<DIM3; j++)
			{
				for(i=0;i<=p;i++) 
					delete[] fmmTreeNode->Lmn[j][i];

				delete [] fmmTreeNode->Lmn[j];
			}
			delete[] fmmTreeNode->Lmn;
		}
	}

	for(i=0; i<8; i++)
		if(fmmTreeNode->child[i])
			freeMoments(fmmTreeNode->child[i]);
}

void CMultipole::addNeighbor(FMMTREENODE* fmmTreeNode1,FMMTREENODE* fmmTreeNode2)
{
	if(fmmTreeNode2==0)	return;

	CUBENODELIST* tempt;
	if(fmmTreeNode1->listNeighbors==0)
	{
		fmmTreeNode1->listNeighbors=new CUBENODELIST;
		fmmTreeNode1->listNeighbors->cubeNode=fmmTreeNode2;
		fmmTreeNode1->listNeighbors->next=0;
	}
	else
	{
		tempt=fmmTreeNode1->listNeighbors;//push process
		fmmTreeNode1->listNeighbors=new CUBENODELIST;
		fmmTreeNode1->listNeighbors->cubeNode=fmmTreeNode2;
		fmmTreeNode1->listNeighbors->next=tempt;
	}
}

void CMultipole::addInteraction(FMMTREENODE* fmmTreeNode1,FMMTREENODE* fmmTreeNode2)
{
	if(fmmTreeNode2==0)	return;

	CUBENODELIST* tempt;
	if(fmmTreeNode1->listInteract==0)
	{
		fmmTreeNode1->listInteract=new CUBENODELIST;
		fmmTreeNode1->listInteract->cubeNode=fmmTreeNode2;
		fmmTreeNode1->listInteract->next=0;
	}
	else
	{
		tempt=fmmTreeNode1->listInteract;//push process
		fmmTreeNode1->listInteract=new CUBENODELIST;
		fmmTreeNode1->listInteract->cubeNode=fmmTreeNode2;
		fmmTreeNode1->listInteract->next=tempt;
	}
}

int CMultipole::countInteraction(FMMTREENODE* fmmTreeNode)
{
	int cnt=0;

	if(fmmTreeNode==0 || fmmTreeNode->listInteract==0)	
		return 0;

	CUBENODELIST* tempt = fmmTreeNode->listInteract;
	
	while(tempt)
	{
		cnt += tempt->cubeNode->tlNode;
		tempt = tempt->next;
	}

	return cnt;
}

int CMultipole::isNeighbor(FMMTREENODE* fmmTreeNode1,FMMTREENODE* fmmTreeNode2)
{
#if 0
	int i, j;
	double x[2][8], y[2][8], z[2][8], r[2], dis;
	D3POINT ct[2];
	bool bNeig;

	ct[0] = fmmTreeNode1->cube.center;
	ct[1] = fmmTreeNode2->cube.center;
	r[0] = fmmTreeNode1->cube.halfD;
	r[1] = fmmTreeNode2->cube.halfD;

	for (i=0;i<2; i++)
	{
		x[i][0] = ct[i].x-r[i];
		y[i][0] = ct[i].y-r[i];
		z[i][0] = ct[i].z-r[i];
		x[i][1] = ct[i].x+r[i];
		y[i][1] = ct[i].y-r[i];
		z[i][1] = ct[i].z-r[i];
		x[i][2] = ct[i].x+r[i];
		y[i][2] = ct[i].y+r[i];
		z[i][2] = ct[i].z-r[i];
		x[i][3] = ct[i].x-r[i];
		y[i][3] = ct[i].y+r[i];
		z[i][3] = ct[i].z-r[i];
		x[i][4] = ct[i].x-r[i];
		y[i][4] = ct[i].y-r[i];
		z[i][4] = ct[i].z+r[i];
		x[i][5] = ct[i].x+r[i];
		y[i][5] = ct[i].y-r[i];
		z[i][5] = ct[i].z+r[i];
		x[i][6] = ct[i].x+r[i];
		y[i][6] = ct[i].y+r[i];
		z[i][6] = ct[i].z+r[i];
		x[i][7] = ct[i].x-r[i];
		y[i][7] = ct[i].y+r[i];
		z[i][7] = ct[i].z+r[i];
	}
	
	bNeig = false;
	for (i=0; i<8; i++)
	{
		for (j=0; j<8; j++)
		{
			dis = sqrt((x[0][i]-x[1][j])*(x[0][i]-x[1][j])+
						(y[0][i]-y[1][j])*(y[0][i]-y[1][j])
						+(z[0][i]-z[1][j])*(z[0][i]-z[1][j]));
			if (dis<1.0*10e-6)
			{
				bNeig = true;
				break;
			}
		}
	}

	if(bNeig)
		return 1;
	else
		return 0;
#else
	double distance;
	distance=sqrt((fmmTreeNode1->cube.center.x-fmmTreeNode2->cube.center.x)
				 *(fmmTreeNode1->cube.center.x-fmmTreeNode2->cube.center.x)
				 +(fmmTreeNode1->cube.center.y-fmmTreeNode2->cube.center.y)
				 *(fmmTreeNode1->cube.center.y-fmmTreeNode2->cube.center.y)
				 +(fmmTreeNode1->cube.center.z-fmmTreeNode2->cube.center.z)
				 *(fmmTreeNode1->cube.center.z-fmmTreeNode2->cube.center.z));
	if(distance<=(3.6*fmmTreeNode1->cube.halfD))
		return 1;
	else
		return 0;
#endif
/*	
	double sqDist;
	sqDist=(fmmTreeNode1->cube.center.x-fmmTreeNode2->cube.center.x)
		  *(fmmTreeNode1->cube.center.x-fmmTreeNode2->cube.center.x)
		  +(fmmTreeNode1->cube.center.y-fmmTreeNode2->cube.center.y)
		  *(fmmTreeNode1->cube.center.y-fmmTreeNode2->cube.center.y)
		  +(fmmTreeNode1->cube.center.z-fmmTreeNode2->cube.center.z)
		  *(fmmTreeNode1->cube.center.z-fmmTreeNode2->cube.center.z);
	if(sqDist <= 13*fmmTreeNode1->cube.halfD*fmmTreeNode1->cube.halfD)
		return 1;
	else
		return 0;
*/
}

int CMultipole::createNeighbor(FMMTREENODE* fmmTreeNode)
{
	int i;

	if(fmmTreeNode==0)
	{
		ErrDisplay("空节点");
		return 1;
	}

	CUBENODELIST* tempt;
	if(fmmTreeNode->father)
	{
		tempt=fmmTreeNode->father->listNeighbors;
		while(tempt)
		{
			if(tempt->cubeNode->tlNode <= MAX_nodes)  //****叶子节点****
				addNeighbor(fmmTreeNode,tempt->cubeNode);
			else
			{
				for(i=0; i<8; i++)
				{
					if(tempt->cubeNode->child[i])
					{
						if(fmmTreeNode->level!=tempt->cubeNode->child[i]->level)
						{
							ErrDisplay("两个节点层次不一样");
							return 1;
						}
						if(isNeighbor(fmmTreeNode,tempt->cubeNode->child[i]))
							addNeighbor(fmmTreeNode,tempt->cubeNode->child[i]);
					}
				}
			}
			tempt=tempt->next;
		}
	}

	for(i=0; i<8; i++)
	{
		if(fmmTreeNode->child[i])
		{
			for(int j=0; j<8; j++)
				if(j != i && fmmTreeNode->child[j])
					addNeighbor(fmmTreeNode->child[i],fmmTreeNode->child[j]);
		}
	}

	for(i=0; i<8; i++)
	{
		if(fmmTreeNode->child[i])
			if(createNeighbor(fmmTreeNode->child[i])) 
				return 1;
	}

	return 0;
}

int CMultipole::createInteraction(FMMTREENODE* fmmTreeNode)
{
	int i;

	if(fmmTreeNode==0)
	{
		ErrDisplay("空节点");
		return 1;
	}

	if(fmmTreeNode->level>=3)
	{
		CUBENODELIST* temp;
		temp=fmmTreeNode->father->listNeighbors;
		while(temp)
		{
			for(i=0; i<8; i++)
			{
				if(temp->cubeNode->child[i])
				{
					if(!isNeighbor(fmmTreeNode,temp->cubeNode->child[i]))
						addInteraction(fmmTreeNode,temp->cubeNode->child[i]);
				}
			}
			temp=temp->next;
		}
	}

	for(i=0; i<8; i++)
	{
		if(fmmTreeNode->child[i])
			if(createInteraction(fmmTreeNode->child[i])) return 1;
	}

	return 0;
}

//inline 
void CMultipole::getExpansionR(D3POINT pDt0, D3POINT pDt)
{
	int i, j;

	complex<double> temp;
	temp.real(pDt.x-pDt0.x);
	temp.imag(pDt.y-pDt0.y);
	double x3=pDt.z-pDt0.z;
	double rr=temp.real()*temp.real()+temp.imag()*temp.imag()+x3*x3;

	RSmn[0][0]=1.0;
//Nishida
	for(i=1; i<=p; i++)
	for(j=0; j<=i; j++)
	{
		if(i == j) 
			RSmn[i][i]=temp*RSmn[i-1][i-1]/(i*2.0);
		else if(i == j+1) 
			RSmn[i][j]=x3*RSmn[j][j];
		else
			RSmn[i][j]=((2*i-1)*x3*RSmn[i-1][j]-rr*RSmn[i-2][j])/((i*i-j*j)*1.0);
	}
/*
//Yoshida
	for(i=1; i<=p; i++) RSmn[i][i]=temp*RSmn[i-1][i-1]/(i*2.0);
	for(j=0; j<p; j++)
	{
		RSmn[j+1][j]=x3*RSmn[j][j];
		for(i=j+1; i<p; i++) RSmn[i+1][j]=((2*i+1)*x3*RSmn[i][j]-rr*RSmn[i-1][j])/(i+j+1)/(i-j+1);
	}
*/
/*
				FILE *stream;
				stream = fopen( "RSmn", "w");
				if( stream != NULL )
				{
					fprintf(stream,"Dt0=%13.4g%13.4g%13.4g", pDt0.x, pDt0.y, pDt0.z);
					fprintf(stream,"     Dtx=%13.4g%13.4g%13.4g\n", pDt.x, pDt.y, pDt.z);
					fprintf(stream,"X:\n");
					for(i=0; i<=2*p; i++)
					{
						fprintf(stream,"%8d", i);
						for(j=0; j<=i; j++)	fprintf(stream,"%13.4g", RSmn[i][j].x);
						fprintf(stream,"\n");
					}
					fprintf(stream,"Y:\n");
					for(i=0; i<=2*p; i++)
					{
						fprintf(stream,"%8d", i);
						for(j=0; j<=i; j++)	fprintf(stream,"%13.4g", RSmn[i][j].y);
						fprintf(stream,"\n");
					}
				}
				fclose(stream);
*/
}

//inline 
void CMultipole::getExpansionR_Local(D3POINT pDt0, D3POINT pDt, complex<double> RSmn0[21][21])
{
	int i, j;

	complex<double> temp;
	temp.real(pDt.x-pDt0.x);
	temp.imag(pDt.y-pDt0.y);
	double x3=pDt.z-pDt0.z;
	double rr=temp.real()*temp.real()+temp.imag()*temp.imag()+x3*x3;

	RSmn0[0][0]=1.0;
#if 0
//Nishida
	for(i=1; i<=p; i++)
	{
		for(j=0; j<=i; j++)
		{
			if(i == j) 
				RSmn0[i][i]=temp*RSmn0[i-1][i-1]/(i*2.0);
			else if(i == j+1) 
				RSmn0[i][j]=x3*RSmn0[j][j];
			else
				RSmn0[i][j]=((2*i-1)*x3*RSmn0[i-1][j]-rr*RSmn0[i-2][j])/((i*i-j*j)*1.0);
		}
	}
#else
	//Yoshida
	for(i=1; i<=2*p; i++) RSmn0[i][i]=temp*RSmn0[i-1][i-1]/(i*2.0);
	for(j=0; j<2*p; j++)
	{
		RSmn0[j+1][j]=x3*RSmn0[j][j];
		for(i=j+1; i<2*p; i++) RSmn0[i+1][j]=((2*i+1)*x3*RSmn0[i][j]-rr*RSmn0[i-1][j])/((i+j+1)*(i-j+1)*1.0);
	}
#endif
}

inline void CMultipole::getExpansionS(D3POINT pDt0, D3POINT pDt)
{
	int i, j;

	complex<double> temp;
	temp.real(pDt.x-pDt0.x);
	temp.imag(pDt.y-pDt0.y);
	double x3=pDt.z-pDt0.z;
	double rr=temp.real()*temp.real()+temp.imag()*temp.imag()+x3*x3;

	RSmn[0][0]=1.0/sqrt(rr);

//Nishida
	for(i=1; i<=p*2; i++)
	{
		for(j=0; j<=i; j++)
		{
			if(i == j) 
				RSmn[i][i]=(i*2-1)*1.0*temp*RSmn[i-1][i-1]/rr;
			else if(i == j+1) 
				RSmn[i][j]=(2*j+1)*x3*RSmn[j][j]/rr;
			else 
				RSmn[i][j]=((2*i-1)*x3*RSmn[i-1][j]-((i-1)*(i-1)-j*j)*1.0*RSmn[i-2][j])/rr;
		}
	}
/*
//Yoshida
	for(i=1; i<=2*p; i++) RSmn[i][i]=temp*RSmn[i-1][i-1]*(2*i-1.0)/rr;
	for(j=0; j<2*p; j++)
	{
		RSmn[j+1][j]=(2*j+1)*x3*RSmn[j][j]/rr;
		for(i=j+1; i<2*p; i++) RSmn[i+1][j]=((2*i+1)*x3*RSmn[i][j]-(i+j)*(i-j)*RSmn[i-1][j])/rr;
	}
*/
/*//
				FILE *stream;
				stream = fopen( "RSmn", "w");
				if( stream != NULL )
				{
					fprintf(stream,"Dt0=%13.4g%13.4g%13.4g", pDt0.x, pDt0.y, pDt0.z);
					fprintf(stream,"     Dtx=%13.4g%13.4g%13.4g\n", pDt.x, pDt.y, pDt.z);
					fprintf(stream,"X:\n");
					for(i=0; i<=2*p; i++)
					{
						fprintf(stream,"%8d", i);
						for(j=0; j<=i; j++)	fprintf(stream,"%13.4g", RSmn[i][j].x);
						fprintf(stream,"\n");
					}
					fprintf(stream,"Y:\n");
					for(i=0; i<=2*p; i++)
					{
						fprintf(stream,"%8d", i);
						for(j=0; j<=i; j++)	fprintf(stream,"%13.4g", RSmn[i][j].y);
						fprintf(stream,"\n");
					}
				}
				fclose(stream);
*///
}

inline void CMultipole::getExpansionS_Local(D3POINT pDt0, D3POINT pDt, complex<double> RSmn0[21][21])
{
	int i, j;

	complex<double> temp;
	temp.real(pDt.x-pDt0.x);
	temp.imag(pDt.y-pDt0.y);
	double x3=pDt.z-pDt0.z;
	double rr=temp.real()*temp.real()+temp.imag()*temp.imag()+x3*x3;

	RSmn0[0][0]=1.0/sqrt(rr);

#if 1
	//Nishida
	for(i=1; i<=p*2; i++)
	{
		for(j=0; j<=i; j++)
		{
			if(i == j) 
				RSmn0[i][i]=(i*2-1)*1.0*temp*RSmn0[i-1][i-1]/rr;
			else if(i == j+1) 
				RSmn0[i][j]=(2*j+1)*x3*RSmn0[j][j]/rr;
			else 
				RSmn0[i][j]=((2*i-1)*x3*RSmn0[i-1][j]-((i-1)*(i-1)-j*j)*1.0*RSmn0[i-2][j])/rr;
		}
	}
#else
	//Yoshida
	for(i=1; i<=2*p; i++) RSmn0[i][i]=temp*RSmn0[i-1][i-1]*(2*i-1.0)/rr;
	for(j=0; j<2*p; j++)
	{
		RSmn0[j+1][j]=(2*j+1)*x3*RSmn0[j][j]/rr;
		for(i=j+1; i<2*p; i++) RSmn0[i+1][j]=((2*i+1)*x3*RSmn0[i][j]-(i+j)*(i-j)*1.0*RSmn0[i-1][j])/rr;
	}
#endif
}

void CMultipole::PrintRSmn(complex<double> RSmn0[21][21], char* filename)
{
	int i, j;
	
	FILE *fout = NULL;
	if ((fout =fopen(filename,"w")) == NULL )
	{
		spdlog::info("Error: can not open file %s!\n");
		return;
	}

	for (i=0; i<=p*2; i++)
	{
		for (j=0; j<=i; j++)
		{
			if(std::imag(RSmn0[i][j]) >= 0)
				fprintf(fout, "%18.16e+%25.16ei ", std::real(RSmn0[i][j]), std::imag(RSmn0[i][j]));
			else
				fprintf(fout, "%18.16e%25.16ei ", std::real(RSmn0[i][j]), std::imag(RSmn0[i][j]));
		}
		fprintf(fout, "\n");
	}
	fclose(fout);
	fout = NULL;
}

//记录叶子节点包含的所有局部单元（所属单元及邻居单元）
int CMultipole::setLeafSequence()
{
    //arryLeafPoint:指向所有叶子节点的指针
	long i, j;
	NODELIST* tempNodeList;
	CUBENODELIST* tempNeibList;

	max_nodes=0;
	min_nodes=MAX_nodes;
	max_level=0;
	min_level=10000000;
	long totalNodes=0L;
	for(i=0; i<tlLeaf; i++)
	{
		totalNodes+=arryLeafPoint[i]->tlNode;
		if(max_nodes < arryLeafPoint[i]->tlNode) max_nodes=arryLeafPoint[i]->tlNode;
		if(min_nodes > arryLeafPoint[i]->tlNode) min_nodes=arryLeafPoint[i]->tlNode;
		if(max_level < arryLeafPoint[i]->level) max_level=arryLeafPoint[i]->level;
		if(min_level > arryLeafPoint[i]->level) min_level=arryLeafPoint[i]->level;
	}
	if(totalNodes != tlNode)
	{
		ErrDisplay("totalNodes != tlNode");
		return 1;
	}
	tempSequence=new long[tlNode];

	max_nearNodes=0;
	min_nearNodes=MAX_nodes*30;//27
	for(i=0; i<tlLeaf; i++)
	{
		//check leaf level
		if(arryLeafPoint[i]->level<3)
		{
			ErrDisplay("arryLeafPoint[i]->level<3");
			return 1;
		}
		//set sequence
		tempNodeList=arryLeafPoint[i]->listNodes;
		tempNum=0;
		while(tempNodeList)
		{
			tempSequence[tempNum]=tempNodeList->nodeNo_;
			tempNum++;
			tempNodeList=tempNodeList->next;
		}
		tempNeibList=arryLeafPoint[i]->listNeighbors;
		while(tempNeibList)
		{
			addSequence(tempNeibList->cubeNode);
			tempNeibList=tempNeibList->next;
		}
		arryLeafPoint[i]->tlLocalNode=tempNum;

		arryLeafPoint[i]->sequence = new int[arryLeafPoint[i]->tlLocalNode];
		for(j=0; j<arryLeafPoint[i]->tlLocalNode; j++) 
			arryLeafPoint[i]->sequence[j] =tempSequence[j];

		if(max_nearNodes < tempNum) max_nearNodes=tempNum;
		if(min_nearNodes > tempNum) min_nearNodes=tempNum;
	}

#if 1
	printf("min_level: %d, max_level:%d\nmin_nearnods:%d, max_nearnodes:%d\n", min_level, max_level, min_nearNodes, max_nearNodes);
#endif

	delete tempSequence;

	return 0;
}

void CMultipole::addSequence(FMMTREENODE* fmmTreeNode)
{
	int i;

	if(fmmTreeNode->tlNode <= max_nodes)//叶子节点
	{
		NODELIST* tempNodeList=fmmTreeNode->listNodes;
		while(tempNodeList)
		{
			tempSequence[tempNum]=tempNodeList->nodeNo_;
			tempNum++;
			tempNodeList=tempNodeList->next;
		}
	}
	else
	{
		for(i=0; i<8; i++)
		{
			if(fmmTreeNode->child[i])
				addSequence(fmmTreeNode->child[i]);
		}
	}
}

void CMultipole::allocatePreCOND()
{
	int i, j;

	for(i=0; i<tlLeaf; i++)
	{
		arryLeafPoint[i]->M = new double*[arryLeafPoint[i]->tlNode];
		for(j=0; j<arryLeafPoint[i]->tlNode; j++) 
			arryLeafPoint[i]->M[j] = new double[arryLeafPoint[i]->tlNode];
	}
}

void CMultipole::freePreCOND()
{
	int i, j;

	for(i=0; i<tlLeaf; i++)
	{
		for(j=0; j<arryLeafPoint[i]->tlNode; j++) delete[] arryLeafPoint[i]->M[j];
		delete[] arryLeafPoint[i]->M;
	}
}
/*
void CMultipole::allocateLeafCoefs_constant()
{
	long k;
	int i, j;

	for(k=0; k<tlLeaf; k++)
	{
		arryLeafPoint[k]->localU = new double*[arryLeafPoint[k]->tlNode];
		for(j=0; j<arryLeafPoint[i]->tlNode; j++) 
			arryLeafPoint[k]->localU[j] =new double[arryLeafPoint[k]->tlLocalNode];

		arryLeafPoint[k]->localQ = NULL;

		arryLeafPoint[k]->localQ = new double*[arryLeafPoint[k]->tlNode];
		for(j=0; j<arryLeafPoint[k]->tlNode; j++) 
			arryLeafPoint[i]->localQ[k] =new double[arryLeafPoint[k]->tlLocalNode];

		arryLeafPoint[k]->INT_RmnU = new complex**[arryLeafPoint[k]->tlNode];
		for(i=0; i<arryLeafPoint[k]->tlNode; i++)
		{
			arryLeafPoint[k]->INT_RmnU[i]=new complex*[p+1];
			for(j=0;j<=p;j++) arryLeafPoint[k]->INT_RmnU[i][j]=new complex[j+1];
		}

		arryLeafPoint[k]->INT_RmnQ = NULL;

		arryLeafPoint[k]->INT_RmnQ= new complex**[arryLeafPoint[k]->tlNode];
		for(i=0; i<arryLeafPoint[k]->tlNode; i++)
		{
			arryLeafPoint[k]->INT_RmnQ[i]=new complex*[p+1];
			for(j=0;j<=p;j++) arryLeafPoint[k]->INT_RmnQ[i][j]=new complex[j+1];
		}

	}
}

void CMultipole::freeLeafCoefs()
{
	long k;
	int i, j;

	for(k=0; k<tlLeaf; k++)
	{
		if(arryLeafPoint[k]->localU)
		{
			for(j=0; j<arryLeafPoint[k]->tlNode; j++)
				delete[] arryLeafPoint[k]->localU[j];
			delete[] arryLeafPoint[k]->localU;
		}
		if(arryLeafPoint[k]->localQ)
		{
			for(j=0; j<arryLeafPoint[k]->tlNode; j++)
				delete[] arryLeafPoint[k]->localQ[j];
			delete[] arryLeafPoint[k]->localQ;
		}

		if(arryLeafPoint[k]->INT_RmnU)
		{
			for(i=0; i<arryLeafPoint[k]->tlNode; i++)
			{
				for(j=0; j<=p; j++)
				{
					delete[] arryLeafPoint[k]->INT_RmnU[i][j];
				}
				delete[] arryLeafPoint[k]->INT_RmnU[i];
			}
			delete[] arryLeafPoint[k]->INT_RmnU;
		}

		if(arryLeafPoint[k]->INT_RmnQ)
		{
			for(i=0; i<arryLeafPoint[k]->tlNode; i++)
			{
				for(j=0; j<=p; j++)
				{
					delete[] arryLeafPoint[k]->INT_RmnQ[i][j];
				}
				delete[] arryLeafPoint[k]->INT_RmnQ[i];
			}
			delete[] arryLeafPoint[k]->INT_RmnQ;
		}
	}
}
*/
void CMultipole::initLeafs(double* xx)
{
    //arryLeafPoint:指向所有叶子节点的指针
	long i, j, k, l;
	D3POINT scen;
/*
	for(k=0; k<tlLeaf; k++)
	{
		for(i=0; i<arryLeafPoint[k]->tlNode; i++)
		{
			scen.x = pDomain->m_pElmCen[arryLeafPoint[k]->sequence[i]].coord[0];
			scen.y = pDomain->m_pElmCen[arryLeafPoint[k]->sequence[i]].coord[1];
			scen.z = pDomain->m_pElmCen[arryLeafPoint[k]->sequence[i]].coord[2];

			getExpansionR(arryLeafPoint[k]->cube.center, scen);
			for(j=0; j<=p; j++)
			for(l=0; l<=j; l++)
				arryLeafPoint[k]->Mmn[j][l]+=RSmn[j][l]*xx[arryLeafPoint[k]->sequence[i]];
		}
	}*/
}

int CMultipole::initALeaf_PotentialAndFlux(int k, complex<double> *uu,complex<double> *qq, complex<double>*** INT_RmnU, complex<double>*** INT_RmnQ)//计算叶子树节点的多级矩远端系数--常量场
{
	//该函数使用前提是：树结构系数清零
	long i, j, l,n,m, idx, kk;

	if(!arryLeafPoint[k]->Mmn)
	{
		spdlog::info("Error: please allocate memory for moment first!\n");
		return BLM_FAILED;
	}
	for(i=0; i<arryLeafPoint[k]->tlNode; i++)//定义Mmn[j][-l]=Mmn[j][l+j]
	{
		idx = arryLeafPoint[k]->sequence[i];

		for(kk=0; kk<DIM3; kk++)
		{
			for(j=0; j<=p; j++)
			{
				for(l=0; l<=j; l++)
				{
	#if 0
					arryLeafPoint[k]->Mmn[kk][j][l]+=conj(pDomain->INT_RmnU[i][j][l])*real(qq[idx*DIM3+kk]);
					arryLeafPoint[k]->Mmn[kk][j][l]-=pDomain->INT_RmnQ[i][j][l]*real(uu[idx*DIM3+kk]);
	#else
					arryLeafPoint[k]->Mmn[kk][j][l]+=/*pDomain->*/INT_RmnU[i][j][l]*real(qq[idx*DIM3+kk]);
					arryLeafPoint[k]->Mmn[kk][j][l]-=/*pDomain->*/INT_RmnQ[i][j][l]*real(uu[idx*DIM3+kk]);
	#endif
	#if 0
					if(l!=0)
					{
						if(l%2)
						{						   
							arryLeafPoint[k]->Mmn[kk][j][l+j]-=conj(pDomain->INT_RmnU[i][j][l])*qq[idx*DIM3+kk];
							arryLeafPoint[k]->Mmn[kk][j][l+j]+=conj(pDomain->INT_RmnQ[i][j][l])*uu[idx*DIM3+kk];

						}
						else
						{
							arryLeafPoint[k]->Mmn[kk][j][l+j]+=conj(pDomain->INT_RmnU[i][j][l])*qq[idx*DIM3+kk];
							arryLeafPoint[k]->Mmn[kk][j][l+j]-=conj(pDomain->INT_RmnQ[i][j][l])*uu[idx*DIM3+kk];//求刘轶军公式3.53
						}
					}
	#endif
				}
			}	//each term
		}
	}	//each node

	return BLM_SUCCESS;
}

/*    this subroutine is also right
void CMultipole::SonToFatherTrans(FMMTREENODE* fmmTreeNode)
{
	if(fmmTreeNode==0||fmmTreeNode->father==0||fmmTreeNode->father->level<2)
		return;

	int i, j, n, m;

	getExpansionR(fmmTreeNode->father->cube.center, fmmTreeNode->cube.center);
	for(i=0; i<=p; i++)
	for(j=0; j<=i; j++)
	{
		for(n=0; n<=i; n++)
		{
			complex tempMmn;
			for(m=0; m<=n; m++)
			{
				if(abs(m-j) > i-n) continue;
				if(j-m<0)
				{
					if((m-j) % 2)
					{
						tempMmn=-1.0*conjugate(fmmTreeNode->Mmn[i-n][m-j]);
					}
					else tempMmn=conjugate(fmmTreeNode->Mmn[i-n][m-j]);
				}
				else tempMmn=fmmTreeNode->Mmn[i-n][j-m];

				fmmTreeNode->father->Mmn[i][j]+=RSmn[n][m]*tempMmn;
			}
			for(m=1; m<=n; m++)
			{
				if( m+j > i-n ) continue;
				if(m % 2)
					fmmTreeNode->father->Mmn[i][j]-=conjugate(RSmn[n][m])*fmmTreeNode->Mmn[i-n][j+m];
				else
					fmmTreeNode->father->Mmn[i][j]+=conjugate(RSmn[n][m])*fmmTreeNode->Mmn[i-n][j+m];
			}
		}
	}
}
*/
void CMultipole::SonToFatherTrans(FMMTREENODE* fmmTreeNode)		//M2M
{
	if(fmmTreeNode==0||fmmTreeNode->father==0||fmmTreeNode->father->level<2)
		return;

	int i, j, n, m, k;
	complex<double> tempMmn, tempRmn;

	getExpansionR(fmmTreeNode->father->cube.center, fmmTreeNode->cube.center);
	for(k=0; k<DIM3; k++)
	{
		for(i=0; i<=p; i++)
		for(j=0; j<=i; j++)
		{
			for(n=0; n<=i; n++)
			{
				for(m=-n; m<=n; m++)
				{
					if(abs(j-m) > i-n) continue;

					if(j-m<0)
					{
						if((m-j) % 2)
						{
							tempMmn=-1.0*conj(fmmTreeNode->Mmn[k][i-n][m-j]);	//R[n][-m]=R[n][m]*(-1)的m次方
						}
						else tempMmn=conj(fmmTreeNode->Mmn[k][i-n][m-j]);
					}
					else tempMmn=fmmTreeNode->Mmn[k][i-n][j-m];

					if(m<0)
					{
						if((-m) % 2)
							tempRmn=-1.0*conj(RSmn[n][-m]);
						else
							tempRmn=conj(RSmn[n][-m]);
					}
					else tempRmn=RSmn[n][m];

					fmmTreeNode->father->Mmn[k][i][j]+=tempRmn*tempMmn;
				}
			}
		}
	}
}

void CMultipole::SonToFatherTrans_para(FMMTREENODE* fmmTreeNode)		//M2M
{
	if(fmmTreeNode==0||fmmTreeNode->father==0||fmmTreeNode->father->level<2)
		return;

	complex<double> RSmn0[21][21];
	int i, j, n, m, k;
	complex<double> tempMmn, tempRmn;

	getExpansionR_Local(fmmTreeNode->father->cube.center, fmmTreeNode->cube.center, RSmn0);
	for(k=0; k<DIM3; k++)
	{
		for(i=0; i<=p; i++)
			for(j=0; j<=i; j++)
			{
				for(n=0; n<=i; n++)
				{
					for(m=-n; m<=n; m++)
					{
						if(abs(j-m) > i-n) continue;

						if(j-m<0)
						{
							if((m-j) % 2)
							{
								tempMmn=-1.0*conj(fmmTreeNode->Mmn[k][i-n][m-j]);	//R[n][-m]=R[n][m]*(-1)的m次方
							}
							else tempMmn=conj(fmmTreeNode->Mmn[k][i-n][m-j]);
						}
						else tempMmn=fmmTreeNode->Mmn[k][i-n][j-m];

						if(m<0)
						{
							if((-m) % 2)
								tempRmn=-1.0*conj(RSmn0[n][-m]);
							else
								tempRmn=conj(RSmn0[n][-m]);
						}
						else tempRmn=RSmn0[n][m];

						#pragma omp critical
						fmmTreeNode->father->Mmn[k][i][j]+=tempRmn*tempMmn;
					}
				}
			}
	}
}

void CMultipole::UpWardPass(FMMTREENODE* fmmTreeNode)
{
	int i;

	if(fmmTreeNode==0||fmmTreeNode->tlNode<=max_nodes)
		return;

	for(i=0; i<8; i++)
		UpWardPass(fmmTreeNode->child[i]);
	
	if(fmmTreeNode->level<2)
		return;

	for(i=0; i<8; i++)
	{
		if(fmmTreeNode->child[i])
			SonToFatherTrans(fmmTreeNode->child[i]);
	}
}

void CMultipole::FarToNearField(FMMTREENODE* fmmTreeNode1,FMMTREENODE* fmmTreeNode2)
{
	if(fmmTreeNode1==0 || fmmTreeNode2==0 || fmmTreeNode1->level<3 || fmmTreeNode2->level<3)
		return;

	//	if(isNeighbor(fmmTreeNode1,fmmTreeNode2)) return;

	int i, j, n, m, k;
	complex<double> tempMmn;
	complex<double> tempSmn;

	getExpansionS(fmmTreeNode2->cube.center, fmmTreeNode1->cube.center);

	for(k=0; k<DIM3; k++)
	{
	for(i=0; i<=p; i++)
		for(j=0; j<=i; j++)
		{
			for(n=0; n<=p; n++)
			{
				for(m=-n; m<=n; m++)
				{
					if(j+m < 0)
					{
						if(-(j+m) % 2)
							tempSmn=-1.0*RSmn[i+n][-(j+m)];
						else
							tempSmn=RSmn[i+n][-(j+m)];
					}
					else
						tempSmn=conj(RSmn[i+n][j+m]);

					if(m < 0)
					{
						if(-m % 2)
							tempMmn=-1.0*conj(fmmTreeNode2->Mmn[k][n][-m]);
						else
							tempMmn=conj(fmmTreeNode2->Mmn[k][n][-m]);
					}
					else
						tempMmn=fmmTreeNode2->Mmn[k][n][m];

					if(i % 2)
						fmmTreeNode1->Lmn[k][i][j]-=tempSmn*tempMmn;
					else
						fmmTreeNode1->Lmn[k][i][j]+=tempSmn*tempMmn;
				}
			}
		}
	}
}

#if 1
void CMultipole::FarToNearField_para(FMMTREENODE* fmmTreeNode1,FMMTREENODE* fmmTreeNode2)
{
	if(fmmTreeNode1==0 || fmmTreeNode2==0 || fmmTreeNode1->level<3 || fmmTreeNode2->level<3)
		return;

	//	if(isNeighbor(fmmTreeNode1,fmmTreeNode2)) return;

	int i, j, n, m, k;
	complex<double> tempMmn;
	complex<double> tempSmn;
	complex<double> RSmn0[21][21];	//suppose p<=10

	getExpansionS_Local(fmmTreeNode2->cube.center, fmmTreeNode1->cube.center, RSmn0);

	for(k=0; k<DIM3; k++)
	{
		for(i=0; i<=p; i++)
			for(j=0; j<=i; j++)
			{
				for(n=0; n<=p; n++)
				{
					for(m=-n; m<=n; m++)
					{
						if(j+m < 0)
						{
							if(-(j+m) % 2)
								tempSmn=-1.0*RSmn0[i+n][-(j+m)];
							else
								tempSmn=RSmn0[i+n][-(j+m)];
						}
						else
							tempSmn=conj(RSmn0[i+n][j+m]);

						if(m < 0)
						{
							if(-m % 2)
								tempMmn=-1.0*conj(fmmTreeNode2->Mmn[k][n][-m]);
							else
								tempMmn=conj(fmmTreeNode2->Mmn[k][n][-m]);
						}
						else
							tempMmn=fmmTreeNode2->Mmn[k][n][m];

						//#pragma omp critical
						{
						if(i % 2)
							fmmTreeNode1->Lmn[k][i][j]-=tempSmn*tempMmn;
						else
							fmmTreeNode1->Lmn[k][i][j]+=tempSmn*tempMmn;
						}
					}
				}
			}
	}
}
#else 
void CMultipole::FarToNearField_para(FMMTREENODE* fmmTreeNode1,FMMTREENODE* fmmTreeNode2)
{
	if(fmmTreeNode1==0 || fmmTreeNode2==0 || fmmTreeNode1->level<3 || fmmTreeNode2->level<3)
		return;

	//	if(isNeighbor(fmmTreeNode1,fmmTreeNode2)) return;

	int i, j, n, m, k;
	complex<double> tempMmn;
	complex<double> tempSmn;
	complex<double> RSmn0[21][21];	//suppose p<=10

	getExpansionS_Local(fmmTreeNode2->cube.center, fmmTreeNode1->cube.center, RSmn0);

	for(k=0; k<DIM3; k++)
	{
		for(i=0; i<=p; i++)
			for(j=0; j<=i; j++)
			{
				for(n=0; n<=p; n++)
				{
					for(m=-n; m<=n; m++)
					{
						if(j+m < 0)
						{
							if(-(j+m) % 2)
								tempSmn=-1.0*conj(RSmn0[i+n][-(j+m)]);
							else
								tempSmn=conj(RSmn0[i+n][-(j+m)]);
						}
						else
							tempSmn=RSmn0[i+n][j+m];

						if(m < 0)
						{
							if(-m % 2)
								tempMmn=-1.0*conj(fmmTreeNode2->Mmn[k][n][-m]);
							else
								tempMmn=conj(fmmTreeNode2->Mmn[k][n][-m]);
						}
						else
							tempMmn=fmmTreeNode2->Mmn[k][n][m];

						//#pragma omp critical
						{
							if(i % 2)
								fmmTreeNode1->Lmn[k][i][j]-=conj(tempSmn)*tempMmn;
							else
								fmmTreeNode1->Lmn[k][i][j]+=conj(tempSmn)*tempMmn;
						}
					}
				}
			}
	}
}
#endif

void CMultipole::UpWardPass_para(FMMTREENODE* fmmTreeNode)
{
	int i;

	if(fmmTreeNode==0||fmmTreeNode->tlNode<=max_nodes)
		return;

	for(i=0; i<8; i++)
		UpWardPass_para(fmmTreeNode->child[i]);

	if(fmmTreeNode->level<2)
		return;

	for(i=0; i<8; i++)
	{
		if(fmmTreeNode->child[i])
			SonToFatherTrans_para(fmmTreeNode->child[i]);
	}
}


void CMultipole::FarToNear(FMMTREENODE* fmmTreeNode1,FMMTREENODE* fmmTreeNode2)	//M2L
{
	int i, j, n, m, k;
	complex<double> tempMmn;
	complex<double> tempSmn;

	getExpansionS(fmmTreeNode2->cube.center, fmmTreeNode1->cube.center);
	for (k=0; k<DIM3; k++)
	{
		for(i=0; i<=p; i++)
		{
			for(j=-i; j<=i; j++)
			{
				for(n=0; n<=p; n++)
				{
					for(m=-n; m<=n; m++)
					{
						if(abs(j+m) > i+n) continue;
						if(j+m < 0)
						{
							if(-(j+m) % 2)
								tempSmn=-1.0*RSmn[i+n][-(j+m)];
							else
								tempSmn=RSmn[i+n][-(j+m)];
						}
						else
							tempSmn=conj(RSmn[i+n][j+m]);

						if(m>=0)
							tempMmn=fmmTreeNode2->Mmn[k][n][m];
						else
							tempMmn=fmmTreeNode2->Mmn[k][n][-m+n];

						if(j>=0)
						{
							if(i%2)
								fmmTreeNode1->Lmn[k][i][j]-=tempSmn*tempMmn;
							else
								fmmTreeNode1->Lmn[k][i][j]+=tempSmn*tempMmn;
						}
						else
						{
							if(i%2)
								fmmTreeNode1->Lmn[k][i][-j+i]-=tempSmn*tempMmn;
							else
								fmmTreeNode1->Lmn[k][i][-j+i]+=tempSmn*tempMmn;
						}
					}
				}
			}
		}
	}
}

void CMultipole::FarToNear1(FMMTREENODE* fmmTreeNode1,FMMTREENODE* fmmTreeNode2)	//M2L
{
	int i, j, n, m, k;
	complex<double> tempMmn, tempSmn;

	getExpansionS(fmmTreeNode2->cube.center, fmmTreeNode1->cube.center);
	for (k=0; k<DIM3; k++)
	{
	for(i=0; i<=p; i++)
		for(j=0; j<=i; j++)
		{
			for(n=0; n<=p; n++)
			{
				for(m=-n; m<=n; m++)
				{
					if(j+m < 0)
					{
						if(-(j+m) % 2)
							tempSmn=-1.0*RSmn[i+n][-(j+m)];
						else
							tempSmn=RSmn[i+n][-(j+m)];
					}
					else
						tempSmn=conj(RSmn[i+n][j+m]);

					if(m < 0)
					{
						if(-m % 2)
							tempMmn=-1.0*conj(fmmTreeNode2->Mmn[k][n][-m]);
						else
							tempMmn=conj(fmmTreeNode2->Mmn[k][n][-m]);
					}
					else
						tempMmn=fmmTreeNode2->Mmn[k][n][m];

					if(i % 2)
						fmmTreeNode1->Lmn[k][i][j]-=tempSmn*tempMmn;
					else
						fmmTreeNode1->Lmn[k][i][j]+=tempSmn*tempMmn;

					if(fabs(std::real(fmmTreeNode1->Lmn[k][i][j])) > 100)
						int aaaa = 0.1;
				}
			}
		}
	}
}


int CMultipole::isFarField(D3POINT* pDt, FMMTREENODE* fmmTreeNode)
{
	if(fmmTreeNode==0 || fmmTreeNode->level<2) return 0;

	double rr=(pDt->x-fmmTreeNode->cube.center.x)*(pDt->x-fmmTreeNode->cube.center.x)
			 +(pDt->y-fmmTreeNode->cube.center.y)*(pDt->y-fmmTreeNode->cube.center.y)
			 +(pDt->z-fmmTreeNode->cube.center.z)*(pDt->z-fmmTreeNode->cube.center.z);
	if(rr >= fmmTreeNode->cube.halfD*fmmTreeNode->cube.halfD*12)	return 1;
	else return 0;
}

void CMultipole::FatherToSonTrans(FMMTREENODE* fmmTreeNode)		//L2L
{
	if(fmmTreeNode==0||fmmTreeNode->father==0||fmmTreeNode->father->level<3)
		return;

	int i, j, n, m, k;
	complex<double> tempRmn, tempLmn;

	getExpansionR(fmmTreeNode->father->cube.center, fmmTreeNode->cube.center);
	for(k=0; k<DIM3; k++)
	{
		for(i=0; i<=p; i++)
		for(j=0; j<=i; j++)
		{
			for(n=i; n<=p; n++)
			{
				for(m=-n; m<=n; m++)
				{
					if(abs(m-j) > n-i) continue;

					if(m-j < 0)
					{
						if((j-m) % 2)
							tempRmn=-1.0*conj(RSmn[n-i][j-m]);
						else
							tempRmn=conj(RSmn[n-i][j-m]);
					}
					else
						tempRmn=RSmn[n-i][m-j];

					if(m < 0)
					{
						if(-m % 2)
							tempLmn=-1.0*conj(fmmTreeNode->father->Lmn[k][n][-m]);
						else
							tempLmn=conj(fmmTreeNode->father->Lmn[k][n][-m]);
					}
					else
						tempLmn=fmmTreeNode->father->Lmn[k][n][m];

					fmmTreeNode->Lmn[k][i][j]+=tempRmn*tempLmn;
				}
			}
		}
	}
}

void CMultipole::FatherToSonTrans_para(FMMTREENODE* fmmTreeNode)		//L2L
{
	if(fmmTreeNode==0||fmmTreeNode->father==0||fmmTreeNode->father->level<3)
		return;

	int i, j, n, m, k;
	complex<double> tempRmn, tempLmn;
	complex<double> RSmn0[21][21];

	getExpansionR_Local(fmmTreeNode->father->cube.center, fmmTreeNode->cube.center, RSmn0);
	for(k=0; k<DIM3; k++)
	{
		for(i=0; i<=p; i++)
			for(j=0; j<=i; j++)
			{
				for(n=i; n<=p; n++)
				{
					for(m=-n; m<=n; m++)
					{
						if(abs(m-j) > n-i) continue;

						if(m-j < 0)
						{
							if((j-m) % 2)
								tempRmn=-1.0*conj(RSmn0[n-i][j-m]);
							else
								tempRmn=conj(RSmn0[n-i][j-m]);
						}
						else
							tempRmn=RSmn0[n-i][m-j];

						if(m < 0)
						{
							if(-m % 2)
								tempLmn=-1.0*conj(fmmTreeNode->father->Lmn[k][n][-m]);
							else
								tempLmn=conj(fmmTreeNode->father->Lmn[k][n][-m]);
						}
						else
							tempLmn=fmmTreeNode->father->Lmn[k][n][m];

						fmmTreeNode->Lmn[k][i][j]+=tempRmn*tempLmn;
					}
				}
			}
	}
}

void CMultipole::GetOneByPIRSmn(D3POINT pNode1,D3POINT pNode2)		//L2L
{
	int i, j, n, m;
	complex<double> tempRmn, tempLmn;

	getExpansionR(pNode1, pNode2);

	complex<double> w=ONEby4PI;

	for(n=0; n<=p; n++)
		for(m=0; m<=n; m++)
		{
			RSmn[n][m]*=w;
		}
}

void CMultipole::DownWardPass(FMMTREENODE* fmmTreeNode)
{
	int i;

	if(fmmTreeNode==0) return;

	if(fmmTreeNode->level>=3)
	{
		CUBENODELIST* temp;
		temp=fmmTreeNode->listInteract;
		while(temp)
		{
			FarToNearField(fmmTreeNode,temp->cubeNode);
			temp=temp->next;
		}

		if(fmmTreeNode->level>=4)
			FatherToSonTrans(fmmTreeNode);
	}

	for(i=0; i<8; i++)
	{
		if(fmmTreeNode->child[i])
			DownWardPass(fmmTreeNode->child[i]);
	}
}

void CMultipole::DownWardPass_para(FMMTREENODE* fmmTreeNode)
{
	int i;

	if(fmmTreeNode==0) return;

	if(fmmTreeNode->level>=3)
	{
		CUBENODELIST* temp;
		temp=fmmTreeNode->listInteract;
		while(temp)
		{
			FarToNearField_para(fmmTreeNode,temp->cubeNode);
			temp=temp->next;
		}

		if(fmmTreeNode->level>=4)
			FatherToSonTrans_para(fmmTreeNode);
	}

	for(i=0; i<8; i++)
	{
		if(fmmTreeNode->child[i])
			DownWardPass_para(fmmTreeNode->child[i]);
	}
	/*
	//for testing
	if(fmmTreeNode->cube.center.x == 1.125 && fmmTreeNode->cube.center.y == 1.125 && fmmTreeNode->cube.center.z == 1.125)
	{
		PrintMmn(fmmTreeNode);
		PrintLmn(fmmTreeNode);
	}
	if(fmmTreeNode->cube.center.x == 0.8125 && fmmTreeNode->cube.center.y == 0.8125 && fmmTreeNode->cube.center.z == 0.8125)
	{
		PrintMmn(fmmTreeNode);
		PrintLmn(fmmTreeNode);
	}
	if(fmmTreeNode->cube.center.x == 0.96875 && fmmTreeNode->cube.center.y == 0.65625 && fmmTreeNode->cube.center.z == 0.65625)
	{
		PrintMmn(fmmTreeNode);
		PrintLmn(fmmTreeNode);
	}
	//for testing
	*/
	
}

void CMultipole::ComputeTransform()
{
	int i;
	//complex<double> RSmn0[21][21];	//suppose p<=10
	complex<double> RSmn0[41][41];	//suppose p<=10
#if 0
	UpWardPass(fmmTreeRoot);         //向上集成系数
	DownWardPass(fmmTreeRoot);		 //向下集成系数
#else
	#ifndef DEBUG
#endif
// 
// 	#pragma omp parallel for
// 	for (i=0; i<8; i++)
// 		UpWardPass_para(fmmTreeRoot->child[i]);
	UpWardPass(fmmTreeRoot);         //向上集成系数

#ifndef DEBUG
#pragma omp parallel for
#endif
	for(i=0; i<8; i++)
		DownWardPass_para(fmmTreeRoot->child[i]);
#endif
}

int CMultipole::GetPtLeaf(Node* pNod, int *k)
{
	int i;
	*k = -1;
	for (i=0; i<tlLeaf; i++)
	{
		if(isInCube(pNod, &(arryLeafPoint[i]->cube)))
		{
			*k = i;
			break;
		}
	}

	return BLM_SUCCESS;
}

FMMTREENODE* CMultipole::GetPtTreeNode(Node* pNod, FMMTREENODE* fmmTreeNode)
{
	int i;
	bool bin = false;
	FMMTREENODE * fn;

	//for testing
// 	if(fmmTreeNode->level == 5)
// 		return fmmTreeNode;

	for(i=0; i<8; i++)
	{
		if(fmmTreeNode->child[i] && isInCube(pNod, &(fmmTreeNode->child[i]->cube)))
		{
			fn = GetPtTreeNode(pNod, fmmTreeNode->child[i]);
			bin = true;
			break;
		}
	}

	if(bin)
		return fn;
	else
		return fmmTreeNode;
}

int CMultipole::GetResultPt(double *pta, double *val, int ilyer)
{
	int i, j, k, n, m, l;
	D3POINT pt;
	double /*pta[3],*/ vneig[3] = {0.0}, vchild = 0.0;
	double w=ONEby4PI;
	bool bleaf = false;
	/*
	complex<double>** RSmn0;
	RSmn0=new  complex<double>*[2*p+1];
	for(i=0;i<=2*p;i++) RSmn0[i]=new  complex<double>[i+1];*/
	complex<double> RSmn0[21][21];	//suppose p<=10
	//complex<double> RSmn0[41][41];	//suppose p<=20
	int pp = 2*10;

	pt.init(pta[0], pta[1], pta[2]);
// 	pta[0] = pNod->coord[0];
// 	pta[1] = pNod->coord[1];
// 	pta[2] = pNod->coord[2];

	Node pNod;
	pNod.coord[0] = pta[0];
	pNod.coord[1] = pta[1];
	pNod.coord[2] = pta[2];

	complex <double>temp[3];
	for (i=0; i<DIM3; i++)
	{
		temp[i].real(0.0);
		temp[i].imag(0.0);
	}

#if 0
	//策略0
	GetPtLeaf(&pNod, &k);
	if(k<0)
	{
		spdlog::info("Error: no leaf contains the node!\n");
		return BLM_FAILED;
	}
	getExpansionR(pt, arryLeafPoint[k]->cube.center);
	for(n=0; n<=p; n++)
	{
		for(m=0; m<=n; m++)
			temp+=RSmn[n][m]*arryLeafPoint[k]->Lmn[n][m];

		for(m=1; m<=n; m++)
			temp+=conj(RSmn[n][m])*conj(arryLeafPoint[k]->Lmn[n][m]);
	}
#endif

#if 0
	//策略1
	//bool bleaf = false;
	FMMTREENODE* fnode = GetPtTreeNode(&pNod, fmmTreeRoot);
	bleaf = isLeafNode(fnode);
	//spdlog::info("Node level: {} {}\n", fnode->level, bleaf);
	for(k=0; k<tlLeaf; k++)
	{
		if(bleaf && (fnode == arryLeafPoint[k] || isNeighbor(fnode, arryLeafPoint[k], 1)))
			continue;

		if(!bleaf && isChild(fnode, arryLeafPoint[k]))
			continue;

		getExpansionS_Local(arryLeafPoint[k]->cube.center, pt, RSmn0);
		//getExpansionS(pt, arryLeafPoint[k]->cube.center);
		for(l=0; l<DIM3; l++)
		{
			for(n=0; n<=p; n++)
			{
				for(m=0; m<=n; m++)
// 					if(m%2)	//odd number
// 						temp[l]+=(-1.0)*conj(RSmn[n][m])*arryLeafPoint[k]->Mmn[l][n][m];
// 					else
						temp[l]+=conj(RSmn0[n][m])*arryLeafPoint[k]->Mmn[l][n][m];

				for(m=1; m<=n; m++)
					if(m%2)	//odd number
						temp[l]+=(-1.0)*RSmn0[n][m]*conj(arryLeafPoint[k]->Mmn[l][n][m]);
					else
						temp[l]+=RSmn0[n][m]*conj(arryLeafPoint[k]->Mmn[l][n][m]);
			}
		}
	}
	if(bleaf)
		getNeigValue(/*fmmTreeRoot*/fnode, pta, vneig);

	if(!bleaf)
		getChildValue(fnode, pta, vneig);
#endif

#if 0
	//策略2
	//bool bleaf = false;
	FMMTREENODE* fnode = GetPtTreeNode(&pNod, fmmTreeRoot);
	spdlog::info("Node level: {}\n", fnode->level);
	bleaf = isLeafNode(fnode);
	for(k=0; k<tlLeaf; k++)
	{
		//calculate M2L
		if(bleaf && (fnode == arryLeafPoint[k] || isNeighbor(fnode, arryLeafPoint[k], 1)))
			continue;

		if(!bleaf && isChild(fnode, arryLeafPoint[k]))
			continue;

		FarToNear1(fnode, arryLeafPoint[k]);
		//FarToNear(fnode, arryLeafPoint[k]);
	}

	getExpansionR(fnode->cube.center, pt);
	for (l=0; l<DIM3; l++)
	{
		for(n=0; n<=p; n++)	//Lmn[n][m+n]=Lmn[n][-m]
		{
			for(m=0; m<=n; m++)
				temp[l]+=RSmn[n][m]*fnode->Lmn[l][n][m];
#if 1
			for(m=1; m<=n; m++)
				temp[l]+=conj(RSmn[n][m])*conj(fnode->Lmn[l][n][m]);
#else
			for(m=1; m<=n; m++)
			{
				if(m%2)
					temp[l]+=(-1.0)*conj(RSmn[n][m])*fnode->Lmn[l][n][m+n];
				else
					temp[l]+=conj(RSmn[n][m])*fnode->Lmn[l][n][m+n];
			}
#endif
		}
	}
	if(bleaf)
		getNeigValue(/*fmmTreeRoot*/fnode, pta, vneig);

	if(!bleaf)
		getChildValue(fnode, pta, vneig);
#endif

#if 1
	//策略3
	FMMTREENODE* fnode = GetPtTreeNode(&pNod, fmmTreeRoot);

	//int interCnt = countInteraction(fnode);

	//for testing
	//spdlog::info("Node level: {}\n", fnode->level);

	getExpansionR_Local(fnode->cube.center, pt, RSmn0);

	//for testing
// 	char filename[256];
// 	memset(filename, 0, sizeof(filename));
// 	sprintf(filename, "RSmn_%d_%d.out", fnode->level, fnode->tlNode);
// 	PrintRSmn(RSmn0, filename);
	//for testing

	if (fnode->level <= 2)
	{
		spdlog::info("Warning: the level of the located tree node is not larger than 2!\n");
	}

	if(fnode->level > 2)
	{
		for(l=0; l<DIM3; l++)
		{
			for(n=0; n<=p; n++)
			{
				for(m=0; m<=n; m++)
					temp[l]+=RSmn0[n][m]*fnode->Lmn[l][n][m];

				for(m=1; m<=n; m++)
					temp[l]+=conj(RSmn0[n][m])*conj(fnode->Lmn[l][n][m]);
			}
		}
		getNeigValue(/*fmmTreeRoot*/fnode, pta, vneig);
	}
	else
	{
		bleaf = isLeafNode(fnode);
		for(k=0; k<tlLeaf; k++)
		{
			if(bleaf && (fnode == arryLeafPoint[k] || isNeighbor(fnode, arryLeafPoint[k], 1)))
				continue;

			if(!bleaf && isChild(fnode, arryLeafPoint[k]))
				continue;

			for (n=0; n<=pp; n++)
			{
				for(m=0; m<=pp; m++)
				{
					RSmn0[n][m].real(0.0);
					RSmn0[n][m].imag(0.0);
				}
			}

			getExpansionS_Local(arryLeafPoint[k]->cube.center, pt, RSmn0);
			for(l=0; l<DIM3; l++)
			{
				for(n=0; n<=p; n++)
				{
					for(m=0; m<=n; m++)
							temp[l]+=conj(RSmn0[n][m])*arryLeafPoint[k]->Mmn[l][n][m];

					for(m=1; m<=n; m++)
						if(m%2)	//odd number
							temp[l]+=(-1.0)*RSmn0[n][m]*conj(arryLeafPoint[k]->Mmn[l][n][m]);
						else
							temp[l]+=RSmn0[n][m]*conj(arryLeafPoint[k]->Mmn[l][n][m]);
				}
			}
		}
		if(bleaf)
			getNeigValue(/*fmmTreeRoot*/fnode, pta, vneig);

		if(!bleaf)
			getChildValue(fnode, pta, vneig);
	}
#endif

	for(l=0; l<DIM3; l++)
		val[l] = std::real(temp[l])*w+vneig[l];

	return BLM_SUCCESS;
}

void CMultipole::getNodNeigs(FMMTREENODE* fmmTreeNode, int *nNeig, FMMTREENODE** NodNeigs)
{
	*nNeig = 0;

	NodNeigs[*nNeig] = fmmTreeNode;
	(*nNeig) += 1;

	CUBENODELIST* neighbor = fmmTreeNode->listNeighbors;
	while(neighbor)
	{
		NodNeigs[*nNeig] = neighbor->cubeNode;
		(*nNeig) += 1;

		if(*nNeig>60)
		{
			spdlog::info("Error: the number of neighbors exceeds 60!\n");
			//exit(0);
		}

		neighbor = neighbor->next;
	}
}

void CMultipole::getNeigValue(FMMTREENODE* fmmTreeNode, double *pt, double *val)
{
	int i, j, npt, nelm, eidx, nNeig;
	double f[3], g[3], *pBc, *pU, *pNorm;
	int *elms = NULL;
	elms= new int[tlNode];

	FMMTREENODE * nodeNeigs[60];
	getNodNeigs(fmmTreeNode, &nNeig, nodeNeigs);

	getTreeNodeElms(nNeig, nodeNeigs, &npt, &m_pPts, &nelm, &elms/*m_pElms*/);

	//spdlog::info("neig: {}\n", nelm);

#if 0
	pBc = new double[(DIM3+1)*nelm];
	pU = new double[DIM3*nelm];
	pNorm = new double[DIM3*nelm];

	for (i=0; i</*nelm*/tlNode; i++)
	{
		eidx = el2g[i];

		pBc[(DIM3+1)*i + 0] = pDomain->m_pBC[(DIM3+1)*eidx + 0];
		for (j=0; j<DIM3; j++)
		{
			pBc[(DIM3+1)*i + (j+1)] = pDomain->m_pBC[(DIM3+1)*eidx + (j+1)];
			pU[DIM3*i + j] = pDomain->m_pUv[DIM3*eidx + j];
			pNorm[DIM3*i + j] = pDomain->m_pNorm[DIM3*eidx + j];
		}
	}
	m_pPotentialBEM->InitSurfElmInfo(nelm, m_pPts, m_pElms);
	m_pPotentialBEM->FieldEvalBEM(1, pt, f, g, nelm, elms, pBc, pU, true);
#endif
	m_pPotentialBEM->FieldEvalBEM(1, pt, f, g, nelm, elms, pDomain->m_pBC, pDomain->m_pUv, true);
// 
// 	delete []pBc;
// 	delete []pU;
// 	delete []pNorm;

	val[0] = f[0];
	val[1] = f[1];
	val[2] = f[2];

	if(elms){ delete []elms; elms = NULL; }
}


void CMultipole::getChildValue(FMMTREENODE* fmmTreeNode, double *pt, double *val)
{
	int i, j, npt, nelm, eidx, nNeig;
	double f[3], g[3], *pBc, *pU, *pNorm;
	int *elms = NULL;
	elms= new int[tlNode];

	FMMTREENODE * nodeNeigs[1];
	nNeig = 1;
	nodeNeigs[0] = fmmTreeNode;

	getTreeNodeElms(nNeig, nodeNeigs, &npt, &m_pPts, &nelm, &elms/*m_pElms*/);
#if 0
	pBc = new double[(DIM3+1)*nelm];
	pU = new double[DIM3*nelm];
	pNorm = new double[DIM3*nelm];

	for (i=0; i<nelm; i++)
	{
		eidx = el2g[i];

		pBc[(DIM3+1)*i + 0] = pDomain->m_pBC[(DIM3+1)*eidx + 0];
		for (j=0; j<DIM3; j++)
		{
			pBc[(DIM3+1)*i + (j+1)] = pDomain->m_pBC[(DIM3+1)*eidx + (j+1)];
			pU[DIM3*i + j] = pDomain->m_pUv[DIM3*eidx + j];
			pNorm[DIM3*i + j] = pDomain->m_pNorm[DIM3*eidx + j];
		}
	}

	m_pPotentialBEM->InitSurfElmInfo(nelm, m_pPts, m_pElms);
	m_pPotentialBEM->FieldEvalBEM(1, pt, f, g, nelm, npt, m_pPts, pBc, m_pElms, pNorm, pU, true);

	delete []pBc;
	delete []pU;
	delete []pNorm;
#endif
	m_pPotentialBEM->FieldEvalBEM(1, pt, f, g, nelm, elms, pDomain->m_pBC, pDomain->m_pUv, true);

	val[0] = f[0];
	val[1] = f[1];
	val[2] = f[2];
	if(elms){ delete []elms; elms = NULL; }
}

bool CMultipole::isLeafNode(FMMTREENODE* fmmTreeNode)
{
	bool bLeaf = true;
	for (int i=0; i<8; i++)
		bLeaf = bLeaf && (fmmTreeNode->child[i] == 0);

	return bLeaf;
}

bool CMultipole::isNeighbor(FMMTREENODE* fmmTreeNode1,FMMTREENODE* fmmTreeNode2, bool caled)
{
	CUBENODELIST* neighbor = fmmTreeNode1->listNeighbors;
	while(neighbor)
	{
		if(fmmTreeNode2 == neighbor->cubeNode)
			return true;

		neighbor = neighbor->next;
	}
	return false;
}

bool CMultipole::isChild(FMMTREENODE* fmmTreeNode1,FMMTREENODE* fmmTreeNode2)
{
	//note that fmmTreeNode1 is not a leaf node
	int i;
	bool bchild = false;

	if(fmmTreeNode1 == NULL)
		return false;

	for (i=0; i<8; i++)
	{
		if(fmmTreeNode2 == fmmTreeNode1->child[i])
			return true;
	}

	for (i=0; i<8; i++)
		bchild |= isChild(fmmTreeNode1->child[i], fmmTreeNode2);

	return bchild;
}

///////////////////////////////get main diagonal entry for mxQ///////////////////////////////
//建立矩阵与向量相乘A*X--常量场
int CMultipole::build_AX_constantU(double* x0,double* bt)
{
	//x0:A*X中的X
	//bt:待求向量bt＝AX

	long i, j, k, m, n;

	resetFmmTree(fmmTreeRoot);   //树结构系数清零
    initLeafs(x0);               //将初始向量x0集成到叶子节点
	UpWardPass(fmmTreeRoot);         //向上集成系数
    DownWardPass(fmmTreeRoot);		 //向下集成系数

	for(i=0; i<tlNode; i++) bt[i]=0.0;
#if 0
	for(k=0; k<tlLeaf; k++)
	{
		for(i=0; i<arryLeafPoint[k]->tlNode; i++)
		{
			fread(pDomain->localU[i], sizeof(double), arryLeafPoint[k]->tlLocalNode, pDomain->pBaseU);
			for(n=0; n<=p; n++)
			for(m=0; m<=n; m++)
			{
				fread(&(pDomain->INT_RmnU[i][n][m].x), sizeof(double), 1, pDomain->pBaseRmnU);
				fread(&(pDomain->INT_RmnU[i][n][m].y), sizeof(double), 1, pDomain->pBaseRmnU);
			}

			for(j=0; j<arryLeafPoint[k]->tlLocalNode; j++)
				bt[arryLeafPoint[k]->sequence[i]]+=pDomain->localU[i][j]*x0[arryLeafPoint[k]->sequence[j]];

			complex temp;
			temp.x=0.0;
			temp.y=0.0;
			for(n=0; n<=p; n++)
			{
				for(m=0; m<=n; m++)
					temp+=pDomain->INT_RmnU[i][n][m]*arryLeafPoint[k]->Lmn[n][m];

				for(m=1; m<=n; m++)
					temp+=conjugate(pDomain->INT_RmnU[i][n][m])*conjugate(arryLeafPoint[k]->Lmn[n][m]);
			}
		    if(fabs(im(temp)/re(temp)) > 0.0001)
			{
				ErrDisplay("Wrong for far field");
				return 1;
			}

			bt[arryLeafPoint[k]->sequence[i]]+=re(temp);
		}
	}
	rewind(pDomain->pBaseU);
	rewind(pDomain->pBaseRmnU);
#endif

	return 0;
}

//建立矩阵与向量相乘A*X--主对角元出外
int CMultipole::build_AX_constantQ(double* bt)
{
	//bt:待求向量bt＝AX

	long i, j, k, m, n;

	for(i=0; i<tlNode; i++) bt[i]=0.0;
#if 0
	for(k=0; k<tlLeaf; k++)
	{
		for(i=0; i<arryLeafPoint[k]->tlNode; i++)
		{
			fread(pDomain->localQ[i], sizeof(double), arryLeafPoint[k]->tlLocalNode, pDomain->pBaseQ);
			for(n=0; n<=p; n++)
			for(m=0; m<=n; m++)
			{
				fread(&(pDomain->INT_RmnQ[i][n][m].x), sizeof(double), 1, pDomain->pBaseRmnQ);
				fread(&(pDomain->INT_RmnQ[i][n][m].y), sizeof(double), 1, pDomain->pBaseRmnQ);
			}

			for(j=0; j<arryLeafPoint[k]->tlLocalNode; j++)
			{
				if(i==j) continue;
				bt[arryLeafPoint[k]->sequence[i]]+=pDomain->localQ[i][j]*x[arryLeafPoint[k]->sequence[j]];
			}

			complex temp;
			temp.x=0.0;
			temp.y=0.0;
			for(n=0; n<=p; n++)
			{
				for(m=0; m<=n; m++)
					temp+=pDomain->INT_RmnQ[i][n][m]*arryLeafPoint[k]->Lmn[n][m];

				for(m=1; m<=n; m++)
					temp+=conjugate(pDomain->INT_RmnQ[i][n][m])*conjugate(arryLeafPoint[k]->Lmn[n][m]);
			}
		    if(fabs(im(temp)/re(temp))  > 0.0001)
			{
				ErrDisplay("Wrong for far field");
				return 1;
			}

			bt[arryLeafPoint[k]->sequence[i]]+=re(temp);
		}
	}
#endif
	return 0;
}

//求解快速算法的预处理矩阵--常量场
int CMultipole::build_preCOND_constant()
{
	long i, k;
/*
	//No preconditioner
	for(k=0; k<tlLeaf; k++)
	{
		for(i=0; i<arryLeafPoint[k]->tlNode; i++)
		for(j=0; j<arryLeafPoint[k]->tlNode; j++)
		{
			if(i == j) arryLeafPoint[k]->M[i][j]=1.0;
			else arryLeafPoint[k]->M[i][j]=0.0;
		}
	}
*/
#if 0
	//Blocked diagonal preconditioner
	for(k=0; k<tlLeaf; k++)
	{
		for(i=0; i<arryLeafPoint[k]->tlNode; i++)
		{
			fread(arryLeafPoint[k]->M[i], sizeof(double), arryLeafPoint[k]->tlNode, pDomain->pBaseU);
			fseek(pDomain->pBaseU, (arryLeafPoint[k]->tlLocalNode-arryLeafPoint[k]->tlNode)*sizeof(double), SEEK_CUR);
		}
/*	
				FILE *stream;
				stream = fopen( "mxLocalM", "w");
				if( stream != NULL )
				{
					fprintf(stream,"leaf(k)=%d; tlNode=%d\n", k, arryLeafPoint[k]->tlNode);
					for(i=0; i<arryLeafPoint[k]->tlNode; i++)
					{
						for(j=0; j<arryLeafPoint[k]->tlNode; j++)
							fprintf(stream,"%10.4g  ", arryLeafPoint[k]->M[i][j]);
						fprintf(stream,"\n");
					}
				}
				fclose(stream);
*/
		if(inv0(arryLeafPoint[k]->M, arryLeafPoint[k]->tlNode, 1.0e-15)) 
		{
			ErrDisplay("Biuld preconditioner for constant field fail");
			return 1;
		}
	}
#endif
	return 0;
}

inline void CMultipole::deleteTemps(double* r, double* vh, double* u, double* u2, double* q, 
							double* h, double* f, double** V, double** QT, double** R, double** W)
{
	int k;
	delete[] r;	    delete[] vh;    delete[] u;
	delete[] u2;    delete[] q;	    delete[] h;      delete[] f;

	for(k=0; k<tlUnknown; k++) delete[] V[k];
	delete[] V;
	for(k=0; k<IT_inner+1; k++) delete[] QT[k];
	delete[] QT;
	for(k=0; k<IT_inner; k++) delete[] R[k];
	delete[] R;
	for(k=0; k<tlUnknown; k++) delete[] W[k];
	delete[] W;
}

inline void CMultipole::preconditioning_constant(double* v, double* v0)
{
	long k, l, m;
/*
	for(k=0; k<tlNode; k++) v[k]=v0[k];
*/
	for(k=0; k<tlNode; k++) v[k]=0.0;
	for(k=0; k<tlLeaf; k++)
	{
		for(l=0; l<arryLeafPoint[k]->tlNode; l++)
		{
			for(m=0; m<arryLeafPoint[k]->tlNode; m++)
				v[arryLeafPoint[k]->sequence[l]]
					+=arryLeafPoint[k]->M[l][m]*v0[arryLeafPoint[k]->sequence[m]];
		}
	}//v=M*v0
}

#if 0
int CMultipole::fmm_GMRES_constant()
{
	char str[50];
	long i,j,k,l;

	long n=tlUnknown;

	double tol_B;                   
	double nrmB;                     
	double* r=new double[n];        
	double* vh=new double[n];      
    double* u=new double[n];
	double* u2=new double[n];
	double* q=new double[IT_inner];
	double nrmR;                      

	double fj_1;
	int flag;
//    0 GMRES converged to the desired tolerance TOL within MAXIT iterations.
//    1 GMRES iterated MAXIT times but did not converge.
//    2 preconditioner M was ill-conditioned.
//    3 GMRES stagnated (two consecutive iterates were the same).

	double phibar;
	double rt,c,s,temp;

	double** V=new double*[n];
	for(i=0; i<n; i++) V[i]=new double[IT_inner+1];

	double* h=new double[IT_inner+1];

	double** QT=new double*[IT_inner+1];
	for(i=0; i<IT_inner+1; i++) QT[i]=new double[IT_inner+1];

	double** R=new double*[IT_inner];
	for(i=0; i<IT_inner; i++) R[i]=new double[IT_inner];

	double* f=new double[IT_inner]; 
//	double* residual=new double[IT_inner];  //temporary

	double** W=new double*[n];  
	for(i=0; i<n; i++) W[i]=new double[IT_inner];

	nrmB = norm0(b,n);                     
	if(nrmB==0.0)
	{
		for(i=0; i<n; i++) x[i]=0.0;
		flag=0;
		it_inner=0;
		it_outer=0;
		deleteTemps(r, vh, u, u2, q, h, f, V, QT, R, W);
		return 0;
	}

	preconditioning_constant(x, b);//x0 = M \ b, initial guess
//	for(i=0; i<n; i++) x0[i]=1.0;
//	for(i=0; i<n; i++) x[i]=x0[i];

	flag=1;
	tol_B=tolerance*nrmB;

	if(build_AX_constantU(x, r))
	{
		deleteTemps(r, vh, u, u2, q, h, f, V, QT, R, W);
		sprintf(str, "Current iteration: it_inner= 0, it_outer= 0");
		ErrDisplay(str);
		ErrDisplay("build_AX_constant fail");
		return 1;
	}
	for(i=0; i<n; i++) r[i]=b[i]-r[i];//r=b-A*x

	nrmR=norm0(r,n);
	if(nrmR<=tol_B)
	{
		flag=0;
		it_inner=0;
		it_outer=0;
		deleteTemps(r, vh, u, u2, q, h, f, V, QT, R, W);
		return 0;
	}

	////
	for(i=1; i<=IT_outer; i++)
	{		
		for(j=0; j<IT_inner+1; j++)
		{
			for(k=0; k<IT_inner+1; k++) QT[j][k]=0.0;
		}

		for(j=0; j<IT_inner; j++)
		{
			for(k=0; k<IT_inner; k++) R[j][k]=0.0;
		}

		preconditioning_constant(vh, r); 

		h[0]=norm0(vh,n); 
		
		for(j=0; j<n; j++) V[j][0]=vh[j]/h[0]; //V(:,1) = vh / h(1)
     
		QT[0][0]=1.0;

		phibar=h[0];

		fj_1=0.0;
		for(j=0; j<IT_inner; j++)  
		{
			FILE *logStream;
			logStream = fopen( pDomain->pStructure->pDoc->csLogFile, "a+");
			if( logStream != NULL )
			{
				CTime tCurrent;
				CString csTime;
				tCurrent = CTime::GetCurrentTime();
				csTime = tCurrent.Format( "%H:%M:%S, on %A, %B %d, %Y" );
				fprintf(logStream,"----Outer loop: %4d;   Inner loop: %4d;   Time: %s\n", i, j, csTime);
				fprintf(logStream,"----tol_B= %10.4g; nrmR= %10.4g; fj_1= %10.4g\n\n", tol_B, nrmR, fj_1);
				fclose(logStream);
			}
			for(k=0; k<n; k++) u[k]=V[k][j]; 
			
			if(build_AX_constantU(u, u2))//////////u2=A*v(:,j),可以利用u
			{
				deleteTemps(r, vh, u, u2, q, h, f, V, QT, R, W);
				sprintf(str, "Current iteration: it_inner=%4d, it_outer=%4d", j+1, i);
				ErrDisplay(str);
				ErrDisplay("build_AX_constant fail");
				return 1;
			}

			preconditioning_constant(u, u2);//u=M*u2

			for(k=0; k<=j; k++)
			{
				h[k]=0.0;
				for(l=0; l<n; l++) h[k]+=V[l][k]*u[l]; 
					
				for(l=0; l<n; l++) u[l]-=h[k]*V[l][k]; 
			}

			h[j+1]=norm0(u,n);

			for(k=0; k<n; k++) V[k][j+1]=u[k]/h[j+1]; 
				
			for(k=0; k<=j; k++)
			{
				R[k][j]=0.0;
				for(l=0; l<=j; l++) R[k][j]+=QT[k][l]*h[l];
			} 

			rt=R[j][j];

			if(h[j+1]==0.0)
			{
				c=1.0;
				s=0.0;
			}
			else if(fabs(h[j+1])>fabs(rt))
			{
				temp=rt/h[j+1];
                s=1.0/sqrt(1.0+temp*temp);
                c=-temp*s;
			}
			else
			{
				temp=h[j+1]/ rt;
                c=1.0/sqrt(1.0+temp*temp);
                s=-temp*c;
			}

			R[j][j]=c*rt-s*h[j+1];

			for(k=0; k<=j; k++) q[k]=QT[j][k];

			for(k=0; k<=j; k++) QT[j][k]=c*q[k];

			for(k=0; k<=j; k++) QT[j+1][k]=s*q[k];

			QT[j][j+1]=-s;
            QT[j+1][j+1]=c;
            f[j]=c*phibar;
            phibar=s*phibar;

			if(f[j]==0.0) 
 			{
 				deleteTemps(r, vh, u, u2, q, h, f, V, QT, R, W);
				sprintf(str, "Current iteration: it_inner=%4d, it_outer=%4d", j+1, i);
				ErrDisplay(str);
				ErrDisplay("GMRES stagnated (f[j]=0)");
				return 1;
			}
			if(fabs((f[j]-fj_1)/f[j]) <=0.001) 
			{
 				deleteTemps(r, vh, u, u2, q, h, f, V, QT, R, W);
				sprintf(str, "Current iteration: it_inner=%4d, it_outer=%4d", j+1, i);
				ErrDisplay(str);
				ErrDisplay("GMRES stagnated (two consecutive iterates were the same)");
				return 1;
			}
			fj_1=f[j];

			if(j==0)
			{
				for(k=0; k<n; k++) W[k][j]=V[k][j]/R[j][j];
			}
			else
			{
				for(k=0; k<n; k++)
				{
					W[k][j]=V[k][j]/R[j][j];
					for(l=0; l<=j-1; l++) W[k][j]-=W[k][l]*R[l][j]/R[j][j];
				}
			}
			for(k=0; k<n; k++) x[k]+=f[j]*W[k][j];

			if(build_AX_constantU(x, vh))////////
			{
				deleteTemps(r, vh, u, u2, q, h, f, V, QT, R, W);
				ErrDisplay("build_AX_constant fail");
				sprintf(str, "Current iteration: it_inner=%4d, it_outer=%4d", j+1, i);
				ErrDisplay(str);
				return 1;
			}
			for(k=0; k<n; k++) vh[k]=b[k]-vh[k];

			nrmR=norm0(vh,n);

			if(nrmR<=tol_B)
			{
				flag = 0;
				it_inner=j+1;
				it_outer=i;
				deleteTemps(r, vh, u, u2, q, h, f, V, QT, R, W);
				return 0;
			}
		}

		if(i==IT_outer)
		{
			deleteTemps(r, vh, u, u2, q, h, f, V, QT, R, W);
			sprintf(str, "Current iteration: it_inner=%4d, it_outer=%4d", j+1, i);
			ErrDisplay(str);
			ErrDisplay("GMRES iterated MAXIT times but did not converge");
			return 1;
		}

		if(flag==1)
		{
			if(build_AX_constantU(x, r))////////
			{
				it_inner=j;
				it_outer=i;
				deleteTemps(r, vh, u, u2, q, h, f, V, QT, R, W);
				sprintf(str, "Current iteration: it_inner=%4d, it_outer=%4d", it_inner, it_outer);
				ErrDisplay(str);
				ErrDisplay("build_AX_constant fail");
				return 1;
			}
			for(k=0; k<n; k++) r[k]=b[k]-r[k];//r0=b-A*x0
		}
	}

	deleteTemps(r, vh, u, u2, q, h, f, V, QT, R, W);

	return 0;
}
#endif
//////////////////////////////real problem/////////////////////////////////////////////
#if 0
//建立右端向量B
void CMultipole::build_B()
{
	pDomain->getRightHandVector();
	b=pDomain->vtF.vt;
}

//建立矩阵与向量相乘A*X
int CMultipole::build_AX()
{
	//vtF:A*X中的X
	//vtU:待求向量vtU＝mxU*X
	//vtQ:待求向量vtQ＝mxQ*X

	long i, j, k, m, n;

	resetFmmTree(fmmTreeRoot);   //树结构系数清零
    initLeafs(pDomain->vtF.vt);  //将初始向量x0集成到叶子节点
	UpWardPass(fmmTreeRoot);         //向上集成系数
    DownWardPass(fmmTreeRoot);		 //向下集成系数

	for(i=0; i<tlNode; i++)
	{
		pDomain->vtU.vt[i]=0.0;
		pDomain->vtQ.vt[i]=0.0;
	}

	if((pDomain->pBaseU = fopen( pDomain->strBaseU, "rb")) == NULL ||
	   (pDomain->pBaseQ = fopen( pDomain->strBaseQ, "rb")) == NULL ||
	   (pDomain->pBaseRmnU = fopen( pDomain->strBaseRmnU, "rb")) == NULL ||
	   (pDomain->pBaseRmnQ = fopen( pDomain->strBaseRmnQ, "rb")) == NULL)
	{
		ErrDisplay("Cannot open database files");
		ErrDisplay("build_AX fail");
		return 1;
	}

	for(k=0; k<tlLeaf; k++)
	{
		for(i=0; i<arryLeafPoint[k]->tlNode; i++)
		{
			complex tempU, tempQ;
			if(pDomain->NodeArray[arryLeafPoint[k]->sequence[i]]->pFace->BC.known==potential)
			{
				fread(pDomain->localU[i], sizeof(double), arryLeafPoint[k]->tlLocalNode, pDomain->pBaseU);
				fseek(pDomain->pBaseQ, arryLeafPoint[k]->tlLocalNode*sizeof(double), SEEK_CUR);
				for(j=0; j<arryLeafPoint[k]->tlLocalNode; j++)
					pDomain->vtU.vt[arryLeafPoint[k]->sequence[i]]
						+=pDomain->localU[i][j]*pDomain->vtF.vt[arryLeafPoint[k]->sequence[j]];

				for(n=0; n<=p; n++)
				for(m=0; m<=n; m++)
				{
					fread(&(pDomain->INT_RmnU[i][n][m].x), sizeof(double), 1, pDomain->pBaseRmnU);
					fread(&(pDomain->INT_RmnU[i][n][m].y), sizeof(double), 1, pDomain->pBaseRmnU);
				}
				fseek(pDomain->pBaseRmnQ, (p+1)*(p+2)*sizeof(double), SEEK_CUR);

				tempU.x=0.0; tempU.y=0.0;
				for(n=0; n<=p; n++)
				{
					for(m=0; m<=n; m++)
						tempU+=pDomain->INT_RmnU[i][n][m]*arryLeafPoint[k]->Lmn[n][m];

					for(m=1; m<=n; m++)
						tempU+=conjugate(pDomain->INT_RmnU[i][n][m])*conjugate(arryLeafPoint[k]->Lmn[n][m]);
				}
				if(im(tempU)/re(tempU) > 0.0001)
				{
					ErrDisplay("Wrong for far field");
					return 1;
				}
				pDomain->vtU.vt[arryLeafPoint[k]->sequence[i]]+=re(tempU);
			}

			else if(pDomain->NodeArray[arryLeafPoint[k]->sequence[i]]->pFace->BC.known==flux)
			{
				fread(pDomain->localQ[i], sizeof(double), arryLeafPoint[k]->tlLocalNode, pDomain->pBaseQ);
				pDomain->localQ[i][i]=pDomain->vtDiagonalQ[arryLeafPoint[k]->sequence[i]];
				fseek(pDomain->pBaseU, arryLeafPoint[k]->tlLocalNode*sizeof(double), SEEK_CUR);
				for(j=0; j<arryLeafPoint[k]->tlLocalNode; j++)
					pDomain->vtQ.vt[arryLeafPoint[k]->sequence[i]]
						+=pDomain->localQ[i][j]*pDomain->vtF.vt[arryLeafPoint[k]->sequence[j]];

				for(n=0; n<=p; n++)
				for(m=0; m<=n; m++)
				{
					fread(&(pDomain->INT_RmnQ[i][n][m].x), sizeof(double), 1, pDomain->pBaseRmnQ);
					fread(&(pDomain->INT_RmnQ[i][n][m].y), sizeof(double), 1, pDomain->pBaseRmnQ);
				}
				fseek(pDomain->pBaseRmnU, (p+1)*(p+2)*sizeof(double), SEEK_CUR);

				tempQ.x=0.0; tempQ.y=0.0;
				for(n=0; n<=p; n++)
				{
					for(m=0; m<=n; m++)
						tempQ+=pDomain->INT_RmnQ[i][n][m]*arryLeafPoint[k]->Lmn[n][m];

					for(m=1; m<=n; m++)
						tempQ+=conjugate(pDomain->INT_RmnQ[i][n][m])*conjugate(arryLeafPoint[k]->Lmn[n][m]);
				}
				if(im(tempQ)/re(tempQ) > 0.0000001)
				{
					ErrDisplay("Wrong for far field");
					return 1;
				}
				pDomain->vtQ.vt[arryLeafPoint[k]->sequence[i]]+=re(tempQ);
			}
			
			else //pDomain->NodeArray[arryLeafPoint[k]->sequence[i]]->pFace->BC.known==none
			{
				fread(pDomain->localU[i], sizeof(double), arryLeafPoint[k]->tlLocalNode, pDomain->pBaseU);
				fread(pDomain->localQ[i], sizeof(double), arryLeafPoint[k]->tlLocalNode, pDomain->pBaseQ);
				pDomain->localQ[i][i]=pDomain->vtDiagonalQ[arryLeafPoint[k]->sequence[i]];
				for(j=0; j<arryLeafPoint[k]->tlLocalNode; j++)
				{
					pDomain->vtU.vt[arryLeafPoint[k]->sequence[i]]
						+=pDomain->localU[i][j]*pDomain->vtF.vt[arryLeafPoint[k]->sequence[j]];
					pDomain->vtQ.vt[arryLeafPoint[k]->sequence[i]]
						+=pDomain->localQ[i][j]*pDomain->vtF.vt[arryLeafPoint[k]->sequence[j]];
				}

				for(n=0; n<=p; n++)
				for(m=0; m<=n; m++)
				{
					fread(&(pDomain->INT_RmnU[i][n][m].x), sizeof(double), 1, pDomain->pBaseRmnU);
					fread(&(pDomain->INT_RmnU[i][n][m].y), sizeof(double), 1, pDomain->pBaseRmnU);
					fread(&(pDomain->INT_RmnQ[i][n][m].x), sizeof(double), 1, pDomain->pBaseRmnQ);
					fread(&(pDomain->INT_RmnQ[i][n][m].y), sizeof(double), 1, pDomain->pBaseRmnQ);
				}

				tempU.x=0.0; tempU.y=0.0;
				tempQ.x=0.0; tempQ.y=0.0;
				for(n=0; n<=p; n++)
				{
					for(m=0; m<=n; m++)
					{
						tempU+=pDomain->INT_RmnU[i][n][m]*arryLeafPoint[k]->Lmn[n][m];
						tempQ+=pDomain->INT_RmnQ[i][n][m]*arryLeafPoint[k]->Lmn[n][m];
					}

					for(m=1; m<=n; m++)
					{
						tempU+=conjugate(pDomain->INT_RmnU[i][n][m])*conjugate(arryLeafPoint[k]->Lmn[n][m]);
						tempQ+=conjugate(pDomain->INT_RmnQ[i][n][m])*conjugate(arryLeafPoint[k]->Lmn[n][m]);
					}
				}
				if(im(tempU)/re(tempU) > 0.0000001 || im(tempQ)/re(tempQ) > 0.0000001)
				{
					ErrDisplay("Wrong for far field");
					return 1;
				}
				pDomain->vtU.vt[arryLeafPoint[k]->sequence[i]]+=re(tempU);
				pDomain->vtQ.vt[arryLeafPoint[k]->sequence[i]]+=re(tempQ);
			}
		}
	}
	fclose(pDomain->pBaseU); pDomain->pBaseU=0;
	fclose(pDomain->pBaseQ); pDomain->pBaseQ=0;
	fclose(pDomain->pBaseRmnU); pDomain->pBaseRmnU=0;
	fclose(pDomain->pBaseRmnQ); pDomain->pBaseRmnQ=0;

	return 0;
}
#endif

#if 0
//建立矩阵与向量相乘A*X; called by CDomain::solveUnknown()
int CMultipole::build_AX_2()
{
	//vtF:A*X中的X
	//vtU:待求向量vtU＝mxU*X
	//vtQ:待求向量vtQ＝mxQ*X

	long i, j, k, m, n;

	resetFmmTree(fmmTreeRoot);   //树结构系数清零
    initLeafs(pDomain->vtF.vt);  //将初始向量x0集成到叶子节点
	UpWardPass(fmmTreeRoot);         //向上集成系数
    DownWardPass(fmmTreeRoot);		 //向下集成系数

	for(i=0; i<tlNode; i++)
	{
		pDomain->vtU.vt[i]=0.0;
		pDomain->vtQ.vt[i]=0.0;
	}

	if((pDomain->pBaseU = fopen( pDomain->strBaseU, "rb")) == NULL ||
	   (pDomain->pBaseQ = fopen( pDomain->strBaseQ, "rb")) == NULL ||
	   (pDomain->pBaseRmnU = fopen( pDomain->strBaseRmnU, "rb")) == NULL ||
	   (pDomain->pBaseRmnQ = fopen( pDomain->strBaseRmnQ, "rb")) == NULL)
	{
		ErrDisplay("Cannot open database files");
		ErrDisplay("build_AX fail");
		return 1;
	}

	for(k=0; k<tlLeaf; k++)
	{
		for(i=0; i<arryLeafPoint[k]->tlNode; i++)
		{
			complex tempU, tempQ;
			if(pDomain->NodeArray[arryLeafPoint[k]->sequence[i]]->pFace->BC.known==flux)
			{
				fread(pDomain->localU[i], sizeof(double), arryLeafPoint[k]->tlLocalNode, pDomain->pBaseU);
				fseek(pDomain->pBaseQ, arryLeafPoint[k]->tlLocalNode*sizeof(double), SEEK_CUR);
				for(j=0; j<arryLeafPoint[k]->tlLocalNode; j++)
					pDomain->vtU.vt[arryLeafPoint[k]->sequence[i]]
						+=pDomain->localU[i][j]*pDomain->vtF.vt[arryLeafPoint[k]->sequence[j]];

				for(n=0; n<=p; n++)
				for(m=0; m<=n; m++)
				{
					fread(&(pDomain->INT_RmnU[i][n][m].x), sizeof(double), 1, pDomain->pBaseRmnU);
					fread(&(pDomain->INT_RmnU[i][n][m].y), sizeof(double), 1, pDomain->pBaseRmnU);
				}
				fseek(pDomain->pBaseRmnQ, (p+1)*(p+2)*sizeof(double), SEEK_CUR);

				tempU.x=0.0; tempU.y=0.0;
				for(n=0; n<=p; n++)
				{
					for(m=0; m<=n; m++)
						tempU+=pDomain->INT_RmnU[i][n][m]*arryLeafPoint[k]->Lmn[n][m];

					for(m=1; m<=n; m++)
						tempU+=conjugate(pDomain->INT_RmnU[i][n][m])*conjugate(arryLeafPoint[k]->Lmn[n][m]);
				}
				if(im(tempU)/re(tempU) > 0.0001)
				{
					ErrDisplay("Wrong for far field");
					return 1;
				}
				pDomain->vtU.vt[arryLeafPoint[k]->sequence[i]]+=re(tempU);
			}

			else if(pDomain->NodeArray[arryLeafPoint[k]->sequence[i]]->pFace->BC.known==potential)
			{
				fread(pDomain->localQ[i], sizeof(double), arryLeafPoint[k]->tlLocalNode, pDomain->pBaseQ);
				pDomain->localQ[i][i]=pDomain->vtDiagonalQ[arryLeafPoint[k]->sequence[i]];
				fseek(pDomain->pBaseU, arryLeafPoint[k]->tlLocalNode*sizeof(double), SEEK_CUR);
				for(j=0; j<arryLeafPoint[k]->tlLocalNode; j++)
					pDomain->vtQ.vt[arryLeafPoint[k]->sequence[i]]
						+=pDomain->localQ[i][j]*pDomain->vtF.vt[arryLeafPoint[k]->sequence[j]];

				for(n=0; n<=p; n++)
				for(m=0; m<=n; m++)
				{
					fread(&(pDomain->INT_RmnQ[i][n][m].x), sizeof(double), 1, pDomain->pBaseRmnQ);
					fread(&(pDomain->INT_RmnQ[i][n][m].y), sizeof(double), 1, pDomain->pBaseRmnQ);
				}
				fseek(pDomain->pBaseRmnU, (p+1)*(p+2)*sizeof(double), SEEK_CUR);

				tempQ.x=0.0; tempQ.y=0.0;
				for(n=0; n<=p; n++)
				{
					for(m=0; m<=n; m++)
						tempQ+=pDomain->INT_RmnQ[i][n][m]*arryLeafPoint[k]->Lmn[n][m];

					for(m=1; m<=n; m++)
						tempQ+=conjugate(pDomain->INT_RmnQ[i][n][m])*conjugate(arryLeafPoint[k]->Lmn[n][m]);
				}
				if(im(tempQ)/re(tempQ) > 0.0000001)
				{
					ErrDisplay("Wrong for far field");
					return 1;
				}
				pDomain->vtQ.vt[arryLeafPoint[k]->sequence[i]]+=re(tempQ);
			}
			
			else //pDomain->NodeArray[arryLeafPoint[k]->sequence[i]]->pFace->BC.known==none
			{
				fread(pDomain->localU[i], sizeof(double), arryLeafPoint[k]->tlLocalNode, pDomain->pBaseU);
				fread(pDomain->localQ[i], sizeof(double), arryLeafPoint[k]->tlLocalNode, pDomain->pBaseQ);
				pDomain->localQ[i][i]=pDomain->vtDiagonalQ[arryLeafPoint[k]->sequence[i]];
				for(j=0; j<arryLeafPoint[k]->tlLocalNode; j++)
				{
					pDomain->vtU.vt[arryLeafPoint[k]->sequence[i]]
						+=pDomain->localU[i][j]*pDomain->vtF.vt[arryLeafPoint[k]->sequence[j]];
					pDomain->vtQ.vt[arryLeafPoint[k]->sequence[i]]
						+=pDomain->localQ[i][j]*pDomain->vtF.vt[arryLeafPoint[k]->sequence[j]];
				}

				for(n=0; n<=p; n++)
				for(m=0; m<=n; m++)
				{
					fread(&(pDomain->INT_RmnU[i][n][m].x), sizeof(double), 1, pDomain->pBaseRmnU);
					fread(&(pDomain->INT_RmnU[i][n][m].y), sizeof(double), 1, pDomain->pBaseRmnU);
					fread(&(pDomain->INT_RmnQ[i][n][m].x), sizeof(double), 1, pDomain->pBaseRmnQ);
					fread(&(pDomain->INT_RmnQ[i][n][m].y), sizeof(double), 1, pDomain->pBaseRmnQ);
				}

				tempU.x=0.0; tempU.y=0.0;
				tempQ.x=0.0; tempQ.y=0.0;
				for(n=0; n<=p; n++)
				{
					for(m=0; m<=n; m++)
					{
						tempU+=pDomain->INT_RmnU[i][n][m]*arryLeafPoint[k]->Lmn[n][m];
						tempQ+=pDomain->INT_RmnQ[i][n][m]*arryLeafPoint[k]->Lmn[n][m];
					}

					for(m=1; m<=n; m++)
					{
						tempU+=conjugate(pDomain->INT_RmnU[i][n][m])*conjugate(arryLeafPoint[k]->Lmn[n][m]);
						tempQ+=conjugate(pDomain->INT_RmnQ[i][n][m])*conjugate(arryLeafPoint[k]->Lmn[n][m]);
					}
				}
				if(im(tempU)/re(tempU) > 0.0000001 || im(tempQ)/re(tempQ) > 0.0000001)
				{
					ErrDisplay("Wrong for far field");
					return 1;
				}
				pDomain->vtU.vt[arryLeafPoint[k]->sequence[i]]+=re(tempU);
				pDomain->vtQ.vt[arryLeafPoint[k]->sequence[i]]+=re(tempQ);
			}
		}
	}
	fclose(pDomain->pBaseU); pDomain->pBaseU=0;
	fclose(pDomain->pBaseQ); pDomain->pBaseQ=0;
	fclose(pDomain->pBaseRmnU); pDomain->pBaseRmnU=0;
	fclose(pDomain->pBaseRmnQ); pDomain->pBaseRmnQ=0;

	return 0;
}
#endif 

void CMultipole::PrintMmn(FMMTREENODE* fmmTreeRoot)
{
	int i, j, k, l, m, n;
	FILE *fout = NULL;
	char filename[256];
	memset(filename, 0, sizeof(filename));
	sprintf(filename, "Mmn_%d_%d.out",fmmTreeRoot->level,fmmTreeRoot->tlNode);

	if ((fout =fopen(filename,"w")) == NULL )
	{
		spdlog::info("Error: can not open file %s!\n");
		return;
	}

	for(n=0; n<=p; n++)
	{
		for(m=0; m<=n; m++)
		{
			if(std::imag(fmmTreeRoot->Mmn[0][n][m]) >= 0)
				fprintf(fout, "%10.8lf+%10.8lf  ", std::real(fmmTreeRoot->Mmn[0][n][m]), std::imag(fmmTreeRoot->Mmn[0][n][m]));
			else
				fprintf(fout, "%10.8lf%10.8lf  ", std::real(fmmTreeRoot->Mmn[0][n][m]), std::imag(fmmTreeRoot->Mmn[0][n][m]));
		}
		fprintf(fout, "\n");
	}
	fprintf(fout, "\n");

	fclose(fout);
	fout = NULL;
}

void CMultipole::PrintLmn(FMMTREENODE* fmmTreeRoot)
{
	int i, j, k, l, m, n;
	char filename[256];
	memset(filename, 0, sizeof(filename));
	sprintf(filename, "Lmn_%d_%d.out",fmmTreeRoot->level,fmmTreeRoot->tlNode);

	FILE *fout = NULL;
	if ((fout =fopen(filename,"w")) == NULL )
	{
		spdlog::info("Error: can not open file %s!\n");
		return;
	}

	for(n=0; n<=p; n++)
	{
		for(m=0; m<=n; m++)
		{
			if(std::imag(fmmTreeRoot->Lmn[0][n][m]) >= 0)
				fprintf(fout, "%10.8lf+%10.8lf  ", std::real(fmmTreeRoot->Lmn[0][n][m]), std::imag(fmmTreeRoot->Lmn[0][n][m]));
			else
				fprintf(fout, "%10.8lf%10.8lf  ", std::real(fmmTreeRoot->Lmn[0][n][m]), std::imag(fmmTreeRoot->Lmn[0][n][m]));
		}
		fprintf(fout, "\n");
	}
	fprintf(fout, "\n");

	fclose(fout);
	fout = NULL;
}
