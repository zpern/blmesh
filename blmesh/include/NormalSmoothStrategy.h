#pragma once
#include "symmetry.h"
#include "BLNode.h"
#include "BLMesh_define.h"
#include "blmeshapi.h"
#include <atomic>
#include <queue>
class NormalSmoothStrategy
{
public:
	NormalSmoothStrategy();
	NormalSmoothStrategy(BLNode **node, Node* pnodes,int num_front) ;
	void virtual SmoothNormal() ;
	/*二阶牛顿差值光滑化*/
	void horsmooth(BLNode *blNod);
	void virtual SmoothNormalOnce(BLNode *blNod, int id) = 0;
	void GetEnegy() ;
	virtual ~NormalSmoothStrategy();
	void SetBlNode(BLNode *node);
	void SetSmoothTimes(int s_times);
	void SetNumberofFrontNode(int num);
	void SetSymm(BLNode* blNod,int id);
	void SetFaceidSP(const std::map<int, TiGER::SymmetryPlane>& sp);
protected :
	std::map<int, TiGER::SymmetryPlane> faceid_to_sp;
	BLNode * * node_array;
	std::atomic_ulong energy;
	int smooth_times;
	//BLNode *blNod;
	Node* pNodes;
	int nFrtNods;
	vector<BLNode*> tmp;
	vector<BLVector> myvec;
};

