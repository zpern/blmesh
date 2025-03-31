#ifndef __BLMesh_BLFront_h__
#define	__BLMesh_BLFront_h__

#include "BLMesh_define.h"
#include "vector.h"
#include "BLNode.h"
#include "Octree.h"
#include "DynamicArray.h"
#include "LightMemoryPool.h"
#include <array>
struct PyramidSet{
	int expose_idx[3];
	int exposed_layer[3];
	int node_idx[3];
};
class BLNode;
class BLFrontList;
using namespace OCT;
class BLFront
#ifdef USE_MEMORY_POOL
	:public LightWightMomeryPool<BLFront>
#endif
{
	friend class BLFrontList;

public:
	BLFront(void);
	~BLFront(void);

public:

#ifdef  MEMORY_DEBUG
	
	static std::set<BLFront*> my_vec;
#endif


	BLVector GetNormal() {
		return m_normal;
	};
	void clearTopo() {
		for (int i = 0; i < 3; i++)
		{
			m_pBLNods[i] = nullptr;
		}
		m_nNeiFront = 0;
	}
	std::array<BLNode*,3> GetNodes() {
		std::array<BLNode*, 3> ans;

		for (int i = 0; i < 3; i++)
		{
			ans[i] = m_pBLNods[i];
		}
		return ans;
	}
	void GetNodes(int* nNods, BLNode** blNods) {
		int i = 0;
		*nNods = 3;

		for (i = 0; i < *nNods; i++)
		{
			blNods[i] = m_pBLNods[i];
		}
	}
	void GetNeigbourFronts(int *nblFront, BLFront** blFronts);
	void AddNeigbourFronts(BLFront* blFront);
	void RmvNeigbourFronts(BLFront* blFront);

	void SetLayerNum(int iLayer){ m_ilayNum = iLayer; }
	int GetLayerNum(){ return m_ilayNum; }
	

	void SetUpperFront(BLFront* blFront){ m_pUpperFront = blFront; }
	inline BLFront*& GetUpperFront(){ return m_pUpperFront; }

	void SetLowerFront(BLFront* blFront) {
		if (blFront == nullptr) {
			m_pLowerFront = nullptr;
			return;
		}
	m_pLowerFront = blFront; }
	BLFront* GetLowerFront(){ return m_pLowerFront; }


	void SetFrontSize(double size){ m_dSize = size; }
	double GetFrontSize(){ return m_dSize; }

	void SetSqrtFrontSize(double size) { m_dSqrtSize = size; }
	double GetSqrtFrontSize() { return m_dSqrtSize; }


	void SetMinFrontSize(double size) { m_dminSize = size; }
	double GetMinFrontSize() { return m_dminSize; }

	void AddBLFrontNods(int i, BLNode* blNod){ m_pBLNods[i] = blNod; }

	 void ReplaceFrontNods(BLNode* blNod, BLNode* blNodNew) ;
	 bool IncludeNode(BLNode* blNod) ;
	 void SetNormal(BLVector normal);
	 void CalNormal(MBLNode* pNodes);
	inline void SetElmIdx(int idx) { this->eidx = idx; }
	inline int GetElmIdx() { return eidx; }
	inline void SetTriIdx(int idx) {this->tridx = idx;}

	

	inline TreeNode* GetOuterNode() { return  outer_node; }
	inline void SetOuterNode(TreeNode* node) { outer_node = node; }

	void SetSurfaceElmIdx(int idx) { this->surface_element_index_ = idx; }
	int GetSurfaceElmIdx() { return surface_element_index_; }


	inline int GetTriIdx(){return tridx;}

	void RmAllNode() {
		for (int i = 0; i < MAX_FRONT_NODES; i++)
			m_pBLNods[i] = nullptr;
	}
	 void RmvNeigbourNode(BLNode * blFront);
	 void SetSTriIdx(const int &nidx, const int &idx, const int &itri) {
		 stridx[nidx * 2 + idx] = itri;
	 }
	 int GetSTriIdx(const int &nidx, const int &idx) {
		 return stridx[nidx * 2 + idx];
	 }
	 int GetThirdNodIdx(BLNode* blNod1, BLNode* blNod2);
	void SetSymm(int idx) 
	{
		if(idx >= 0)
		{
			m_bSymm = true; 
			m_iSymmidx = idx; 
		}
	}

	int GetSymm() 
	{ 
		int idx = -1;
		if(m_bSymm)
			idx = m_iSymmidx;

		return idx;
	}
	inline bool IsSymm() {return m_bSymm;}
	inline bool IsSimple() { return isSimple; }
	inline void SetiSSimple(bool b) { isSimple=b; }
	inline bool IsPyramid() {return m_bPyramid;}
	void SetPyramidFlag(bool bPyramid) {m_bPyramid = bPyramid;}

	

public:
	PyramidSet ps;
	
	BLNode* m_pBLNods[MAX_FRONT_NODES];
	BLFront* m_pUpperFront;
	BLFront* m_pLowerFront;
	std::array<int, 6> conn;	//connection idx
	int is_prism_valid;  
	bool interact_with_other; 
	bool m_bPyramid;	//indicate whether there are attached pyramid triangles
	bool isVirtual;  //indicate whether the front is 
	bool m_bSymm;
	bool isSimple;  //simplenode
	bool is_boundary_;//表征是否为前沿

public:
	std::array<BLFront*,3> m_arrNeigFronts;
	short m_nNeiFront;
protected:
	
	short m_ilayNum;
	
	int eidx;
	int surface_element_index_;
	int tridx;	//for intersecting test: front triangle
	std::array<int,12> stridx;	//for intersecting test: side triangles
	int* m_iPyramidTri; //for intersecting test: pyramid triangles
	double m_dSize;
	double m_dminSize;
	double m_dSqrtSize;
	BLFront* m_pPrevFront;
	BLFront* m_pNxtFront;
	BLVector m_normal;
	TreeNode* outer_node;//用于求交快速判断在哪个node内，依据局部性
	
	int m_iSymmidx;
};

#endif
