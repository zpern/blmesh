#ifndef __BLMesh_BLNode_h__
#define __BLMesh_BLNode_h__

#include "BLDomain.h"
#include "BLFront.h"
#include "BLMesh_define.h"
#include "BLVector.h"
#include "LightMemoryPool.h"
#include "PotentialBEM.h"
#include "common.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <list>
#include <set>
#include <stdio.h>
#include <unordered_set>
#include <vector>
enum NODETYPE {
    CONVEX = 0,
    CONCAVE,
    PLAIN
};
class BLFront;

#define _CONFIG_FILE

typedef struct {
    int nfield;
    double *pt;
    double *fv;
    double *g;
    int nelm;
    int npt;
    double *bpt;
    double *bc;
    int *belm;
    double *norm;
    double *u;
} FieldArg;

class BLNode
#ifdef USE_MEMORY_POOL
    : public LightWightMomeryPool<BLNode>
#endif
{
public:
    BLNode(int iNod = -1)
        : delete_count(0)
        , m_BLLowerNode(nullptr)
        , m_BLUpperNode(nullptr)
        , m_dfixHightRtio(0)
        , m_iNod(iNod)
        , iso_stop(false)
        , expect_height_(-1)
        , distance_ratio(1.0)
        , m_pBLUpperNode(nullptr)
        , m_pBLLowerNode(nullptr)
        , m_pBLDescentNode(nullptr)
        , m_bStopPropagate(false)
        , adjacent(false)
        , isVirtualP(false)
        , m_Normal(0.0, 0.0, 0.0)
        , m_h0(cf.step_len)
        , bsys(false)
        , m_bBdry(false)
        , m_nNormContri(1)
        , m_iSymAxis(-1)
        , m_dHightRtio(0)
        , if_need_smooth(true)
    {
        m_vecNeighFronts.reserve(10);
#ifdef MEMORY_DEBUG
        my_vec.insert(this);
#endif
        per_node = nullptr;
    }

    ~BLNode(void);

public:
    /*
     * @brife 删除掉一个周围front，如果删除完了，则删除本指针（类似智能指针）
     * @author yhf
     */
    void AddDeleteCount();

    /// 设置周期性点
public:
    // add by yhf
    BLNode *per_node;
    enum PerType {
        Forward = 0,
        Reverse = 1
    } per_type;
    BLNode::PerType getPerNodeType() { return per_type; }
    BLNode *getPerNode() { return per_node; }
#ifdef MEMORY_DEBUG
    int getDelCount() { return delete_count; }
    static std::set<BLNode *> my_vec;
#endif
    void CalConcaveConvex(MBLNode *pNode);
    BLVector GetCoord(MBLNode *pNode);
    int GetNodIdx() { return m_iNod; }
    void SetNodIdx(int iNod)
    {
        assert(iNod >= 0);
        m_iNod = iNod;
    }

    double GetDistanceRatio() { return distance_ratio; }
    void SetDistanceRatio(double val) { distance_ratio = val; }
    void SetLayerNum(int iLayer) { m_ilayNum = iLayer; }
    int GetLayerNum() { return m_ilayNum; }

    inline bool HasUpperNode()
    {
        // return m_pBLUpperNode != nullptr;
        return m_BLUpperNode != nullptr;
    }

    inline bool HasLowerNode()
    {
        // return m_pBLLowerNode != nullptr;
        return m_BLLowerNode != nullptr;
    }

    // need to fix
    inline BLNode *GetUpperNode()
    {
        // return m_pBLLowerNode;
        return m_BLUpperNode;
    }

    // need to fix
    inline BLNode *GetLowerNode()
    {
        // return m_pBLLowerNode;
        return m_BLLowerNode;
    }

    // delete
    int UpperNodesCount() { return -1; }
    // delete
    int LowerNodesCount() { return -1; }

    void SetUpperNode(BLNode *blNode) { m_pBLUpperNode = blNode; }

    void SetLowerNode(BLNode *blNode) { m_pBLLowerNode = blNode; }

    void ReplaceUpperNode(BLNode *blNode, BLNode *blNodNew) {}

    void ReplacetLowerNode(BLNode *blNode, BLNode *blNodNew) {}

    void AddUpperNode(BLNode *blNode) { m_BLUpperNode = blNode; }

    void AddLowerNode(BLNode *blNode) { m_BLLowerNode = blNode; }
#if 0
	//left-right
	bool HasLeftNode()
	{
		return m_pBLLeftNode != nullptr;
	}

	bool HasRightNode()
	{
		return m_pBLRightNode != nullptr;
	}

	BLNode* GetLeftNode()
	{
		return m_pBLLeftNode;
	}

	BLNode* GetRightNode()
	{
		return m_pBLRightNode;
	}

	void SetLeftNode(BLNode* blNode)
	{
		m_pBLLeftNode = blNode;
	}

	void SetRightNode(BLNode* blNode)
	{
		m_pBLRightNode = blNode;
	}
#endif

    void AddNeigFronts(BLFront *blFront) { m_vecNeighFronts.push_back(blFront); }

    int NeighFrontCount() { return m_vecNeighFronts.size(); }
    std::vector<BLFront *> &GetNeigFronts() { return m_vecNeighFronts; }
    void GetNeigFronts(BLFront **blFronts, int *nBLFront)
    {
        int i = 0;

        for (i = 0; i < m_vecNeighFronts.size(); i++) {
            blFronts[i] = m_vecNeighFronts[i];
        }
        *nBLFront = (int)m_vecNeighFronts.size();
    }
    bool checkValid()
    {
        if (delete_count == m_vecNeighFronts.size()) {
            return true;
        }
        return false;
    }
    void RmvNeigFront(BLFront *blFront)
    {
        std::vector<BLFront *>::iterator it;
        it = find(m_vecNeighFronts.begin(), m_vecNeighFronts.end(), blFront);

        if (it != m_vecNeighFronts.end()) {
            m_vecNeighFronts.erase(it);
        }
    }

    vector<BLNode *> &GetNeigNods();
    // void GetNeigNods(BLNode** blNods, int *nBLNod, bool flag = false);

    void RmvUpperNod(BLNode *blNod)
    {
        if (m_BLUpperNode == blNod) {
            m_BLUpperNode = nullptr;
        }
    }

    void RmvLowerNod(BLNode *blNod)
    {
        if (m_pBLLowerNode == blNod) {
            m_pBLLowerNode = nullptr;
        }
    }

    // Retrieve the calculated normal
    const BLVector &GetNormal() const;
    void SetNormal(BLVector vec);

    void GetRealNeigNods(BLNode **blNods, int *nBLNod, bool flag = false);
    void GetVirtualNeigNods(BLNode **blNods, int *nBLNod, MBLNode *pnode, bool flag = false);

    // To calculate the normal by averaging the normals of surrounding fronts
    BLVector GetNormal(MBLNode *pNodes);

    BLVector GetNormal(MBLNode *pNodes, int flag);

    bool GetVirtualFlag() { return isVirtualP; };
    void SetVirtualFlag(bool isv) { isVirtualP = isv; }
    double GetBeitaVisu(BLVector normal);
    double GetBeitaVisu();
    void SetBeitaVisu(double val);

    /*迭代求解，代价最大*/
    BLVector GetTheMostnNormalNormal(MBLNode *pNodes);
    /*几何中心方法求解法向，比较中庸*/
    BLVector GetCenterNormal(MBLNode *pNodes);

    /*简单平均*/
    BLVector GetSimpleNormal(MBLNode *pNodes);

    /*简单分类平均*/
    BLVector GetNaiveNormal(MBLNode *pNodes);

    /*求圆心的方法求法向，应该是理论最优解*/
    BLVector GetCirculeDotNormal(MBLNode *pNodes);

    BLVector GetGeometryNormal(MBLNode *pNodes);

    // To calculate the normal
    /*
     * @pPotBEM: the object of PotentialBEM
     * @pNodes: the coordinates of grid points
     * @bpt: the coordinates of boundary nodes
     * @npt: the number of boundary nodes
     * @belm: the connectivity of boundary elements
     * @nelm: the number of boundary elements
     * @bc: the boundary conditions & values
     * @u: the boundary values(potential value or flux value)
     */

    BLVector GetHeight();
    BLVector GetHeight(double aver);

    // add by yhf
    void SetNodeType(NODETYPE nt) { nodet = nt; }
    NODETYPE GetNodeType() { return nodet; }

    void SetHeightLength(double val) { expect_height_ = val; }
    double GetHeightLength();

    void SetStopFlag(bool status = true);
    bool GetStopFlag() { return m_bStopPropagate; }

    void SetBSys(bool flag, int isym = -1)
    {
        bsys = flag;
        m_iSymAxis = isym;
    }
    bool GetBSys() { return bsys; }
    int GetSymAxis() { return m_iSymAxis; }

    inline int GetDecentID() { return desentID; }
    inline void SetDecentID(int val) { desentID = val; }
    void SetBdryPt(bool flag) { m_bBdry = flag; }
    bool GetBdryPt() { return m_bBdry; }

    void SetNormCntribt(int nCntri) { m_nNormContri = nCntri; }
    bool IsComplexNode() { return m_nNormContri > 3; }

    void SetDescentNode(BLNode *blNod) { m_pBLDescentNode = blNod; }
    BLNode *GetDescentNod() { return m_pBLDescentNode ? m_pBLDescentNode : this; }

    // set and get added increase ratio
    void SetHightRatio(double vrtio) { m_dHightRtio = vrtio; }
    double GetHightRatio() { return m_dHightRtio; }

    void SetFixedHightRatio(double vrtio)
    {
        if (std::isnan(vrtio) || std::isinf(vrtio)) {
            cout << "123";
        }
        m_dfixHightRtio = vrtio;
    }
    double GetFixedHightRatio() { return m_dfixHightRtio; }

    void setPerNodeType(PerType type)
    {
        per_type = type;
        if (GetLowerNode()) {
            if (GetLowerNode()->per_type != per_type) {}
        }
    }
    void setPerNode(BLNode *blnode, PerType isforward)
    {

        per_node = blnode;
        per_type = isforward;

        if (GetLowerNode()) {
            if (GetLowerNode()->per_type != per_type) {}
        }
    }

public:
    bool adjacent;
    double adjacent_value;

private:
    double distance_ratio;

    enum NODETYPE nodet;
    short m_ilayNum;

    bool m_bStopPropagate;
    bool isVirtualP;

    // double m_uValue;
    // double m_averU;
    // double m_uVpre;
    // double m_uVprt;
    // double m_gValue;

    double expect_height_;

    BLVector m_Normal;
    BLNode *m_pBLUpperNode;
    BLNode *m_pBLLowerNode;
    BLNode *m_BLUpperNode;
    BLNode *m_BLLowerNode;

    BLNode *m_pBLDescentNode;

    std::vector<BLFront *> m_vecNeighFronts;
    std::vector<BLNode *> m_vecNeighNodes;

    /****************/
    // for testing vector field
    // double m_uVx;
    // double m_uVy;
    // double m_uVz;
    // double m_uV;
    // double m_gX[MAX_NORMAL_COMPONENT];
    // double m_gY[MAX_NORMAL_COMPONENT];
    /****************/

    // index of mesh node
    int m_iNod;
    bool bsys;
    bool m_bBdry;

    int desentID;

    short delete_count;
    short m_nNormContri; // m_nNormContri>2: complex node
    short m_iSymAxis;    // x: 0, y: 1, z: 2, default: -1;

    double m_dHightRtio;
    double m_dfixHightRtio;

public:
    double ave_length;
    int m_symmplane_id;
    double m_h0;
    int respect_layer;
    double respect_height;
    double respect_ratio;
    double recommend_height;
    double beitaVisu;
    bool iso_stop;
    bool if_need_smooth;
};

#endif
