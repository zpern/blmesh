/**
 * @file BLFront.cpp
 * @author yhf (hfye@zju.edu.cn)
 * @brief define the basic front
 * @version 1.6
 * @date 2020-06-13
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#include <set>
#include <spdlog/spdlog.h> 
 #include "BLFront.h"
#ifdef  MEMORY_DEBUG
std::set<BLFront*> BLFront::my_vec = std::set<BLFront*>();
#endif
#ifdef USE_MEMORY_POOL
INIT(BLFront);
#endif
BLFront::BLFront(void) :
	interact_with_other(false),
	m_pPrevFront(nullptr),
	m_pNxtFront(nullptr),
	is_boundary_(true),
	m_pUpperFront(nullptr),
	m_pLowerFront(nullptr),
	m_bSymm(false),
	m_iSymmidx(-1),
	m_bPyramid(false),
	isVirtual(false),
	m_iPyramidTri(nullptr),
	outer_node(nullptr),
	isSimple(true),
	is_prism_valid(1),
	stridx({-1}),
	m_nNeiFront(0),
	m_dSize(-1){
#ifdef  MEMORY_DEBUG
	my_vec.insert(this);
#endif
}
void BLFront::RmvNeigbourNode(BLNode* blFront) {
	for (auto &i : m_pBLNods) {
		if (i == blFront)
			i = nullptr;
	}
}
int BLFront::GetThirdNodIdx(BLNode* blNod1, BLNode* blNod2)
{
	int i;

	for (i = 0; i <3; i++)
	{
		if (blNod1 != m_pBLNods[i] && blNod2 != m_pBLNods[i])
			return i;
	}

	return -1;
}
void BLFront::ReplaceFrontNods(BLNode* blNod, BLNode* blNodNew)
{
	int i;

	for (i = 0; i < 3; i++)
	{
		if (blNod == m_pBLNods[i])
			m_pBLNods[i] = blNodNew;
	}
}
bool BLFront::IncludeNode(BLNode* blNod)
{
	int i;

	for (i = 0; i <3; i++)
	{
		if (blNod == m_pBLNods[i])
			return true;
	}

	return false;
}
BLFront::~BLFront(void)
{	
#ifdef  MEMORY_DEBUG
	my_vec.erase(this);
#endif
	if (m_iPyramidTri) {
		delete[]m_iPyramidTri; m_iPyramidTri = nullptr;
	}
	if (m_pNxtFront) {
		m_pNxtFront->m_pPrevFront = nullptr;
	}
	if (m_pPrevFront) {
		m_pPrevFront->m_pNxtFront = nullptr;
	}
	for (int i = 0; i < 3; i++)
		if (m_pBLNods[i])
			m_pBLNods[i]->AddDeleteCount();
	//cout << "Destory Front" << endl;
	m_nNeiFront = 0;
}





void BLFront::SetNormal(BLVector normal)
{
	m_normal = normal;
}
void BLFront::GetNeigbourFronts(int *nblFront, BLFront** blFronts)
{
	int i;
	*nblFront = m_nNeiFront;
	for (i = 0; i < *nblFront; i++)
	{
		blFronts[i] = m_arrNeigFronts[i];
	}
}

void BLFront::AddNeigbourFronts(BLFront* blFront)
{
	if (m_nNeiFront == 3)
		throw std::logic_error("non-manifold input surface mesh! first node="+std::to_string((blFront->GetLowerFront()->m_pBLNods[0]->GetNodIdx())));
	m_arrNeigFronts[m_nNeiFront++]=blFront;
}

void BLFront::RmvNeigbourFronts(BLFront* blFront)
{
	int count = 0;
	for (int i = 0; i < m_nNeiFront; i++) {
		if (m_arrNeigFronts[i] != blFront) {
			m_arrNeigFronts[count++] = m_arrNeigFronts[i];
		}
	}
	m_nNeiFront = count;
}
void BLFront::CalNormal(MBLNode* pNodes)
{
	int iNod1 = m_pBLNods[0]->GetNodIdx();
	int iNod2 = m_pBLNods[1]->GetNodIdx();
	int iNod3 = m_pBLNods[2]->GetNodIdx();
	BLVector vectmp1(pNodes[iNod2].coord[0] - pNodes[iNod1].coord[0], pNodes[iNod2].coord[1] - pNodes[iNod1].coord[1], pNodes[iNod2].coord[2] - pNodes[iNod1].coord[2]);
	BLVector vectmp2(pNodes[iNod3].coord[0] - pNodes[iNod1].coord[0], pNodes[iNod3].coord[1] - pNodes[iNod1].coord[1], pNodes[iNod3].coord[2] - pNodes[iNod1].coord[2]);
	m_normal = (vectmp2 ^ vectmp1).normalized();
}
