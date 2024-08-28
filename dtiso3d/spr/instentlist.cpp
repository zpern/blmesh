#include "instentlist.h"
#include <spdlog/spdlog.h> 
 #include <stdio.h>
#include <assert.h>
#include <memory.h>
#include <stdlib.h>

/* 将其组织成哈希表 */
InstEntList::InstEntList()
{
	m_pNodeHash = NULL;
	m_nHashMaxIdx = 0;
	m_nInstEntArrSize = 1; /* valid index starts from 1 */
	m_nValidInstEnts = 0;
	
	m_bLocked = false;
	m_nLockedEntArrSize = 0;
}

InstEntList::~InstEntList()
{
	if (m_pNodeHash) {
		delete[] m_pNodeHash;
		m_pNodeHash = nullptr;
	}
}

/* 初始化 */
int InstEntList::initialise(int maxIdx)
{
	//memset(m_pNodeHash, 0, sizeof(int)*m_nHashMaxIdx);
	if (m_nHashMaxIdx < maxIdx)
	{
		int newMaxIdx = maxIdx * 2; //#zyj 添加一个额外容积
		if (m_pNodeHash)
			delete[] m_pNodeHash;
		m_pNodeHash = new int[newMaxIdx]();
		m_nHashMaxIdx = newMaxIdx;
		//memset(&m_pNodeHash, 0, MAX_INST_ENT_SIZE * sizeof(InstEntityLnk));
	}
	//#zyj
	else {
		if (addedNode.size() > maxIdx / 100)
			memset(m_pNodeHash, 0, maxIdx * sizeof(int));
		else
			for (const auto& i : addedNode) {
				m_pNodeHash[i] = 0;
			}
	}
	addedNode.clear();

	m_nInstEntArrSize = 1; /* valid index starts from 1 */
	m_nValidInstEnts = 0;
	
	
	return 1;
}

/* 初始化 */
int InstEntList::resetDirtNodes()
{
	int i, dirtNodeId, dirtNodeSize = m_vecDirtNodes.size();
	
	for (i = 0; i < dirtNodeSize; i++)
	{
		dirtNodeId = m_vecDirtNodes[i];
		assert(dirtNodeId >= 0 && dirtNodeId < m_nHashMaxIdx);
		m_pNodeHash[dirtNodeId] = 0;
	}
	m_nInstEntArrSize = 1; /* valid index starts from 1 */
	m_nValidInstEnts = 0;
	m_vecDirtNodes.clear();
	//memset(&m_arrInstEntLnks[0], 0, MAX_INST_ENT_SIZE*sizeof(InstEntityLnk));
	
	return 1;
}
/* 查询、增加和减少 */
int InstEntList::findEdge(int i1, int i2)
{
	int iKey, ii2;
	int iEntIt;
	InstEntityLnk *pEntLnk = NULL;

	iKey = i1 < i2 ? i1 : i2;
	ii2  = i1 < i2 ? i2 : i1;

	assert(iKey < m_nHashMaxIdx);
	iEntIt = m_pNodeHash[iKey];
	while (iEntIt > 0)
	{
		pEntLnk = &m_arrInstEntLnks[iEntIt];
		assert(pEntLnk->instEnt.i1 == iKey);
		if (pEntLnk->instEnt.type == IET_EDG && pEntLnk->instEnt.i2 == ii2)
			return iEntIt;
		iEntIt = pEntLnk->next;
	}

	return 0;
}

int InstEntList::findFace(int i1, int i2, int i3)
{
	int iKey, ii2, ii3;
	int iEntIt;
	InstEntityLnk *pEntLnk = NULL;

	if (i1 < i2 && i1 < i3)
	{
		iKey = i1;
		ii2  = i2;
		ii3  = i3;
	}
	else if (i2 < i1 && i2 < i3)
	{
		iKey = i2;
		ii2  = i3;
		ii3  = i1;
	}
	else
	{
		iKey = i3;
		ii2  = i1;
		ii3  = i2;
	}

	assert(iKey < m_nHashMaxIdx);
	iEntIt = m_pNodeHash[iKey];
	while (iEntIt > 0)
	{
		pEntLnk = &m_arrInstEntLnks[iEntIt];
		if (pEntLnk->instEnt.type == IET_FAC && (
			(pEntLnk->instEnt.i2 == ii2 && pEntLnk->instEnt.i3 == ii3) ||
			(pEntLnk->instEnt.i2 == ii3 && pEntLnk->instEnt.i3 == ii2)))
			return iEntIt;
		iEntIt = pEntLnk->next;
	}

	return 0;
}

int InstEntList::insertEdge(int i1, int i2)
{
	int iKey, ii2;
	int iEnt = m_nInstEntArrSize;
	int iEntHead;

	if (iEnt >= MAX_INST_ENT_SIZE)
	{
		printf("Increase MAX_INST_ENT_SIZE(%d), please.\n ", MAX_INST_ENT_SIZE);
		exit(1);
	}

	assert(iEnt < MAX_INST_ENT_SIZE);
	iKey = i1 < i2 ? i1 : i2;
	ii2  = i1 < i2 ? i2 : i1;
	assert(iKey < m_nHashMaxIdx);

	memset(&m_arrInstEntLnks[iEnt], 0, sizeof(InstEntityLnk));
	m_arrInstEntLnks[iEnt].instEnt.type = IET_EDG;
	m_arrInstEntLnks[iEnt].instEnt.i1 = iKey;
	m_arrInstEntLnks[iEnt].instEnt.i2 = ii2;
	m_arrInstEntLnks[iEnt].instEnt.idx = iEnt;
	

	/* 加入到哈希表 */
	iEntHead = m_pNodeHash[iKey];
	m_pNodeHash[iKey] = iEnt;
	addedNode.push_back(iKey);
	if (iEntHead > 0)
		m_arrInstEntLnks[iEnt].next = iEntHead;
	m_nInstEntArrSize++;
	m_nValidInstEnts++;

	m_vecDirtNodes.push_back(iKey); /* add a dirted node */
	return iEnt;
}

int InstEntList::insertFace(int i1, int i2, int i3)
{
	int iKey, ii2, ii3;
	int iEnt = m_nInstEntArrSize;
	int iEntHead;

	if (iEnt >= MAX_INST_ENT_SIZE)
	{
		printf("Increase MAX_SPR_POLY_SIZE(%d), please. insertFace\n", MAX_INST_ENT_SIZE);
		exit(1);
	}
	assert(iEnt < MAX_INST_ENT_SIZE);
	if (i1 < i2 && i1 < i3)
	{
		iKey = i1;
		ii2  = i2;
		ii3  = i3;
	}
	else if (i2 < i1 && i2 < i3)
	{
		iKey = i2;
		ii2  = i3;
		ii3  = i1;
	}
	else
	{
		iKey = i3;
		ii2  = i1;
		ii3  = i2;
	}
	assert(iKey < m_nHashMaxIdx);

	memset(&m_arrInstEntLnks[iEnt], 0, sizeof(InstEntityLnk));
	m_arrInstEntLnks[iEnt].instEnt.type = IET_FAC;
	m_arrInstEntLnks[iEnt].instEnt.i1 = iKey;
	m_arrInstEntLnks[iEnt].instEnt.i2 = ii2;
	m_arrInstEntLnks[iEnt].instEnt.i3 = ii3;
	m_arrInstEntLnks[iEnt].instEnt.idx = iEnt;
	
	/* 加入到哈希表 */
	iEntHead = m_pNodeHash[iKey];
	m_pNodeHash[iKey] = iEnt;
	addedNode.push_back(iKey);
	if (iEntHead > 0)
		m_arrInstEntLnks[iEnt].next = iEntHead;
	m_nInstEntArrSize++;
	m_nValidInstEnts++;

	m_vecDirtNodes.push_back(iKey); /* add a dirted node */
	return iEnt;
}

int InstEntList::removeEdge(int i1, int i2)
{
	int iKey, ii2;
	int iEnt = m_nInstEntArrSize;
	int iEntIt, iEntPrev;
	InstEntityLnk *pEntLnk = NULL;

	iKey = i1 < i2 ? i1 : i2;
	ii2  = i1 < i2 ? i2 : i1;

	assert(iKey < m_nHashMaxIdx);
	iEntIt = m_pNodeHash[iKey];
	iEntPrev = 0;
	while (iEntIt > 0)
	{
		pEntLnk = &m_arrInstEntLnks[iEntIt];
		if (pEntLnk->instEnt.type == IET_EDG && pEntLnk->instEnt.i2 == ii2)
		{
			pEntLnk->instEnt.flag |= IET_INVALID_FLAG; /* 设置为无效 */
			
			/* 从哈希表中移走 */
			if (iEntPrev > 0)
				m_arrInstEntLnks[iEntPrev].next = pEntLnk->next;
			else
				m_pNodeHash[iKey] = pEntLnk->next;
			
			m_nValidInstEnts--;
			return iEntIt;
		}
		iEntPrev = iEntIt;
		iEntIt = pEntLnk->next;
	}
	
	return 0;
}

int InstEntList::removeFace(int i1, int i2, int i3)
{
	int iKey, ii2, ii3;
	int iEnt = m_nInstEntArrSize;
	int iEntIt, iEntPrev;
	InstEntityLnk *pEntLnk = NULL;

	if (i1 < i2 && i1 < i3)
	{
		iKey = i1;
		ii2  = i2;
		ii3  = i3;
	}
	else if (i2 < i1 && i2 < i3)
	{
		iKey = i2;
		ii2  = i3;
		ii3  = i1;
	}
	else
	{
		iKey = i3;
		ii2  = i1;
		ii3  = i2;
	}

	assert(iKey < m_nHashMaxIdx);
	iEntIt = m_pNodeHash[iKey];
	iEntPrev = 0;
	while (iEntIt > 0)
	{
		pEntLnk = &m_arrInstEntLnks[iEntIt];
		if (pEntLnk->instEnt.type == IET_FAC && (
			(pEntLnk->instEnt.i2 == ii2 && pEntLnk->instEnt.i3 == ii3) ||
			(pEntLnk->instEnt.i2 == ii3 && pEntLnk->instEnt.i3 == ii2)))
		{
			pEntLnk->instEnt.flag |= IET_INVALID_FLAG; /* 设置为无效 */
			
			/* 从哈希表中移走 */
			if (iEntPrev > 0)
				m_arrInstEntLnks[iEntPrev].next = pEntLnk->next;
			else
				m_pNodeHash[iKey] = pEntLnk->next;

			m_nValidInstEnts--;
			return iEntIt;
		}
		iEntPrev = iEntIt;
		iEntIt = pEntLnk->next;
	}
	
	return 0;
}

int InstEntList::removeEnt(int iEnt)
{
	int iKey;
	int iEntIt, iEntPrev;
	InstEntityLnk *pEntLnk = NULL, *pEntLnkIt = NULL;

	pEntLnk = &m_arrInstEntLnks[iEnt];
	iKey = pEntLnk->instEnt.i1;
	iEntIt = m_pNodeHash[iKey];
	iEntPrev = 0;
	while (iEntIt > 0)
	{
		pEntLnkIt = &m_arrInstEntLnks[iEntIt];
		if (iEntIt == iEnt)
		{
			pEntLnkIt->instEnt.flag |= IET_INVALID_FLAG; /* 设置为无效 */
			
			/* 从哈希表中移走 */
			if (iEntPrev > 0)
				m_arrInstEntLnks[iEntPrev].next = pEntLnkIt->next;
			else
				m_pNodeHash[iKey] = pEntLnkIt->next;

			m_nValidInstEnts--;

			return 1;
		}
		iEntPrev = iEntIt;
		iEntIt =pEntLnkIt->next;
	}

	return 0;
}

/* 遍历所有InstEntities */
int InstEntList::getFirstEnt()
{
	int iEntIt = 0;
	int m_nIterSize = m_bLocked ? m_nLockedEntArrSize : m_nInstEntArrSize;

	if (m_nValidInstEnts > 0)
	{
		while (iEntIt < m_nIterSize - 1 && !isValid(++iEntIt));
	}
	return iEntIt;
}

int InstEntList::getNextEnt(int iEnt)
{
	int iEntIt = iEnt;
	int m_nIterSize = m_bLocked ? m_nLockedEntArrSize : m_nInstEntArrSize;

	assert(iEntIt < m_nIterSize);
	do 
	{
		iEntIt++;
	}
	while (iEntIt < m_nIterSize && !isValid(iEntIt));
	
	return iEntIt == m_nIterSize ? 0 : iEntIt;
}

bool InstEntList::isLocked()
{
	return m_bLocked;
}

void InstEntList::lock()
{
	m_bLocked = true;
	m_nLockedEntArrSize = m_nInstEntArrSize;
}

void InstEntList::unlock()
{
	m_bLocked = false;
	m_nLockedEntArrSize = 0;
}

bool InstEntList::isEmpty()
{
	return m_nValidInstEnts <= 0;
}

int InstEntList::getSize()
{
	return m_nValidInstEnts;
}

/* 状态标志 */
void InstEntList::setRecoverable(int iEnt, bool flag)
{
	flag ? m_arrInstEntLnks[iEnt].instEnt.flag &= ~IET_IRRECOVERABLE_FLAG :
		   m_arrInstEntLnks[iEnt].instEnt.flag |=  IET_IRRECOVERABLE_FLAG;
}

void InstEntList::setValid(int iEnt, bool flag)
{
	flag ? m_arrInstEntLnks[iEnt].instEnt.flag &= ~IET_INVALID_FLAG :
		   m_arrInstEntLnks[iEnt].instEnt.flag |=  IET_INVALID_FLAG;
}

void InstEntList::setActive(int iEnt, bool flag)
{
	flag ? m_arrInstEntLnks[iEnt].instEnt.flag &= ~IET_INACTIVE_FLAG :
		   m_arrInstEntLnks[iEnt].instEnt.flag |=  IET_INACTIVE_FLAG;
}

void InstEntList::setNewBorn(int iEnt, bool flag)
{
	flag ? m_arrInstEntLnks[iEnt].instEnt.flag |= IET_NEWBORN_FLAG :
		   m_arrInstEntLnks[iEnt].instEnt.flag &=  ~IET_NEWBORN_FLAG;
}

void InstEntList::setNoEdit(int iEnt, bool flag)
{
	flag ? m_arrInstEntLnks[iEnt].instEnt.flag |= IET_NOEDIT_FLAG :
		   m_arrInstEntLnks[iEnt].instEnt.flag &=  ~IET_NOEDIT_FLAG;
}

bool InstEntList::isRecoverable(int iEnt)
{
	return m_arrInstEntLnks[iEnt].instEnt.type > IET_NUL && !(m_arrInstEntLnks[iEnt].instEnt.flag & IET_IRRECOVERABLE_FLAG);
}

bool InstEntList::isValid(int iEnt)
{
	return m_arrInstEntLnks[iEnt].instEnt.type > IET_NUL && !(m_arrInstEntLnks[iEnt].instEnt.flag & IET_INVALID_FLAG);
}

bool InstEntList::isActive(int iEnt)
{
	return m_arrInstEntLnks[iEnt].instEnt.type > IET_NUL && !(m_arrInstEntLnks[iEnt].instEnt.flag & IET_INACTIVE_FLAG);
}

bool InstEntList::isNewBorn(int iEnt)
{
	return m_arrInstEntLnks[iEnt].instEnt.type > IET_NUL && (m_arrInstEntLnks[iEnt].instEnt.flag & IET_NEWBORN_FLAG);
}

bool InstEntList::isNoEdit(int iEnt)
{
	return m_arrInstEntLnks[iEnt].instEnt.type > IET_NUL && (m_arrInstEntLnks[iEnt].instEnt.flag & IET_NOEDIT_FLAG);
}

int InstEntList::counterIncre(int iEnt)
{
#if 0
	if (m_arrInstEntLnks[iEnt].instEnt.counter > 3)
		spdlog::info("Counter > 3; Should not be here.\n");
#endif
	return ++m_arrInstEntLnks[iEnt].instEnt.counter;
}
int InstEntList::counterDecre(int iEnt)
{
#if 0
	if (m_arrInstEntLnks[iEnt].instEnt.counter < 1)
		spdlog::info("Counter < 0; Should not be here.\n");
#endif
	return --m_arrInstEntLnks[iEnt].instEnt.counter;
}
int InstEntList::entCounter(int iEnt)
{
#if 0
	if (m_arrInstEntLnks[iEnt].instEnt.counter > 3)
		spdlog::info("Counter > 3; Should not be here.\n");
#endif
	return m_arrInstEntLnks[iEnt].instEnt.counter;
}

void InstEntList::setEntCounter(int iEnt, int entCnt)
{
	m_arrInstEntLnks[iEnt].instEnt.counter = entCnt;
}

InstEntType InstEntList::entType(int iEnt)
{
	return m_arrInstEntLnks[iEnt].instEnt.type;
}

int InstEntList::nodeIndices(int iEnt, int indices[])
{
	int type = entType(iEnt);
	switch (type)
	{
	case IET_FAC:
		indices[2] = m_arrInstEntLnks[iEnt].instEnt.i3;
	case IET_EDG:
		indices[1] = m_arrInstEntLnks[iEnt].instEnt.i2;
	case IET_NOD:
		indices[0] = m_arrInstEntLnks[iEnt].instEnt.i1;
		break;
	}

	return type;
}

