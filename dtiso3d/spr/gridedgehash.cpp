/* ----------------------------------------------------------------------------
 * ----------------------------------------------------------------------------
 *
 * 多学科应用模拟的赋能环境
 * Enabling Environment for Multi-displinary Application Simulations
 *
 * 陈建军 中国 浙江大学工程与科学计算研究中心
 * 版权所有	  2009年10月30日
 * Chen Jianjun  Center for Engineering & Scientific Computation,
 * Zhejiang University, P. R. China
 * Copyright reserved, Oct. 30, 2009
 * 
 * 联系方式 (For further information, please conctact)
 *   电话 (Tel)：+86-571-87953166
 *   传真 (Fax)：+86-571-87953167
 *   邮箱 (Mail)：chenjj@zju.edu.cn
 *
 * 文件名称 (File Name)：gridedgehash.h
 * 初始版本 (Initial Version): V1.0
 * 功能介绍 (Function Introduction：
 *     实现边表类GridEdgeHash (哈希表)，用于记录顶点之间的邻接关系
 *     Implement class GridEdgeHash, to record the neighboring info. between grid nodes
 *
 * -----------------------------修改记录 (Revision Record)------------------------
 * 修改者 (Revisor):
 * 修改日期 (Revision Date):
 * 当前版本 (Current Version):
 * 修改介绍 (Revision Introduction):
 * ------------------------------------------------------------------------------
 * ------------------------------------------------------------------------------*/
#pragma warning(disable:4786)

#include <spdlog/spdlog.h> 
 #include "gridedgehash.h"
#include "celltopology.h"
#include <assert.h>
#include <utility>
#include <stdlib.h>

/* -----------------------------------------------------------------------------
 * 实现边表类GridEdgeHash (哈希表)，用于记录顶点之间的邻接关系
 * ----------------------------------------------------------------------------*/
GridEdgeHash::~GridEdgeHash()
{
	release();
}

/* ---------------------------------------------------------------------------
 * 根据输入网格信息初始化哈希表
 * ---------------------------------------------------------------------------*/
int GridEdgeHash::initialize(
	std::vector<int> vecETopus,	/* 单元拓扑 */
	std::vector<int> vecENodes	/* 单元索引 */
	)
{
	int i = 0, j = 0, nodeS = 0;
	const unsigned int *localNodeIdx = NULL;
	int globNodeIdx[2], i1, i2;
	int etopu, nodeCount, edgeCount, addCount;
	int nRet = 0;

	for (i = 0; i < vecETopus.size(); i++)
	{
		etopu = vecETopus[i];
		nodeCount = EEMAS::TopologyInfo::corners((EEMAS::EntityTopology)etopu);
		if (EEMAS::TopologyInfo::dimension((EEMAS::EntityTopology)etopu) >= 2)
		{
			edgeCount = EEMAS::TopologyInfo::edges((EEMAS::EntityTopology)etopu);
			for (j = 0; j < edgeCount; j++)
			{
				EEMAS::TopologyInfo::edge_vertices((EEMAS::EntityTopology)etopu, j, &localNodeIdx);
				assert(localNodeIdx != NULL);

				globNodeIdx[0] = vecENodes[localNodeIdx[0] + nodeS];
				globNodeIdx[1] = vecENodes[localNodeIdx[1] + nodeS];

				if (m_eHashType == MIN_KEY || m_eHashType == MAX_KEY)
				{
					i1 = globNodeIdx[0] > globNodeIdx[1] ? globNodeIdx[1] : globNodeIdx[0];
					i2 = globNodeIdx[0] > globNodeIdx[1] ? globNodeIdx[0] : globNodeIdx[1];
					
					nRet = m_eHashType == MIN_KEY ? 
						add_edge(i1, i2, i, &addCount) : 
						add_edge(i2, i1, i, &addCount);
				}
				else
				{
					assert(m_eHashType == DOUB_KEY);
					
					nRet = add_edge(globNodeIdx[0], globNodeIdx[1], i, &addCount);
					if (nRet == 0 || nRet == 1)
					{
						nRet = add_edge(globNodeIdx[1], globNodeIdx[0], i, &addCount);
					}
				}

				if (nRet != 0 && nRet != 1)
					goto FAIL;
			}
		}
		else if (etopu == EEMAS::LINE)
		{
			assert(nodeCount == 2);

			globNodeIdx[0] = vecENodes[nodeS];
			globNodeIdx[1] = vecENodes[nodeS + 1];

			if (m_eHashType == MIN_KEY || m_eHashType == MAX_KEY)
			{
				i1 = globNodeIdx[0] > globNodeIdx[1] ? globNodeIdx[1] : globNodeIdx[0];
				i2 = globNodeIdx[0] > globNodeIdx[1] ? globNodeIdx[0] : globNodeIdx[1];
				
				nRet = m_eHashType == MIN_KEY ? 
					add_edge(i1, i2, i, &addCount) : 
					add_edge(i2, i1, i, &addCount);
			}
			else
			{
				assert(m_eHashType == DOUB_KEY);
				
				nRet = add_edge(globNodeIdx[0], globNodeIdx[1], i, &addCount);
				if (nRet == 0 || nRet == 1)
				{
					nRet = add_edge(globNodeIdx[1], globNodeIdx[0], i, &addCount);
				}
			}

			if (nRet != 0 && nRet != 1)
				goto FAIL;
		}
		else
		{
			assert(nodeCount == 1);
		}

		nodeS += nodeCount;
	}

	goto END;

FAIL:
	release();

END:
	return nRet;
}

/* ---------------------------------------------------------------------------
 * 根据哈希表查询节点的所有邻接节点
 * ---------------------------------------------------------------------------*/
int GridEdgeHash::find_neighbors(int nodeKey, std::vector<int>& vecNeigNodes)
{
	GridEdgeHash::iterator beg = lower_bound(nodeKey);
	GridEdgeHash::iterator end = upper_bound(nodeKey);
	GEdgeHashNode *hNode = NULL;

	while (beg != end)
	{
		hNode = beg->second;
		if (hNode)
		{
			vecNeigNodes.push_back(hNode->node2);
		}
		++beg;
	}
	
	/* 如果类型是MIN_KEY/MAX_KEY，上述算法将遗漏一些边 */

	return m_eHashType == DOUB_KEY ? 0 : 1;	/* EEMAS_S_OK / EEMAS_S_FALSE */
}

/* ---------------------------------------------------------------------------
 * 根据哈希表查询包含节点的所有单元
 * ---------------------------------------------------------------------------*/
int GridEdgeHash::find_node_elems(int nodeKey, std::vector<int>& vecNeigElems)
{
	GridEdgeHash::iterator beg, end;
	GEdgeHashNode *hNode = NULL;
	int i = 0;

	if (m_bElemRecorded)
	{
		beg = lower_bound(nodeKey);
		end = upper_bound(nodeKey);

		while (beg != end)
		{
			hNode = beg->second;
			if (hNode)
			{
				for (i = 0; i < hNode->vecElems.size(); i++)
				{
					vecNeigElems.push_back(hNode->vecElems[i]);
				}
			}
			++beg;
		}

		/* 如果类型是MIN_KEY/MAX_KEY，上述算法将遗漏一些边 */
		return m_eHashType == DOUB_KEY ? 0 : 1;	/* EEMAS_S_OK / EEMAS_S_FALSE */
	}

	return 7;	/* EEMAS_INVALID_CALL */
}

/* ---------------------------------------------------------------------------
 * 根据哈希表查询包含边的所有单元
 * ---------------------------------------------------------------------------*/
int GridEdgeHash::find_edge_elems(int node1, int node2, std::vector<int>& vecNeigElems)
{
	GridEdgeHash::iterator beg, end;
	GEdgeHashNode *hNode = NULL;
	int i = 0, i1, i2, nodeKey, nodeOth;

	if (m_bElemRecorded)
	{
		i1 = node1 > node2 ? node2 : node1;
		i2 = node1 > node2 ? node1 : node2;
		nodeKey = m_eHashType == MAX_KEY ? i2 : i1;
		nodeOth = m_eHashType == MAX_KEY ? i1 : i2;

		beg = lower_bound(nodeKey);
		end = upper_bound(nodeKey);

		while (beg != end)
		{
			hNode = beg->second;
			if (hNode && hNode->node2 == nodeOth)
			{
				for (i = 0; i < hNode->vecElems.size(); i++)
				{
					vecNeigElems.push_back(hNode->vecElems[i]);
				}
			}
			++beg;
		}

		return 0;
	}

	return 7;	/* EEMAS_INVALID_CALL */
}


/* ---------------------------------------------------------------------------
 * 受保护函数
 * ---------------------------------------------------------------------------*/
void GridEdgeHash::release()	/* 释放资源 */
{
	int i = 0;
	GEdgeHashNode *nodeIter = NULL, *nodeIterAft = NULL;
	GridEdgeHash::iterator it = begin();

	while (it != end())
	{
		nodeIter = it->second;
		if (nodeIter)
		{
			delete nodeIter;
		}
		++it;
	}
	clear();
}

int GridEdgeHash::add_edge(
	int nodeKey,		/* 关键点 */
	int node2,			/* 第二节点 */
	int elem,			/* 单元 */
	int *addCount		/* 增加的类型：
							0. 第一次增加该关键字，
							1. 第一次增加该边
						   >1. 再次增加该边 */
	)
{
	GridEdgeHash::iterator beg = lower_bound(nodeKey);
	GridEdgeHash::iterator end = upper_bound(nodeKey);
	GEdgeHashNode *hNode = NULL;

	*addCount = 0;
	if (beg != end)
	{
		*addCount = 1;
		while (beg != end)
		{
			hNode = beg->second;
			if (hNode && hNode->node2 == node2)
			{/* 再次增加该边 */
				if (m_bElemRecorded)
				{
					hNode->vecElems.push_back(elem);
					*addCount = hNode->vecElems.size();
// #ifdef _DEBUG
// 					if (*addCount > 2)
// 					{
// 						printf("More than two elements share the edge (%d, %d).\n", nodeKey, node2);
// 					}
// #endif
				}	
				else
					*addCount = 2;

				break;
			}
			++beg;
		}
	}
	
	if (*addCount <= 1)
	{/* 增加一条边 */
		if (!(hNode = new GEdgeHashNode))
			return -1;	/* EEMAS_OUT_OF_MEMORY */

		hNode->node2 = node2;
		insert(std::make_pair(nodeKey, hNode));
		if (m_bElemRecorded)
			hNode->vecElems.push_back(elem);
	}

	return 0;
}

