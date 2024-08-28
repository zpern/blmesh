#include <spdlog/spdlog.h> 
 #include "skirtpolylist.h"

#if 1
SkirtPolyNode* SkirtPolyNode::getSkirtPolyNode(int edgeNode[2], int sideFace[2], float minQual, SkirtPolyNode* next)
{
	SkirtPolyNode* polynode = new SkirtPolyNode(edgeNode, sideFace, minQual);
	polynode->next = next;

	return polynode;
}

// setter / getter functions
void SkirtPolyNode::setQual(float minQual)
{
	this->minQual = minQual;
}

void SkirtPolyNode::getQual(float *minQual)
{
	*minQual = this->minQual;
}


void SkirtPolyNode::setEdge(int edgeNode[2])
{
	this->edgeNode[0] = edgeNode[0];
	this->edgeNode[1] = edgeNode[1];
}

void SkirtPolyNode::getEdge(int edgeNode[2])
{
	edgeNode[0] = this->edgeNode[0];
	edgeNode[1] = this->edgeNode[1];
}


void SkirtPolyNode::setFace(int sideFace[2])
{
	this->sideFace[0] = sideFace[0];
	this->sideFace[1] = sideFace[1];
}

void SkirtPolyNode::getFace(int sideFace[2])
{
	sideFace[0] = this->sideFace[0];
	sideFace[1] = this->sideFace[1];
}

void SkirtPolyNode::insertAfter(SkirtPolyNode* polynode)
{
	polynode->next = this->next;
	this->next = polynode;
}

SkirtPolyNode* SkirtPolyNode::removeAfter()
{
	SkirtPolyNode* plyNode = this->next;
	if(this->next == nullptr)
		return nullptr;
	this->next = plyNode->next;
	return plyNode;
}

/* print info */
void SkirtPolyNode::printinfo()
{
#ifdef _FILE_OUT
#else
	spdlog::info("\t\t********** Shell Information **********\n");
	spdlog::info("\t\tEdge Index ({} {})\n", this->edgeNode[0], this->edgeNode[1]);
	//spdlog::info("\t\tElement Index: {}\n", this->eleIdx);
	spdlog::info("\t\tQuality: min(%10.5lf)\n", this->minQual);
	//spdlog::info("\t\tSkirt Polygon (");
	//for (int i=0; i<this->nSkirtNodes; i++)
	//	spdlog::info(" {}", this->skirtPoly[i]);
	spdlog::info(" )\n\n");
#endif
}

// HashPolyNode

void HashPolyNode::insertAfter(HashPolyNode* node)
{
	node->next = this->next;
	this->next = node;
}

HashPolyNode* HashPolyNode::removeAfter()
{
	HashPolyNode* hashNode = this->next;
	if(this->next == nullptr)
		return nullptr;
	this->next = hashNode->next;
	return hashNode;
}

/************************************************************************/
/*  按质量大小从小到大依次插入到链表中                                  */
/************************************************************************/
SkirtPolyNode * HashPolyNode::addSkirtPoly(int edgeNode[2], int sideFace[2], float minQual)
{

	int key1, key2;
	SkirtPolyNode *cur, *pre;
	key1 = edgeNode[0];
	key2 = edgeNode[1];

	if(this->ptIdx != key1)
		return nullptr;

	SkirtPolyNode *notetmp, *newNode;
	notetmp = new SkirtPolyNode;
	newNode = notetmp->getSkirtPolyNode(edgeNode, sideFace, minQual);
	
	pre = cur = this->pSkirtPolyList;
	if(cur == nullptr)
		this->pSkirtPolyList = newNode;

	while(cur != nullptr)
	{
		if(cur->minQual > minQual)
		{
			if( pre != cur)
				pre->insertAfter(newNode);
			else	//表头节点
			{
				this->pSkirtPolyList = newNode;
				newNode->next = cur;
			}

			break;
		}
		else
		{
			pre = cur;
			cur = cur->next;
		}
	}
	if(cur == nullptr && pre != nullptr)
		pre->insertAfter(newNode);

	return newNode;
}

/************************************************************************/
/*                                                                      */
/************************************************************************/
void HashPolyNode::deleteSkirtPoly(int edgeNode[2], int sideFace[2])
{
	int key1, key2;
	SkirtPolyNode *cur, *pre, *temp;
	key1 = edgeNode[0];
	key2 = edgeNode[1];

	if(this->ptIdx != key1)
		return ;

	pre = cur = this->pSkirtPolyList;
	while(cur != nullptr)
	{
		if(cur->edgeNode[1] == key2 && cur->sideFace[0] == sideFace[0] && cur->sideFace[1] == sideFace[1])
		{
			if( pre != cur)
			{
				temp = pre->removeAfter();
				if(temp != nullptr)
				{
					//delete temp;
					//temp = nullptr;
					temp->setDeleted(true);
				}
			}
			else	//表头
			{
				this->pSkirtPolyList = cur->next;
				//delete cur;
				//cur = nullptr;
				cur->setDeleted(true);
			}

			this->nPolyCnt--;

			break;
		}
		else
		{
			pre = cur;
			cur = cur->next;
		}
	}
}

/************************************************************************/
/*                                                                      */
/************************************************************************/
bool HashPolyNode::updateSkirtPoly(int edgeNode[2], int sideFace[2], float minQual)
{
	int key1, key2;
	SkirtPolyNode *cur, *pre, *temp;
	key1 = edgeNode[0];
	key2 = edgeNode[1];

	if(this->ptIdx != key1)
		return false;

	pre = cur = this->pSkirtPolyList;
	if(cur == nullptr)
	{
		spdlog::info("update failed: no such skirtpoly.\n");
		return false;
	}

	while(cur != nullptr)
	{
		if(cur->edgeNode[1] == key2 && cur->sideFace[0] == sideFace[0] && cur->sideFace[1] == sideFace[1])
		{
			cur->setEdge(edgeNode);
			cur->setQual(minQual);
			resortSkirtPoly(cur);

			break;
		}
		else
		{
			pre = cur;
			cur = cur->next;
		}
	}

	return true;
}

void HashPolyNode::resortSkirtPoly(SkirtPolyNode* polynode)
{
	SkirtPolyNode *cur, *pre, *temp;

	pre = cur = this->pSkirtPolyList;
	while (cur != nullptr)
	{
		if(cur == polynode)
			break;
		pre = cur;
		cur = cur->next;
	}
	if(cur == nullptr)
		return;
	if(pre != cur)
		temp = pre;
	else	//polynode为表头节点
		temp = nullptr;

	pre = cur = this->pSkirtPolyList;
	while(cur != nullptr)
	{
		if(cur->minQual > polynode->minQual && cur != polynode)
		{
			if(temp != nullptr)
			{
				temp->removeAfter();
			}
			else	//polynode为表头节点
			{
				this->pSkirtPolyList = polynode->next;
			}
			
			if(pre != polynode)
				pre->insertAfter(polynode);
			else
				this->pSkirtPolyList->insertAfter(polynode);
			break;
		}
		else
		{
			pre = cur;
			cur = cur->next;
		}
	}

	if(cur == nullptr && pre != nullptr && pre != polynode)
	{
		if(temp != nullptr)
			temp->removeAfter();
		else	//polynode为表头节点
			this->pSkirtPolyList = polynode->next;
		
			pre->insertAfter(polynode);
	}
}

/************************************************************************/
/*                                                                      */
/************************************************************************/
bool HashPolyNode::isPolyExist(int edgeNode[2], int sideFace[2], SkirtPolyNode ** pNode)
{
	int key1, key2;
	key1 = edgeNode[0];
	key2 = edgeNode[1];

	if(this->ptIdx != key1)
		return false;
	*pNode = nullptr;
	SkirtPolyNode* cur = pSkirtPolyList;
	while(cur != nullptr)
	{
		if(cur->edgeNode[1] == key2 && cur->sideFace[0] == sideFace[0] && cur->sideFace[1] == sideFace[1])
		{
			*pNode = cur;
			return true;
		}
		else
			cur = cur->next;
	}

	return false;
}

/************************************************************************/
/*                                                                      */
/************************************************************************/
SkirtPolyNode* HashPolyNode::getSkirtPoly(int edgeNode[2], int sideFace[2])
{
	int key1, key2;
	key1 = edgeNode[0];
	key2 = edgeNode[1];
	
	if(this->ptIdx != key1)
		return nullptr;
	
	SkirtPolyNode* cur = pSkirtPolyList;
	while(cur != nullptr)
	{
		if(cur->edgeNode[1] == key2 && cur->sideFace[0] == sideFace[0] && cur->sideFace[1] == sideFace[1])
			return cur;
		else
			cur = cur->next;
	}
	
	return nullptr;
}

/* print info */
void HashPolyNode::printinfo()
{
#ifdef _FILE_OUT
#else
	spdlog::info("\t********** Hash Node Information **********\n");
	spdlog::info("\tHash Key: {}\n", this->ptIdx);
	spdlog::info("\tShell Number: {}\n", this->nPolyCnt);
#ifdef _OUTPUT_LEVEL_III
	SkirtPolyNode* cur = pSkirtPolyList;
	while(cur != nullptr)
	{
		cur->printinfo();
		cur = cur->next;
	}
#endif
	spdlog::info("\n");
#endif
}

//HashPolyList
/************************************************************************/
/*                                                                      */
/************************************************************************/
SkirtPolyNode * HashPolyList::addSkirtPoly(int edgeNode[2], int sideFace[2], float minQual)
{
	int key, p;
	if(edgeNode[0] > edgeNode[1])
	{
		p = edgeNode[0]; edgeNode[0] = edgeNode[1]; edgeNode[1] = p;
	}
	if(sideFace[0] > sideFace[1])
	{
		p = sideFace[0]; sideFace[0] = sideFace[1]; sideFace[1] = p;
	}
	key = edgeNode[0];
	SkirtPolyNode * newNode = nullptr;
	HashPolyNode *cur, *pre, *temp;

	pre = cur = this->first;
	while (cur != nullptr)
	{
		if (cur->ptIdx != key)
		{
			pre = cur;
			cur = cur->next;
		}
		else
		{
			//cur->addSkirtPoly(edgeNode,eleIdx,skirtPolyNode, npt, minQual, maxQual, aveQual);
			break;
		}
	}

	if(cur == nullptr)		//链表中不存在该节点
	{
		temp = new HashPolyNode(key);
		if(this->first == nullptr)
			this->first = temp;
		else
		{
			pre->insertAfter(temp);
		}
		cur = temp;
	}

	newNode = cur->addSkirtPoly(edgeNode, sideFace, minQual);
	if(newNode)
		cur->nPolyCnt++;

	return newNode;
}

void HashPolyList::deleteSkirtPoly(int edgeNode[2], int sideFace[2])
{
	int key, p;

	if(edgeNode[0] > edgeNode[1])
	{
		p = edgeNode[0]; edgeNode[0] = edgeNode[1]; edgeNode[1] = p;
	}
	if (sideFace[0] > sideFace[1])
	{
		p = sideFace[0]; sideFace[0] = sideFace[1]; sideFace[1] = p;
	}

	key = edgeNode[0];

	HashPolyNode* cur = this->first;
	while (cur != nullptr)
	{
		if (cur->ptIdx != key)
			cur = cur->next;
		else
		{
			cur->deleteSkirtPoly(edgeNode, sideFace);
			break;
		}
	}
}

void HashPolyList::updateSkirtPoly(int edgeNode[2], int sideFace[2], float minQual)
{
	int key, p;

	if(edgeNode[0] > edgeNode[1])
	{
		p = edgeNode[0]; edgeNode[0] = edgeNode[1]; edgeNode[1] = p;
	}
	if (sideFace[0] > sideFace[1])
	{
		p = sideFace[0]; sideFace[0] = sideFace[1]; sideFace[1] = p;
	}

	key = edgeNode[0];
	HashPolyNode* cur = this->first;
	while (cur != nullptr)
	{
		if (cur->ptIdx != key)
			cur = cur->next;
		else
		{
			cur->updateSkirtPoly(edgeNode, sideFace, minQual);
			break;
		}
	}
}

bool HashPolyList::isPolyExist(int edgeNode[2], int sideFace[2], SkirtPolyNode ** pNode)
{
	int key, p;

	if(edgeNode[0] > edgeNode[1])
	{
		p = edgeNode[0]; edgeNode[0] = edgeNode[1]; edgeNode[1] = p;
	}
	if (sideFace[0] > sideFace[1])
	{
		p = sideFace[0]; sideFace[0] = sideFace[1]; sideFace[1] = p;
	}

	key = edgeNode[0];
	HashPolyNode* cur = this->first;
	while (cur != nullptr)
	{
		if (cur->ptIdx != key)
			cur = cur->next;
		else
		{
			return cur->isPolyExist(edgeNode, sideFace, pNode);
		}
	}

	return false;
}

void HashPolyList::getSkirtPoly(int edgeNode[2], int sideFace[2], float *minQual)
{

}

SkirtPolyNode* HashPolyList::getSkirtPoly(int edgeNode[2], int sideFace[2])
{
	int key, p;

	if(edgeNode[0] > edgeNode[1])
	{
		p = edgeNode[0]; edgeNode[0] = edgeNode[1]; edgeNode[1] = p;
	}
	if (sideFace[0] > sideFace[1])
	{
		p = sideFace[0]; sideFace[0] = sideFace[1]; sideFace[1] = p;
	}
	key = edgeNode[0];
	HashPolyNode* cur = this->first;
	while (cur != nullptr)
	{
		if (cur->ptIdx != key)
			cur = cur->next;
		else
		{
			return cur->getSkirtPoly(edgeNode, sideFace);
		}
	}
	return nullptr;
}

void HashPolyList::removeZeroSkirtPoly()
{
	
}

/* print info */
void HashPolyList::printinfo()
{
#ifdef _FILE_OUT
#else
	spdlog::info("********** Hash List Information **********\n");
	HashPolyNode* cur = first;
	while(cur != nullptr)
	{
		cur->printinfo();
		cur = cur->next;
	}
	spdlog::info("\n");
#endif
}

#else

SkirtPolyNode* SkirtPolyNode::getSkirtPolyNode(int edgeNode[2], float minQual, SkirtPolyNode* next)
{
	SkirtPolyNode* polynode = new SkirtPolyNode(edgeNode, minQual);
	polynode->next = next;

	return polynode;
}

// setter / getter functions
void SkirtPolyNode::setQual(float minQual)
{
	this->minQual = minQual;
}

void SkirtPolyNode::getQual(float *minQual)
{
	*minQual = this->minQual;
}


void SkirtPolyNode::setEdge(int edgeNode[2])
{
	this->edgeNode[0] = edgeNode[0];
	this->edgeNode[1] = edgeNode[1];
}

void SkirtPolyNode::getEdge(int edgeNode[2])
{
	edgeNode[0] = this->edgeNode[0];
	edgeNode[1] = this->edgeNode[1];
}

void SkirtPolyNode::insertAfter(SkirtPolyNode* polynode)
{
	polynode->next = this->next;
	this->next = polynode;
}

SkirtPolyNode* SkirtPolyNode::removeAfter()
{
	SkirtPolyNode* plyNode = this->next;
	if(this->next == nullptr)
		return nullptr;
	this->next = plyNode->next;
	return plyNode;
}

/* print info */
void SkirtPolyNode::printinfo()
{
#ifdef _FILE_OUT
#else
	spdlog::info("\t\t********** Shell Information **********\n");
	spdlog::info("\t\tEdge Index ({} {})\n", this->edgeNode[0], this->edgeNode[1]);
	//spdlog::info("\t\tElement Index: {}\n", this->eleIdx);
	spdlog::info("\t\tQuality: min(%10.5lf)\n", this->minQual);
	//spdlog::info("\t\tSkirt Polygon (");
	//for (int i=0; i<this->nSkirtNodes; i++)
	//	spdlog::info(" {}", this->skirtPoly[i]);
	spdlog::info(" )\n\n");
#endif
}

// HashPolyNode

void HashPolyNode::insertAfter(HashPolyNode* node)
{
	node->next = this->next;
	this->next = node;
}

HashPolyNode* HashPolyNode::removeAfter()
{
	HashPolyNode* hashNode = this->next;
	if(this->next == nullptr)
		return nullptr;
	this->next = hashNode->next;
	return hashNode;
}

/************************************************************************/
/*  按质量大小从小到大依次插入到链表中                                  */
/************************************************************************/
SkirtPolyNode * HashPolyNode::addSkirtPoly(int edgeNode[2], float minQual)
{

	int key1, key2;
	SkirtPolyNode *cur, *pre;
	key1 = edgeNode[0];
	key2 = edgeNode[1];

	if(this->ptIdx != key1)
		return nullptr;

	SkirtPolyNode *notetmp, *newNode;
	notetmp = new SkirtPolyNode;
	newNode = notetmp->getSkirtPolyNode(edgeNode, minQual);

	pre = cur = this->pSkirtPolyList;
	if(cur == nullptr)
		this->pSkirtPolyList = newNode;

	while(cur != nullptr)
	{
		if(cur->minQual > minQual)
		{
			if( pre != cur)
				pre->insertAfter(newNode);
			else	//表头节点
			{
				this->pSkirtPolyList = newNode;
				newNode->next = cur;
			}

			break;
		}
		else
		{
			pre = cur;
			cur = cur->next;
		}
	}
	if(cur == nullptr && pre != nullptr)
		pre->insertAfter(newNode);

	return newNode;
}

/************************************************************************/
/*                                                                      */
/************************************************************************/
void HashPolyNode::deleteSkirtPoly(int edgeNode[2])
{
	int key1, key2;
	SkirtPolyNode *cur, *pre, *temp;
	key1 = edgeNode[0];
	key2 = edgeNode[1];

	if(this->ptIdx != key1)
		return ;

	pre = cur = this->pSkirtPolyList;
	while(cur != nullptr)
	{
		if(cur->edgeNode[1] == key2)
		{
			if( pre != cur)
			{
				temp = pre->removeAfter();
				if(temp != nullptr)
				{
					//delete temp;
					//temp = nullptr;
					temp->setDeleted(true);
				}
			}
			else	//表头
			{
				this->pSkirtPolyList = cur->next;
				//delete cur;
				//cur = nullptr;
				cur->setDeleted(true);
			}

			this->nPolyCnt--;

			break;
		}
		else
		{
			pre = cur;
			cur = cur->next;
		}
	}
}

/************************************************************************/
/*                                                                      */
/************************************************************************/
bool HashPolyNode::updateSkirtPoly(int edgeNode[2], float minQual)
{
	int key1, key2;
	SkirtPolyNode *cur, *pre, *temp;
	key1 = edgeNode[0];
	key2 = edgeNode[1];

	if(this->ptIdx != key1)
		return false;

	pre = cur = this->pSkirtPolyList;
	if(cur == nullptr)
	{
		spdlog::info("update failed: no such skirtpoly.\n");
		return false;
	}

	while(cur != nullptr)
	{
		if(cur->edgeNode[1] == key2)
		{
			cur->setEdge(edgeNode);
			cur->setQual(minQual);
			resortSkirtPoly(cur);

			break;
		}
		else
		{
			pre = cur;
			cur = cur->next;
		}
	}

	return true;
}

void HashPolyNode::resortSkirtPoly(SkirtPolyNode* polynode)
{
	SkirtPolyNode *cur, *pre, *temp;

	pre = cur = this->pSkirtPolyList;
	while (cur != nullptr)
	{
		if(cur == polynode)
			break;
		pre = cur;
		cur = cur->next;
	}
	if(cur == nullptr)
		return;
	if(pre != cur)
		temp = pre;
	else	//polynode为表头节点
		temp = nullptr;

	pre = cur = this->pSkirtPolyList;
	while(cur != nullptr)
	{
		if(cur->minQual > polynode->minQual && cur != polynode)
		{
			if(temp != nullptr)
			{
				temp->removeAfter();
			}
			else	//polynode为表头节点
			{
				this->pSkirtPolyList = polynode->next;
			}

			if(pre != polynode)
				pre->insertAfter(polynode);
			else
				this->pSkirtPolyList->insertAfter(polynode);
			break;
		}
		else
		{
			pre = cur;
			cur = cur->next;
		}
	}

	if(cur == nullptr && pre != nullptr && pre != polynode)
	{
		if(temp != nullptr)
			temp->removeAfter();
		else	//polynode为表头节点
			this->pSkirtPolyList = polynode->next;

		pre->insertAfter(polynode);
	}
}

/************************************************************************/
/*                                                                      */
/************************************************************************/
bool HashPolyNode::isPolyExist(int edgeNode[2])
{
	int key1, key2;
	key1 = edgeNode[0];
	key2 = edgeNode[1];

	if(this->ptIdx != key1)
		return false;

	SkirtPolyNode* cur = pSkirtPolyList;
	while(cur != nullptr)
	{
		if(cur->edgeNode[1] == key2)
			return true;
		else
			cur = cur->next;
	}

	return false;
}

/************************************************************************/
/*                                                                      */
/************************************************************************/
SkirtPolyNode* HashPolyNode::getSkirtPoly(int edgeNode[2])
{
	int key1, key2;
	key1 = edgeNode[0];
	key2 = edgeNode[1];

	if(this->ptIdx != key1)
		return nullptr;

	SkirtPolyNode* cur = pSkirtPolyList;
	while(cur != nullptr)
	{
		if(cur->edgeNode[1] == key2)
			return cur;
		else
			cur = cur->next;
	}

	return nullptr;
}

/* print info */
void HashPolyNode::printinfo()
{
#ifdef _FILE_OUT
#else
	spdlog::info("\t********** Hash Node Information **********\n");
	spdlog::info("\tHash Key: {}\n", this->ptIdx);
	spdlog::info("\tShell Number: {}\n", this->nPolyCnt);
#ifdef _OUTPUT_LEVEL_III
	SkirtPolyNode* cur = pSkirtPolyList;
	while(cur != nullptr)
	{
		cur->printinfo();
		cur = cur->next;
	}
#endif
	spdlog::info("\n");
#endif
}

//HashPolyList
/************************************************************************/
/*                                                                      */
/************************************************************************/
SkirtPolyNode * HashPolyList::addSkirtPoly(int edgeNode[2], float minQual)
{
	int key, p;
	if(edgeNode[0] > edgeNode[1])
	{
		p = edgeNode[0]; edgeNode[0] = edgeNode[1]; edgeNode[1] = p;
	}
	key = edgeNode[0];
	SkirtPolyNode * newNode = nullptr;
	HashPolyNode *cur, *pre, *temp;

	pre = cur = this->first;
	while (cur != nullptr)
	{
		if (cur->ptIdx != key)
		{
			pre = cur;
			cur = cur->next;
		}
		else
		{
			//cur->addSkirtPoly(edgeNode,eleIdx,skirtPolyNode, npt, minQual, maxQual, aveQual);
			break;
		}
	}

	if(cur == nullptr)		//链表中不存在该节点
	{
		temp = new HashPolyNode(key);
		if(this->first == nullptr)
			this->first = temp;
		else
		{
			pre->insertAfter(temp);
		}
		cur = temp;
	}

	newNode = cur->addSkirtPoly(edgeNode, minQual);
	if(newNode)
		cur->nPolyCnt++;

	return newNode;
}

void HashPolyList::deleteSkirtPoly(int edgeNode[2])
{
	int key, p;

	if(edgeNode[0] > edgeNode[1])
	{
		p = edgeNode[0]; edgeNode[0] = edgeNode[1]; edgeNode[1] = p;
	}

	key = edgeNode[0];

	HashPolyNode* cur = this->first;
	while (cur != nullptr)
	{
		if (cur->ptIdx != key)
			cur = cur->next;
		else
		{
			cur->deleteSkirtPoly(edgeNode);
			break;
		}
	}
}

void HashPolyList::updateSkirtPoly(int edgeNode[2], float minQual)
{
	int key, p;

	if(edgeNode[0] > edgeNode[1])
	{
		p = edgeNode[0]; edgeNode[0] = edgeNode[1]; edgeNode[1] = p;
	}

	key = edgeNode[0];
	HashPolyNode* cur = this->first;
	while (cur != nullptr)
	{
		if (cur->ptIdx != key)
			cur = cur->next;
		else
		{
			cur->updateSkirtPoly(edgeNode, minQual);
			break;
		}
	}
}

bool HashPolyList::isPolyExist(int edgeNode[2])
{
	int key, p;

	if(edgeNode[0] > edgeNode[1])
	{
		p = edgeNode[0]; edgeNode[0] = edgeNode[1]; edgeNode[1] = p;
	}

	key = edgeNode[0];
	HashPolyNode* cur = this->first;
	while (cur != nullptr)
	{
		if (cur->ptIdx != key)
			cur = cur->next;
		else
		{
			return cur->isPolyExist(edgeNode);
		}
	}

	return false;
}

void HashPolyList::getSkirtPoly(int edgeNode[2], float *minQual)
{

}

SkirtPolyNode* HashPolyList::getSkirtPoly(int edgeNode[2])
{
	int key, p;

	if(edgeNode[0] > edgeNode[1])
	{
		p = edgeNode[0]; edgeNode[0] = edgeNode[1]; edgeNode[1] = p;
	}
	key = edgeNode[0];
	HashPolyNode* cur = this->first;
	while (cur != nullptr)
	{
		if (cur->ptIdx != key)
			cur = cur->next;
		else
		{
			return cur->getSkirtPoly(edgeNode);
		}
	}
	return nullptr;
}

void HashPolyList::removeZeroSkirtPoly()
{

}

/* print info */
void HashPolyList::printinfo()
{
#ifdef _FILE_OUT
#else
	spdlog::info("********** Hash List Information **********\n");
	HashPolyNode* cur = first;
	while(cur != nullptr)
	{
		cur->printinfo();
		cur = cur->next;
	}
	spdlog::info("\n");
#endif
}
#endif

//SkirtPolyHeap

/************************************************************************/
/* 含有n个节点的堆的高度为lgn，故插入一个节点的时间复杂度为O（lgn）     */
/************************************************************************/
int SkirtPolyHeap::insertSkirtPoly(SkirtPolyNode *pPolyNode)
{
	if(nHeapSize == nMaxHeapSize)
		allocHeap();

	m_pSkirtPolyHeap[nHeapSize] = pPolyNode;
	minHeapifyUp(nHeapSize);
	nHeapSize++;

	return 0;
}

SkirtPolyNode *SkirtPolyHeap::getMinSkirtPoly()
{
	SkirtPolyNode *pNode;
	if(isEmpty())
		return 0;
	
	pNode = m_pSkirtPolyHeap[0];
	//m_pSkirtPolyHeap[0] = m_pSkirtPolyHeap[nHeapSize-1];
	//nHeapSize--;

	//minHeapifyDown(0);
	
	return pNode;
}

void SkirtPolyHeap::rmMinSkirtPoly()
{
	if(isEmpty())
		return ;

	m_pSkirtPolyHeap[0] = m_pSkirtPolyHeap[nHeapSize-1];
	nHeapSize--;

	minHeapifyDown(0);
}

void SkirtPolyHeap::minHeapifyDown(int beg)
{
	int i, j;
	SkirtPolyNode *temp;

	i = beg;
	j = 2*i+1;	//左子女编号
	temp = m_pSkirtPolyHeap[i];
	
	while(j < nHeapSize)
	{
		if (j < nHeapSize-1 && m_pSkirtPolyHeap[j]->getMinQual() > m_pSkirtPolyHeap[j+1]->getMinQual())
			j++;
		if(temp->getMinQual() < m_pSkirtPolyHeap[j]->getMinQual())
			break;
		else
		{
			m_pSkirtPolyHeap[i] = m_pSkirtPolyHeap[j];
			i = j;
			j = 2*j+1;
		}
	}

	m_pSkirtPolyHeap[i] = temp;
}

void SkirtPolyHeap::minHeapifyUp(int beg)
{
	int i, j;
	SkirtPolyNode *temp;

	j = beg;
	i = (j-1)/2;	//父亲节点
	temp = m_pSkirtPolyHeap[j];
	while(j>0)
	{
		if(m_pSkirtPolyHeap[i]->getMinQual() < temp->getMinQual())
			break;
		else
		{
			m_pSkirtPolyHeap[j] = m_pSkirtPolyHeap[i];
			j = i;
			i = (i-1)/2;
		}
	}
	m_pSkirtPolyHeap[j] = temp;
}

void SkirtPolyHeap::printfinfo()
{
	SkirtPolyNode *temp;
	temp = getMinSkirtPoly();
	do 
	{
		int edge[2], eleidx;
		temp->getEdge(edge);
		spdlog::info("edge ({} {})\n", edge[0], edge[1]);
		spdlog::info("value: %lf\n", temp->getMinQual());

		temp = getMinSkirtPoly();
	} while (temp);
}