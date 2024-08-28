#pragma once
#include "NormalSmoothStrategy.h"
class IterationNormalSmoothStrategy:public NormalSmoothStrategy
{
public:
	IterationNormalSmoothStrategy(BLNode **node, Node* pnodes, int num_front);
	void virtual SmoothNormal() ;
	~IterationNormalSmoothStrategy();
};

