#pragma once
#ifndef SIMPLE_NORMAL_SMOOTH_STRATEGY_H_
#define SIMPLE_NORMAL_SMOOTH_STRATEGY_H_
#include "NormalSmoothStrategy.h"
class SimpleNormalSmoothStrategy:public NormalSmoothStrategy
{
public:
	SimpleNormalSmoothStrategy(BLNode **node, MBLNode* pnodes,int num_front);
	void virtual SmoothNormalOnce(BLNode *blNod, int id);
	
	~SimpleNormalSmoothStrategy();
};
#endif //!NORMAL_SMOOTH_STRATEGY_H_

