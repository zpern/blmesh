#pragma once
#ifndef WZ_NORMAL_SMOOTH_STRATEGY_H_
#define WZ_NORMAL_SMOOTH_STRATEGY_H_

#include <functional>

#include "NormalSmoothStrategy.h"
#include "sphericalvoronoigeometry.h"
#include "BLVector.h"
#include "BLMesh_define.h"



namespace WZSMOOTHING {










	class WangZhiNormalSmoothStategy :public NormalSmoothStrategy {
	public:
		WangZhiNormalSmoothStategy(BLNode **node, MBLNode* pnodes, int num_front);
		void classification();
		void virtual SmoothNormalOnce(BLNode *blNod, int id);

		~WangZhiNormalSmoothStategy();
		set<int> types_;
	};
}
#endif
