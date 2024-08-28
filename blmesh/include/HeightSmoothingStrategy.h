#pragma once
#ifndef HEIGHT_NORMAL_SMOOTH_STRATEGY_H_
#define HEIGHT_NORMAL_SMOOTH_STRATEGY_H_

#include <functional>

#include "NormalSmoothStrategy.h"
#include "sphericalvoronoigeometry.h"
#include "BLVector.h"
#include "BLMesh_define.h"



namespace HEIGHTSMOOTHING {










	class HeightSmoothStategy :public NormalSmoothStrategy {
	public:
		HeightSmoothStategy(BLNode **node, Node* pnodes, int num_front);
		void virtual SmoothNormalOnce(BLNode *blNod, int id);

		~HeightSmoothStategy();
	};
}
#endif
