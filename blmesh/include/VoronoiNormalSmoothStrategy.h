#pragma once
#ifndef VORONOI_NORMAL_SMOOTH_STRATEGY_H_
#define VORONOI_NORMAL_SMOOTH_STRATEGY_H_

#include <functional>

#include "NormalSmoothStrategy.h"
#include "sphericalvoronoigeometry.h"
#include "BLVector.h"
#include "BLMesh_define.h"



namespace VORONOISMOOTHING {










	class VoronoiNormalSmoothStategy :public NormalSmoothStrategy {
	public:
		VoronoiNormalSmoothStategy(BLNode **node, Node* pnodes, int num_front);
		void virtual SmoothNormalOnce(BLNode *blNod, int id);

		~VoronoiNormalSmoothStategy();
	};
}
#endif
