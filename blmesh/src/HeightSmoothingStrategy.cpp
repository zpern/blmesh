#include <spdlog/spdlog.h> 
 #include "../include/HeightSmoothingStrategy.h"

HEIGHTSMOOTHING::HeightSmoothStategy::HeightSmoothStategy(BLNode ** node, MBLNode * pnodes, int num_front)
{
}

void HEIGHTSMOOTHING::HeightSmoothStategy::SmoothNormalOnce(BLNode * blNod, int id)
{


	blNod->if_need_smooth = false;
	constexpr double thredhold = 25;
	double min_front_size = std::numeric_limits<double>::max();


	int nNeigFrts, nNeigNods;
	BLVector centor(pNodes[id].coord);;

	const std::vector<BLNode*>& blNodes = blNod->GetNeigNods();
	nNeigNods = blNodes.size();
	if (nNeigNods == 0)
		return;

	std::vector<BLVector> neighbour_point_normals;
	std::vector<BLVector> neighbout_front_noamls;

	const std::vector <BLFront*>& blfronts = blNod->GetNeigFronts();
	// 周围点数应该等于周围前沿数
	for (int i = 0; i < nNeigNods; i++) {
		//	min_front_size = min(min_front_size, blfronts[i]->GetMinFrontSize());
		neighbour_point_normals.push_back(blNodes[i]->GetNormal());
		neighbout_front_noamls.push_back(blfronts[i]->GetNormal());
	}

	std::vector<BLVector> neighbourpoints;
	std::vector<double> neighbourheight;
	for (int i = 0; i < nNeigNods; i++) {
		neighbourpoints.push_back(blNodes[i]->GetCoord(pNodes));
		neighbourheight.push_back(blNodes[i]->GetHeightLength());
	}

	std::vector<double> ideal_height;
	for (int i = 0; i < nNeigNods; i++) {
		double H=((neighbourpoints[i]+ neighbourheight[i]* neighbour_point_normals[i]- centor)*(blNod->GetNormal()+ neighbour_point_normals[i])) \
			/ ((blNod->GetNormal() + neighbour_point_normals[i])*neighbour_point_normals[i]);
		ideal_height.push_back(H);
	}

	double H_max = *std::max_element(ideal_height.begin(),ideal_height.end());
	double H_min = *std::min_element(ideal_height.begin(), ideal_height.end());
	double H_init = blNod->GetHeightLength()/(1 + blNod->GetHightRatio());






}

HEIGHTSMOOTHING::HeightSmoothStategy::~HeightSmoothStategy()
{
}
