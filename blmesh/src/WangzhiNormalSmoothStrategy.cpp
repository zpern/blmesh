#include <spdlog/spdlog.h> 
 #include "../include/WangzhiNormalSmoothStrategy.h"

namespace WZSMOOTHING {
	WangZhiNormalSmoothStategy::WangZhiNormalSmoothStategy(BLNode **node, MBLNode* pnodes, int num_front) :NormalSmoothStrategy(node, pnodes, num_front)
	{

		classification();

	}

	void WangZhiNormalSmoothStategy::classification()
	{
	//	types_.resize(nFrtNods, 4);
		for (int i = 0; i < nFrtNods; i++) {
		//	int type = 0;
			auto front=node_array[i]->GetNeigFronts();
			for (int j = 0; j < front.size(); j++) {
				for (int k = j+1; k < front.size(); k++) {
					if (front[j]->GetNormal().normalized()*front[k]->GetNormal().normalized() < cos(PI / 3)) {
					//	type = 1;
						types_.insert(node_array[i]->GetNodIdx());


						
					}
				}
			}
		
		}
		

		for (int i = 0; i < nFrtNods; i++) {
			//	int type = 0;
		
			

		}


	}

	void WangZhiNormalSmoothStategy::SmoothNormalOnce(BLNode * blNod, int iq)
	{
		blNod->if_need_smooth = false;
		if (types_.find(blNod->GetNodIdx()) != types_.end())
			return;
		
		constexpr double thredhold = 25;
		double min_front_size = std::numeric_limits<double>::max();


		int nNeigFrts, nNeigNods;
		int id = blNod->GetNodIdx();
		BLVector centor(pNodes[id].coord);;

		const std::vector<BLNode*>& blNodes = blNod->GetNeigNods();
		nNeigNods = blNodes.size();
		if (nNeigNods == 0)
			return;


		BLVector sum(0, 0, 0);
		double height = 0;
		int count = 0;
		const std::vector <BLFront*>& blfronts = blNod->GetNeigFronts();
		// 周围点数应该等于周围前沿数
		for (int i = 0; i < nNeigNods; i++) {
			min_front_size = min(min_front_size, blfronts[i]->GetMinFrontSize());
		}
		//frontsize /= nNeigFrts;
		BLVector desentNode(pNodes[blNod->GetDecentID()].coord);
		double height1 = (desentNode - centor).magnitude();

		/*高度宽度比**/
		double ratio = height1 / min_front_size;
		static double lengths[200];
		BLVector norm0(0,0,0);
		//std::vector<double> lengths;
		double length_sum = 0;
		int type = 4;
		for (int i = 0; i < nNeigNods; i++) {
			int idx = blNodes[i]->GetNodIdx();
			if (types_.find(blNodes[i]->GetNodIdx()) != types_.end()) {
				type = 3;
			}
			BLVector coord1(pNodes[idx].coord[0], pNodes[idx].coord[1], pNodes[idx].coord[2]);
			double length = (coord1 - centor).magnitude();
			lengths[i] = length;
			norm0 =norm0+ (10 / length)* blNodes[i]->GetNormal();

			length_sum += 1/length;
		}
		norm0.normalize();
		length_sum /= nNeigNods;

		if (type == 3) {
			
		}
		else {
		
		}





		int t = 0;

		if(blNod->GetLowerNode())
		while (blNod->GetLowerNode()->GetNormal()*(norm0)  < cos(PI/2/9)) {
			if (t > 10) {
				norm0 = blNod->GetLowerNode()->GetNormal();
				break;
			}
			norm0 = (norm0 + 0.7*blNod->GetLowerNode()->GetNormal());
			norm0.normalize();
			t++;
		}


		t = 0;
		while (blNod->GetBeitaVisu(norm0) * 180 / PI < 10) {
			if (t > 10) {
				return;
			}
			norm0 = (norm0 + 0.7*blNod->GetNormal());
			norm0.normalize();
			t++;
		}

		if (blNod->GetBSys())
		{
			norm0[blNod->GetSymAxis()] = 0;
		}

		if (norm0*  blNod->GetNormal() < 0.9985) {
			for (int i = 0; i < nNeigNods; i++) {
				blNodes[i]->if_need_smooth = true;

			}
			blNod->if_need_smooth = true;

		}
	//	myvec[iq] = norm0;
			blNod->SetNormal(norm0);	
	}

	WangZhiNormalSmoothStategy::~WangZhiNormalSmoothStategy()
	{

	}
}
