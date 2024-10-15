#include <spdlog/spdlog.h> 
 #include "../include/VoronoiNormalSmoothStrategy.h"
namespace VORONOISMOOTHING {
	VoronoiNormalSmoothStategy::VoronoiNormalSmoothStategy(BLNode **node, MBLNode* pnodes, int num_front) :NormalSmoothStrategy(node, pnodes, num_front)
	{


	}

	void VoronoiNormalSmoothStategy::SmoothNormalOnce(BLNode * blNod,int iq)
	{
		//if (blNod->GetVirtualFlag())
	//	return;
		blNod->if_need_smooth = false;
		constexpr double thredhold = 25;
		double min_front_size = std::numeric_limits<double>::max();


		int nNeigFrts, nNeigNods;
		int id = blNod->GetNodIdx();
		BLVector centor(pNodes[id].coord);;

		const std::vector<BLNode*>& blNodes = blNod->GetNeigNods();
		nNeigNods = blNodes.size();
		if (nNeigNods == 0)
			return;

		std::vector<BLVector> neighbour_point_normals;
		std::vector<BLVector> neighbout_front_normals;
		std::vector<BLVector> neighbout_upper_normals;

		const std::vector <BLFront*>& blfronts = blNod->GetNeigFronts();
		// 周围点数应该等于周围前沿数
		for (int i = 0; i < nNeigNods; i++) {
			//	min_front_size = min(min_front_size, blfronts[i]->GetMinFrontSize());
			neighbour_point_normals.push_back(blNodes[i]->GetNormal());
			neighbout_front_normals.push_back(blfronts[i]->GetNormal());

			auto nodes = blfronts[i]->GetNodes();
			std::array<BLVector, 3> blvs;
			for (int j = 0; j < 3; j++) {
				blvs[j]=nodes[j]->GetCoord(pNodes);
				if (nodes[j] != blNod) {
					blvs[j] = blvs[j] + 1.3*nodes[j]->GetHeight();
				}
			}
			BLVector front_normal = (blvs[2] - blvs[0]) ^ (blvs[1] - blvs[0]);
			if (front_normal * blNod->GetNormal() < 0)
				front_normal = -front_normal;
			neighbout_upper_normals.push_back(front_normal.normalized());
			
		}


		/// 过滤差质量单元
		BLVector norm0(0,0,0);
		double min_dot_product=1;
		for (int i = 0; i < nNeigNods; i++) {
			for (int j = 0; j < nNeigNods; j++) {
				min_dot_product=min(min_dot_product,neighbour_point_normals[i]*neighbout_front_normals[j]);
			}
		}
		//std::cout << min_dot_product<<" ";
		if (min_dot_product > 0.95) {
			for (int i = 0; i < nNeigNods; i++) {
				norm0 += neighbour_point_normals[i];
			}
			norm0.normalize();
		}
		else {


			std::vector<VoronoiPoint> nnormals;


			double all_length = 0;
			double all_height = 0;
			for (int i = 0; i < nNeigNods; i++) {
				int idx = blNodes[i]->GetNodIdx();
				BLVector coord1(pNodes[idx].coord[0], pNodes[idx].coord[1], pNodes[idx].coord[2]);
				double length = (coord1 - centor).magnitude();
				all_length += length;
				double height = blNodes[i]->GetHeightLength();
				all_height += height;
			}



			std::vector<double> cost(nNeigNods);
			double ave_cost = 0;
			for (int i = 0; i < nNeigNods; i++) {
				int idx = blNodes[i]->GetNodIdx();

				int n_front_nodes;
				BLNode* front_nodes[3];
				blfronts[i]->GetNodes(&n_front_nodes, front_nodes);;
				BLVector pt[3];
				for (int j = 0; j < 3; j++) {
					pt[j] = front_nodes[j]->GetCoord(pNodes);
				}

				double aera = BLVector::crossProduct((pt[1] - pt[0]), (pt[2] - pt[0])).magnitude();
				double cir_length = (pt[0] - pt[1]).magnitude() + (pt[2] - pt[1]).magnitude() + (pt[0] - pt[2]).magnitude();
				double ratio = aera / cir_length / cir_length * 18 / sqrt(3);
				cost[i] = ratio;
				ave_cost += ratio;
			}
			ave_cost /= nNeigNods;



			double ave_up_cost = 0;
			for (int i = 0; i < nNeigNods; i++) {
				int idx = blNodes[i]->GetNodIdx();

				int n_front_nodes;
				BLNode* front_nodes[3];
				blfronts[i]->GetNodes(&n_front_nodes, front_nodes);;
				BLVector pt[3];
				for (int j = 0; j < 3; j++) {
					pt[j] = front_nodes[j]->GetCoord(pNodes) + front_nodes[j]->GetNormal() * front_nodes[j]->GetHeightLength();
				}

				double aera = BLVector::crossProduct((pt[1] - pt[0]), (pt[2] - pt[0])).magnitude();
				double cir_length = (pt[0] - pt[1]).magnitude() + (pt[2] - pt[1]).magnitude() + (pt[0] - pt[2]).magnitude();
				double ratio = aera / cir_length / cir_length * 18 / sqrt(3);
				if (ratio > cost[i]) {
					cost[i] = 1;
				}
				else
				{

					//	cost[i] = ratio;// -cost[i] + ratio;
				}

				ave_up_cost += cost[i];
			}
			ave_up_cost /= nNeigNods;



			for (int i = 0; i < nNeigNods; i++) {
				int idx = blNodes[i]->GetNodIdx();

				int n_front_nodes;
				BLNode* front_nodes[3];
				blfronts[i]->GetNodes(&n_front_nodes, front_nodes);;





				BLVector coord1(pNodes[idx].coord[0], pNodes[idx].coord[1], pNodes[idx].coord[2]);
				VoronoiPoint p;
				double length = (coord1 - centor).magnitude();
				double height = blNodes[i]->GetHeightLength();
				p.pos = front_nodes[0]->GetNormal() + front_nodes[1]->GetNormal() + front_nodes[2]->GetNormal() - blNod->GetNormal();
				p.pos.normalize();
				p.func = std::shared_ptr<ExpWightFunc<2>>(new ExpWightFunc<2>());
				//	double t = (height / all_length / 10);
			//		double supposed_length = t * all_length + (1-t) * length;
				if (blNod->GetLowerNode()) {
					//		BLVector p1 = coord1 + blNodes[i]->GetNormal()*blNodes[i]->GetHeightLength();
					//		BLVector p2 =(centor + blNod->GetLowerNode()->GetNormal()*blNod->GetHeightLength());

					ave_up_cost = 1;
					double cost_all = ave_up_cost / cost[i] + cost[i] / ave_up_cost;
					//		cout << cost <<" "<< actual_length<<" "<< supposed_length<< endl;
							//	p.func->setExpoent(ratio2Exponent(height/length));
					p.func->setExpoent(2 / cost_all);// 0.3 + 0.7*(sqrt(length / all_length)));
					//	 p.func->setExpoent((sqrt(length / all_length)));
				}
				else
					p.func->setExpoent(1);// 0.3 + 0.7*(sqrt(length / all_length)));



				nnormals.push_back(p);

			}







			if (blNod->GetLowerNode()) {
				VoronoiPoint p;
				p.pos = blNod->GetLowerNode()->GetNormal();
				p.func = std::shared_ptr<ExpWightFunc<1>>(new ExpWightFunc<1>());
				p.func->setExpoent(200 * all_height / all_length);
				p.func->setExpoent(1);
				nnormals.push_back(p);


				for (const auto i : neighbout_upper_normals) {
					VoronoiPoint p;
					p.pos = i;
					p.func = std::shared_ptr<ExpWightFunc<1>>(new ExpWightFunc<1>());
					//p.func->setExpoent(2);
					p.func->setExpoent(200 * all_height / all_length);
					//std::cout << 200 * all_height / all_length << " ";
					nnormals.push_back(p);
				}
			}
			else {
				for (auto i : neighbout_front_normals) {
					VoronoiPoint p;
					p.pos = i;
					p.func = std::shared_ptr<ExpWightFunc<1>>(new ExpWightFunc<1>());
					p.func->setExpoent(0.01);
					p.func->setExpoent(1);
					nnormals.push_back(p);
				}
			}

			std::vector<BLVector> ans;
			BLVector b1(10, 0, 0);
			BLVector b2(10, 0, 0);
			VoronoiPoint tp1, tp2;

			std::map<pair<int, int>, pair<double, BLVector>> cost_map;





			double max_cost = -100;
			for (int i = 0; i < nnormals.size(); i++) {
				for (int j = i + 1; j < nnormals.size(); j++) {
					double cost = 0;
					auto n = CentorOfTwoPoint(nnormals[i], nnormals[j], cost);
					if (n.x > 9)
						continue;
#if 1
					if (cost > max_cost) {
						max_cost = cost;
						b1 = n;
						tp1 = nnormals[i];
						tp2 = nnormals[j];
						/*	if (ans.size())
								ans.back() = n;
							else
								ans.push_back(n);*/

					}
					//#else
					ans.push_back(n);
#endif

					cost_map[pair<int, int>(std::min(i, j), std::max(i, j))].first = nnormals[i].func->operator()(PointGeodesicLength(nnormals[i].pos, n));
					cost_map[pair<int, int>(std::min(i, j), std::max(i, j))].second = n;// nnormals[i].func->operator()(PointGeodesicLength(nnormals[i].pos, ans.back()));

				}
			}

			max_cost = -100;
			for (int i = 0; i < nnormals.size(); i++) {
				for (int j = i + 1; j < nnormals.size(); j++) {
					for (int k = j + 1; k < nnormals.size(); k++) {

						if (nnormals[i].func->operator()(PointGeodesicLength(cost_map[std::pair<int, int>(j, k)].second, nnormals[i].pos)) > cost_map[std::pair<int, int>(j, k)].first)
							if (nnormals[j].func->operator()(PointGeodesicLength(cost_map[std::pair<int, int>(i, k)].second, nnormals[j].pos)) > cost_map[std::pair<int, int>(i, k)].first)
								if (nnormals[k].func->operator()(PointGeodesicLength(cost_map[std::pair<int, int>(i, j)].second, nnormals[k].pos)) > cost_map[std::pair<int, int>(i, j)].first) {
									double cost = 0;
									auto n = (CentorOfThreePoint(nnormals[i], nnormals[j], nnormals[k], cost));
									if (cost > max_cost) {
										max_cost = cost;
										b2 = n;
									}
									ans.push_back(n);

								}
					}
				}
			}


			//for (int i = 0; i < neighbour_point_normals.size(); i++) {
			//	for (int j = i + 1; j < neighbour_point_normals.size(); j++) {
			//		for (int k = j + 1; k < neighbour_point_normals.size(); k++) {
			//			ans.push_back(CentorOfThreePoint(neighbour_point_normals[i], neighbour_point_normals[j], neighbour_point_normals[k]));
			//		}
			//	}
			//}



			//for (int i = 0; i < neighbout_front_noamls.size(); i++) {
			//	for (int j = i + 1; j < neighbout_front_noamls.size(); j++) {
			//		for (int k = j + 1; k < neighbout_front_noamls.size(); k++) {
			//			ans.push_back(CentorOfThreePoint(neighbout_front_noamls[i], neighbout_front_noamls[j], neighbout_front_noamls[k]));
			//		}
			//	}
			//}
		//	ans.clear();
			//auto ansp = ans;
			//ansp.clear();
			//ans.push_back(b1);
			//ans.push_back(b2);


			
			double ncost = 1000;
			for (auto i : ans) {
				double min_cost = -1e20;
				for (auto j : nnormals) {
					if (min_cost < j.func->operator()(PointGeodesicLength(i, j.pos))) {
						min_cost = j.func->operator()(PointGeodesicLength(i, j.pos));
					}
				}
				if (min_cost < ncost && min_cost != -1e20) {
					ncost = min_cost;
					norm0 = i;
				}
			}



			//	cout <<  ncost << " "<< nnormals.size()<<endl;

				//BLVector sum(0, 0, 0);
				//double height = 0;
				//int count = 0;
				//const std::vector <BLFront*>& blfronts = blNod->GetNeigFronts();
				//// 周围点数应该等于周围前沿数
				//for (int i = 0; i < nNeigNods; i++) {
				//	min_front_size = min(min_front_size, blfronts[i]->GetMinFrontSize());
				//}
				////frontsize /= nNeigFrts;
				//BLVector desentNode(pNodes[blNod->GetDecentID()].coord);
				//double height1 = (desentNode - centor).magnitude();

				///*高度宽度比*/
				//double ratio = height1 / min_front_size;
				//static double lengths[200];
				////std::vector<double> lengths;
				//double length_sum = 0;
				//for (int i = 0; i < nNeigNods; i++) {
				//	int idx = blNodes[i]->GetNodIdx();
				//	BLVector coord1(pNodes[idx].coord[0], pNodes[idx].coord[1], pNodes[idx].coord[2]);
				//	double length = (coord1 - centor).magnitude2();
				//	lengths[i] = length;
				//	length_sum += length;
				//}
				//length_sum /= nNeigNods;

				//for (int i = 0; i < nNeigNods; i++) {
				//	count++;
				//	height += blNodes[i]->GetHightRatio();
				//	int id = blNodes[i]->GetNodIdx();

				//	double sp = blNodes[i]->GetBeitaVisu() * 2 / (PI);

				//	int idx = blNodes[i]->GetNodIdx();
				//	BLVector coord1(pNodes[idx].coord[0], pNodes[idx].coord[1], pNodes[idx].coord[2]);
				//	double length = lengths[i] / length_sum;
				//	double angle = (coord1 - centor).normalized()*blNodes[i]->GetNormal();

				//	sum += blNodes[i]->GetNormal() * (atan(0.01*pow(length, 5 * angle) / (sp*sp)));


				//}
				//sum.normalize();
				//if (blNod->GetLowerNode()) {
				//	sum += blNod->GetLowerNode()->GetNormal();
				//	sum.normalize();
				//}
				//double r_borm;
				//r_borm = (pow(1.7, ratio) - 1) * 15;
				//BLVector norm0 = (r_borm *sum + blNod->GetNormal());
				//norm0.normalize();
				//int t = 0;

				//constexpr double output_thredhold = thredhold * 0.8;

				//while (blNod->GetBeitaVisu(norm0) * 180 / PI < output_thredhold) {
				//	if (t > 10) {
				//		return;
				//	}
				//	norm0 = (norm0 + 0.7*blNod->GetNormal());
				//	norm0.normalize();
				//	t++;
				//}
			//if ((norm0 - b1).magnitude2()>1e-3) {
			//	double min_cost2 = -1e20;
			//	for (auto j : nnormals) {
			//		if (min_cost2 < j.func->operator()(PointGeodesicLength(norm0, j.pos))) {
			//			min_cost2 = j.func->operator()(PointGeodesicLength(norm0, j.pos));
			//		}
			//	}

			//	double min_cost1 = -1e20;
			//	for (auto j : nnormals) {
			//		if (min_cost1 < j.func->operator()(PointGeodesicLength(b1, j.pos))) {
			//			min_cost1 = j.func->operator()(PointGeodesicLength(b1, j.pos));
			//		}
			//	}



			//	cout <<"single"<< std::max(tp1.func->operator()(PointGeodesicLength(b1, tp1.pos)), tp2.func->operator()(PointGeodesicLength(b1, tp2.pos)));
			//	//	cout << endl;
			//	cout << ":" << min_cost1 << "             all";
			//cout<<	std::max(tp1.func->operator()(PointGeodesicLength(norm0, tp1.pos)), tp2.func->operator()(PointGeodesicLength(norm0, tp2.pos)));
			//cout << ":" << min_cost2 << endl; ;

			//}
		}
		if (blNod->GetBSys())
		{
			norm0[blNod->GetSymAxis()] = 0;
		}

		if (norm0*  blNod->GetNormal() < 0.999) {
			for (int i = 0; i < nNeigNods; i++) {
				blNodes[i]->if_need_smooth = true;

			}
			blNod->if_need_smooth = true;

		}

		
		myvec[iq] = norm0;
	//	blNod->SetNormal(norm0);
	}

	VoronoiNormalSmoothStategy::~VoronoiNormalSmoothStategy()
	{

	}
}
