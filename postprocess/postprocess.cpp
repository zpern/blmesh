#include "postprocess.h"
#include "../dtiso3d/geom/geom_func.h"
#include <spdlog/spdlog.h> 
 #include <map>
#include <set>
#include <algorithm>
#include "../common/singleton_terminate.h"
#ifndef PI
#define PI 3.14159265358979323846
#endif
double Quality(std::vector<BLVector> coordinate) {
	double c[4][3];
	for (int j = 0; j < 4; j++) {
		c[j][0] = coordinate[j].x;
		c[j][1] = coordinate[j].y;
		c[j][2] = coordinate[j].z;
	}
	int edge[6][2] = { {0,1},{1,2},{2,0},{0,3},{1,3},{2,3} };
	double volume = GEOM_FUNC::orient3d(c[1], c[0], c[2], c[3]);
	double length = 0;
	for (int i = 0; i < 6; i++) {
		double l = 0;
		for (int j = 0; j < 3; j++) {
			l += (pow(c[edge[i][0]][j] - c[edge[i][1]][j], 2));
		}
		length += sqrt(l);
	}
	length /= 6;
	return volume / pow(length, 3);
}
bool MeshOptimize(std::vector<BLVector> &coordinate, std::vector<std::vector<int>> &connector, std::set<int>& boundary_point)
{
	std::vector<bool> point_map(coordinate.size(),false);// true: optimize ; false: do nothing 
	std::vector<pair<double,int>> point_quality; // just for selceting the worst point
	for (auto i : connector) {
		if (i.size() == 4) {
			//calculate quality
			vector<BLVector> normal(4);
			//minsin
			double minsin = 0;
			for (int j = 0; j < 4; j++) {
				normal[j] = (coordinate[i[j]] - coordinate[i[(j + 1) % 4]]) ^ (coordinate[i[j]] - coordinate[i[(j + 3) % 4]]);
				normal[j].normalize();
			}

			for (int j = 0; j < 4; j++) {
				minsin = std::max(minsin, abs(normal[j] * normal[(j + 1) % 4]));
			}
			if (minsin > cos(PI / 180 * 5)) {
			//	cout << "Low quality tetra found! Quality=" << minsin << endl;
				for (int j = 0; j < 4; j++) {
					point_quality.push_back(pair<double, int>(minsin,i[j]));					
				}
			}
		}
	}
	if (checkterminate()) {
		return false;
	}
	std::sort(point_quality.rbegin(), point_quality.rend());
	std::size_t MAX_POINT = 200;
	for (std::size_t i = 0; i < std::min(point_quality.size(), MAX_POINT); i++) {
		if(boundary_point.find(point_quality[i].second)== boundary_point.end())
			point_map[point_quality[i].second] = true;
	}
	int split_prism[3][4] = { {1,0,2,5}, {3,4,5,0},{5,4,1,0} };
	int split_pyramid[4][4] = { {0,1,2,4}, {0,2,3,4}, {0,1,3,4}, {1,2,3,4} };
	std::vector<vector<int>> extra_ele;
	for (auto i : connector) {
		if (checkterminate()) {
			return false;
		}
		bool need_split = false;
		if (i.size() == 6) {//prism
			for (int j = 0; j < 6; j++) {
				if (point_map[i[j]])
					need_split = true;
			}
			if (need_split) {
				for (int k = 0; k < 3; k++) {
					vector<int> e;
					for (int l = 0; l < 4; l++) {
						e.push_back(i[split_prism[k][l]]);
					}
					extra_ele.push_back(e);
				}
			}
		}
		if (i.size() == 5) {//pyramid
			for (int j = 0; j < 5; j++) {
				if (point_map[i[j]])
					need_split = true;
			}
			if (need_split) {
				for (int k = 0; k < 4; k++) {
					vector<int> e;
					for (int l = 0; l < 4; l++) {
						e.push_back(i[split_pyramid[k][l]]);
					}
					extra_ele.push_back(e);
				}
			}
		}
	}
	auto quality = [coordinate](std::vector<int> tetra) {
		double c[4][3];
		for (int j = 0; j < 4; j++) {
			c[j][0] = coordinate[tetra[j]].x;
			c[j][1] = coordinate[tetra[j]].y;
			c[j][2] = coordinate[tetra[j]].z;
		}
		int edge[6][2] = { {0,1},{1,2},{2,0},{0,3},{1,3},{2,3} };
		double volume = GEOM_FUNC::orient3d(c[0], c[1], c[2], c[3]);
		double length=0;
		for (int i = 0; i < 6; i++) {
			for (int j = 0; j < 3; j++) {
				length += sqrt(pow(c[edge[i][0]][j]-c[edge[i][1]][j],2));
			}
		}
		length /= 6;
		return volume/pow(length,3);
	};
	if (checkterminate()) {
		return false;
	}
	//averagy edge length
	auto edgeLength = [coordinate](std::vector<int> tetra) {
		double c[4][3];
		for (int j = 0; j < 4; j++) {
			c[j][0] = coordinate[tetra[j]].x;
			c[j][1] = coordinate[tetra[j]].y;
			c[j][2] = coordinate[tetra[j]].z;
		}
		int edge[6][2] = { {0,1},{1,2},{2,0},{0,3},{1,3},{2,3} };
		double length = 0;
		for (int i = 0; i < 6; i++) {
			double l = 0;
			for (int j = 0; j < 3; j++) {
				l += (pow(c[edge[i][0]][j] - c[edge[i][1]][j], 2));
			}
			length += sqrt(l);
		}
		length /= 6;
		return length;
	};
	std::map<int, vector<vector<int>>> spheres;
	std::map<int,pair<double, double>> quality_compare;


	for (auto j : connector) {
		if (j.size() == 4)
			for (int k : j) {
				if (point_map[k])
					spheres[k].push_back(j);
			}
	}


	for (auto j : extra_ele) {
		for (int k : j) {
			if(point_map[k])
				spheres[k].push_back(j);
		}
	}


	for (int p = 0; p < point_map.size(); p++) {
		if (point_map[p]) {
			//find sphere
			double min_volume = std::numeric_limits<double>::max();
			for (auto &i : spheres[p]) {
				min_volume = std::min(min_volume, Quality({ coordinate[i[0]],coordinate[i[1]],coordinate[i[2]],coordinate[i[3]] }));
			}
			cout.precision(15);
			quality_compare[p].first = min_volume;
			//cout << min_volume << endl;
		}
	}

	const int MAX_OPT_LOOP = 3;
	int loop = MAX_OPT_LOOP;
	while (loop--) {
		if (checkterminate()) {
			return false;
		}
		for (auto bad_node_sphere : spheres) {
			if (checkterminate()) {
				return false;
			}
			int idx = bad_node_sphere.first;
			auto sphere = bad_node_sphere.second;
			// find worst quality
			vector<int> worst;
			double min_cost = 1e20;
			for (auto i : sphere) {
				double cost = Quality({ coordinate[i[0]],coordinate[i[1]],coordinate[i[2]],coordinate[i[3]] });
				if (cost < min_cost) {
					min_cost = cost;
					worst = i;
				}
			}
			BLVector face_normal;
			for (int k = 0; k < 4; k++) {
				if (worst[k] == idx) {
					BLVector nodes[3];
					for (int j = 0; j < 3; j++) {
						nodes[j] = coordinate[worst[(k + j + 1) % 3]];
					}
					face_normal = -(nodes[1] - nodes[0]) ^ (nodes[2] - nodes[1]);
					face_normal.normalize();
				}
			}
			
			

			double step = edgeLength(worst)*0.01;
			BLVector old_coordinate = coordinate[idx];


			double start_step = 0 ;
			double end_step = step * 10;
			double old_quality = min_cost;


			while (end_step - start_step > step * 1e-5) {
				double current_step = end_step / 2 + start_step / 2;
				BLVector new_coord = current_step * face_normal + coordinate[idx];

				std::vector<std::vector<BLVector>> new_sphere;
				for (auto i : sphere) {
					std::vector<BLVector> c;
					for (int j = 0; j < 4; j++) {
						if (i[j] == idx) {
							c.push_back(new_coord);
						}
						else
							c.push_back(coordinate[i[j]]);
					}
					new_sphere.push_back(c);
				}
				bool is_too_big = false;
				std::vector<int> current_worst;
				double minq = 1e20;
				for (int i = 0; i < new_sphere.size(); i++) {
					double q=Quality(new_sphere[i]);
					if (q < 0)
						is_too_big = true;
					if (minq > q) {
						minq = q;
						current_worst = sphere[i];
					}
				}
				//if (current_worst != worst)
				//	toobig = true;
				if (minq > old_quality) {
					old_quality = minq;
				}
				else
					is_too_big = true;

				if (is_too_big)
					end_step = current_step;
				else
					start_step = current_step;
			}
			double mid_step = start_step;
			coordinate[idx] = coordinate[idx]+ mid_step*face_normal;
		}
	}


	for (int p = 0; p < point_map.size(); p++) {
		if (point_map[p]) {
			//find sphere
			double min_volume = 1e20;
			for (auto &i : spheres[p]) {				
				min_volume = std::min(min_volume, Quality({ coordinate[i[0]],coordinate[i[1]],coordinate[i[2]],coordinate[i[3]] }));
			}
			cout.precision(15);
			quality_compare[p].second = min_volume;
			//cout << min_volume << endl;
		}
	}
	double worst_before =10;
	double worst_after = 10;

	for (auto i : quality_compare) {
		worst_before =std::min(worst_before, i.second.first);
		worst_after = std::min(worst_after,  i.second.second);
	}
	spdlog::info("Optimization point size = {}",quality_compare.size());
	spdlog::info("{}====>{}",worst_before, worst_after);

	if ( worst_after > worst_before * 1.3 )
		return true;
	return false;

}
