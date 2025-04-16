#include "MNormalMesh.h"
#include "spdlog/spdlog.h"
#include "geometryfunction.h"
#include "../include/MeshEvaluation.h"
#include <fstream>
#include <exception>
#include <sstream>
#include <map>
#include <algorithm>
#include <stack>
#include <array>
#include <queue>
#include <unordered_map>
#include <unordered_set>
using std::stringstream;
using std::ifstream;
using std::ofstream;
using std::map;
using std::array;
void MNormalMesh::SetBehavior(ChamferBehavior& behavior)
{
	behavior_ = behavior;
}

void MNormalMesh::ReadPlsBuf(std::string f,
                             std::vector<std::array<double, 3>>& points) {
        std::istringstream fin(f);

        int number_of_point, number_of_element;
        string line;
        getline(fin, line);
        stringstream ss(line);
        ss >> number_of_element >> number_of_point;

        number_of_node = number_of_point;
        number_of_origin_node = number_of_node;

        number_of_triangles = number_of_element;
        number_of_origin_triangles = number_of_element;

        node_array.resize(number_of_point);
        vector<vector<std::array<int, 3>>> graph;
        graph.resize(number_of_point);
        point_normals.resize(number_of_point);
        coordinate.resize(number_of_point);
        real_node_id_.resize(number_of_point);
        for (int i = 0; i < number_of_point; i++) {
        //  fin >> line >> node_array[i].coordinate.x >>
       //       node_array[i].coordinate.y >> node_array[i].coordinate.z;
        node_array[i].coordinate.x = points[i][0];
        node_array[i].coordinate.y = points[i][1];
        node_array[i].coordinate.z = points[i][2];
          for (int k = 0; k < 3; k++)
            coordinate[i][k] = node_array[i].coordinate[k];
          node_array[i].node_id_ = i;
          real_node_id_[i] = i;  // of course they are themself
        }
        connector.resize(number_of_element);
        attribute.resize(number_of_element);
        for (int i = 0; i < number_of_element; i++) {
          fin >> line >> connector[i][0] >> connector[i][1] >>
              connector[i][2] >> attribute[i];

          for (int k = 0; k < 3; k++) connector[i][k]--;  // start from 0
          for (int k = 0; k < 3; k++)
            graph[connector[i][k]].push_back(connector[i]);
          for (int k = 0; k < 3; k++)
            node_array[connector[i][k]].neighbour_front_index_.push_back(i);
        }

        for (int i = 0; i < number_of_point; i++) {
          for (auto &j : graph[i]) {
            if (j[0] == i) {
              j[0] = j[1];
              j[1] = j[2];
            } else if (j[1] == i) {
              j[1] = j[0];
              j[0] = j[2];
            }
          }
          int start = graph[i][0][0];
          int pos = start;
          while (true) {
            for (auto &j : graph[i]) {
              if (j[0] == pos) {
                pos = j[1];
                node_array[i].neighbour_node_.push_back(node_array.begin() +
                                                        j[1]);
                break;
              }
            }
            if (pos == start) {
              break;
            }
          }
        }
        CaculateFrontNormal();
        CalculateNodeNormal();
        return;
}
void MNormalMesh::ReadPls(std::string filename)
{

	ifstream fin(filename);
	if (!fin.is_open() || filename == "") {
		fin = ifstream(behavior_.input_filename);
		if (!fin.is_open()) {
			throw std::runtime_error("no such file");
		}
	}
	int number_of_point, number_of_element;
	string line;
	getline(fin, line);
	stringstream ss(line);
	ss >> number_of_element >> number_of_point;

	number_of_node = number_of_point;
	number_of_origin_node = number_of_node;

	number_of_triangles = number_of_element;
	number_of_origin_triangles = number_of_element;


	node_array.resize(number_of_point);
	vector<vector<std::array<int, 3>>> graph;
	graph.resize(number_of_point);
	point_normals.resize(number_of_point);
	coordinate.resize(number_of_point);
	real_node_id_.resize(number_of_point);
	for (int i = 0; i < number_of_point; i++) {
		fin >> line >> node_array[i].coordinate.x >> node_array[i].coordinate.y >> node_array[i].coordinate.z;
		for (int k = 0; k < 3; k++)
			coordinate[i][k] = node_array[i].coordinate[k];
		node_array[i].node_id_ = i;
		real_node_id_[i] = i;//of course they are themself
	}
	connector.resize(number_of_element);
	attribute.resize(number_of_element);
	for (int i = 0; i < number_of_element; i++) {

		fin >> line >> connector[i][0] >> connector[i][1] >> connector[i][2] >> attribute[i];

		
		for (int k = 0; k < 3; k++)
			connector[i][k]--;//start from 0
		for (int k = 0; k < 3; k++)
			graph[connector[i][k]].push_back(connector[i]);
		for (int k = 0; k < 3; k++)
			node_array[connector[i][k]].neighbour_front_index_.push_back(i);
	}

	for (int i = 0; i < number_of_point; i++) {
		
		for (auto &j : graph[i]) {
			if (j[0] == i) {
				j[0] = j[1];
				j[1] = j[2];
			}
			else if (j[1] == i) {
				j[1] = j[0];
				j[0] = j[2];
			}
		}
		int start=graph[i][0][0];
		int pos = start;
		while (true) {
			for (auto &j : graph[i]) {				
				if (j[0] == pos) {
					pos = j[1];
					node_array[i].neighbour_node_.push_back(node_array.begin() + j[1]);
					break;
				}
			}
			if (pos == start) {
				break;
			}
		}
	}
	CaculateFrontNormal();
	CalculateNodeNormal();
	return;
}

void MNormalMesh::WriteNorm()
{
	string filename = behavior_.getOutputName(".norm");

	ofstream fout(filename);	
	fout.precision(13);
	for (int i = 0; i < point_normals.size(); i++) {
		fout << i << " " << point_normals[i][0] << " " << point_normals[i][1] << " " << point_normals[i][2] << endl;
	}
	
	fout.close();
	spdlog::info("Writing {0}",filename);
}

void MNormalMesh::WriteMem(std::string& f,
                           std::vector<std::array<double, 3>>& points) {

	WriteMem(f, points, POINT_OFFSET);
}

void MNormalMesh::WriteVol(std::vector<std::array<double, 3>>& v, std::vector<std::vector<int>>& f, int& lower_num, double len,int& add_point_num)
{
	std::map<std::array<double, 3>, int> coord_to_id;
	std::vector<std::array<int, 3>> lower_ids(connector.size());
	for (int i = 0; i < connector.size(); i++) {
		int count = 0;
		int lowerid_1=-1,lower_id2=-1;
		std::array<int, 3>& ids= lower_ids[i];
		for (int k = 0; k < 3; k++) {
			std::array<double, 3> v;
			v[0] = coordinate[connector[i][k]].x;
			v[1] = coordinate[connector[i][k]].y;
			v[2] = coordinate[connector[i][k]].z;
			ids[k] = -1;
			if (coord_to_id.find(v) == coord_to_id.end()) {
				ids[k] = coord_to_id.size();
				coord_to_id[v] = ids[k];
			}
			else
				ids[k] = coord_to_id[v];
		}
	}
	lower_num = coord_to_id.size();
	v.resize(lower_num + coordinate.size());
	for (auto i : coord_to_id) {
		v[i.second] = i.first;
	}
	for (int i = 0; i < coordinate.size(); i++) {

		v[lower_num+i]={ coordinate[i].x + len * point_normals[i].x,
			coordinate[i].y + len * point_normals[i].y, coordinate[i].z + len * point_normals[i].z };
	}

	add_point_num=0;
	for (int i = 0; i < connector.size(); i++) {
		// 四面体情况
		if (lower_ids[i][0] == lower_ids[i][1] && lower_ids[i][2] == lower_ids[i][1]) {
			f.push_back({ lower_ids[i][0],connector[i][0]+lower_num,connector[i][1] + lower_num ,connector[i][2] + lower_num });
		}
		else if (lower_ids[i][0] != lower_ids[i][1] && lower_ids[i][2] != lower_ids[i][1]&& lower_ids[i][2] != lower_ids[i][0]) {
			f.push_back({ lower_ids[i][0] ,lower_ids[i][1]  ,lower_ids[i][2]
					,connector[i][0] + lower_num,connector[i][1] + lower_num ,connector[i][2] + lower_num });
		}
		else {/// 复杂情况
			int k1, k2,k3;
			for (int j = 0; j < 3; j++) {
				if (lower_ids[i][j] == lower_ids[i][(j + 1) % 3]) {
					k1 = j;
					k2 = (j + 1) % 3;
					k3 = (j + 2) % 3;
					break;
				}
			}
			f.push_back({ lower_ids[i][k3] ,lower_ids[i][k1] ,connector[i][k3] + lower_num,connector[i][k1] + lower_num ,connector[i][2] + lower_num,(int)v.size()});
			f.push_back({ lower_ids[i][k2] ,lower_ids[i][k3] ,connector[i][k2] + lower_num,connector[i][k3] + lower_num ,connector[i][2] + lower_num,(int)v.size() });
			f.push_back({ lower_ids[i][k2] + lower_num,lower_ids[i][k1]+ lower_num ,connector[i][k3] + lower_num,(int)v.size() });
			std::array<double, 3> ncoord{0,0,0};
			for (int k = 0; k < 3; k++) {
				for (int j = 0; j < 3; j++) {
					ncoord[k] += coordinate[connector[i][j]][k];
					ncoord[k] += coordinate[connector[i][j]][k] + point_normals[connector[i][j]][k] * len;
				}
			}
			for (int k = 0; k < 3; k++) {
				ncoord[k] /= 6;
			}
			add_point_num++;
			v.push_back(ncoord);			
		}
	}
}

void MNormalMesh::WriteMem(std::string& f, std::vector<std::array<double, 3>>& points, double len)
{

	points.resize(coordinate.size());
	for (int k = 0; k < coordinate.size(); k++) {
		for (int i = 0; i < 3; i++) {
			points[k][i] = coordinate[k][i] +
				len * point_normals[k][i];
		}
	}

	std::ostringstream ss;
	ss << connector.size() << " " << coordinate.size() << " "
		<< "0 0 0 0"
		<< std::endl;
	for (int i = 0; i < connector.size(); i++) {
		ss << i + 1 << " " << connector[i][0] + 1 << " "
			<< connector[i][1] + 1 << " " << connector[i][2] + 1
			<< " " << attribute[i] << std::endl;
	}
	f = ss.str();

}

void MNormalMesh::WriteVtk()
{
	string filename = behavior_.getOutputName(".vtk");
	int i, j, npt, nelm, sidx;
	FILE *fout = nullptr;

	npt = coordinate.size();
	//nelm = m_nElems;
	nelm = connector.size();

	fout = fopen(filename.data(), "w");
	if (fout == nullptr)
	{
		printf("Can not open file %s\n", filename.data());
	}

	//vtk
	::fprintf(fout, "# vtk DataFile Version 2.0\n");
	::fprintf(fout, "boundary layer mesh\n");
	::fprintf(fout, "ASCII\n");
	::fprintf(fout, "DATASET UNSTRUCTURED_GRID\n");
	::fprintf(fout, "POINTS %d double\n", npt);

	for (i = 0; i < npt; i++)
	{
		::fprintf(fout, "%f %f %f\n", coordinate[i][0]+POINT_OFFSET*point_normals[i][0], coordinate[i][1]+POINT_OFFSET*point_normals[i][1], coordinate[i][2]+POINT_OFFSET*point_normals[i][2]);
	}

	//sidx = m_nSurfElems;
	sidx = 0;

	//sidx = 2914;
	//nelm = sidx + 100;
	::fprintf(fout, "CELLS %d %d\n", nelm - sidx, (nelm - sidx) * 4);
	for (i = 0; i < nelm - sidx; i++)
		//for (i=sidx; i<nelm; i++)
	{
		::fprintf(fout, "3 %d %d %d\n", connector[i][0],
			connector[i][1], connector[i][2]);
	}

	::fprintf(fout, "CELL_TYPES %d\n", nelm - sidx);

	for (i = 0; i < nelm - sidx; i++)
		::fprintf(fout, "%d\n", 5);

#if 1
	//for testing
	::fprintf(fout, "POINT_DATA %d\n", npt);
	 		//::fprintf(fout, "SCALARS fixed int\n");
	 		//::fprintf(fout, "LOOKUP_TABLE default\n");
	 		//for (i=0; i<nelm; i++)
	 		//{
	 		//	::fprintf(fout, "%d\n", attribute[i]);
	 		//}
	::fprintf(fout, "NORMALS node_normals double\n");
	for (i = 0; i < npt; i++)
	{
		
		::fprintf(fout, "%lf %lf %lf\n", point_normals[i].x, point_normals[i].y, point_normals[i].z);
	}
#endif
#if 0
	//for testing
	::fprintf(fout, "CELL_DATA %d\n", nelm);
	// 		::fprintf(fout, "SCALARS fixed double\n");
	// 		::fprintf(fout, "LOOKUP_TABLE default\n");
	// 		for (i=0; i<npt; i++)
	// 		{
	// 			::fprintf(fout, "%lf\n", m_pNodes[i].uvalue);
	// 		}
	::fprintf(fout, "NORMALS face_normals double\n");
	for (i = 0; i < nelm; i++)
	{
		BLVector norm = m_pElems[i].norm;
		::fprintf(fout, "%lf %lf %lf\n", norm.x, norm.y, norm.z);
	}
#endif

	fclose(fout);
	fout = nullptr;
	spdlog::info("Writing {0}",filename);
}

void MNormalMesh::WritePls(string filename)
{
	if (filename == "")
		filename = behavior_.getOutputName();
	ofstream fout(filename);
	fout.precision(13);
	fout << connector.size() << " " << coordinate.size() << " 0 0 0 0" << endl;

	for (int i = 0; i < coordinate.size(); i++) {
		fout << i+1 << " " << coordinate[i][0]+ POINT_OFFSET*point_normals[i][0] << " " << coordinate[i][1] + POINT_OFFSET * point_normals[i][1] << " " << coordinate[i][2] + POINT_OFFSET * point_normals[i][2] << endl;
	}
	for (int i = 0; i < connector.size(); i++) {
		fout << i + 1 << " " << connector[i][0]+1 << " " << connector[i][1] + 1 << " " << connector[i][2] + 1 <<" "<<attribute[i]<< endl;
	}
	fout.close();
	spdlog::info("Writing {0}",filename);
}

void MNormalMesh::CalculateMultiNormal()
{
	clock_t start_time = clock();
	MeshEvaluator::GetSingleton().SetDefaultParameter(1e-5, 0.1);
	vector<vector<ComplexNode>::iterator> need_to_split;

	//////////////////////////////////////////////////////////////////////////
	/* instead of performing multiple normal calculation to every node, we only
		do the splitting operation on node with poor quality */
	//////////////////////////////////////////////////////////////////////////
	double d = 0;
	for (auto i = node_array.begin(); i != node_array.end(); i++) {
		d = max(d, i->original_skewness_);
		if (i->original_skewness_ > SKEWNESS_THREADHOLD)// naca 0.14459
			need_to_split.push_back(i);
	}
	//spdlog::info("//////////////////////////////////");
	spdlog::info("The maximum skewness of single normal method = {0}", d);
	//spdlog::info("//////////////////////////////////");

	int number_of_complex_nodes = 0;
	for (auto i : need_to_split) {
		number_of_complex_nodes++;
		i->SplitNode();	
	}
	/*
	for (auto i : need_to_split) {
		if (!i->getFinalMesh().triangle_lists_.empty()) {
			cout << i->getFinalMesh().max_skewness_ << endl;
		}
		else
			cout << i->original_skewness_ << endl;
	}
	*/
	clock_t end_time = clock();
	spdlog::info("");
	spdlog::info("");
	spdlog::info("/////////////////////////////////");


	spdlog::info("The number of Complex node={0}.", number_of_complex_nodes);
	spdlog::info("The time comsuption of splitting={:03.2f} sec.", (end_time - start_time)*1.0/CLOCKS_PER_SEC);
}

void MNormalMesh::PrintVisibilityConeInfo()
{
	
}


void MNormalMesh::BuildTopo()
{
	auto max_attribute = max_element(attribute.begin(), attribute.end());
	int new_attribute = *max_attribute + 1;
	map<int, map<int, vector<int>>> complex_edges;

	for (auto i = node_array.begin(); i != node_array.end(); i++)
	{
		auto& local_mesh = i->getFinalMesh();
		if (local_mesh.getExtraPointCount())
		{
			std::vector<int> point_index_map(local_mesh.getValidPointCount());
			std::map<int, int> valid_map;

			point_index_map[0] = i->node_id_;
			for (int j = 0; j < local_mesh.getExtraPointCount() - 1; j++)
			{
				point_index_map[j + 1] = coordinate.size();
				coordinate.push_back(i->coordinate);
				point_normals.push_back(BLVector());
				real_node_id_.push_back(i->node_id_);
			}
			real_node_id_[i->node_id_] = i->node_id_;

			int count = -1;
			for (int k = 0; k < local_mesh.virtual_point_lists_.size(); k++)
			{
				if (local_mesh.valid.find(k) == local_mesh.valid.end() ||
					local_mesh.virtual_point_lists_[k].isFarNode())
					continue;
				count++;
				valid_map[k] = count;
				local_mesh.virtual_point_lists_[k].setGlobalIndex(point_index_map[count]);

				if ((!local_mesh.virtual_point_lists_[k].getGlobalNeighbourTriIndex().empty()) ||
					(local_mesh.valid.find(ONE_INSERT) != local_mesh.valid.end() &&
						k == local_mesh.virtual_point_lists_.size() - 1))
				{
					coordinate[point_index_map[count]] = coordinate[point_index_map[count]];
					point_normals[point_index_map[count]] = local_mesh.virtual_point_lists_[k].getCoord();
					if (local_mesh.virtual_point_lists_[k].getCoord().magnitude2() < 0.9 ||
						local_mesh.virtual_point_lists_[k].getCoord().magnitude2() > 1.1)
						throw std::logic_error("zero vector found!");
					for (auto j : local_mesh.virtual_point_lists_[k].getGlobalNeighbourTriIndex())
					{
						for (int l = 0; l < 3; l++)
							if (connector[j][l] == i->node_id_)
							{
								connector[j][l] = point_index_map[count];
							}
					}
				}
			}
		}
	}
	// mesh stitching

	struct pos
	{
		int original_index;
		BLVector left;
		BLVector right;
	};
	std::map<std::array<int, 2>, pos> complex_edges_pair;
	map<std::array<int, 2>, int> global_far_node_map; // store global point idx for each virtual edge
	for (auto i = node_array.begin(); i != node_array.end(); i++)
	{
		auto& local_mesh = i->getFinalMesh();
		if (local_mesh.getExtraPointCount())
		{
			for (int k1 = 0; k1 < local_mesh.boundary_edges_.size(); k1++)
			{
				int k = local_mesh.boundary_edges_[k1][0];
				if (local_mesh.valid.find(k) == local_mesh.valid.end() ||
					!local_mesh.virtual_point_lists_[k].isFarNode()) // far point
					continue;
				int index = local_mesh.virtual_point_lists_[k].getGlobalIndex();
				if (node_array[index].getFinalMesh().getExtraPointCount())
				{
					complex_edges_pair[std::array<int, 2>{i->node_id_, index}] = pos();
					// for (auto l : node_array[index].getFinalMesh().boundary_edges_) {
					//	// TODO should not be the first edge, should sort
					//	if
					//(node_array[index].getFinalMesh().virtual_point_lists_[l[1]].isFarNode()&&node_array[index].getFinalMesh().virtual_point_lists_[l[1]].getGlobalIndex()
					//== i->node_id_) 		global_far_node_map[std::array<int, 2>{i->node_id_, index}] =
					//node_array[index].getFinalMesh().virtual_point_lists_[l[0]].getGlobalIndex();
					// }
				}
				else
					global_far_node_map[std::array<int, 2>{i->node_id_, index}] =
					index; // do nothing ,just record the index to map
			}
		}
	}

	// enumate all the possible item
	static vector<vector<int>> connection[2][2] = { vector<vector<int>>() };
	static bool create = false;
	if (!create)
	{
		create = true;
		connection[0][0].push_back(vector<int>{-1, 1});
		connection[0][0].push_back(vector<int>{1, -1});

		connection[0][1].push_back(vector<int>{1, -2, -1});
		connection[0][1].push_back(vector<int>{-2, 1, -1});
		connection[0][1].push_back(vector<int>{-2, -1, 1});

		connection[1][0].push_back(vector<int>{1, 2, -1});
		connection[1][0].push_back(vector<int>{1, -1, 2});
		connection[1][0].push_back(vector<int>{-1, 1, 2});

		connection[1][1].push_back(vector<int>{-2, -1, 1, 2});
		connection[1][1].push_back(vector<int>{1, -2, -1, 2});
		connection[1][1].push_back(vector<int>{-2, 1, -1, 2});
		connection[1][1].push_back(vector<int>{1, 2, -2, -1});
		connection[1][1].push_back(vector<int>{1, -2, 2, -1});
		connection[1][1].push_back(vector<int>{-2, 1, 2, -1});
	}

	for (auto edge : complex_edges_pair)
	{
		int s = edge.first[0];
		int e = edge.first[1];
		if (s > e)
			continue;
		auto& meshs = node_array[s].getFinalMesh();
		auto& meshe = node_array[e].getFinalMesh();
		// Here we should care about oritation

		vector<pair<int, int>>
			active_triangles_left; // fist is triangle index, second is 0-3 point index in the triangle
		for (int i = 0; i < node_array[s].getFinalMesh().triangle_lists_.size(); i++)
		{
			for (int j = 0; j < 3; j++)
			{
				if (node_array[s]
					.getFinalMesh()
					.virtual_point_lists_[node_array[s].getFinalMesh().triangle_lists_[i].point_index_[j]]
					.getGlobalIndex() == e)
				{
					active_triangles_left.push_back(pair<int, int>{i, j});
				}
			}
		}
		vector<pair<int, int>>
			active_triangles_right; // fist is triangle index, second is 0-3 point index in the triangle
		for (int i = 0; i < node_array[e].getFinalMesh().triangle_lists_.size(); i++)
		{
			for (int j = 0; j < 3; j++)
			{
				if (node_array[e]
					.getFinalMesh()
					.virtual_point_lists_[node_array[e].getFinalMesh().triangle_lists_[i].point_index_[j]]
					.getGlobalIndex() == s)
				{
					active_triangles_right.push_back(pair<int, int>{i, j});
				}
			}
		}
		// Sort the active triangles, TODO ::max triangles = 2, do not support more than 2 triangles
		if (active_triangles_left.size() == 2)
			if (meshs.triangle_lists_[active_triangles_left[0].first]
				.point_index_[(active_triangles_left[0].second + 1) % 3] ==
				meshs.triangle_lists_[active_triangles_left[1].first]
				.point_index_[(active_triangles_left[1].second + 2) % 3])
			{
				swap(active_triangles_left[0], active_triangles_left[1]);
			}
		if (active_triangles_right.size() == 2)
			if (meshe.triangle_lists_[active_triangles_right[0].first]
				.point_index_[(active_triangles_right[0].second + 1) % 3] ==
				meshe.triangle_lists_[active_triangles_right[1].first]
				.point_index_[(active_triangles_right[1].second + 2) % 3])
			{
				swap(active_triangles_right[0], active_triangles_right[1]);
			}
		if ((active_triangles_left.size() > 2) || active_triangles_right.size() > 2)
			throw std::runtime_error("wake up and write code for more than 2!");

		auto combination = connection[active_triangles_left.size() - 1][active_triangles_right.size() - 1];

		vector<int> best_combination;
		double min_cost = 100;
		// find best combination
		for (auto c : combination)
		{
			// TODO : there should consider the real facts
			vector<BLVector> normals;
			for (int i = 0; i < c.size(); i++)
			{
				if (c[i] > 0)
				{
					auto conn = meshs.triangle_lists_[active_triangles_left[c[i] - 1].first].point_index_;
					BLVector n1 =
						meshs.virtual_point_lists_[conn[(active_triangles_left[c[i] - 1].second + 1) % 3]].getCoord() -
						meshs.virtual_point_lists_[conn[(active_triangles_left[c[i] - 1].second + 0) % 3]].getCoord();
					BLVector n2 =
						meshs.virtual_point_lists_[conn[(active_triangles_left[c[i] - 1].second + 2) % 3]].getCoord() -
						meshs.virtual_point_lists_[conn[(active_triangles_left[c[i] - 1].second + 1) % 3]].getCoord();
					normals.push_back((n1 ^ n2).normalized());
				}
				else
				{
					auto conn = meshe.triangle_lists_[active_triangles_right[-c[i] - 1].first].point_index_;
					BLVector n1 =
						meshe.virtual_point_lists_[conn[(active_triangles_right[-c[i] - 1].second + 1) % 3]]
						.getCoord() -
						meshe.virtual_point_lists_[conn[(active_triangles_right[-c[i] - 1].second + 0) % 3]].getCoord();
					BLVector n2 =
						meshe.virtual_point_lists_[conn[(active_triangles_right[-c[i] - 1].second + 2) % 3]]
						.getCoord() -
						meshe.virtual_point_lists_[conn[(active_triangles_right[-c[i] - 1].second + 1) % 3]].getCoord();
					normals.push_back((n1 ^ n2).normalized());
				}
			}
			double cost = 0;
			for (int i = 0; i < c.size() - 1; i++)
			{
				cost += abs((normals[i] * normals[i + 1]) - 1);
			}
			if (min_cost > cost)
			{
				min_cost = cost;
				best_combination = c;
			}
		}

		// find api point
		std::queue<int> api_point_left;
		std::stack<int> api_point_right;

		for (int i = 0; i < active_triangles_left.size(); i++)
		{
			auto conn = meshs.triangle_lists_[active_triangles_left[i].first].point_index_;
			api_point_left.push(
				meshs.virtual_point_lists_[conn[(active_triangles_left[i].second + 1) % 3]].getGlobalIndex());
		}
		api_point_left.push(
			meshs
			.virtual_point_lists_
			[meshs.triangle_lists_[active_triangles_left[active_triangles_left.size() - 1].first]
			.point_index_[(active_triangles_left[active_triangles_left.size() - 1].second + 2) % 3]]
			.getGlobalIndex());
		for (int i = 0; i < active_triangles_right.size(); i++)
		{
			auto conn = meshe.triangle_lists_[active_triangles_right[i].first].point_index_;
			api_point_right.push(
				meshe.virtual_point_lists_[conn[(active_triangles_right[i].second + 1) % 3]].getGlobalIndex());
		}
		api_point_right.push(
			meshe
			.virtual_point_lists_
			[meshe.triangle_lists_[active_triangles_right[active_triangles_right.size() - 1].first]
			.point_index_[(active_triangles_right[active_triangles_right.size() - 1].second + 2) % 3]]
			.getGlobalIndex());

		pair<int, int> last_changed{ 0, 0 };
		for (int i = 0; i < best_combination.size(); i++)
		{
			std::array<int, 3> new_tri;
			if (best_combination[i] > 0)
			{
				int index = best_combination[i] - 1;
				while (last_changed.second)
				{
					last_changed.second--;
					api_point_right.pop();
				}
				last_changed.first++;
				int api_point_index = api_point_right.top();
				for (int j = 0; j < 3; j++)
				{
					if (meshs
						.virtual_point_lists_[meshs.triangle_lists_[active_triangles_left[index].first]
						.point_index_[j]]
						.getGlobalIndex() == e)
					{
						new_tri[j] = api_point_index;
					}
					else
					{
						new_tri[j] = meshs
							.virtual_point_lists_[meshs.triangle_lists_[active_triangles_left[index].first]
							.point_index_[j]]
							.getGlobalIndex();
					}
				}
				meshs.triangle_lists_[active_triangles_left[index].first].added_flag = true;
			}
			else
			{
				int index = -best_combination[i] - 1;
				while (last_changed.first)
				{
					last_changed.first--;
					api_point_left.pop();
				}
				last_changed.second++;
				int api_point_index = api_point_left.front();
				for (int j = 0; j < 3; j++)
				{
					if (meshe
						.virtual_point_lists_[meshe.triangle_lists_[active_triangles_right[index].first]
						.point_index_[j]]
						.getGlobalIndex() == s)
					{
						new_tri[j] = api_point_index;
					}
					else
					{
						new_tri[j] =
							meshe
							.virtual_point_lists_[meshe.triangle_lists_[active_triangles_right[index].first]
							.point_index_[j]]
							.getGlobalIndex();
					}
				}
				meshe.triangle_lists_[active_triangles_right[index].first].added_flag = true;
			}

			connector.push_back(new_tri);
			attribute.push_back(new_attribute);
		}
	}

	// add inner new triangles without other complex node
	for (auto i = node_array.begin(); i != node_array.end(); i++)
	{
		auto& local_mesh = i->getFinalMesh();
		if (local_mesh.getExtraPointCount())
		{
			for (auto k : local_mesh.triangle_lists_)
			{
				if (k.added_flag)
					continue;
				std::array<int, 3> new_tri;
				for (int j = 0; j < 3; j++)
				{
					new_tri[j] = local_mesh.virtual_point_lists_[k.point_index_[j]].getGlobalIndex();
				}
				connector.push_back(new_tri);
				attribute.push_back(new_attribute);
			}
		}
	}

	std::unordered_set<int> used_points;
	for (const auto& tri : connector)
	{
		for (int idx : tri)
		{
			used_points.insert(idx);
		}
	}

	std::map<int, std::vector<int>> real_to_virtual;
	for (int i = 0; i < real_node_id_.size(); ++i)
	{
		real_to_virtual[real_node_id_[i]].push_back(i);
	}
	for (const auto& pair : real_to_virtual)
	{
		std::vector<int> used_group;
		for (int idx : pair.second)
		{
			if (used_points.count(idx))
			{
				used_group.push_back(idx);
			}
		}
		if (used_group.size() > 1)
		{
			split_node_groups.push_back(used_group);
		}
	}
	std::cout << "wohaoshua";
}

void MNormalMesh::GenerateFirstLayer(double step_len)
{
	spdlog::info("Try to generate the first layer of mesh");
	// Generate new nodes;

	int node_count = coordinate.size();
	auto extra = decltype(coordinate)();
	std::vector<std::vector<int>> elements;
	// if a triangle connected more than 2 virtual node, we called the triangle as virutal triangles, and all node in vritual triangles should be taken into consideration.
	std::map<std::set<int>,vector<int>> edges;
	std::vector<int> real_tags = real_node_id_;
	int it = 0;
	std::set<int> ordinary_triangles;
	for (auto i : connector) {
		int virtual_count = 0;
		std::set<int> reals;
		for (int j = 0; j < 3; j++)
			reals.insert(real_node_id_[i[j]]);
		if (reals.size() < 3) {
			edges[reals].push_back(it);
		}
		else {
			ordinary_triangles.insert(it);
		}
		it++;
	}
	

	for (auto i : edges) {
		std::vector<int> target_virtual_triangle_list = i.second;
		int real_count = i.first.size();
		std::vector<int> second_upper;
		int first_node = *i.first.begin();
		int second_node = *i.first.rbegin();

		bool reverse = false;
		switch (real_count)
		{
		case 1:
			//degenerated tetra
			for(auto j:i.second)
				elements.push_back(vector<int>{connector[j][0]+node_count, connector[j][1] + node_count, connector[j][2] + node_count, real_node_id_[connector[j][0]]});
			break;
		case 2:
			//degenerated into two tetra
			for (auto j : i.second) {
				second_upper.clear();
				reverse = false;
				int lower_count = 0;// the numer of node's real node in lower of i.first
				for (int k = 0; k < 3; k++)
					if (real_node_id_[connector[j][k]] == first_node) {
						lower_count++;
						if (k == 1) {
							reverse = true;
						}
					}
					else
						second_upper.push_back(connector[j][k]+node_count);
				if (reverse&&second_upper.size()==2)
					std::swap(second_upper[0], second_upper[1]);

				first_node = *i.first.begin();
				second_node = *i.first.rbegin();

				if (lower_count == 1) {	
					elements.push_back(vector<int>{connector[j][0] + node_count, connector[j][2] + node_count, connector[j][1] + node_count, first_node});
					elements.push_back(vector<int>{second_node, first_node, second_upper[0],second_upper[1]});
				}
				else {
					elements.push_back(vector<int>{connector[j][0] + node_count, connector[j][2] + node_count, connector[j][1] + node_count, first_node});
				}
			}
			break;
		default:
			break;
		}
	}
	

	
	for (auto i : ordinary_triangles) {
		int add_extra_point = -1;
		std::set<int> virtual_count;
		for (int j = 0; j < 3; j++) {
			virtual_count.insert(connector[i][j]);
		}
		std::vector<int> node_id;
		for (auto k : virtual_count) {
			node_id.push_back(k);
		}
		for (int j = 0; j < node_id.size(); j++) {
			for (int l = j + 1; l < node_id.size(); l++) {
				std::set<int> edge{node_id[j],node_id[l]};
				if (edges.find(edge) != edges.end()) {
					if (add_extra_point < 0) {
						add_extra_point = node_count * 2 + extra.size();
						BLVector centor;
						for (int k = 0; k < 3; k++) {
							centor =centor+ coordinate[node_id[k]];
						}
						for (int k = 0; k < 3; k++) {
							centor =centor+ coordinate[node_id[k]] + step_len * point_normals[node_id[k]];
						}
						centor = centor / 6;
						extra.push_back(centor);
						
					}
					int first_node = node_id[j];
					int second_node = node_id[l];

					for (int k = 0; k < 3; k++) {
						if (connector[i][k] == first_node && connector[i][(k + 1) % 3] == second_node) {
							swap(first_node, second_node);
							break;
						}
					}
					elements.push_back(std::vector<int>{first_node,second_node,second_node+node_count,add_extra_point});
					elements.push_back(std::vector<int>{first_node,second_node+node_count,first_node+node_count,add_extra_point});
				}
			}
		}
		if (add_extra_point >= 0) {
			for (int j = 0; j < node_id.size(); j++) {
				for (int l = j + 1; l < node_id.size(); l++) {
					std::set<int> edge{ node_id[j],node_id[l] };
					if (edges.find(edge) == edges.end()) {
						int first_node = node_id[j];
						int second_node = node_id[l];
						for (int k = 0; k < 3; k++) {
							if (connector[i][k] == first_node && connector[i][(k + 1) % 3] == second_node) {
								swap(first_node, second_node);
								break;
							}
						}
						elements.push_back(std::vector<int>{first_node, second_node, second_node + node_count,first_node+node_count, add_extra_point});
					}
				}
			}
		}
		else {
			elements.push_back(std::vector<int>{connector[i][0], connector[i][2], connector[i][1], connector[i][0]+node_count, connector[i][2] + node_count, connector[i][1]+node_count});
		}
	}
	


	// In fact, we can classicfy the node into two set by graph-coloring algorithm, only need to take care of "triangle constrain",
	// so 0-black, means decompostion start in 0-layer, 1-white(1-layer), 2-grey(either) for solving triangle conflict.
	std::vector<int> color(real_tags.size(),2);

	// First of all, we need to build a reliable graph between nodes. 
	std::vector<std::vector<int>> graph;
	

	//save vtk file
	string filename = behavior_.getOutputName("_test_vol.vtk");
	int i, j, npt, nelm, sidx;
	FILE *fout = nullptr;

	npt = coordinate.size();
	//nelm = m_nElems;
	nelm = elements.size();

	fout = fopen(filename.data(), "w");
	if (fout == nullptr)
	{
		printf("Can not open file %s\n", filename.data());
	}

	//vtk
	::fprintf(fout, "# vtk DataFile Version 2.0\n");
	::fprintf(fout, "boundary layer mesh\n");
	::fprintf(fout, "ASCII\n");
	::fprintf(fout, "DATASET UNSTRUCTURED_GRID\n");
	::fprintf(fout, "POINTS %d double\n", 2*npt+extra.size());

	for (i = 0; i < npt; i++)
	{
		::fprintf(fout, "%lf %lf %lf\n", coordinate[i][0], coordinate[i][1], coordinate[i][2]);
	}

	for (i = 0; i < npt; i++)
	{
			::fprintf(fout, "%lf %lf %lf\n", coordinate[i][0] + step_len * point_normals[i][0], coordinate[i][1] + step_len * point_normals[i][1], coordinate[i][2] + step_len * point_normals[i][2]);
	}
	for (auto i : extra) {
	
			::fprintf(fout, "%lf %lf %lf\n", i[0],i[1],i[2]);
	}
	

	//sidx = m_nSurfElems;
	sidx = 0;

	//sidx = 2914;
	//nelm = sidx + 100;
	int number_count = 0;
	for (auto i : elements) {
		number_count += i.size() + 1;
	}
	::fprintf(fout, "CELLS %d %d\n", nelm, number_count);
	for (i = 0; i < nelm - sidx; i++)
		//for (i=sidx; i<nelm; i++)
	{
		if(elements[i].size()==4)
			::fprintf(fout, "4 %d %d %d %d\n", elements[i][0],elements[i][1],elements[i][2],elements[i][3]);
		if (elements[i].size() == 5)
			::fprintf(fout, "5 %d %d %d %d %d\n", elements[i][0], elements[i][1], elements[i][2], elements[i][3],elements[i][4]);
		if (elements[i].size() == 6)
			::fprintf(fout, "6 %d %d %d %d %d %d\n", elements[i][0], elements[i][1], elements[i][2], elements[i][3], elements[i][4], elements[i][5]);
	}

	::fprintf(fout, "CELL_TYPES %d\n", nelm);

	for (i = 0; i < nelm; i++) {
		if (elements[i].size() == 4)
			::fprintf(fout, "10\n");
		if (elements[i].size() == 5)
			::fprintf(fout, "14\n");
		if (elements[i].size() == 6)
			::fprintf(fout, "13\n");
	}

	::fclose(fout);
	fout = nullptr;
	spdlog::info("Writing {0}", filename);
}



void MNormalMesh::CaculateFrontNormal()
{
	for (auto &node : node_array) {
		for (int i = 0; i < node.neighbour_node_.size(); i++) {
			int j = (i + 1) % node.neighbour_node_.size();
			BLVector e1 = node.neighbour_node_[j]->coordinate - node.coordinate;
			BLVector e2 = node.neighbour_node_[i]->coordinate - node.neighbour_node_[j]->coordinate;
			BLVector normal = (e2^e1).normalized();
			node.neighbour_front_direction_.push_back(normal);
			
		}
		vector<int> ans;
		ans.swap(node.neighbour_front_index_);
		for (int i = 0; i < node.neighbour_node_.size(); i++) {
			int j = (i + 1) % node.neighbour_node_.size();
			std::array<int, 3> neighbour{ node.neighbour_node_[i]->node_id_, node.neighbour_node_[j]->node_id_,node.node_id_};
			std::sort(neighbour.begin(), neighbour.end());
			for (auto k : ans) {
				std::array<int, 3> my_order{connector[k][0],connector[k][1] ,connector[k][2] };
				std::sort(my_order.begin(), my_order.end());
				if (my_order == neighbour)
					node.neighbour_front_index_.push_back(k);
			}
		}
	}
}

void MNormalMesh::CalculateNodeNormal()
{
	for (auto &node : node_array) {
		if (node.neighbour_front_direction_.empty())
			throw std::runtime_error("please fill the front normal first");
		node.single_normal_ = GEEOMETRY_FUNCTION::getMostNormal(node.neighbour_front_direction_).normalized();
		node.visible_angle_ = node.CaculateVisableAngle(node.single_normal_);
		node.original_skewness_ = (acos(node.visible_angle_) )  / (PI/2);
		if (node.visible_angle_ < 0)
			node.original_skewness_ = 1;

		point_normals[node.node_id_] = node.single_normal_;
	}
}
