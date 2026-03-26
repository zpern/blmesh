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
#include <unordered_set>
#include <omp.h>       // OpenMP 并行支持
using std::stringstream;
using std::ifstream;
using std::ofstream;
using std::map;
using std::array;
//#define check_intersection
#define check_intersection2

void MNormalMesh::ReadPlsBuf(std::string f,std::vector<std::array<double, 3>>& points) {
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
            node_array[i].coordinate.x = points[i][0];
            node_array[i].coordinate.y = points[i][1];
            node_array[i].coordinate.z = points[i][2];
            for (int k = 0; k < 3; k++) {
                coordinate[i][k] = node_array[i].coordinate[k];
            }
            node_array[i].node_id_ = i;
            real_node_id_[i] = i; // of course they are themself
        }

        connector.resize(number_of_element);
        attribute.resize(number_of_element);

        for (int i = 0; i < number_of_element; i++) {
          fin >> line >> connector[i][0] >> connector[i][1] >> connector[i][2] >> attribute[i];
          for (int k = 0; k < 3; k++) graph[connector[i][k]].push_back(connector[i]);
          for (int k = 0; k < 3; k++) node_array[connector[i][k]].neighbour_front_index_.push_back(i);
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

            // 寻找起始点
            unordered_map<int, int> count_map;
            int start = graph[i][0][0];
            int pos = start;

            // 统计 j[0] 出现的次数
             for (const auto &j : graph[i]) {
                 count_map[j[0]]++;
             }

             // 统计j[1]与j[0]出现的次数
             for (const auto& j : graph[i]) {
                 if (count_map.find(j[1]) != count_map.end()) {
                     count_map[j[1]]++;
                 }
             }

             for (const auto& entry : count_map) {
                 if (entry.second == 1) {
                     avoid_spliteNode.push_back(i);
                     start = entry.first;
                     pos = start;
                     break;
                 }
             }

             while (true) {
                 bool reach_end = true;
                 for (auto& j : graph[i]) {
                     if (j[0] == pos) {
                         pos = j[1];
                         node_array[i].neighbour_node_.push_back(node_array.begin() + j[1]);
                         reach_end = false;
                         break;
                     }
                 }
                 if (pos == start || reach_end) {
                     break;
                 }
             }
        }
        CaculateFrontNormal();
        CalculateNodeNormal();

}

void MNormalMesh::CaculateFrontNormal()
{
    for (auto& node : node_array) {
        for (int i = 0; i < node.neighbour_node_.size(); i++) {
            int j = (i + 1) % node.neighbour_node_.size();
            BLVector e1 = node.neighbour_node_[j]->coordinate - node.coordinate;
            BLVector e2 = node.neighbour_node_[i]->coordinate - node.neighbour_node_[j]->coordinate;
            BLVector normal = (e2 ^ e1).normalized();
            node.neighbour_front_direction_.push_back(normal);

        }
        vector<int> ans;
        ans.swap(node.neighbour_front_index_);
        for (int i = 0; i < node.neighbour_node_.size(); i++) {
            int j = (i + 1) % node.neighbour_node_.size();
            std::array<int, 3> neighbour{ node.neighbour_node_[i]->node_id_, node.neighbour_node_[j]->node_id_,node.node_id_ };
            std::sort(neighbour.begin(), neighbour.end());
            for (auto k : ans) {
                std::array<int, 3> my_order{ connector[k][0],connector[k][1] ,connector[k][2] };
                std::sort(my_order.begin(), my_order.end());
                if (my_order == neighbour)
                    node.neighbour_front_index_.push_back(k);
            }
        }
    }
}
void MNormalMesh::CalculateNodeNormal()
{
    for (auto& node : node_array) {
        if (node.neighbour_front_direction_.empty())
            throw std::runtime_error("please fill the front normal first");
        node.single_normal_ = GEEOMETRY_FUNCTION::getMostNormal(node.neighbour_front_direction_).normalized();
        node.visible_angle_ = node.CaculateVisableAngle(node.single_normal_);
        node.original_skewness_ = (acos(node.visible_angle_)) / (PI / 2);
        if (node.visible_angle_ < 0)
            node.original_skewness_ = 1;

        point_normals[node.node_id_] = node.single_normal_;
    }
}
void MNormalMesh::FixedLength()
{
    for (int i = 0; i < node_array.size(); i++) {
        auto& node = node_array[i];
        if (node.neighbour_node_.empty()) {
            continue;
        }
        double length = 0;
        for (auto neighbour_node : node.neighbour_node_) {
            length += (neighbour_node->coordinate - node.coordinate).magnitude();
        }
        length /= node.neighbour_node_.size();

        std::array<double, 3> key = { node.coordinate.x, node.coordinate.y, node.coordinate.z };

        auto it = point_to_length.find(key);
        if (it != point_to_length.end() && length * 0.8 < it->second) {
            node.highRatio = length * 0.3 / it->second;
        }
    }
    using NodeIter = decltype(node_array.begin());
    for (int i = 0; i < node_array.size(); i++) {
        std::queue<NodeIter> q;
        q.push(node_array.begin() + i);

        while (!q.empty()) {
            NodeIter cur_it = q.front();
            q.pop();

            auto &cur = *cur_it;

            for (auto nb_it : cur.neighbour_node_) {
                auto &nb = *nb_it;

                if (nb.highRatio > cur.highRatio + 0.1) {
                    nb.highRatio = cur.highRatio + 0.1;
                    q.push(nb_it);
                }
            }
        }
    }
    for (int i = 0; i < node_array.size(); i++){
        auto node = node_array[i];
        std::array<double, 3> key = { node.coordinate.x, node.coordinate.y, node.coordinate.z };
        point_to_length[key] *= node.highRatio; 
    }
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
        // 如果节点在 avoid_spliteNode 中，则跳过该节点
        if (std::find(avoid_spliteNode.begin(), avoid_spliteNode.end(), i->node_id_) !=
            avoid_spliteNode.end()) {
            continue; // 跳过当前节点
        }
        if (exist_prism) {
            if (i->original_skewness_ > 1) { // naca 0.14459
                need_to_split.push_back(i);
            }
        } 
		else {
            if (i->original_skewness_ > SKEWNESS_THREADHOLD) { // naca 0.14459
                need_to_split.push_back(i);
            }
        }
           
	}
	//spdlog::info("//////////////////////////////////");
	spdlog::info("The maximum skewness of single normal method = {0}", d);
	//spdlog::info("//////////////////////////////////");

	int number_of_complex_nodes = 0;
	for (auto i : need_to_split) {
		number_of_complex_nodes++;
		i->SplitNode();	
	}

	clock_t end_time = clock();
	spdlog::info("");
	spdlog::info("");
	spdlog::info("/////////////////////////////////");


	spdlog::info("The number of Complex node={0}.", number_of_complex_nodes);
	spdlog::info("The time comsuption of splitting={:03.2f} sec.", (end_time - start_time)*1.0/CLOCKS_PER_SEC);
}
void MNormalMesh::SmoothNormalsSimple(int iterations)
{
    std::vector<BLVector> new_normals(coordinate.size());

    for (int it = 0; it < iterations; ++it)
    {
        for (size_t idx = 0; idx < coordinate.size(); ++idx)
        {
            BLVector avg_normal = point_normals[idx]; // 当前点法向
            int count = 1;

            // 遍历连接到当前点的所有三角形
            for (auto& tri : connector)
            {
                for (int j = 0; j < 3; ++j)
                {
                    if (tri[j] == idx)
                    {
                        // 累加三角形其他顶点的法向
                        for (int k = 0; k < 3; ++k)
                        {
                            if (k != j)
                            {
                                avg_normal += 0.1*point_normals[tri[k]];
                                count++;
                            }
                        }
                    }
                }
            }

            avg_normal = avg_normal / count;
            avg_normal.normalize();
            new_normals[idx] = avg_normal;
        }

        // 更新法向
        point_normals = new_normals;
    }
}

void MNormalMesh::BuildTopo(int faceCount)
{
    int new_attribute = faceCount + 1;
    std::map<int, std::map<int, std::vector<int>>> complex_edges;

    // ==============================
    // Step 1: 初始化所有局部节点
    // ==============================
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
                        {
                            if (connector[j][l] == i->node_id_)
                            {
                                connector[j][l] = point_index_map[count];
                            }
                        }
                    }
                }
            }
        }
    }

    // ==============================
    // Step 2: mesh stitching
    // ==============================
    struct pos
    {
        int original_index;
        BLVector left;
        BLVector right;
    };

    std::map<std::array<int, 2>, pos> complex_edges_pair;
    std::map<std::array<int, 2>, int> global_far_node_map; // store global point idx for each virtual edge

    for (auto i = node_array.begin(); i != node_array.end(); i++) {
        auto &local_mesh = i->getFinalMesh();

        if (local_mesh.getExtraPointCount()) {
            for (int k1 = 0; k1 < local_mesh.boundary_edges_.size(); k1++) {
                int k = local_mesh.boundary_edges_[k1][0];

                if (local_mesh.valid.find(k) == local_mesh.valid.end() ||
                    !local_mesh.virtual_point_lists_[k].isFarNode()) { // far point
                    continue;
                }

                int index = local_mesh.virtual_point_lists_[k].getGlobalIndex();

                // --- 统一顺序，小节点在前，大节点在后 ---
                int u = std::min(i->node_id_, index);
                int v = std::max(i->node_id_, index);

                if (node_array[index].getFinalMesh().getExtraPointCount()) {
                    complex_edges_pair[{u, v}] = pos();  // 复杂边
                } else {
                    global_far_node_map[{u, v}] = index; // 普通边
                }
            }
        }
    }


    // ==============================
    // Step 3: 枚举所有可能连接模式
    // ==============================
    static std::vector<std::vector<int>> connection[2][2] = { std::vector<std::vector<int>>() };
    static bool create = false;

    if (!create)
    {
        create = true;

        connection[0][0].push_back({ -1, 1 });
        connection[0][0].push_back({ 1, -1 });

        connection[0][1].push_back({ 1, -2, -1 });
        connection[0][1].push_back({ -2, 1, -1 });
        connection[0][1].push_back({ -2, -1, 1 });

        connection[1][0].push_back({ 1, 2, -1 });
        connection[1][0].push_back({ 1, -1, 2 });
        connection[1][0].push_back({ -1, 1, 2 });

        connection[1][1].push_back({ -2, -1, 1, 2 });
        connection[1][1].push_back({ 1, -2, -1, 2 });
        connection[1][1].push_back({ -2, 1, -1, 2 });
        connection[1][1].push_back({ 1, 2, -2, -1 });
        connection[1][1].push_back({ 1, -2, 2, -1 });
        connection[1][1].push_back({ -2, 1, 2, -1 });
    }

    // ==============================
    // Step 4: 遍历所有复杂边对
    // ==============================
    for (auto edge : complex_edges_pair) {
        int s = edge.first[0];
        int e = edge.first[1];
        auto &meshs = node_array[s].getFinalMesh();
        auto &meshe = node_array[e].getFinalMesh();

        if (s == 18555) {
            std::cout << "xy";
        }
        // --- Collect left active triangles ---
        std::vector<std::pair<int, int>> active_triangles_left;
        for (int i = 0; i < meshs.triangle_lists_.size(); i++) {
            for (int j = 0; j < 3; j++) {
                if (meshs.virtual_point_lists_[meshs.triangle_lists_[i].point_index_[j]]
                        .getGlobalIndex() == e) {
                    active_triangles_left.push_back({i, j});
                }
            }
        }

        // --- Collect right active triangles ---
        std::vector<std::pair<int, int>> active_triangles_right;
        for (int i = 0; i < meshe.triangle_lists_.size(); i++) {
            for (int j = 0; j < 3; j++) {
                if (meshe.virtual_point_lists_[meshe.triangle_lists_[i].point_index_[j]]
                        .getGlobalIndex() == s) {
                    active_triangles_right.push_back({i, j});
                }
            }
        }

        // --- Sort triangles if 2 on each side ---
        if (active_triangles_left.size() == 2) {
            if (meshs.triangle_lists_[active_triangles_left[0].first]
                    .point_index_[(active_triangles_left[0].second + 1) % 3] ==
                meshs.triangle_lists_[active_triangles_left[1].first]
                    .point_index_[(active_triangles_left[1].second + 2) % 3]) {
                std::swap(active_triangles_left[0], active_triangles_left[1]);
            }
        }

        if (active_triangles_right.size() == 2) {
            if (meshe.triangle_lists_[active_triangles_right[0].first]
                    .point_index_[(active_triangles_right[0].second + 1) % 3] ==
                meshe.triangle_lists_[active_triangles_right[1].first]
                    .point_index_[(active_triangles_right[1].second + 2) % 3]) {
                std::swap(active_triangles_right[0], active_triangles_right[1]);
            }
        }

        if ((active_triangles_left.size() > 2) || (active_triangles_right.size() > 2)) {
            throw std::runtime_error("wake up and write code for more than 2!");
        }

        if (active_triangles_left.size() == 0 || active_triangles_right.size() == 0) {
            struct BLVectorHash {
                std::size_t operator()(const BLVector &v) const
                {
                    auto h1 = std::hash<double>{}(v.x);
                    auto h2 = std::hash<double>{}(v.y);
                    auto h3 = std::hash<double>{}(v.z);

                    std::size_t seed = h1;
                    seed ^= h2 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
                    seed ^= h3 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
                    return seed;
                }
            };
            auto collectRepeated = [&](VirtualSphereMesh &mesh) {
                std::unordered_map<BLVector, std::unordered_set<int>, BLVectorHash> coord2global;

                for (const auto &tri : mesh.triangle_lists_) {
                    for (int v : tri.point_index_) {
                        int globalIdx = mesh.virtual_point_lists_[v].getGlobalIndex();
                        const BLVector &point = coordinate[globalIdx];

                        coord2global[point].insert(globalIdx);
                    }
                }

                std::unordered_set<int> repeated_set;
                for (const auto &[point, ids] : coord2global) {
                    if (ids.size() > 1) {
                        repeated_set.insert(ids.begin(), ids.end());
                    }
                }

                return std::vector<int>(repeated_set.begin(), repeated_set.end());
            };

            std::vector<int> spoints = collectRepeated(meshs);
            std::vector<int> epoints = collectRepeated(meshe);

            std::vector<std::pair<int, int>> results;
            for (const auto &arr : connector) {
                int foundS = -1;
                int foundE = -1;

                // 查找 s 中的任意一个值
                for (int x : arr) {
                    if (std::find(spoints.begin(), spoints.end(), x) != spoints.end()) {
                        foundS = x;
                        break;
                    }
                }

                // 查找 e 中的任意一个值
                for (int x : arr) {
                    if (std::find(epoints.begin(), epoints.end(), x) != epoints.end()) {
                        foundE = x;
                        break;
                    }
                }

                if (foundS != -1 && foundE != -1) {
                    results.emplace_back(foundS, foundE);
                }
            }


            if (active_triangles_right.empty()) {
                for (auto &tri_info : active_triangles_left) {
                    std::array<int, 3> new_tri;
                    auto &tri = meshs.triangle_lists_[tri_info.first];
                    for (int j = 0; j < 3; j++) {
                        if (meshs.virtual_point_lists_[tri.point_index_[j]].getGlobalIndex() == e) {
                            new_tri[j] = results[0].second;
                        } else {
                            new_tri[j] =
                                meshs.virtual_point_lists_[tri.point_index_[j]].getGlobalIndex();
                        }
                    }
                    tri.added_flag = true;
                    connector.push_back(new_tri);
                    attribute.push_back(new_attribute);
                }
            } 
            else if (active_triangles_left.empty()) {

                for (auto &tri_info : active_triangles_right) {
                    std::array<int, 3> new_tri;
                    auto &tri = meshe.triangle_lists_[tri_info.first];
                    for (int j = 0; j < 3; j++) {
                        if (meshe.virtual_point_lists_[tri.point_index_[j]].getGlobalIndex() == s) {
                            new_tri[j] = results[0].first;
                        } else {
                            new_tri[j] =
                                meshe.virtual_point_lists_[tri.point_index_[j]].getGlobalIndex();
                        }
                    }
                    tri.added_flag = true;
                    connector.push_back(new_tri);
                    attribute.push_back(new_attribute);
                }
            }
        } 
        else {

            // --- Determine combination ---
            auto combination =
                connection[active_triangles_left.size() - 1][active_triangles_right.size() - 1];
            std::vector<int> best_combination;
            double min_cost = 100;

            // ==============================
            // Step 5: Evaluate best combination
            // ==============================
            for (auto c : combination) {
                std::vector<BLVector> normals;

                for (int i = 0; i < c.size(); i++) {
                    if (c[i] > 0) {
                        auto conn = meshs.triangle_lists_[active_triangles_left[c[i] - 1].first]
                                        .point_index_;
                        BLVector n1 =
                            meshs
                                .virtual_point_lists_
                                    [conn[(active_triangles_left[c[i] - 1].second + 1) % 3]]
                                .getCoord() -
                            meshs
                                .virtual_point_lists_
                                    [conn[(active_triangles_left[c[i] - 1].second + 0) % 3]]
                                .getCoord();
                        BLVector n2 =
                            meshs
                                .virtual_point_lists_
                                    [conn[(active_triangles_left[c[i] - 1].second + 2) % 3]]
                                .getCoord() -
                            meshs
                                .virtual_point_lists_
                                    [conn[(active_triangles_left[c[i] - 1].second + 1) % 3]]
                                .getCoord();
                        normals.push_back((n1 ^ n2).normalized());
                    } else {
                        auto conn = meshe.triangle_lists_[active_triangles_right[-c[i] - 1].first]
                                        .point_index_;
                        BLVector n1 =
                            meshe
                                .virtual_point_lists_
                                    [conn[(active_triangles_right[-c[i] - 1].second + 1) % 3]]
                                .getCoord() -
                            meshe
                                .virtual_point_lists_
                                    [conn[(active_triangles_right[-c[i] - 1].second + 0) % 3]]
                                .getCoord();
                        BLVector n2 =
                            meshe
                                .virtual_point_lists_
                                    [conn[(active_triangles_right[-c[i] - 1].second + 2) % 3]]
                                .getCoord() -
                            meshe
                                .virtual_point_lists_
                                    [conn[(active_triangles_right[-c[i] - 1].second + 1) % 3]]
                                .getCoord();
                        normals.push_back((n1 ^ n2).normalized());
                    }
                }

                double cost = 0;
                for (int i = 0; i < c.size() - 1; i++) {
                    cost += std::abs((normals[i] * normals[i + 1]) - 1);
                }

                if (min_cost > cost) {
                    min_cost = cost;
                    best_combination = c;
                }
            }

            // ==============================
            // Step 6: Build new triangles
            // ==============================
            std::queue<int> api_point_left;
            std::stack<int> api_point_right;

            for (int i = 0; i < active_triangles_left.size(); i++) {
                auto conn = meshs.triangle_lists_[active_triangles_left[i].first].point_index_;
                api_point_left.push(
                    meshs.virtual_point_lists_[conn[(active_triangles_left[i].second + 1) % 3]]
                        .getGlobalIndex());
            }

            api_point_left.push(
                meshs
                    .virtual_point_lists_
                        [meshs.triangle_lists_[active_triangles_left.back().first]
                             .point_index_[(active_triangles_left.back().second + 2) % 3]]
                    .getGlobalIndex());

            for (int i = 0; i < active_triangles_right.size(); i++) {
                auto conn = meshe.triangle_lists_[active_triangles_right[i].first].point_index_;
                api_point_right.push(
                    meshe.virtual_point_lists_[conn[(active_triangles_right[i].second + 1) % 3]]
                        .getGlobalIndex());
            }

            api_point_right.push(
                meshe
                    .virtual_point_lists_
                        [meshe.triangle_lists_[active_triangles_right.back().first]
                             .point_index_[(active_triangles_right.back().second + 2) % 3]]
                    .getGlobalIndex());

            std::pair<int, int> last_changed{0, 0};

            for (int i = 0; i < best_combination.size(); i++) {
                std::array<int, 3> new_tri;

                if (best_combination[i] > 0) {
                    int index = best_combination[i] - 1;

                    while (last_changed.second) {
                        last_changed.second--;
                        api_point_right.pop();
                    }

                    last_changed.first++;
                    int api_point_index = api_point_right.top();

                    for (int j = 0; j < 3; j++) {
                        if (meshs
                                .virtual_point_lists_
                                    [meshs.triangle_lists_[active_triangles_left[index].first]
                                         .point_index_[j]]
                                .getGlobalIndex() == e) {
                            new_tri[j] = api_point_index;
                        } else {
                            new_tri[j] =
                                meshs
                                    .virtual_point_lists_
                                        [meshs.triangle_lists_[active_triangles_left[index].first]
                                             .point_index_[j]]
                                    .getGlobalIndex();
                        }
                    }

                    meshs.triangle_lists_[active_triangles_left[index].first].added_flag = true;
                } else {
                    int index = -best_combination[i] - 1;

                    while (last_changed.first) {
                        last_changed.first--;
                        api_point_left.pop();
                    }

                    last_changed.second++;
                    int api_point_index = api_point_left.front();

                    for (int j = 0; j < 3; j++) {
                        if (meshe
                                .virtual_point_lists_
                                    [meshe.triangle_lists_[active_triangles_right[index].first]
                                         .point_index_[j]]
                                .getGlobalIndex() == s) {
                            new_tri[j] = api_point_index;
                        } else {
                            new_tri[j] =
                                meshe
                                    .virtual_point_lists_
                                        [meshe.triangle_lists_[active_triangles_right[index].first]
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
    }
    // ==============================
    // Step 7: Add remaining inner triangles
    // ==============================
    for (auto i = node_array.begin(); i != node_array.end(); i++)
    {
        auto& local_mesh = i->getFinalMesh();

        if (local_mesh.getExtraPointCount())
        {
            for (auto& k : local_mesh.triangle_lists_)
            {
                if (k.added_flag)
                    continue;

                std::array<int, 3> new_tri;
                for (int j = 0; j < 3; j++) {
                    new_tri[j] = local_mesh.virtual_point_lists_[k.point_index_[j]].getGlobalIndex();
                }
                connector.push_back(new_tri);
                attribute.push_back(new_attribute);
            }
        }
    }
}


void MNormalMesh::pre_WriteVol(std::vector<std::array<double, 3>> &v,std::vector<std::vector<int>> &f,int &lower_num,int &add_point_num)
{

    std::map<std::array<double, 3>, int> coord_to_id; 
    std::set<int> duplicate_lower_ids;
    std::vector<std::array<int, 3>> lower_ids(connector.size());

    // 建立编号映射
    int current_id = 0;
    for (int i = 0; i < connector.size(); ++i) {
        std::array<int, 3> &ids = lower_ids[i];
        for (int k = 0; k < 3; ++k) {
            int coord_idx = connector[i][k];
            std::array<double, 3> pt = {coordinate[coord_idx].x,
                                       coordinate[coord_idx].y,
                                       coordinate[coord_idx].z};

            // 建立 coord_to_id 映射
            auto it = coord_to_id.find(pt);
            if (it == coord_to_id.end()) {
                coord_to_id[pt] = current_id;
                ids[k] = current_id++;
            } else {
                ids[k] = it->second;
            }
        }
    }

    for (int i = 0; i < lower_ids.size(); ++i) {
        const auto &ids = lower_ids[i];
        const auto &conn = connector[i];  // 原始坐标索引

        if (ids[0] == ids[1]) duplicate_lower_ids.insert(conn[0]);
        if (ids[0] == ids[2]) duplicate_lower_ids.insert(conn[0]);
        if (ids[1] == ids[2]) duplicate_lower_ids.insert(conn[1]);
    }


    lower_num = coord_to_id.size();

    auto idx = [&](int base, int layer = 0) {
        if (layer == -1) {
            return static_cast<int>(base);
        }
        return static_cast<int>(base + lower_num + layer * coordinate.size());
    };
    auto has_duplicate = [&](int i1, int i2, int i3, int i4) {
        std::set<std::array<double, 3>> s = {v[i1], v[i2], v[i3], v[i4]};
        return s.size() == 4; // 如果 size < 4，说明有重复
    };


    add_point_num = (number_of_layer-1) * coordinate.size();
    v.resize(lower_num + number_of_layer * coordinate.size());

    for (auto i : coord_to_id) {
        v[i.second] = i.first;
    }
    length.assign(coordinate.size(),0);
    for (int i = 0; i < coordinate.size(); i++) {
        length[i] = point_to_length[{coordinate[i].x, coordinate[i].y, coordinate[i].z}];
    }

        // 抽出去的相交检测 + 步长回缩
    if (!IntersectionCheck(lower_ids, v, lower_num, length)) {
        multiplySuccess = false;
        std::cout << " multiply intersection" << std::endl;
        return;
    }
    std::cout << "finish intersection" << std::endl;

    for (int i = 0; i < coordinate.size(); i++) {
        if (duplicate_lower_ids.find(i) == duplicate_lower_ids.end()) {
            length[i] = 0;
        }
    }
    for (int j = 1; j <= number_of_layer; j++) {
        for (int i = 0; i < coordinate.size(); i++) {
            v[idx(i, j - 1)] = {coordinate[i].x + j * length[i] * point_normals[i].x,
                                coordinate[i].y + j * length[i] * point_normals[i].y,
                                coordinate[i].z + j * length[i] * point_normals[i].z};
        }
    }

    // 切割方案
    // 第一层
    for (int i = 0; i < connector.size(); i++) {
        // 三法向分裂情况
        if (lower_ids[i][0] == lower_ids[i][1] && lower_ids[i][2] == lower_ids[i][1]) {
            f.push_back({lower_ids[i][0],
                            idx(connector[i][0]),
                            idx(connector[i][1]),
                            idx(connector[i][2])});
        }
        // 无法向分裂情况
        else if (lower_ids[i][0] != lower_ids[i][1] && lower_ids[i][2] != lower_ids[i][1] &&
                    lower_ids[i][2] != lower_ids[i][0]) {

            int k1, k2, k3;
            for (int j = 0; j < 3; j++) {
                if (lower_ids[i][j] < lower_ids[i][(j + 1) % 3] &&
                    lower_ids[i][j] < lower_ids[i][(j + 2) % 3]) {
                    k1 = j;
                    k2 = (lower_ids[i][(j + 2) % 3] < lower_ids[i][(j + 1) % 3]) ? (j + 2) % 3
                                                                                    : (j + 1) % 3;
                    k3 = (lower_ids[i][(j + 2) % 3] > lower_ids[i][(j + 1) % 3]) ? (j + 2) % 3
                                                                                    : (j + 1) % 3;
                    break;
                }
            }

            int count = 0;
            if (duplicate_lower_ids.find(connector[i][k1]) == duplicate_lower_ids.end()) count++;
            if (duplicate_lower_ids.find(connector[i][k2]) == duplicate_lower_ids.end()) count++;
            if (duplicate_lower_ids.find(connector[i][k3]) == duplicate_lower_ids.end()) count++;
            if (count == 3) {
                continue;
            } else {
                if(has_duplicate(lower_ids[i][k1],lower_ids[i][k2],idx(connector[i][k2]),idx(connector[i][k3]))) 
                    f.push_back({lower_ids[i][k1],lower_ids[i][k2],idx(connector[i][k2]),idx(connector[i][k3])});

                if(has_duplicate(lower_ids[i][k1],idx(connector[i][k1]),idx(connector[i][k2]),idx(connector[i][k3])))
                    f.push_back({lower_ids[i][k1],idx(connector[i][k1]),idx(connector[i][k2]),idx(connector[i][k3])});

                if(has_duplicate(lower_ids[i][k1], lower_ids[i][k2], lower_ids[i][k3], idx(connector[i][k3])))
                    f.push_back({lower_ids[i][k1], lower_ids[i][k2], lower_ids[i][k3], idx(connector[i][k3])});
            }
        }
        // 法向二分情况
        else {
            int k1, k2, k3;
            for (int j = 0; j < 3; j++) {
                if (lower_ids[i][j] == lower_ids[i][(j + 1) % 3]) {
                    k1 = j;
                    k2 = (j + 1) % 3;
                    k3 = (j + 2) % 3;
                    break;
                }
            }

            if (lower_ids[i][k1] < lower_ids[i][k3]) {
                f.push_back({lower_ids[i][k1],
                                idx(connector[i][k1]),
                                idx(connector[i][k2]),
                                idx(connector[i][k3])});
            } else {
                f.push_back({lower_ids[i][k3],
                                idx(connector[i][k1]),
                                idx(connector[i][k2]),
                                idx(connector[i][k3])});

                f.push_back({lower_ids[i][k1],
                                idx(connector[i][k1]),
                                idx(connector[i][k2]),
                                lower_ids[i][k3]});
            }
        }
    }
    return;
    }
//void MNormalMesh::WriteVol(std::vector<std::array<double, 3>> &v,std::vector<std::vector<int>> &f,int &lower_num,int &add_point_num)
//{
//    std::map<std::array<double, 3>, int> coord_to_id;
//    std::vector<std::array<int, 3>> lower_ids(connector.size());
//    for (int i = 0; i < connector.size(); i++) {
//        int count = 0;
//        int lowerid_1 = -1, lower_id2 = -1;
//        std::array<int, 3> &ids = lower_ids[i];
//        for (int k = 0; k < 3; k++) {
//            std::array<double, 3> v;
//            v[0] = coordinate[connector[i][k]].x;
//            v[1] = coordinate[connector[i][k]].y;
//            v[2] = coordinate[connector[i][k]].z;
//            ids[k] = -1;
//            if (coord_to_id.find(v) == coord_to_id.end()) {
//                ids[k] = coord_to_id.size();
//                coord_to_id[v] = ids[k];
//            } else {
//                ids[k] = coord_to_id[v];
//            }
//        }
//    }
//    lower_num = coord_to_id.size();
//
//    auto idx = [&](int base, int layer = 0) {
//        if (layer == -1) {
//            return static_cast<int>(base);
//        }
//        return static_cast<int>(base + lower_num + layer * coordinate.size());
//    };
//
//    add_point_num = (number_of_layer-1) * coordinate.size();
//    v.resize(lower_num + number_of_layer * coordinate.size());
//    for (auto i : coord_to_id) {
//        v[i.second] = i.first;
//    }
//
//    length.assign(coordinate.size(), step_of_length);
//    if (point_to_length.size()) {
//        for (int i = 0; i < coordinate.size(); i++) {
//            length[i] = point_to_length[{coordinate[i].x, coordinate[i].y, coordinate[i].z}];
//        }
//    }
//
//     //Intersection
//    if (true) {
//            std::set<int> record_point;
//            int iter_count = 0;      // 加循环计数
//            const int MAX_ITER = 20; // 最大迭代次数
//            bool zero_step_retry_done = false;
//            do {
//                IntersecChecker checker_;
//                BoundingBox box({std::numeric_limits<double>::max(),
//                                 std::numeric_limits<double>::max(),
//                                 std::numeric_limits<double>::max(),
//                                 std::numeric_limits<double>::lowest(),
//                                 std::numeric_limits<double>::lowest(),
//                                 std::numeric_limits<double>::lowest()});
//                for (int i = 0; i < coordinate.size(); i++) {
//                    box[0] = std::min(box[0], coordinate[i].x);
//                    box[1] = std::min(box[1], coordinate[i].y);
//                    box[2] = std::min(box[2], coordinate[i].z);
//                    box[3] = std::max(box[3], coordinate[i].x);
//                    box[4] = std::max(box[4], coordinate[i].y);
//                    box[5] = std::max(box[5], coordinate[i].z);
//                }
//                checker_.init(box);
//
//                // inner
//                int first_id = checker_.addPoint([&]() {
//                    std::vector<BLVector> tmp;
//                    tmp.reserve(lower_num);
//                    for (size_t i = 0; i < lower_num; ++i) {
//                        tmp.emplace_back(v[i][0], v[i][1], v[i][2]);
//                    }
//                    return tmp;
//                }());
//
//                std::vector<std::pair<HexaTag, std::vector<int>>> elements;
//                for (size_t i = 0; i < connector.size(); ++i) {
//                    if (lower_ids[i][0] != lower_ids[i][1] && lower_ids[i][1] != lower_ids[i][2] &&
//                        lower_ids[i][0] != lower_ids[i][2]) {
//                        elements.emplace_back(
//                            HexaTag(i, 0, TRI_BOTTOM),
//                            std::vector<int>{lower_ids[i][0], lower_ids[i][1], lower_ids[i][2]});
//                    }
//                }
//                checker_.addElements(elements);
//
//                // side and top
//                if (iter_count < MAX_ITER) {
//                    // 常规回退：相交点步长缩小
//                    for (auto i : record_point) {
//                        length[i] *= 0.8;
//                    }
//                    iter_count++;
//                    record_point.clear();
//                    // continue;  // 让外层重新跑一遍检测（如果你是 while/for 结构）
//                } 
//                else {
//                    // 已经到最大次数：先做一次“步长置0重试”
//                    if (!zero_step_retry_done) {
//#ifdef _DEBUG
//                        std::cout << "first" << std::endl;
//                        for (auto p : record_point) {
//                            spdlog::error("point id = {}, coord = ({:.16f}, {:.16f}, {:.16f}), "
//                                          "length = {:.16f}",
//                                          p,
//                                          coordinate[p].x,
//                                          coordinate[p].y,
//                                          coordinate[p].z,
//                                          length[p]);
//                        }
//#endif // DEBUG
//                        for (auto i : record_point) {
//                            length[i] *= 0.0;
//                        }
//                        zero_step_retry_done = true;
//                        record_point.clear();
//                        break;
//                    } 
//                    else {
//#ifdef _DEBUG
//                        std::cout << "second" << std::endl;
//                        for (auto p : record_point) {
//                            spdlog::error("point id = {}, coord = ({:.16f}, {:.16f}, {:.16f}), "
//                                          "length = {:.16f}",
//                                          p,
//                                          coordinate[p].x,
//                                          coordinate[p].y,
//                                          coordinate[p].z,
//                                          length[p]);
//                        }
//#endif // DEBUG
//       // 置0重试后仍相交：退出
//                        spdlog::info("Temporarily revert to using a single normal.");
//                        multiplySuccess = false;
//                        return;
//                    }
//                    break;
//                }
//
//                std::vector<BLVector> grown_coordinate;
//                grown_coordinate.resize(length.size());
//                for (int i = 0; i < coordinate.size(); i++) {
//                    grown_coordinate[i] = {
//                        coordinate[i].x + number_of_layer * length[i] * point_normals[i].x,
//                        coordinate[i].y + number_of_layer * length[i] * point_normals[i].y,
//                        coordinate[i].z + number_of_layer * length[i] * point_normals[i].z};
//                }
//                checker_.addPoint(grown_coordinate);
//
//                int number = 0;
//                for (const auto &tri : connector) {
//                    checker_.addElement(HexaTag(number++, 1, TRI_TOP),
//                                        std::vector<int>{idx(tri[0]), idx(tri[1]), idx(tri[2])});
//                }
//                int i = -1;
//                for (const auto &tri : connector) {
//                    i++;
//                    std::vector<int> candidate = {idx(tri[0]), idx(tri[1]), idx(tri[2])};
//                    if (grown_coordinate[tri[0]] == grown_coordinate[tri[1]] ||
//                        grown_coordinate[tri[0]] == grown_coordinate[tri[2]] ||
//                        grown_coordinate[tri[1]] == grown_coordinate[tri[2]]) {
//                        continue;
//                    }
//                    if (checker_.checkIntersect(candidate)) {
//                        record_point.insert(tri[0]);
//                        record_point.insert(tri[1]);
//                        record_point.insert(tri[2]);
//                        continue;
//                    }
//                }
//
//
//            } while (!record_point.empty());
//            std::cout << "finish intersection" << std::endl;
//    }
//    for (int j = 1; j <= number_of_layer; j++) {
//        for (int i = 0; i < coordinate.size(); i++) {
//
//            v[idx(i, j - 1)] = {coordinate[i].x + j * length[i] * point_normals[i].x,
//                                coordinate[i].y + j * length[i] * point_normals[i].y,
//                                coordinate[i].z + j * length[i] * point_normals[i].z};
//        }
//    }
//
//    //预处理后的直接生长
//    if (exist_prism) {
//        for (int i = 0; i < connector.size(); i++) {
//            f.push_back({lower_ids[i][0],
//                         lower_ids[i][1],
//                         lower_ids[i][2],
//                         idx(connector[i][0]),
//                         idx(connector[i][1]),
//                         idx(connector[i][2])});
//        }
//    } 
//    // 切割方案
//    else {
//        // 第一层
//        for (int i = 0; i < connector.size(); i++) {
//            // 三法向分裂情况
//            if (lower_ids[i][0] == lower_ids[i][1] && lower_ids[i][2] == lower_ids[i][1]) {
//                f.push_back({lower_ids[i][0],
//                             idx(connector[i][0]),
//                             idx(connector[i][1]),
//                             idx(connector[i][2])});
//            }
//            // 无法向分裂情况
//            else if (lower_ids[i][0] != lower_ids[i][1] && lower_ids[i][2] != lower_ids[i][1] &&
//                     lower_ids[i][2] != lower_ids[i][0]) {
//                int k1, k2, k3;
//                for (int j = 0; j < 3; j++) {
//                    if (lower_ids[i][j] < lower_ids[i][(j + 1) % 3] &&
//                        lower_ids[i][j] < lower_ids[i][(j + 2) % 3]) {
//                        k1 = j;
//                        k2 = (lower_ids[i][(j + 2) % 3] < lower_ids[i][(j + 1) % 3]) ? (j + 2) % 3
//                                                                                     : (j + 1) % 3;
//                        k3 = (lower_ids[i][(j + 2) % 3] > lower_ids[i][(j + 1) % 3]) ? (j + 2) % 3
//                                                                                     : (j + 1) % 3;
//                        break;
//                    }
//                }
//                f.push_back({lower_ids[i][k1],
//                             lower_ids[i][k2],
//                             idx(connector[i][k2]),
//                             idx(connector[i][k3])});
//                f.push_back({lower_ids[i][k1],
//                             idx(connector[i][k1]),
//                             idx(connector[i][k2]),
//                             idx(connector[i][k3])});
//                f.push_back(
//                    {lower_ids[i][k1], lower_ids[i][k2], lower_ids[i][k3], idx(connector[i][k3])});
//            }
//            // 法向二分情况
//            else {
//                int k1, k2, k3;
//                for (int j = 0; j < 3; j++) {
//                    if (lower_ids[i][j] == lower_ids[i][(j + 1) % 3]) {
//                        k1 = j;
//                        k2 = (j + 1) % 3;
//                        k3 = (j + 2) % 3;
//                        break;
//                    }
//                }
//                if (lower_ids[i][k1] < lower_ids[i][k3]) {
//                    f.push_back({lower_ids[i][k1],
//                                 idx(connector[i][k1]),
//                                 idx(connector[i][k2]),
//                                 idx(connector[i][k3])});
//                } else {
//                    f.push_back({lower_ids[i][k3],
//                                 idx(connector[i][k1]),
//                                 idx(connector[i][k2]),
//                                 idx(connector[i][k3])});
//                    f.push_back({lower_ids[i][k1],
//                                 idx(connector[i][k1]),
//                                 idx(connector[i][k2]),
//                                 lower_ids[i][k3]});
//                }
//            }
//        }
//    }
//
//     //后续层
//    if (number_of_layer > 1) {
//        for (int j = 0; j < number_of_layer-1 ; j++) {
//            for (int i = 0; i < connector.size(); i++) {
//                f.push_back({idx(connector[i][0], j ),
//                                idx(connector[i][1], j ),
//                                idx(connector[i][2], j ),
//                                idx(connector[i][0], j+1 ),
//                                idx(connector[i][1], j+1 ),
//                                idx(connector[i][2], j+1)});
//            }
//        }
//    }
//
//    return;
//    }
void MNormalMesh::WriteMesh(std::string& f, std::vector<std::array<double, 3>>& points, double len)
{

	points.resize(coordinate.size());
	for (int k = 0; k < coordinate.size(); k++) {
		for (int i = 0; i < 3; i++) {
            points[k][i] = coordinate[k][i] + number_of_layer * length[k] * point_normals[k][i];
		}
	}

	std::ostringstream ss;
	ss << connector.size() << " " << coordinate.size() << " "
		<< "0 0 0 0"
		<< std::endl;
	for (int i = 0; i < connector.size(); i++) {
		ss << i + 1 << " " << connector[i][0] << " "
			<< connector[i][1] << " " << connector[i][2]
			<< " " << attribute[i] << std::endl;
	}
	f = ss.str();

}

bool MNormalMesh::IntersectionCheck(
    const std::vector<std::array<int, 3>>& lower_ids,
    const std::vector<std::array<double, 3>>& bottom_points,
    int lower_num,
    std::vector<double>& length)
{
    auto idx = [&](int base, int layer = 0) {
        if (layer == -1) {
            return static_cast<int>(base);
        }
        return static_cast<int>(base + lower_num + layer * coordinate.size());
    };

    auto buildBoundingBox = [&]() {
        BoundingBox box({
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::lowest(),
            std::numeric_limits<double>::lowest(),
            std::numeric_limits<double>::lowest()
        });

        for (int i = 0; i < coordinate.size(); ++i) {
            box[0] = std::min(box[0], coordinate[i].x);
            box[1] = std::min(box[1], coordinate[i].y);
            box[2] = std::min(box[2], coordinate[i].z);
            box[3] = std::max(box[3], coordinate[i].x);
            box[4] = std::max(box[4], coordinate[i].y);
            box[5] = std::max(box[5], coordinate[i].z);
        }
        return box;
    };

    const int kMaxIter = 20;
    bool zero_step_retry_done = false;

    for (int iter = 0; ; ++iter) {
        IntersecChecker checker_;
        checker_.init(buildBoundingBox());

        // 1. 加入底层点
        std::vector<BLVector> base_pts;
        base_pts.reserve(lower_num);
        for (int i = 0; i < lower_num; ++i) {
            base_pts.emplace_back(
                bottom_points[i][0],
                bottom_points[i][1],
                bottom_points[i][2]
            );
        }
        checker_.addPoint(base_pts);

        // 2. 加入底层三角片
        std::vector<std::pair<HexaTag, std::vector<int>>> elements;
        elements.reserve(connector.size());
        for (size_t i = 0; i < connector.size(); ++i) {
            if (lower_ids[i][0] != lower_ids[i][1] &&
                lower_ids[i][1] != lower_ids[i][2] &&
                lower_ids[i][0] != lower_ids[i][2]) {
                elements.emplace_back(
                    HexaTag(i, 0, TRI_BOTTOM),
                    std::vector<int>{lower_ids[i][0], lower_ids[i][1], lower_ids[i][2]}
                );
            }
        }
        checker_.addElements(elements);

        // 3. 生成顶层候选点
        std::vector<BLVector> grown_coordinate(coordinate.size());
        for (int i = 0; i < coordinate.size(); ++i) {
            grown_coordinate[i] = {
                coordinate[i].x + number_of_layer * length[i] * point_normals[i].x,
                coordinate[i].y + number_of_layer * length[i] * point_normals[i].y,
                coordinate[i].z + number_of_layer * length[i] * point_normals[i].z
            };
        }
        checker_.addPoint(grown_coordinate);

        // 4. 加入顶层三角片
        int number = 0;
        for (const auto& tri : connector) {
            checker_.addElement(
                HexaTag(number++, 1, TRI_TOP),
                std::vector<int>{idx(tri[0]), idx(tri[1]), idx(tri[2])}
            );
        }

        // 5. 检测哪些三角片相交，并记录相关点
        std::set<int> record_point;
        for (const auto& tri : connector) {
            std::vector<int> candidate = {idx(tri[0]), idx(tri[1]), idx(tri[2])};
            if (checker_.checkIntersect(candidate)) {
                record_point.insert(tri[0]);
                record_point.insert(tri[1]);
                record_point.insert(tri[2]);
            }
        }

        // 6. 没有相交，结束
        if (record_point.empty()) {
            std::cout << "finish intersection" << std::endl;
            return true;
        }

        // 7. 有相交，执行步长回缩
        if (iter < kMaxIter) {
            for (int p : record_point) {
                length[p] *= 0.8;
            }
        } else if (!zero_step_retry_done) {
            for (int p : record_point) {
                length[p] = 0.0;
            }
            zero_step_retry_done = true;
        } else {
            // 置0重试后仍然相交，明确失败，避免死循环
            std::cout << "intersection resolve failed after zero-step retry" << std::endl;
            return false;
        }
    }
}




void splite_by_faceID(std::vector<std::array<double, 3>> &points,
                      std::vector<std::array<double, 3>> &points_multiply,
                      std::vector<std::array<double, 3>> &points_nonwall,
                      std::string &f,
                      std::string &f_multiply,
                      std::string &f_nonwall,
                      std::vector<int> surfaceID)
{
    std::istringstream iss(f);
    std::ostringstream oss_multiply, oss_nonwall;
    std::string line;

    std::map<int, int> used_ids_multiply;
    std::map<int, int> used_ids_nonwall;

    int pointid_multiply = 0;
    int pointid_nonwall = 0;

    int faceid_multiply = 0;
    int faceid_nonwall = 0;

    std::getline(iss, line); // 用来跳过 header 行

    while (std::getline(iss, line)) {
        std::istringstream linestream(line);
        int faceid, id1, id2, id3, surfaceid;
        linestream >> faceid >> id1 >> id2 >> id3 >> surfaceid;

        std::array<int, 3> ids = {id1, id2, id3};
        std::array<int, 3> face;

        bool is_multiply =
            std::find(surfaceID.begin(), surfaceID.end(), surfaceid) != surfaceID.end();

        auto &used_ids = is_multiply ? used_ids_multiply : used_ids_nonwall;
        auto &out_points = is_multiply ? points_multiply : points_nonwall;
        auto &out_stream = is_multiply ? oss_multiply : oss_nonwall;
        int &pointid = is_multiply ? pointid_multiply : pointid_nonwall;
        int &out_faceid = is_multiply ? faceid_multiply : faceid_nonwall;

        for (int i = 0; i < 3; ++i) {
            if (!used_ids.count(ids[i])) {
                face[i] = pointid;
                used_ids[ids[i]] = pointid++;
                out_points.push_back(points[ids[i]]);
            } else {
                face[i] = used_ids[ids[i]];
            }
        }

        out_stream << ++out_faceid << " " << face[0] << " " << face[1] << " " << face[2] << " "
                   << surfaceid << "\n";
    }
    std::ostringstream header_multiply, header_nonwall;
    header_multiply << faceid_multiply << " " << pointid_multiply << " 0 0 0 0\n";
    header_nonwall << faceid_nonwall << " " << pointid_nonwall << " 0 0 0 0\n";

    f_multiply = header_multiply.str() + oss_multiply.str();
    f_nonwall = header_nonwall.str() + oss_nonwall.str();
}
void combine_by_faceID(std::vector<std::array<double, 3>> &points,
                       std::vector<std::array<double, 3>> points_multiply,
                       std::vector<std::array<double, 3>> points_nonwall,
                       std::string &f,
                       std::string f_multiply,
                       std::string f_nonwall)
{
    std::map<std::array<double, 3>, int> point_to_new_id;
    std::vector<std::array<double, 3>> new_points;
    std::ostringstream oss;

    // 添加点，并记录新编号（去重）
    auto add_point = [&](const std::array<double, 3> &pt) -> int {
        if (point_to_new_id.find(pt) != point_to_new_id.end()) {
            return point_to_new_id[pt];
        }
        int id = new_points.size();
        point_to_new_id[pt] = id;
        new_points.push_back(pt);
        return id;
    };

    // 处理一个面数据字符串，转换面编号
    auto process_faces = [&](const std::string &face_data,
                             const std::vector<std::array<double, 3>> &source_points,
                             int &face_id_counter) {
        std::istringstream iss(face_data);
        std::string line;
        std::getline(iss, line); // 用来跳过 header 行

        while (std::getline(iss, line)) {
            std::istringstream linestream(line);
            int fid, p1, p2, p3, sid;
            linestream >> fid >> p1 >> p2 >> p3 >> sid;

            int new_p1 = add_point(source_points[p1]);
            int new_p2 = add_point(source_points[p2]);
            int new_p3 = add_point(source_points[p3]);

            oss << ++face_id_counter << " " << new_p1 << " " << new_p2 << " " << new_p3 << " "
                << sid << "\n";
        }
    };

    int face_id = 0;
    process_faces(f_multiply, points_multiply, face_id);
    process_faces(f_nonwall, points_nonwall, face_id);

    std::ostringstream header;
    header << face_id << " " << new_points.size() << " 0 0 0 0\n";

    f = header.str() + oss.str();
    points = new_points;
}
