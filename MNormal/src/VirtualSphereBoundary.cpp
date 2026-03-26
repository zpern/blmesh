
#include "..\include\VirtualSphereBoundary.h"
#include "..\include\geometryfunction.h"
#include <set>
#include <map>
#include <algorithm>
#define BAD_QUALITY 0.85
map<set<int>, VirtualSphereMesh> VirtualSphereMesh::best_result = map<set<int>, VirtualSphereMesh>();
void VirtualSphereMesh::removeInvalidNode()
{
	valid.clear();
	for (auto i : boundary_edges_) {
		valid.insert(i[0]);
		valid.insert(i[1]);
	}
	for (auto i : triangle_lists_) {
		valid.insert(i.point_index_[0]);
		valid.insert(i.point_index_[1]);
		valid.insert(i.point_index_[2]);
	}
	valid.insert(-1 - add_node_count);


}
std::vector<VEDGE> VirtualSphereMesh::findBoundaryBetween(int start, int end)
{
	int left_start_index, left_end_index;
	for (int i = 0; i < boundary_edges_.size(); i++) {
		if (start == boundary_edges_[i][1]) {
			left_end_index = i;
		}
		if (end == boundary_edges_[i][0]) {
			left_start_index = i;
		}
	}
	std::vector<VEDGE> leaf_edges;
	for (int ptr = left_start_index; ptr%boundary_edges_.size() != (left_end_index+1) % boundary_edges_.size(); ptr++) {
		leaf_edges.push_back(boundary_edges_[ptr%boundary_edges_.size()]);
	}
	leaf_edges.push_back(VEDGE{ start,end });
	return leaf_edges;
}
VirtualSphereMesh VirtualSphereMesh::mergeMesh(VirtualSphereMesh & b1, VirtualSphereMesh & b2)
{
	VirtualSphereMesh ans=b1;
	// merge points
	for (auto i : b2.triangle_lists_)
		ans.triangle_lists_.push_back(i);
	for (auto i : b2.valid)
		ans.valid.insert(i);


	if (b2.virtual_point_lists_.size() > ans.virtual_point_lists_.size()) {
		ans.virtual_point_lists_ = b2.virtual_point_lists_;
		ans.valid = b2.valid;
	}

	//modify negative bit
	if (ans.valid.find(ONE_INSERT) != ans.valid.end())
		ans.valid.erase(NO_INSERT);
	ans.add_node_count = b1.add_node_count + b2.add_node_count;
	

	//merge boundary edge, by xor
	ans.boundary_edges_.clear();
	for (int j = 0; j < b2.boundary_edges_.size();j++) {
		bool in = false;
		for (int i = 0; i < b1.boundary_edges_.size(); i++) {
			if (b2.boundary_edges_[j][0] == b1.boundary_edges_[i][1] && b2.boundary_edges_[j][1] == b1.boundary_edges_[i][0])
				in = true;
		}
		if (!in)
			ans.boundary_edges_.push_back(b2.boundary_edges_[j]);
	}
	for (int j = 0; j < b1.boundary_edges_.size(); j++) {
		bool in = false;
		for (int i = 0; i < b2.boundary_edges_.size(); i++) {
			if (b2.boundary_edges_[i][0] == b1.boundary_edges_[j][1] && b2.boundary_edges_[i][1] == b1.boundary_edges_[j][0])
				in = true;
		}
		if (!in)
			ans.boundary_edges_.push_back(b1.boundary_edges_[j]);
	}
	return ans;


}
const int VirtualSphereMesh::getExtraPointCount() const
{
	if (valid.size())
		return (valid.size())/2;
	return 0;
}
const int VirtualSphereMesh::getValidPointCount() const
{
	if(valid.size())
	return valid.size()-1;
	return 1;
}
bool VirtualSphereMesh::isPointInBoundary(BLVector pos)
{
	for (auto j : boundary_edges_) {

		bool valid = true;

		
		for (auto i : boundary_edges_) {
			
			if ((!GEEOMETRY_FUNCTION::isEdgeintersected(pos,virtual_point_lists_[j[0]].getCoord(), virtual_point_lists_[i[0]].getCoord(), virtual_point_lists_[i[1]].getCoord()))&&\
				(!GEEOMETRY_FUNCTION::isEdgeintersected(pos, virtual_point_lists_[j[1]].getCoord(), virtual_point_lists_[i[0]].getCoord(), virtual_point_lists_[i[1]].getCoord()))) {
				valid = false;
			}
		}
		if (valid) {
			BLVector v1 = virtual_point_lists_[j[0]].getCoord() - pos;
			BLVector v2 = virtual_point_lists_[j[1]].getCoord() - pos;
			return (v1^v2)*virtual_point_lists_[j[1]].getCoord() >0;
		}
	}

	return false;
}

void VirtualSphereMesh::saveBoundaryVtk(std::string filename)
{
	int nBdry = boundary_edges_.size();

	const char *filename_ = filename.data();
	int i, idx;
	std::set<int> setpnt;
	std::set<int>::iterator sit;
	FILE *fout = nullptr;
	
	fout = fopen(filename_, "w");
	if (fout == nullptr)
	{
		throw(std::string("no authority to write file"));
	}
	//fprintf(fout, "%d %d 0 0 0 0 0\n", setpnt.size(), nBdry);
	fprintf(fout, "# vtk DataFile Version 2.0\nboundary layer mesh\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS %d double\n", virtual_point_lists_.size() +1);
	sit = setpnt.begin();
	i = 0;
	fprintf(fout, "0 0 0\n");
	for(int i=0;i< virtual_point_lists_.size();i++)
	{
		idx = i;
		VPoint* ptr = dynamic_cast<VPoint*>(&virtual_point_lists_[idx]);
		fprintf(fout, "%f %f %f\n", ptr->getCoord()[0], ptr->getCoord()[1], ptr->getCoord()[2]);
	}
	nBdry = boundary_edges_.size();
	fprintf(fout, "CELLS %d %d\n", nBdry, nBdry * 3);

	for (i = 0; i < nBdry; i++)
	{
		fprintf(fout, " 2 %d %d\n", boundary_edges_[i][0]+1, boundary_edges_[i][1]+1);
	}
	fprintf(fout, "CELL_TYPES %d\n", nBdry);
	for (i = 0; i < nBdry; i++)
	{
		fprintf(fout, " 3\n");
	}

	fclose(fout);
	fout = nullptr;
}

void VirtualSphereMesh::saveTriVtk(std::string filename)
{

	int nBdry = triangle_lists_.size();

	const char *filename_ = filename.data();
	int i, idx;
	std::set<int> setpnt;
	std::set<int>::iterator sit;
	FILE *fout = nullptr;

	for (i = 0; i < nBdry; i++)
	{
		setpnt.insert(triangle_lists_[i].point_index_[0]);
		setpnt.insert(triangle_lists_[i].point_index_[1]);
		setpnt.insert(triangle_lists_[i].point_index_[2]);
	}
	map<int, int> ptmap;


	fout = fopen(filename_, "w");
	if (fout == nullptr)
	{
		throw(std::string("no authority to write file"));
	}
	//fprintf(fout, "%d %d 0 0 0 0 0\n", setpnt.size(), nBdry);
	fprintf(fout, "# vtk DataFile Version 2.0\nboundary layer mesh\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS %d double\n", setpnt.size() + 1);
	sit = setpnt.begin();
	i = 0;
	fprintf(fout, "0 0 0\n");
	while (sit != setpnt.end())
	{
		idx = *sit;
		VPoint* ptr = dynamic_cast<VPoint*>(&virtual_point_lists_[idx]);
		fprintf(fout, "%f %f %f\n", ptr->getCoord()[0], ptr->getCoord()[1], ptr->getCoord()[2]);
		ptmap[idx] = i++;
		++sit;
	}

	fprintf(fout, "CELLS %d %d\n", nBdry+ boundary_edges_.size(), nBdry * 4+ boundary_edges_.size()*3);

	for (i = 0; i < nBdry; i++)
	{
		fprintf(fout, " 3 %d %d %d\n", ptmap[triangle_lists_[i].point_index_[0]] + 1, ptmap[triangle_lists_[i].point_index_[1]] + 1, ptmap[triangle_lists_[i].point_index_[2]] + 1);
	}

	for (i = 0; i < boundary_edges_.size(); i++)
	{
		fprintf(fout, " 2 %d %d\n", ptmap[boundary_edges_[i][0]] + 1, ptmap[boundary_edges_[i][1]] + 1);
	}
	fprintf(fout, "CELL_TYPES %d\n", nBdry+ boundary_edges_.size());
	for (i = 0; i < nBdry; i++)
	{
		fprintf(fout, " 5\n");
	}

	for (i = 0; i < boundary_edges_.size(); i++)
	{
		fprintf(fout, " 3\n");
	}
	fclose(fout);
	fout = nullptr;

}

map<int, int> VirtualSphereMesh::getMapNodeAfterFarNode()
{
	map<int, int> ans;



	return ans;
}

void VirtualSphereMesh::addTriangle(std::array<int, 3> tri)
{
	triangle_lists_.push_back(VTriangle());
	triangle_lists_.back().point_index_ = tri;
}

double VirtualSphereMesh::getEdgeSkewness(int i)
{
	double ans = 0;
	if ((!virtual_point_lists_[boundary_edges_[i][0]].isFarNode()) && (!virtual_point_lists_[boundary_edges_[i][1]].isFarNode())) {
		ans = max(ans, abs(ANGLE(virtual_point_lists_[boundary_edges_[i][0]].getCoord(), virtual_point_lists_[boundary_edges_[i][1]].getCoord()) - 90) / 90);
	}
	return ans;
}

// skweness 




double VirtualSphereMesh::getTriSkewness(int i)
{
	return getTriSkewness(triangle_lists_[i].point_index_);
}

double VirtualSphereMesh::getTriSkewness(std::array<int, 3> tri)
{
	double ans = 0;
	BLVector p[3];
	int count = 0;
	for (int j = 0; j < 3; j++) {
		p[j] = virtual_point_lists_[tri[j]].getCoord();
	}
	for (int k = 0; k < 3; k++)
		for (int j = k + 1; j < 3; j++)
			if ((!virtual_point_lists_[tri[j]].isFarNode()) && (!virtual_point_lists_[tri[k]].isFarNode())) {
				ans = max(ans, abs(ANGLE(p[k], p[j]) - 90) / 90);
			}
	for (int k = 0; k < 3; k++)
		if (!virtual_point_lists_[tri[k]].isFarNode()) {
			if ((!virtual_point_lists_[tri[(k + 1) % 3]].isFarNode()) && (!virtual_point_lists_[tri[(k + 2) % 3]].isFarNode())) {
				ans = max(ans, (ANGLE(p[(k + 1) % 3] - p[k], p[(k + 2) % 3] - p[k]) - 60) / 120);
				ans = max(ans, (60 - ANGLE(p[(k + 1) % 3] - p[k], p[(k + 2) % 3] - p[k])) / 60);
			}
			else {
				ans = max(ans, abs((ANGLE(p[(k + 1) % 3] - p[k], p[(k + 2) % 3] - p[k]) - 90)) / 90);
			}
		}
	for (int k = 0; k < 3; k++) {
		if (virtual_point_lists_[tri[k]].isFarNode() && virtual_point_lists_[tri[(k + 1) % 3]].isFarNode()) {
			ans = 1;
		}
	}
	return ans;
}

VirtualSphereMesh::VirtualSphereMesh()
{

}

VirtualSphereMesh::VirtualSphereMesh(const vector<BLVector>& convex_edge, const std::vector<BLVector>& facets_set_normal, vector<vector<int>> edge2face, vector<int> splitter)
{

	max_skewness_ = 0;
	if (convex_edge.size() != facets_set_normal.size())
		throw std::logic_error("impossbile input!");
	for (int i = 0; i < convex_edge.size(); i++) {
		// build point lists
		VPoint point ;
		point.setCoord(convex_edge[i]);
		point.setIndex(virtual_point_lists_.size());
		point.setGlobalIndex(splitter[i]);
		virtual_point_lists_.push_back(point);
		


		VPoint point_sphere  ;
		point_sphere.setCoord(facets_set_normal[i]);
		point_sphere.setIndex(virtual_point_lists_.size());
		point_sphere.setGlobalNeighbourTriIndex(edge2face[i]);
		virtual_point_lists_.push_back(point_sphere);
		
	}

	for (int i = 0; i < convex_edge.size(); i++) {
		// build edge lists
		auto edge_1 = VEDGE();
		auto edge_2 = VEDGE();
		(edge_1)[0] = 2 * i + 0;
		(edge_1)[1] = 2 * i + 1;
		int j = (i + 1) % convex_edge.size();
		(edge_2)[0] = 2 * i + 1;
		(edge_2)[1] = 2 * j + 0;
		boundary_edges_.push_back(edge_1);
		boundary_edges_.push_back(edge_2);
	}
}

bool VirtualSphereMesh::isAddnode()
{
	return add_node_count>0;
}

void VirtualSphereMesh::addNode()
{
	BLVector ans;
	triangle_lists_.clear();
	

	//get valid point
	vector<int> points;
	for (auto i : valid)
		if (i >= 0&&(!virtual_point_lists_[i].isFarNode()))
			points.push_back(i);
	double min_skewness = 1;
	for (int i = 0; i < points.size(); i++) {
		for (int j = i+1; j < points.size(); j++) {
			for (int k = j+1; k < points.size(); k++) {
				BLVector center=GEEOMETRY_FUNCTION::solveCenterPointOfSphereCap(virtual_point_lists_[i].getCoord(), virtual_point_lists_[j].getCoord(), virtual_point_lists_[k].getCoord());	
				if (center.magnitude() < 0.9)
					continue;
				virtual_point_lists_.push_back(VPoint());
				virtual_point_lists_.back().setCoord(center);
				virtual_point_lists_.back().setIndex(virtual_point_lists_.size()-1);
				double max_skewness = 0;
				for (int l = 0; l < boundary_edges_.size(); l++) {
					std::array<int, 3> tri{boundary_edges_[l][0],boundary_edges_[l][1] ,virtual_point_lists_.size() - 1 };
					max_skewness=max(max_skewness,getTriSkewness(tri));
				}
				if (min_skewness > max_skewness) {
					min_skewness = max_skewness;
					ans = center;
				}
				virtual_point_lists_.pop_back();
			}
		}
	}
	

	if (ans.magnitude() < 0.9)
		return;
	int index = virtual_point_lists_.size() ;
	ans.normalize();
	virtual_point_lists_.push_back(VPoint());
	virtual_point_lists_.back().setCoord(ans);
	virtual_point_lists_.back().setIndex(index);

	for (int i = 0; i < boundary_edges_.size(); i++) {
		VTriangle t;
		t.point_index_[0] = boundary_edges_[i][0];
		t.point_index_[1] = boundary_edges_[i][1];
		t.point_index_[2] = index;
		triangle_lists_.push_back(t);
	}
	valid.erase(NO_INSERT);
	valid.insert(ONE_INSERT);

	max_skewness_ = 0.0;
	for (int i = 0; i < triangle_lists_.size(); i++) {
		max_skewness_ = max(max_skewness_, getTriSkewness(i));
	}
	add_node_count++;
}

double VirtualSphereMesh::AddNodeSphereSpr()
{
	if (VirtualSphereMesh::best_result.find(valid) != VirtualSphereMesh::best_result.end()) {
		*this = VirtualSphereMesh::best_result[valid];
		return max_skewness_;
	}
	removeInvalidNode();
	int point_count = getValidPointCount();
	switch (point_count)
	{
	case 1:
		break;
	case 2:
		max_skewness_ = getEdgeSkewness(0);
		break;
	case 3:
		addNode();
		if(isAddnode())
		for(int k=0;k<3;k++)
			max_skewness_ =max(max_skewness_, getTriSkewness(k));
		break;
	default:
		int p1 = boundary_edges_[0][0];
		int p2 = boundary_edges_[0][1];
		VirtualSphereMesh best_one;
		double min_skewness = 1.0;
		for (int i = 0; i < virtual_point_lists_.size(); i++) {
			if (valid.find(virtual_point_lists_[i].getIndex()) == valid.end())
				continue;
			if (virtual_point_lists_[i].getIndex() == p1 || virtual_point_lists_[i].getIndex() == p2)
				continue;
			int p3 = virtual_point_lists_[i].getIndex();
			std::array<int, 3> my_tri{ p1,p2,p3 };

			// create new edges
			VirtualSphereMesh left = *this;
			VirtualSphereMesh right = *this;

		

			left.boundary_edges_ = findBoundaryBetween(p1, p3);
			right.boundary_edges_ = findBoundaryBetween(p3, p2);

			
			VirtualSphereMesh median;
			median.virtual_point_lists_ = this->virtual_point_lists_;
			median.boundary_edges_.resize(3);
			median.boundary_edges_[0] = std::array<int, 2>{p1, p2};
			median.boundary_edges_[1] = std::array<int, 2>{p2, p3};
			median.boundary_edges_[2] = std::array<int, 2>{p3, p1};

			median.SphereSpr();
			median.removeInvalidNode(); 
			double tri_skewness = getTriSkewness(my_tri);


			VirtualSphereMesh left_with_node = left;
			VirtualSphereMesh right_with_node = right;
			double right_skewness = right.SphereSpr();
			double leaf_skewness = left.SphereSpr();
			double left_with_node_skewness = left_with_node.AddNodeSphereSpr();
			double right_with_node_skewness = right_with_node.AddNodeSphereSpr();


			VirtualSphereMesh result[3];
			result[0] = VirtualSphereMesh::mergeMesh(left, median);
			result[0] = VirtualSphereMesh::mergeMesh(result[0], right_with_node);			
			result[0].max_skewness_ = max(leaf_skewness, tri_skewness);
			result[0].max_skewness_ = max(result[0].max_skewness_, right_with_node_skewness);
			result[0].removeInvalidNode();

			result[1] = VirtualSphereMesh::mergeMesh(left_with_node, median);
			result[1] = VirtualSphereMesh::mergeMesh(result[1], right);
			result[1].max_skewness_ = max(right_skewness, tri_skewness);
			result[1].max_skewness_ = max(result[1].max_skewness_, left_with_node_skewness);
			result[1].removeInvalidNode();

			result[2] = *this;
			result[2].addNode();


			for (int k = 0; k < 3; k++) {				
				if(result[k].isAddnode())
				if (min_skewness > result[k].max_skewness_&&result[k].triangle_lists_.size() == point_count) {
					min_skewness = result[k].max_skewness_;
					best_one = result[k];
				}
			}
			
		
		}
		*this = best_one;

		break;
	}

	removeInvalidNode();
	VirtualSphereMesh::best_result[valid] = *this;
	return max_skewness_;
	//*this = VirtualSphereMesh::best_result[valid];

}

double VirtualSphereMesh::SphereSpr()
{
	if (VirtualSphereMesh::best_result.find(valid) != VirtualSphereMesh::best_result.end()) {
		*this = VirtualSphereMesh::best_result[valid];
		return max_skewness_;
	}
	removeInvalidNode();
	int point_count = getValidPointCount();
	switch (point_count)
	{
	case 1:
		break;
	case 2:
		max_skewness_ = getEdgeSkewness(0);
		//if (skewness > 0.85)
		//	cout << "warning! low quality mesh find!" << endl;
		break;
	case 3:
		addTriangle(std::array<int, 3>{boundary_edges_[0][0], boundary_edges_[1][0], boundary_edges_[2][0]});
		max_skewness_ = getTriSkewness(0);
		//if (skewness > BAD_QUALITY)
		//	cout << "warning! low quality mesh find!" << endl;
		break;
	default:
		int p1 = boundary_edges_[0][0];
		int p2 = boundary_edges_[0][1];
		VirtualSphereMesh best_one;
		double min_skewness = 1.0;
		for (int i = 0; i < virtual_point_lists_.size(); i++) {
			if (valid.find(virtual_point_lists_[i].getIndex()) == valid.end())
				continue;
			if (virtual_point_lists_[i].getIndex() == p1 || virtual_point_lists_[i].getIndex() == p2)
				continue;
			int p3 = virtual_point_lists_[i].getIndex();
			std::array<int, 3> my_tri{ p1,p2,p3 };

			// create new edges
			VirtualSphereMesh left = *this;
			VirtualSphereMesh right = *this;


			VirtualSphereMesh median;
			//saveBoundaryVtk("test_bdry1.vtk");

			left.boundary_edges_ = findBoundaryBetween(p1, p3);
			right.boundary_edges_ = findBoundaryBetween(p3, p2);


			
			median.virtual_point_lists_ = this->virtual_point_lists_;
			median.boundary_edges_.resize(3);
			median.boundary_edges_[0] = std::array<int, 2>{p1, p2};
			median.boundary_edges_[1] = std::array<int, 2>{p2, p3};
			median.boundary_edges_[2] = std::array<int, 2>{p3, p1};


			median.SphereSpr();
			median.removeInvalidNode();


			 
			int size = valid.size();

			

			left.removeInvalidNode();
			right.removeInvalidNode();

			//left.saveBoundaryVtk("leftBoundary.vtk");
			//right.saveBoundaryVtk("rightBoundary.vtk");
			//median.saveBoundaryVtk("medianBoundary.vtk");

			double right_skewness = right.SphereSpr();
			double leaf_skewness = left.SphereSpr();

			//left.saveTriVtk("left.vtk");
			//right.saveTriVtk("right.vtk");
			//median.saveTriVtk("median.vtk");

			auto result = VirtualSphereMesh::mergeMesh(left, median);
			result = VirtualSphereMesh::mergeMesh(result, right);

			result.removeInvalidNode();
			double tri_skewness = getTriSkewness(my_tri);
			result.max_skewness_ = max(leaf_skewness, right_skewness);
			result.max_skewness_ = max(result.max_skewness_, tri_skewness);
			if (min_skewness > result.max_skewness_&&result.triangle_lists_.size() == point_count - 2) {
				min_skewness = result.max_skewness_;
				best_one = result;
			}

		}
		*this = best_one;

		break;
	}

	removeInvalidNode();
	VirtualSphereMesh::best_result[valid] = *this;
	return max_skewness_;
	//*this = VirtualSphereMesh::best_result[valid];


	

}
