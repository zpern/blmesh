#pragma once
#ifndef _MNORMAL_MESH_H_
#define _MNORMAL_MESH_H_
#include<vector>
#include<string>
#include <array>
#include "complexnode.h"
#include "./hexa_octree.hpp"   // Uses the previously defined binary tree.
#include "./hexa_tag.h"        // Provides HexaTag and HexaTagHash.
#include "./geometry_check.hpp"  // Assumes it contains the declaration for tri_tri_intersect.
#include "./intersection_check.hpp"

class MNormalMesh {
public:
	MNormalMesh() {}

	void ReadPlsBuf(std::string f, std::vector<std::array<double, 3>>& points);
	void FixedLength();
	void CalculateMultiNormal();
	void SmoothNormalsSimple(int itertion);
	void BuildTopo(int faceCount);

	bool IntersectionCheck(const std::vector<std::array<int, 3>> &lower_ids,
                                       const std::vector<std::array<double, 3>> &bottom_points,
                                       int lower_num,
                                       std::vector<double> &length);
    void pre_WriteVol(std::vector<std::array<double, 3>> &v,std::vector<std::vector<int>> &f,int &lower_num,int &add_point_num);
	void WriteVol(std::vector<std::array<double, 3>>& v, std::vector<std::vector<int>>& f, int& lower_num,int& add_point_num);
	void WriteMesh(std::string& f,std::vector<std::array<double, 3>>& points,double d);


	
    
    

public:
	unsigned int number_of_origin_node;
	unsigned int number_of_node;
	std::vector<BLVector> coordinate;
	std::vector<BLVector> point_normals;
	std::vector<int> real_node_id_;/// mark as virutal tags, id=the coorespoding real node id;

	unsigned int number_of_origin_triangles;
	unsigned int number_of_triangles;
	std::vector<std::array<int, 3>> connector;
	std::vector<int> attribute;

	vector<ComplexNode> node_array;

	int number_of_layer;
	std::map<std::array<double,3>, double> point_to_length;
    std::vector<double> length;
    bool exist_prism =false;
    bool multiplySuccess = true;
    std::vector<int> avoid_spliteNode;
	

protected:
	void CaculateFrontNormal();
	void CalculateNodeNormal();
	
};
void splite_by_faceID(std::vector<std::array<double, 3>>& point, std::vector<std::array<double, 3>>& point_multiply,
                      std::vector<std::array<double, 3>>& point_nonwall, std::string& f, std::string& f_multiply,
                      std::string& f_nonwall, std::vector<int> srufaceID);
void combine_by_faceID(std::vector<std::array<double, 3>>& points, std::vector<std::array<double, 3>> points_multiply,
                       std::vector<std::array<double, 3>> points_nonwall, std::string& f, std::string f_multiply,
                       std::string f_nonwall);
#endif