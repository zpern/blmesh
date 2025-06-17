#pragma once
#ifndef _MNORMAL_MESH_H_
#define _MNORMAL_MESH_H_
#include<vector>
#include<string>
#include <array>
#include "complexnode.h"
#include "ChamferBehavior.h"
class MNormalMesh {
public:
	MNormalMesh() {}

	void SetBehavior(ChamferBehavior& behavior);

	void ReadPlsBuf(std::string f, std::vector<std::array<double, 3>>& points);
	void ReadPls(string filename="");

	void WriteNorm();
	void WriteVol(std::vector<std::array<double, 3>>& v, std::vector<std::vector<int>>& f, int& lower_num,double len,int& add_point_num);
	void WriteMesh(std::string& f,
                           std::vector<std::array<double, 3>>& points,double d);
    void WriteMem(std::string& f,
                      std::vector<std::array<double, 3>>& points);
	void WriteVtk();
	void WritePls(string filename="");

	/**
	 * @brief Calculate multiple normal visible graph of each vertex on surfacemesh, which is the core of mesh Chamfer.  
	 * 
	 */
	void CalculateMultiNormal();
	void PrintVisibilityConeInfo();
	/**
	 * @brief stiching the complex node with ordinary node 
	 * 
	 */
	void BuildTopo(int faceCount);
	
	/**
	* @brief Generate only one layer of boundary layer mesh. 
	* @note This function only can be call after virutal normal generation and virtual topology generation
	*/
	void GenerateFirstLayer(double step_len);
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
protected:
	void CaculateFrontNormal();
	void CalculateNodeNormal();
	


private:
	ChamferBehavior behavior_;
};
void splite_by_faceID(std::vector<std::array<double, 3>>& point, std::vector<std::array<double, 3>>& point_multiply,
                      std::vector<std::array<double, 3>>& point_nonwall, std::string& f, std::string& f_multiply,
                      std::string& f_nonwall, std::vector<int> srufaceID);
void combine_by_faceID(std::vector<std::array<double, 3>>& points, std::vector<std::array<double, 3>> points_multiply,
                       std::vector<std::array<double, 3>> points_nonwall, std::string& f, std::string f_multiply,
                       std::string f_nonwall);
#endif