// Octree.h: interface for the Octree class.
//
//////////////////////////////////////////////////////////////////////

#if !defined _TIGER_OCTREE_H_
#define _TIGER_OCTREE_H_
#include <unordered_map>
#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
#include <array>
#include <vector>
#include <atomic>
#include <set>
#include <cstdlib>
#include <cstdio>
#include "common.h"
#include "OctreeAgent.h"
#include "geom_func.h"
#include "blmeshapi.h"
#include <atomic>
#include <limits>
#include<unordered_set>
#include <map>
#include <fstream>
#include "Eigen/Dense"

#ifndef FLT_EPSILON
#define FLT_EPSILON 1e-20
#endif
namespace OCT {
	class OctreeAgentAbstract;
	class Octree;
	struct Tri {
		Tri(const std::array<double, 6>& b,const int& i) :
			box(b),

			max_value(-1e20),

			id(i),
			valid(true)
		{}
		Tri() :
			box({ -1e20,-1e20,-1e20,1e20,1e20,1e20 }),

			max_value(-1e20),

			id(-1),
			valid(true)
		{}

		double min_value;
		double max_value;


		
		std::array<double, 6> box;
		bool valid = true;
		int id;

		bool operator < (const Tri& t) {
			if (min_value != t.min_value) {
				return min_value < t.min_value;
			}
			return max_value < t.max_value;
		}

	};
	class TreeNode
	{
		friend class Octree;

	public:
		TreeNode();

		Eigen::Vector3d getPCA();

		TreeNode(const OCCUBE& cube, const std::vector<int>& data, TreeNode* father, int layer, int fnormal, OctreeAgent * pOctree);
		virtual ~TreeNode();
		TreeNode* getFather();
		TreeNode*& getChild(int i);
		OCCUBE getCube() { return ocCube; }
		void appendData(vector<int>& data);
		void appendData_unsorted(int i, std::array<double, 6>& box);
		void appendData(int i, std::array<double, 6>& box);
		void rmData_unsorted(int i);
		void rmData(int i);
		void sortElement();
		int ref_id;

		bool changed;//Marked as true when new element inserted, false when check intersection process done
		vector<int> memory_queue;

		

	private:

	public:
		OctreeAgent * pOctreeAgent;

		int num_delete;
		OCCUBE ocCube;
		int layer;
		BLVector target_direction;

		std::set<int> vecData;
		int num_obj;

		std::vector<Tri> order_data;

		
		int most_dispersed_direction;
		TreeNode* pChildNodes[OC_SUBNODES];
		TreeNode* pFather;
		int num_data;
	};

	class Octree
	{
	public:
#ifdef _DEBUG
		int num_inter;
#endif
		inline int  tri_overlap_test_no_box(std::array<int, 3>& triIdx1, std::array<int, 3>& triIdx2, MBLNode *&pNods, int &share_node_num, int &share1, int& share2);
		static int tri_overlap_test(int triIdx1, int triIdx2, MBLNode * pNods, int * ele);
		//int tri_overlap_test(int triIdx1, int triIdx2, Node * pNods, int * ele);
		Octree();
		Octree(OctreeAgent *agent, int dep = 3) :pOctreeAgent(agent), depth(dep) { sorted = true; node_before = NULL; pRoot = NULL; special_id = -1; sort_telerant = false; }
		virtual ~Octree();

		/*
		 * build an octree
		 * @ cube: initial bounding box of objects
		 * @ vecData: index of objects
		 * @ maxObjNum: max objects number of each nodes
		 * @ maxDepth: max depth of tree
		 **/
		int buildOcTree(const OCCUBE& cube, const std::vector<int>& vecData, int maxObjNum);

		void rebuild();

		int createBranch(TreeNode** subNode, const OCCUBE& cube, const std::vector<int>& vecData, int layer, int maxObjNum, int maxDepth);

		void traverseOcTree();
		void printNodeSize();
		void printNodeSize(TreeNode* r, vector<int> & p);

		/*
		 *
		 * @brife 
		 **/
		void printElement(std::string filename);

		bool TriIn(int id);
		int getTriNum();

		void saveOCTreeVectorVTK(std::string filename);
		void saveOCTreeVTK(std::string filename);

		void saveBiggestOctantVectorVTK(std::string filename);

		/*print transfer info to cmd windows**/
		void printTransferCost();
		/*
		 * trim the nodes that does not contain any object
		 * @ return number of trimmed nodes
		 **/
		int trim(TreeNode* node);

		/*
		 * insert an object into the tree
		 * @ data: index of the object
		 **/
		void insert(int data, TreeNode *node, const OCCUBE cube, std::array<double, 6>& box);
		/*
		*@brife insert a lot elements
		**/
		void insert(vector<int>& data, TreeNode* node);
		void insertPreProcess(int data);
		void rmDataPreProcess(int data);
		void sortNode(TreeNode *node);

		/*
		 * insert an object into the tree
		 * @ data: index of the object
		 **/
		void rmData(int data, TreeNode *node, const OCCUBE cube, std::array<double, 6>& box);

		bool chckIntersectPreProcess(int data);

		void destroyNode(TreeNode* p);



		double chckIntersectWithLine(BLVector start, BLVector end, TreeNode *node, const OCCUBE& cube);

		double chckIntersectWithLine(BLVector start, BLVector end);
		/*
		 * query an object
		 * @ data: index of the object
		 **/
		TreeNode* query(int data);

		int getDepth() { return depth; }

		TreeNode* getRootNode() { return pRoot; }

		TreeNode* getNodeBefore() { return node_before; }
		inline void setNodeBefore(TreeNode* t) { node_before = t ? t : pRoot; }
	private:
		bool chckIntersect(int i, TreeNode *node, const OCCUBE cube, std::array<double, 6>& box);
		
	public:
		std::set<int> possible_items;


		static int time_stamp;
		
		int last_intersection_;
		bool sort_telerant;
		bool sorted;

		int special_id;
		OctreeAgent * pOctreeAgent;
#ifdef _DEBUG
		double transfer_time_cost_;
		int transfer_times_;
		int call_times_;
#endif
		void record_stone_position();
		inline int get_stone();
		void sort_octants();
		void find_intersected_triangles();
		bool check_intersection_in_set(int id);
	public:
		//static int time_machine;
	private:
		

		std::vector<TreeNode*> leaf_nodes;
		int time_stone_;
		int most_length_;
		int second_length_;
		int last_length_;
		int max_depth_;
		TreeNode* pRoot;
		TreeNode* node_before;
		int depth;
		std::set<int> vec_data;

	};



	int tri_tri_inter(double * A, double * B, double * C, double * O, double * P, double * Q);
}

#endif // !defined(OCTREE_H__)