
#if !defined(_OCTREEAGENT_H_)
#define _OCTREEAGENT_H_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
#include <array>
#include "common.h"
#include <cmath>
#include <vector>
#include "BLMesh_define.h"
#include "DynamicArray.h"
#include <algorithm>
#include<unordered_set>

class OctreeAgent 
{
public:
	OctreeAgent()=delete;
	//OctreeAgent::OctreeAgent(int* pEle, double* pCoord);
	OctreeAgent(DynamicArray<SearchTriangle>& pEle, MBLNode *pNod);
	virtual ~OctreeAgent();
	

	BLVector getnormal(int data);

	std::array<BLVector,3> getcoord(int data);
	std::array<double, 3> getcoord(int data, int dim);
	void divideData(const OCCUBE& cube, const std::vector<int>& vecData, std::vector<int>(& outData)[8]);
	int patchCube(int patchIdx, const OCCUBE& cube);
	bool LineInCube(BLVector start, BLVector end, const OCCUBE & cube) ;
	inline bool patchTotalInCube(int patchIdx, const OCCUBE & cube);
	bool LineTotalInCube(BLVector start, BLVector end, const OCCUBE& cube);
	bool patchInCube(int patchIdx, const OCCUBE& cube);
	//void setElm(DynamicArray<std::array<int, 3>> &elm) {this->pEle = elm;}


	void setBoxsArray(std::vector<std::array<double, 6>>* boxs, int offset) {
		boxs_ptr = boxs;
	}

	DynamicArray<SearchTriangle>& getElm(){ return pEle; }
	void setNod(MBLNode *nod) {this->pNod = nod;}
	MBLNode *getNod(){ return pNod;}
    void getBox(const int & index, std::array<double, 6>& box);
	
public:
	std::unordered_set<int> mset;
	std::vector<std::array<double, 6>>* boxs_ptr;
	DynamicArray<SearchTriangle>& pEle;

	//double* pCoord;
	MBLNode *pNod;
};

#endif // !defined(_OCTREEAGENT_H_)
