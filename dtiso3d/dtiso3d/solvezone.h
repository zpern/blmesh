#ifndef _SOLVEZONE_H_
#define _SOLVEZONE_H_

#include <vector>
#include <map>

// Zone <=> Subdomain <=> 封闭连通区域
class BlockedArray_Node;
class BlockedArray_Elem;

struct ZoneAttachInfo {
	std::vector<int> tri2Faces;
	std::vector<int> tri2FaceOffsets; //dim = nFace + 1, 方便差值计算数量
	std::vector<int> block2Faces;
	std::vector<int> block2FaceOffsets; //dim = nFace + 1, 方便差值计算数量
	std::vector<int> pecBlocks;
	std::vector<double> blockSize; //未指定则为负数
	std::vector<int> tri2Domain;
	bool valid = false;

	std::map<int, double> getFaceSizes();
};

struct ZoneOutInfo {
	std::vector<int> tri2Domain;
	std::vector<int> tri2Face;
	std::vector<int> tri2FaceOffsets; //dim = nFace + 1, 方便差值计算数量
	std::vector<int> elem2Domains;
	std::vector<std::vector<int>> domain2Blocks;
};

class ZoneSolver //给面网格划分子连通区域，并且确定某个连通区域属于哪些体
{
public:
	ZoneSolver(BlockedArray_Node* pnodes, BlockedArray_Elem* pelems);

	int divideZone(const std::vector<int>& bdTris, const std::vector<std::vector<int>>& trisToBlock); //用输入表面划分网格为不同区域
	int deleteBlock(int block);
	int getShrinkedMesh(std::vector<double>& pts_out, std::vector<int>& elems_out, std::vector<int>& ptshash); //在网格数据中去除没有用到的点
	std::vector<std::vector<int>> getZoneToBlock() {
		return zoneToBlock;
	}
	std::vector<int> getElemToZone() {
		return elemToZone;
	}
protected:
	std::vector<double> pts;
	std::vector<int> elems;
	std::vector<int> elemToZone;
	std::vector<std::vector<int>> zoneToBlock;

	std::vector<int> tris;
	std::vector<std::vector<int>> triToFaces;
};

template <typename T>
std::vector<std::vector<T>> fold2DVector(int d1num, int* d2nums, T* arr) {
	std::vector<std::vector<T>> ret;
	int off = 0;
	for (int i = 0; i < d1num; i++) {
		ret.emplace_back();
		for (int j = 0; j < d2nums[i]; j++) {
			ret.back().push_back(arr[off + j]);
		}
		off += d2nums[i];
	}
	return ret;
}

template <typename T>
int stretch2DVector(const std::vector<std::vector<T>>& vec2, std::vector<int>& d2nums, std::vector<T>& vec1) {
	for (const auto& elem : vec2) {
		d2nums.push_back(elem.size());
		for (const auto& val : elem) {
			vec1.push_back(val);
		}
	}
	return 0;
}

inline int numVec2Offset(const std::vector<int>& nums, std::vector<int>& offs) {
	offs.clear();
	offs.push_back(0);
	int off = 0;
	for (int i = 0; i < nums.size(); i++) {
		off += nums[i];
		offs.push_back(off);
	}
	return 0;
}

#endif //_SOLVEZONE_H_