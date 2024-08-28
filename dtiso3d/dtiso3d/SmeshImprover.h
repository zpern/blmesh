#pragma once
#include <vector>
class SmeshImprover {
public:
	void initSmesh(std::vector<double> pts, std::vector<int> tris, std::vector<double> sizes);
	int setSizeRatio(int val); // 0-1归一化参数，用于控制整体尺寸，val越大尺寸越大
	bool needImprove(); //之后再定是否需要improve的标准
	int doImprove();
	void ouputSmesh(std::vector<double>& pts, std::vector<int>& tris, std::vector<int>& triToOri); //triToOri为每个面到原始面的映射
};