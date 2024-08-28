/*
* @brife ： 适应结构化存储的三棱柱存储仓库的网格数据结构，降低内存消耗
* @author ： hfye@zju.edu.cn
*/
#ifndef  _PRISM_STORE_HOUSE_
#define  _PRISM_STORE_HOUSE_
#include <array> 
#include <vector>
class PrismStorehouse

{
public:
	PrismStorehouse();
	~PrismStorehouse();

	/*
	* @brife 新增一个单元
	*/
	void AddPrism(const int& decent_id,const short& layer,const std::array<int,6>& nconn);
private:
	std::vector<std::vector<int>> point_list_;///存储点以及其后续点
	int* graph_;///存储第一层的拓扑关系
	int num_graph_cells_;///graph 数组的大小
};


#endif
