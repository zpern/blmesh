#pragma once

#include "Eigen/Dense"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace BLMESH
{

namespace zju
{

//! \addtogroup core
//! @{

//! 齐次坐标类型
using vec4d = typename Eigen::Vector4d;

//! 三维点坐标类型
using vec3d = typename Eigen::Vector3d;

//! UV参数区域点坐标类型
using vec2d = typename Eigen::Vector2d;

//! 三角形三个顶点的索引类型
using vec3i = typename Eigen::Vector3i;

//! 一条边两个顶点索引类型
using vec2i = typename Eigen::Vector2i;


/**
 * @brief 这是网格数据的结构体
 * @struct Mesh
 * 
*/
struct Mesh
{
	//! 网格维度
	size_t dim;

	//! 网格点的坐标 n * dim
	Eigen::MatrixXd vertex;

	//! 网格拓扑 n * m, m是网格边的数目
	Eigen::MatrixXi topo;

	//! 网格面标记，记录网格面的属性信息
	Eigen::MatrixXi masks;

	//! 记录网格面的法向
	Eigen::MatrixXd f_normal;
	
	Eigen::MatrixXd v_normal;
};

/**
 * @brief 这是pair类型哈希映射结构体
 * @struct Hashfunc
*/
struct Hashfunc {
	/**
	 * @brief 重载结构体对象函数
	 * 
	 * @param key 
	 * @return size_t 
	*/
	size_t operator() (const std::pair<int,int>& key) const{
		return std::hash<int>()(key.first) ^ std::hash<int>()(key.second);
	}
};

/**
 * @brief 这是pair类型比较结构体
 * @struct Equalfunc
 * 
*/
struct Equalfunc {
	/**
	 * @brief 重载结构体对象函数
	 * 
	 * @param a 
	 * @param b 
	 * @return true 
	 * @return false 
	*/
	bool operator() (const std::pair<int,int>& a, const std::pair<int,int>& b) const{
		return a.first == b.first && a.second == b.second;
	}
};


//! @}


//! \defgroup core core
//! \brief 核心数据定义

//! \defgroup algorithms algorithms
//! \brief 几何算法处理



} // namespace zju

} // namespace VGModel