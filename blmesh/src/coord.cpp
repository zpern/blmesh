/* ----------------------------------------------------------------------------
 * ----------------------------------------------------------------------------
 *
 * 多学科应用模拟的赋能环境
 * Enabling Environment for Multi-displinary Application Simulations
 *
 * 陈建军 中国 浙江大学工程与科学计算研究中心
 * 版权所有	  2007年10月15日
 * Chen Jianjun  Center for Engineering & Scientific Computation,
 * Zhejiang University, P. R. China
 * Copyright reserved, Oct. 15, 2007
 * 
 * 联系方式 (For further information, please conctact)
 *   电话 (Tel)：+86-571-87953165
 *   传真 (Fax)：+86-571-87953167
 *   邮箱 (Mail)：chenjj@zju.edu.cn
 *
 * 文件名称 (File Name)：vector.h
 * 初始版本 (Initial Version): V1.0
 * 功能介绍 (Function Introduction：
 *     定义了一套密度控制机制
 *     Define a set of element spacing controlling scheme.
 * 
 *
 * -----------------------------修改记录 (Revision Record)------------------------
 * 修改者 (Revisor):
 * 修改日期 (Revision Date):
 * 当前版本 (Current Version):
 * 修改介绍 (Revision Introduction):
 * ------------------------------------------------------------------------------
 * ------------------------------------------------------------------------------*/

#include <spdlog/spdlog.h> 
 #include "coord.h"
#include <math.h>


/* *********************************************************************************
 * 定义三维矢量
 * define 3-dimensional vector
 * **********************************************************************************/


/* 矢量求模 calc. magnitude */
double Coord::magnitude() const
{
	return sqrt(x * x + y * y + z * z);
}

/* 归一 normalization */
void Coord::normalize()
{
	double m = magnitude();
	x = x / m;
	y = y / m;
	z = z / m;
}

/* **********************************
 * 重载操作符 overloaded operator 
 * **********************************/

/* 赋值 assignment */
const Coord& Coord::operator = (const Coord& v)
{
	x = v.x;
	y = v.y;
	z = v.z;
	return *this;
}

const Coord& Coord::operator = (double v)
{
	x = y = z = v;
	return *this;
}

/* 四则运算 simple mathamatic operator */
const Coord Coord::operator + (const Coord& v) const
{
	return Coord(x + v.x, y + v.y, z + v.z);
}

const Coord Coord::operator - (const Coord& v) const
{
	return Coord(x - v.x, y - v.y, z - v.z);
}

const Coord Coord::operator * (double s) const
{
	return Coord(x *s, y*s, z *s);
}
 
const Coord Coord::operator / (double s) const
{
	return Coord(x / s, y / s, z / s);
}


const Coord& Coord::operator += (const Coord& v)
{
	x += v.x;
	y += v.y;
	z += v.z;
	return *this;
}

const Coord& Coord::operator -= (const Coord& v)
{
	x -= v.x;
	y -= v.y;
	z -= v.z;
	return *this;
}

const Coord& Coord::operator += (double delta)
{
	x += delta;
	y += delta;
	z += delta;
	return *this;
}

const Coord& Coord::operator -= (double delta)
{
	x -= delta;
	y -= delta;
	z -= delta;
	return *this;
}
bool Coord::operator <= (Coord& v) {
	return (x <= v.x) && (y <= v.y) && (z <= v.z);
}
bool Coord::operator <= (const Coord& v)
{
	return (x <= v.x) && (y <= v.y) && (z <= v.z);
}
bool Coord::operator >= (Coord& v)
{
	return (x >= v.x) && (y >= v.y) && (z >= v.z);
}
bool Coord::operator >= (const Coord& v)
{
	return (x >= v.x) && (y >= v.y) && (z >= v.z);
}

/* 点积 & 叉积 dot & cross */
double Coord::operator * (const Coord& v) const
{
	return x * v.x + y * v.y + z * v.z;
}

const Coord Coord::operator ^ (const Coord& v) const
{
	Coord vr;
    vr.x = y * v.z - v.y * z;
	vr.y = v.x * z - x * v.z;
	vr.z = x * v.y - y * v.x;
	return vr;
}


const Coord operator - (const Coord& v)
{
	Coord vr;
	vr.x = -v.x;
	vr.y = -v.y;
	vr.z = -v.z;
	return vr;
}

const Coord operator * (double scale, const Coord& v)
{
	 return v * scale;
}

