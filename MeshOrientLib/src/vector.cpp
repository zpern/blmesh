
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
 * 修改者 (Revisor): yhf
 * 修改日期 (Revision Date):20190823
 * 当前版本 (Current Version):
 * 修改介绍 (Revision Introduction):
 * ------------------------------------------------------------------------------
 * ------------------------------------------------------------------------------*/

#include <spdlog/spdlog.h> 
 #include "vector.h"
#include <cmath>


 /* *********************************************************************************
  * 定义三维矢量
  * define 3-dimensional vector
  * **********************************************************************************/


  /* 矢量求模 calc. magnitude */
double BLVector::magnitude() const
{
	return sqrt(x * x + y * y + z * z);
}

double BLVector::length() const
{
	return magnitude();
}

bool BLVector::isApproximatelyEqualTo(const BLVector v1)
{
	return (abs(x - v1.x) < length() * (1e-6)) && (abs(y - v1.y) < length() * (1e-6)) && (abs(z - v1.z) < length() * (1e-6));
}

/* 矢量求模 calc. magnitude */
double BLVector::magnitude2() const
{
	return (x * x + y * y +
		z * z);
}

double BLVector::dotProduct(const BLVector v1, const BLVector v2)
{
	return v1 * v2;
}

BLVector BLVector::crossProduct(const BLVector v1, const BLVector v2)
{
	return v1 ^ v2;
}
BLVector BLVector::normalized() {
	double m = magnitude();
	if (m == 0)
	{
		return BLVector(0, 0, 0);
	}
	else {
		return BLVector(x / m, y / m, z / m);
	}
}
/* 归一 normalization */
void BLVector::normalize()
{
	double m = magnitude();
	if (abs(m) < 1e-10)//修改浮点数等于0为绝对值小于误差 yhf20190823
	{
		x = 0.0;
		y = 0.0;
		z = 0.0;
	}
	else {
		x = x / m;
		y = y / m;
		z = z / m;
	}
}

/* **********************************
 * 重载操作符 overloaded operator
 * **********************************/

 /* 赋值 assignment */
const BLVector& BLVector::operator = (const BLVector& v)
{
	x = v.x;
	y = v.y;
	z = v.z;
	return *this;
}

const BLVector& BLVector::operator = (double v)
{
	x = y = z = v;
	return *this;
}

/* 四则运算 simple mathamatic operator */
BLVector BLVector::operator + (const BLVector& v) const
{
	BLVector vr;
	vr.x = x + v.x;
	vr.y = y + v.y;
	vr.z = z + v.z;
	return vr;
}

BLVector BLVector::operator - (const BLVector& v) const
{
	BLVector vr;
	vr.x = x - v.x;
	vr.y = y - v.y;
	vr.z = z - v.z;
	return vr;
}

BLVector BLVector::operator * (double s) const
{
	BLVector vr;
	vr.x = x * s;
	vr.y = y * s;
	vr.z = z * s;
	return vr;
}
BLVector BLVector::operator / (double s) const
{
	BLVector vr;
	vr.x = x / s;
	vr.y = y / s;
	vr.z = z / s;
	return vr;
}


BLVector& BLVector::operator += (const BLVector& v)
{
	x += v.x;
	y += v.y;
	z += v.z;
	return *this;
}
BLVector& BLVector::operator -= (const BLVector& v)
{
	x -= v.x;
	y -= v.y;
	z -= v.z;
	return *this;
}

BLVector& BLVector::operator += (double delta)
{
	x += delta;
	y += delta;
	z += delta;
	return *this;
}
BLVector& BLVector::operator -= (double delta)
{
	x -= delta;
	y -= delta;
	z -= delta;
	return *this;
}


/* 点积 & 叉积 dot & cross */
double BLVector::operator * (const BLVector& v) const
{
	return x * v.x + y * v.y + z * v.z;
}

BLVector BLVector::operator ^ (const BLVector& v) const
{
	BLVector vr;
	vr.x = y * v.z - v.y * z;
	vr.y = v.x * z - x * v.z;
	vr.z = x * v.y - y * v.x;
	return vr;
}


BLVector operator - (const BLVector& v)
{
	BLVector vr;
	vr.x = -v.x;
	vr.y = -v.y;
	vr.z = -v.z;
	return vr;
}
BLVector operator * (double scale, const BLVector& v)
{
	return v * scale;
}

