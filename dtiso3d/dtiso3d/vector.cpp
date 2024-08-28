/* ----------------------------------------------------------------------------
 * ----------------------------------------------------------------------------
 *
 * 三维各向同性Delaunay网格生成器 (版本号：0.3)
 * 3D Isotropic Delaunay Mesh Generation (Version 0.3)
 *
 * 陈建军 中国 浙江大学工程与科学计算研究中心
 * 版权所有	  2005年9月15日
 * Chen Jianjun  Center for Engineering & Scientific Computation,
 * Zhejiang University, P. R. China
 * Copyright reserved, 2005, 10, 26
 * 
 * 联系方式
 *   电话：+86-571-87953165
 *   传真：+86-571-87953167
 *   邮箱：zdchenjj@yahoo.com.cn
 * For further information, please conctact
 *  Tel: +86-571-87953165
 *  Fax: +86-571-87953167
 * Mail: zdchenjj@yahoo.com.cn
 *
 * ------------------------------------------------------------------------------
 * ------------------------------------------------------------------------------*/

#include <spdlog/spdlog.h> 
 #include "vector.h"
#include <math.h>

/* *********************************************************************************
 * 定义三维矢量
 * define 3-dimensional vector
 * **********************************************************************************/

/* 矢量求模 calc. magnitude */
double ISOVector::magnitude() const
{
	return sqrt(x * x + y * y + z * z);
}

/* 归一 normalization */
void ISOVector::normalize()
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
const ISOVector& ISOVector::operator = (const ISOVector& v)
{
	x = v.x;
	y = v.y;
	z = v.z;
	return *this;
}

const ISOVector& ISOVector::operator = (double v)
{
	x = y = z = v;
	return *this;
}

/* 四则运算 simple mathamatic operator */
const ISOVector ISOVector::operator + (const ISOVector& v) const
{
	ISOVector vr;
	vr.x = x + v.x;
	vr.y = y + v.y;
	vr.z = z + v.z;
	return vr;
}

const ISOVector ISOVector::operator - (const ISOVector& v) const
{
	ISOVector vr;
	vr.x = x - v.x;
	vr.y = y - v.y;
	vr.z = z - v.z;
	return vr;
}

const ISOVector ISOVector::operator * (double s) const
{
	ISOVector vr;
	vr.x = x * s;
	vr.y = y * s;
	vr.z = z * s;
	return vr;
}
 
const ISOVector ISOVector::operator / (double s) const
{
	ISOVector vr;
	vr.x = x / s;
	vr.y = y / s;
	vr.z = z / s;
	return vr;
}


const ISOVector& ISOVector::operator += (const ISOVector& v)
{
	x += v.x;
	y += v.y;
	z += v.y;
	return *this;
}

const ISOVector& ISOVector::operator -= (const ISOVector& v)
{
	x -= v.x;
	y -= v.y;
	z -= v.y;
	return *this;
}

const ISOVector& ISOVector::operator += (double delta)
{
	x += delta;
	y += delta;
	z += delta;
	return *this;
}

const ISOVector& ISOVector::operator -= (double delta)
{
	x -= delta;
	y -= delta;
	z -= delta;
	return *this;
}


/* 点积 & 叉积 dot & cross */
double ISOVector::operator * (const ISOVector& v) const
{
	return x * v.x + y * v.y + z * v.z;
}

const ISOVector ISOVector::operator ^ (const ISOVector& v) const
{
	ISOVector vr;
    vr.x = y * v.z - v.y * z;
	vr.y = v.x * z - x * v.z;
	vr.z = x * v.y - y * v.x;
	return vr;
}


const ISOVector operator - (const ISOVector& v)
{
	ISOVector vr;
	vr.x = -v.x;
	vr.y = -v.y;
	vr.z = -v.z;
	return vr;
}

const ISOVector operator * (double scale, const ISOVector& v)
{
	 return v * scale;
}

// Cross products.
// vector1 cross vector2
inline ISOVector operator*(const ISOVector &vector1, const ISOVector &vector2)
{
	//return Vector(vector1) ^= vector2;
}