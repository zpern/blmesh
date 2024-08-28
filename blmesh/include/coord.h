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

#ifndef __coord_h__
#define	__coord_h__


/* *********************************************************************************
 * 定义三维矢量
 * define 3-dimensional vector
 * **********************************************************************************/
class Coord
{
public:
	Coord() { x = y = z = 0.; }
	/*这里应该用成员初始化列表加快速度yhf*/
	Coord(double v1, double v2, double v3):x(v1),y(v2),z(v3) {
	}
	Coord(const Coord& v) {
		x = v.x;  y = v.y; z = v.z;
	}
	
	/* 矢量求模 calc. magnitude */
	double magnitude() const;
	
	/* 归一 normalization */
	void normalize();

	/* **********************************
	 * 重载操作符 overloaded operator 
	 * **********************************/

	/* 赋值 assignment */
	const Coord& operator = (const Coord& v);
	const Coord& operator = (double v);

	/* 四则运算 simple mathamatic operator */
    const Coord operator + (const Coord& v) const;
    const Coord operator - (const Coord& v) const;
    const Coord operator * (double s) const;  
    const Coord operator / (double s) const;
    
    const Coord& operator += (const Coord& v);
    const Coord& operator -= (const Coord& v);
    const Coord& operator += (double delta);
    const Coord& operator -= (double delta);

	/* 关系运算 operator */
    bool operator <= (const Coord& v);
    bool operator >= (const Coord& v);

	bool operator <= (Coord& v);
	bool operator >= (Coord& v);

	/* 点积 & 叉积 dot & cross */
	double operator * (const Coord& v) const;
	const Coord operator ^ (const Coord& V) const;

	double operator [](int &index) {
		switch (index)
		{
			case 0:
				return x;
				break;
			case 1:
				return y;
		default:
			return z;
			break;
		}
	}
	friend const Coord operator - (const Coord& p);
	friend const Coord operator * (double scale, const Coord& v);

public:
	double x, y, z;
};

#endif /* __coord_h__*/