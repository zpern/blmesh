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
 #include "iso3d_utility.h"
#include <math.h>
#include <assert.h>
#include <time.h>
#include "matrix.h"
#include "spr.h"

/*
 * 两点距离平方
 * square of distance between two points
 */
REAL squaDist(MYPOINT p1, MYPOINT p2)
{
	REAL s = 0.0;
	int i;
	for (i = 0; i < DIM; i++)
		s += (p2[i] - p1[i]) * (p2[i] - p1[i]);
	return s;
}

/*
 * 三点面积绝对值
 * absolute value of area of a triangle 
 */
REAL areas(REAL xa, REAL ya, REAL za, 
		   REAL xb, REAL yb, REAL zb,
		   REAL xc, REAL yc, REAL zc)
{
	REAL aa = 0.0;
	REAL pq1, pq2, pq3, pr1, pr2, pr3, aq1, aq2, aq3;

	pq1     =   xb-xa;
	pq2     =   yb-ya;
	pq3     =   zb-za;
	pr1     =   xc-xa;
	pr2     =   yc-ya;
	pr3     =   zc-za;
	aq1     =   pq2*pr3-pq3*pr2;
	aq2     = -(pq1*pr3-pq3*pr1);
	aq3     =   pq1*pr2-pq2*pr1;
	aa      =   0.5*sqrt(aq1*aq1+aq2*aq2+aq3*aq3);
	return aa;
}

/* 
 * 表面右手法向量
 * normal of a face based on the right-hand rule
 */
int facNorm(MYPOINT p1, MYPOINT p2, MYPOINT p3, VECTOR norm)
{
	VECTOR v1, v2;
	int m;

	for (m = 0; m < DIM; m++)
	{
		v1[m] = p3[m] - p2[m];
		v2[m] = p1[m] - p2[m];
	}
	norm[0] = v1[1] * v2[2] - v1[2] * v2[1];
	norm[1] = v1[2] * v2[0] - v1[0] * v2[2];
	norm[2] = v1[0] * v2[1] - v1[1] * v2[0];
	return 1;
}

/*
 * 单位化向量
 * normalize a vector 
 */
int normVect(VECTOR v)
{
	int m;
	REAL mag = 0.0;

	for (m = 0; m < DIM; m++)
		mag += v[m] * v[m];
	mag = sqrt(mag);

	if(mag >= EPS_ABS_ZERO)	
    {
		for (m = 0; m < DIM; m++)
			v[m] /= mag;
	  return 1;
    }
	else
	{
		mag = mag;
	}

	return 0;
}

/*
 * 向量叉乘
 * cross operation of two vectors
 */
int crossVect(VECTOR v1, VECTOR v2, VECTOR *v) /* v = v1^v2 */
{
	assert(v);

	(*v)[0]  = v1[1] * v2[2] - v1[2] * v2[1];
	(*v)[1]  = v1[2] * v2[0] - v1[0] * v2[2];
	(*v)[2]  = v1[0] * v2[1] - v1[1] * v2[0];

	return 1;
}

/*
 * 比较两个向量的方向 compare directions of two vectors
 * -1 未知 undefined
 * 0  反向 of same directions
 * 1  同向 of opposite directions
 */
int compVectDir(VECTOR v1, VECTOR v2)
{
	int m;
	if (normVect(v1) && normVect(v2))
	{
		for (m = 0; m < DIM; m++)
		{
			if (fabs(v1[m]) > EPS_ZERO)
			{
				if ((v1[m] > 0.0 && v2[m] > 0.0) ||
					(v1[m] < 0.0 && v2[m] < 0.0))
					return 1;
				else
					return 0;
			}
		}
	}
	return -1;
}

/* 
 * 获取将给定向量变换为Z轴（0，0，1）所需的变换矩阵(4×4)
 * get the transform matrix to convert vector v to (0, 0, 1) (Z axis)
 */
int vectorToZ(ISOVector n, double m[16])
{
	double Rx[] = 
	{
		1, 0, 0, 0, 0, 1, 0, 0,
		0, 0, 1, 0, 0, 0, 0, 1
	}; /* default value is identical */
	double Ry[] = 
	{
		0, 0, 0, 0, 0, 1, 0, 0, 
		0, 0, 0, 0, 0, 0, 0, 1
	};
	double v, cos_alpha, sin_alpha, cos_beta, sin_beta;

#ifdef _PDMG_TEST
	Vector vtest;
#endif /* _TEST */

	n.normalize();
	v = sqrt(n.y * n.y + n.z * n.z);

	/* v等于0时Rx是单位矩阵, Rx is identical while v is zero */
	if (fabs(v) >= EPS_ZERO)
	{
		cos_alpha = n.z / v;
		sin_alpha = n.y / v;
		Rx[5] = Rx[10] = cos_alpha;
		Rx[6] = sin_alpha;
		Rx[9] = -sin_alpha;
	}

	cos_beta = v;
	sin_beta = -n.x;
	Ry[0] = Ry[10] = cos_beta;
	Ry[2] = -sin_beta;
	Ry[8] = sin_beta;

	matrixMultiply(m, Rx, Ry);

#ifdef _PDMG_TEST
	vtest = vectorTransform(n, m);
	printf("TEST FUNC vectorToZ: (%f, %f, %f)\n", vtest.x, vtest.y, vtest.z);
#endif /* _TEST */

	return 1;
}

/* 
 * p4点对三角面p1p2p3是否可视，即p1p2p3p4的体积是否为正
 * check if p4 is visible to p1p2p3, i.e. the volume of p1p2p3p4 is positive 
 */
bool isVisible(MYPOINT p1, MYPOINT p2, MYPOINT p3, MYPOINT p4)
{
#if 1
	REAL aa, bb, cc, id1, id2, d3, d4, dd, hh, al, ee, ff, gg, ai, aj, ak;
 	REAL s1, s2, s3, s4, s5, s6;
 	REAL te;
 	aa         =  2 * (p2[0] - p1[0]); 
 	bb         =  2 * (p2[1] - p1[1]); 
 	cc         =  2 * (p2[2] - p1[2]); 
 	id1         =  p1[0] * p1[0] + p1[1] * p1[1] + p1[2] * p1[2];
 	id2         =  p2[0] * p2[0] + p2[1] * p2[1] + p2[2] * p2[2];
 	d3         =   p3[0] * p3[0] + p3[1] * p3[1] + p3[2] * p3[2];
 	d4         =   p4[0] * p4[0] + p4[1] * p4[1] + p4[2] * p4[2];
 	dd         =  id2-id1;
 	hh         =  d3-id1;
 	al         =  d4-id1;
 	ee         =  2 * (p3[0] - p1[0]);
 	ff         =  2 * (p3[1] - p1[1]);
 	gg         =  2 * (p3[2] - p1[2]);
 	ai         =  2 * (p4[0] - p1[0]);
 	aj         =  2 * (p4[1] - p1[1]);
 	ak         =  2 * (p4[2] - p1[2]);
 	s1         =  ff*ak-aj*gg;
 	s2         =  ee*ak-ai*gg;
 	s3         =  ee*aj-ai*ff;
 	s4         =  hh*ak-al*gg;
 	s5         =  hh*aj-al*ff;
 	s6         =  ee*al-hh*ai;	
	te         =  aa*s1-bb*s2+cc*s3;
#else
	REAL te = -GEOM_FUNC::orient3d(p1,p2,p3,p4);
#endif 

	return te > 0.0;
}

/* 
 * ia ib ic id 是四面体单元的四个点的点码值，判断四面体单元
 * A->B->C->D是否有体积(volume-positive)
 * ia, ib, ic and id are four node codes of four forming points of Tetrahedra 
 * ABCD, judge if the volume of ABCD is positive
 */
bool isVolumePositive(int ia, int ib, int ic, int id)
{
	if (((ia == 0) && (ib == 1) && (ic == 2) && (id == 3)) ||
		((ia == 0) && (ib == 2) && (ic == 3) && (id == 1)) ||
		((ia == 0) && (ib == 3) && (ic == 1) && (id == 2)) ||
		((ia == 1) && (ib == 2) && (ic == 0) && (id == 3)) ||
		((ia == 1) && (ib == 3) && (ic == 2) && (id == 0)) ||
		((ia == 2) && (ib == 3) && (ic == 0) && (id == 1)) ||
		((ia == 1) && (ib == 0) && (ic == 3) && (id == 2)) ||
		((ia == 2) && (ib == 0) && (ic == 1) && (id == 3)) ||
		((ia == 3) && (ib == 0) && (ic == 2) && (id == 1)) ||
		((ia == 2) && (ib == 1) && (ic == 3) && (id == 0)) ||
		((ia == 3) && (ib == 1) && (ic == 0) && (id == 2)) ||
		((ia == 3) && (ib == 2) && (ic == 1) && (id == 0)))
		return true;
	else
		return false;
		
}

#if 0
int generateRandArray(int *arr, int maxValue)
{
	int i, x, sum1, sum2;

	memset(arr, 0, sizeof(INTEGER)*maxValue);
	srand(time(0)); /* 选取种子，避免伪随机 */

	for (i = 1;i <= maxValue; i++)
	{
		while(arr[x=rand()%maxValue])
			;
		arr[x]=i;
	}

#ifdef _DEBUG
	/* ----------------------------------------------------
	 * 检查随机数对不对 
	 * --------------------------------------------------*/
	for (i = 0, sum1 = 0, sum2 = 0; i < maxValue; i++)
	{
		assert(arr[i] > 0 && arr[i] <= maxValue);
		sum1 += i + 1;
		sum2 += arr[i]; 
	}
	assert(sum1 == sum2);
#endif

	return 1;
}
#else

int generateRandArray(int *arr, int maxValue)
{
	int i, x, sum1, sum2, idx0, idx1, tmp;

//	srand(time(0)); /* 选取种子，避免伪随机 */
	for (i = 0; i < maxValue; i++)
		arr[i] = i+1;

    for (i = 0; i < maxValue; i++)
    {
		idx0 = rand()%(maxValue - i);
        idx1 = maxValue - i - 1;

		tmp = arr[idx0];
        arr[idx0] = arr[idx1];
        arr[idx1] = tmp;
    }

#ifdef _DEBUG
	/* ----------------------------------------------------
	 * 检查随机数对不对 
	 * --------------------------------------------------*/
	for (i = 0, sum1 = 0, sum2 = 0; i < maxValue; i++)
	{
		assert(arr[i] > 0 && arr[i] <= maxValue);
		sum1 += i + 1;
		sum2 += arr[i]; 
	}
	assert(sum1 == sum2);
#endif

	return 1;

}

#endif