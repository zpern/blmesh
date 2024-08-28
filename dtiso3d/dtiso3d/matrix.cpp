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
 #include "matrix.h"
#include <math.h>

/* 
 * 矩阵复制 m = a
 * copy a matrix m = a
 */
void matrixCopy(double *m, const double *a, const int row, const int col)
{
	int i;

	for(i = 0; i <  row * col; i++)
		m[i] = a[i];
}

/* 
 * 矩阵相加 m = a + b
 * add two matrices m = a + b
 */
void matrixAdd(double *m, const double *a, const double *b,
		   const int row, const int col)
{
	int i;

	for(i = 0; i < row * col; i++)
		m[i] = a[i] + b[i];
}

/* 
 * 矩阵相减 m = a - b
 * subtract two matrices m = a - b
 */
void matrixSub(double *m, const double *a, const double *b,
		   const int row, const int col)
{
	int i;

	for(i = 0; i < row * col; i++)
		m[i] = a[i] - b[i];
} 

/* 
 * 矩阵相乘 m = a * b 仅当col1 = row2时有效
 * multiply two matrices m = a * b ( valid only while col1 = row2)
 */
bool matrixMultiply(double *m, const double *a, const double *b, 
			    const int  row1, const int  col1, 
			    const int  row2, const int  col2)
{
	int i, i1, j, k, k1;
	int size = row1 * col2;

	if (col1 != row2) return false;

	for (i = 0, i1 = 0; i < size; i += col2, i1 += col1)
		for (j = 0; j < col2; j++)
		{
			m[i+j] = 0.0;

			for (k = 0, k1 = 0; k < col1; k++, k1 += col2)
				m[i+j] += a[i1+k] * b[k1+j];
		}    

	return true;
}

/* 
 * 矩阵转置 (m = transpose(a))
 * transpose a matrix (m = transpose(a))
 */
void matrixTranspose(double *m, const double *a, const int row, const int col)
{
	int i, ii, j, jj;
	
	for(i = 0, ii = 0; i < row; i++, ii += row)
		for(j = 0, jj = 0; j < col; j++, jj += col)
			m[ii+j] = a[jj+i];
}

/* 
 * 矩阵转置 (m = transpose(m))
 * transpose a matrix (m = transpose(m))
 */
void matrixTranspose(double *m, const int row, const int col)
{
	double *temp = new double[row * col];
	int i, ii, j, jj;
	
	for (i = 0, ii = 0; i < row; i++,ii += row)
		for (j = 0, jj = 0; j < col; j++, jj += col)
			 temp[ii+j] = m[jj+i];
	
	for (i = 0; i < row * col; i++)  
		m[i]=temp[i];	

	delete[] temp;
	return;
}

/* 
 * 矩阵求逆用的 LUD 函数 
 * 输入参数：data(n, n)
 * 输出参数：data(n, n)  index(n, n)
 * LUD function for invert-matrix computation
 * Input: data(n, n)
 * Output: data(n, n) index(n, n)
 */
static bool LUDcmp( double *data, int *index, int n)
{
	const double TINY = 1.0e-20f;
	double *vv = new double [n];
	double temp, big;
	int i, ii, j, jj, k, kk, imax = 0;

	for(i=0, ii=0; i<n; i++, ii+=n)
	{
		big=0.0f;
		for(j=0; j<n; j++)
 		if((temp = double(fabs(data[ii+j])) ) > big)  big=temp;
       
		if( big == 0.0f )
		{
			delete vv;
			return false;	// FALSE
		}
     
		vv[i] = 1.0f/big;
	}

	for(j=0, jj=0; j<n; j++, jj+= n)
	{
		for(i=0, ii=0; i<=j; i++, ii+= n)
		{
			temp = data[ii+j];
			for( int k=0, kk=0; k<i; k++, kk+= n)
			temp -= data[ii+k] * data[kk+j];
			data[ii+j] = temp;
		}
		
		big=0.0f;
		imax=j;
		
		for(i=j+1, ii=i*n; i<n; i++, ii+= n)
		{
			temp = data[ii+j];
			for(k=0, kk=0; k<j; k++, kk+=n)
				temp -= data[ii+k] * data[kk+j];
			data[ii+j] = temp;
			if(( temp = vv[i] * double(fabs(temp)) ) > big)
			{
				big=temp;
				imax= i;
			}
		}
		
		if( j != imax)
		{
			ii = imax * n;
			for(k=0; k<n; k++)
			{
				temp = data[ii+k];
				data[ii+k] = data[jj+k];
				data[jj+k] = temp;
			}
			vv[imax] = vv[j];
		}
		
		index[j] = imax;
		if(data[jj+j] == 0.0f) data[jj+j] = TINY;
		if(j<n)
		{
			temp = 1.0f/data[jj+j];
			for(i=j+1, ii=i*n; i<n; i++, ii+=n)
				data[ii+j] *= temp;
		}
	}
	delete vv;

	return true;	
}

/* 
 * 矩阵求逆用的 LUD 函数 
 * 输入参数：data(n, n)  index(n, n)
 * 输出参数：b(n, n) 
 * LUD function for invert-matrix computation
 * Input: data(n, n)  index(n, n)
 * Output: b(n, n) 
 */
void LUBksb(const double *data, const int *index, double *b, const int n)
{
	int ii, i, j, k = -1, ip;
	double sum;
	for(i=0, ii=0; i<n; i++, ii+= n)
	{
		ip=index[i];
		sum=b[ip];
		b[ip]=b[i];
		if(k != -1)	 for(j=k; j<i; j++) sum -= data[ii+j]*b[j];
		else if(sum) k=i;
		b[i]=sum;
	}
	for(i=n-1, ii=i*n; i>=0; i--, ii-= n)
	{
		sum=b[i];
		for(j=i+1; j<n; j++) sum -= data[ii+j]*b[j];
		b[i]=sum/data[ii+i];
	}
}

/*
 * 矩阵求逆 m = invert(a)
 * true 表示成功，false 表示该逆矩阵不存在
 * invert a matrix m = invert(a)
 * true: successful; false: no inverse matrix
 */
bool matrixInverse(double *m, const double *a, const int n)
{
	int *index = new int[n];
	double *data = new double[n*n];
	int i, j, jj;

	for (i = 0; i < n*n; i++) 
		data[i] = a[i];  

	if (!LUDcmp(data, index, n))
	{
		 delete[] index;
		 delete[] data;
		 return false;	
	}
	
	double *b = new double[n];
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++) 
			b[j]=0.0f;
		b[i] = 1.0f;

		LUBksb(data, index, b, n);

		for (j = 0, jj = 0; j < n; j++, jj += n)  
			m[jj+i] = b[j];
	}


	delete[] b;
	delete[] index;
    delete[] data;

	return true;
}

/*
 * 矩阵求逆 m = invert(m)
 * true 表示成功，false 表示该逆矩阵不存在
 * invert a matrix m = invert(m)
 * true: successful; false: no inverse matrix
 */
bool matrixInverse(double *m, const int n)
{
	int *index = new int[n];
	double *data  = new double[n*n];
	int i, j, jj;

	for (i = 0; i < n*n; i++)  
		data[i] = m[i];  

	if (!LUDcmp(data, index, n))
	{
		delete[] index;
		delete[] data;
		return false;	
	}
	
	double *b = new double[n];
	double *inverse = new double[n*n];
	for (i = 0; i < n; i++)
	{
		for(j = 0; j < n; j++)  
			b[j] = 0.0f;
		b[i] = 1.0f;

		LUBksb(data, index, b, n);

		for (j = 0, jj = 0; j < n; j++, jj += n)  
			inverse[jj+i] = b[j];
	}

	for (i = 0; i < n*n; i++)  
		m[i] = inverse[i];  

	delete[] b;
	delete[] inverse;
	delete[] index;
	delete[] data;

	return true;
}

/*
 * 矩阵归一化 
 * identify a matrix
 */
void matrixIdentity(double *m, const int n)
{
	int i, j;

	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			m[i*n+j] = (i == j) ? 1.0 : 0.0;
}


/* v 经过 m 变换后为 v', 函数返回 v' */
ISOVector vectorTransform(const ISOVector& v, const double m[16])
{
	double  in[4]={ v.x, v.y, v.z, 1.0f };
	double  out[4];
	ISOVector p;
	matrixMultiply(out, in, m, 1, 4, 4, 4);
	p.x=out[0];
	p.y=out[1];
	p.z=out[2];
	return p;
}


