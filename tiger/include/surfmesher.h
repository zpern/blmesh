#pragma once
/** ----------------------------------------------------------------------------
	* TIGER, a Trustable Intelligent GridER
	*
	* TIGER的曲面网格生成API集成接口定义 (版本号：1.0)
	* Definitions of APIs for Surface Mesh Generation (Version 1.0)
	*
	* 陈建军 中国 浙江大学航空航天学院
	* 版权所有	  2020年2月6日
	* Chen Jianjun  School of Aeronautics & Astronautics,
	* Zhejiang University, P. R. China
	* Copyright reserved, Feb 6th, 2020
	*
	* 联系方式
	*   电话：+86-571-87951883
	*   传真：+86-571-87953167
	*   邮箱：chenjj@zju.edu.cn
	* For further information, please conctact
	*  Tel: +86-571-87951883
	*  Fax: +86-571-87953167
	*  Mail: chenjj@zju.edu.cn
	*
	* ------------------------------------------------------------------------------*/
#ifndef __tiger_surfmesher_h__
#define __tiger_surfmesher_h__
#ifdef WIN32
#ifdef API_main_EXPORTS 
#define DECL_VOLTET __declspec(dllexport)
#else
#define DECL_VOLTET __declspec(dllimport)
#endif API_main_EXPORTS
#endif WIN32
	/** -------------------------------------------------------------------------------
		* Given a curve and a user defined sizing function for element spacing, try to
		* discretize the curve into a sequence of a linear segments
		*
		* The inputs include two parts:
		* The first part is the curve defined by Furgeson interpolations. See Chapter 3
		* of the following references for more details.
		*
		* 郑耀，陈建军，著. 非结构网格生成：理论、算法和应用. 科学出版社, 2016. ISBN: 978-7-03-047283-0.
	  * Yao Zheng, Jianjun Chen, Unstructured Mesh Generation: Theories, Algorithms and Applications (in Chinese),
		* the 1st Edition, ISBN 978-7-03-047283-0, Science Press, Beijing, China, March, 2016 (273 Pages).
		*
		* The second part of the input is a user defined sizing function, which is expected to
		* output a sizing value for an arbitrary point inside the problem domain
		*
		* The outputs are a sequence of mesh points defining the resulting linear segments
		*
		* Input Paramerters:
		*	cvPtNum			# of points interpolating the curve.
		*	cvPts				the coordinates of curve points. dim = 3*(cvPtNums[0] + ....+cvPtNums[cvNums-1])
		* sizingFunc	the user specified sizing function
		* Output Paramerters:
		*	mshPtNum		# of mesh points
		*	mshPtUs			the parametric coordinates of mesh points. dim = mshPtNum
		*	mshPtXYZs		the physical coordinates of mesh points. dim = 3 * mshPtNum
		* Returned value:
		*	0				Completely succeed
		*	1				Partially succeed
		*	otherwise		Fail. Error info. will be defined later
		* ------------------------------------------------------------------------------*/
DECL_VOLTET int API_Discretize_Single_Curve(
	int cvPtNum, double cvPts[],
	double(*sizingFunc)(double x, double y, double z),
	int *mshPtNum, double **mshPtUs, double **mshPtXYZs
);


/** -------------------------------------------------------------------------------
	* Given a surface and a user defined sizing function for element spacing, try to
	* discretize the surface into a combination of triangles
	*
	* The inputs include three parts:
	* The first part is the surface defined by Furgeson interpolations. See Chapter 3
	* of the following references for more details.
	*
	* 郑耀，陈建军，著. 非结构网格生成：理论、算法和应用. 科学出版社, 2016. ISBN: 978-7-03-047283-0.
  * Yao Zheng, Jianjun Chen, Unstructured Mesh Generation: Theories, Algorithms and Applications (in Chinese),
	* the 1st Edition, ISBN 978-7-03-047283-0, Science Press, Beijing, China, March, 2016 (273 Pages).
	*
	* Meanwhile, a group of curves trimming the surface are given to define the boundary
	* of the surface.
	*
	* The second part of the input is the discretization of boundary curves (if any).
	*
	* The third part of the input is a user defined sizing function, which is expected to
	* output a sizing value for an arbitrary point inside the problem domain
	*
	* The outputs are a sequence of mesh points defining the resulting triangular mesh &
	* the discretization results of boundary curves (if necessary)
	*
	* Input Paramerters:
	* cvNum					# of boundary curves of the surface
	*	cvPtNums			# of points interpolating a curve. dim = cvNums
	*	cvPts					the coordinates of curve points. dim = 3*(cvPtNums[0] + ....+cvPtNums[cvNum-1])
	* fcPtNums			# of points along U & V directions interpolating the surface. dim = 2
	*	fcPts					the coordinates of curve points. dim = 3*(fcPtNums[0]*fcPtNum[1])
	* cvMshPtNums 	# of mesh points of total boundary curve.
	*	cvMshPts			coordinates of mesh points along each boundary curve (if any).
	*								dim = 3*cvMshPtNums
	* cvElmNums 	# of elements along each boundary curve (if any). dim = cvNum.
	* cvElms		# of points index along each boundary curve (if any). dim = (cvElmNums[0]+1)+ ...+(cvElmNums[cvNum-1]+1)
	* sizingFunc	the user specified sizing function
	*
	* Output Paramerters:
	*	mshPtNum			# of mesh points
	*	mshPtUVs			the parametric coordinates of mesh points. dim = 2*mshPtNum
	*	mshPtXYZs			the physical coordinates of mesh points. dim = 3 * mshPtNum
	*								dim = 3*mshPtNum
	*	sufElmNum			# of surface elements
	*	sufElms				the connections of surface elements. dim = 3*sufElmNum
	*
	* Returned value:
	*	0				Completely succeed
	*	1				Partially succeed
	*	otherwise		Fail. Error info. will be defined later
	* ------------------------------------------------------------------------------*/
DECL_VOLTET int API_Triangulize_Single_Surface(
	int cvNum, /*int cvPtNums[], double cvPts[],*/
	int fcPtNums[], double fcPts[],
	int cvMshPtNums,
	double cvMshPts[],
	int cvElmNums[],
	int cvElms[],
	double(*sizingFunc)(double x, double y, double z),
	int *mshPtNum, double **mshPtUVs, double **mshPtXYZs,
	int *sufElmNum, int **sufElms
);

DECL_VOLTET void APIfreeCV(
	double** PhyResult,
	double** LnResult
);

DECL_VOLTET void APIfreeSF(
	double** uv,
	double** xyz,
	int** topu
);

#endif /* __tiger_surfmesher_h__ */
