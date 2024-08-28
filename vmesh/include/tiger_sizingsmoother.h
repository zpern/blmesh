/** ----------------------------------------------------------------------------
	* TIGER, a Trustable Intelligent GridER
	*
	* TIGER的非结构网格单元尺寸控制函数优化器接口定义 (版本号：1.0)
	* Definitions of APIs for Smoothers of Sizing Functions of Unstructured Grid Generation (Version 1.0) 
	*
	* 陈建军 中国 浙江大学航空航天学院
	* 版权所有	  2020年3月5日
	* Chen Jianjun  School of Aeronautics & Astronautics,
	* Zhejiang University, P. R. China
	* Copyright reserved, Mar 5th, 2020
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
#ifndef __tiger_sizingsmoother_h__
#define __tiger_sizingsmoother_h__

/** -------------------------------------------------------------------------------
	* Given a sizing function defined by a tetrahedral mesh, try to smooth it
	*	NOTE: 
	* 1. no smoothing is done
	* 2. linear function is used in default
	*
	* Input Paramerters:
	*	mshPtNums		# of mesh points
	*	mshPts			the coordinates of mesh points. dim = 3*mshPtNums
	*	volElmNum		# of volume elements
	*	volElms			the connections of volume elements. dim = 4*volElmNum
	* betas				the allowed maximal gradation values dim = volElmNum
	* Input/Output Parameters
	*	ptValues		sizing values at mesh points	
	* Returned value:
	*	0				Completely succeed
	*	1				Partially succeed
	*	otherwise		Fail. Error info. will be defined later
	* ------------------------------------------------------------------------------*/
int API_Smooth_VolmBKG_SF(
	int mshPtNum, double mshPts[], double ptValues[],
	int volElmNum, int volElms[],	double betas[]
	);

/** -------------------------------------------------------------------------------
	* Given a sizing function defined by a triangular mesh, try to create a sizing 
	* function to record it
	*
	*	NOTE: no smoothing is done
	* 1. no smoothing is done
	* 2. linear function is used in default
	*
	* Input Paramerters:
	* bndPtNum			# of boundary points of the tessellation 
	*	bndPts				the coordinates of boundary points. dim = 3*bndPtNum
	*	bndFctNum			# of boundary facets of the tessellation 				
	* bndFcts				the connections of boundary faces. dim = 3*bndFctNum
	* betas				the allowed maximal gradation values dim = bndFctNum
	* Input/Output Parameters
	*	ptValues		sizing values at mesh points	
	* Returned value:
	*	0				Completely succeed
	*	1				Partially succeed
	*	otherwise		Fail. Error info. will be defined later
	* ------------------------------------------------------------------------------*/
int API_Smooth_SurfBKG_SF(
	int bndPtNum, double bndPts[], double ptValues[], 
	int bndFctNum, int bndFcts[], double beltas[]
	);
	
/** -------------------------------------------------------------------------------
	* Given a sizing function defined by a linear mesh, try to create a sizing 
	* function to record it
	*
	*	NOTE: no smoothing is done
	* 1. no smoothing is done
	* 2. linear function is used in default
	*
	* Input Paramerters:
	* bndPtNum			# of boundary points of the tessellation 
	*	bndPts				the coordinates of boundary points. dim = 3*bndPtNum
	*	bndSegNum			# of linear segments of the tessellation 				
	* bndSegs				the connections of linear segements. dim = 2*bndSegNum
	* betas				the allowed maximal gradation values dim = bndSegNum
	* Input/Output Parameters
	*	ptValues		sizing values at mesh points	
	* Returned value:
	*	0				Completely succeed
	*	1				Partially succeed
	*	otherwise		Fail. Error info. will be defined later
	* ------------------------------------------------------------------------------*/
int API_Smooth_CurveBKG_SF(
	int bndPtNum, double bndPts[], double ptValues[], 
	int bndSegNum, int bndSegs[], double beltas
	);
	
#endif /* __tiger_sizingsmoother_h__ */