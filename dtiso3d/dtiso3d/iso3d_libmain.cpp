/* ----------------------------------------------------------------------------
 * ----------------------------------------------------------------------------
 *
 * 三维各向同性Delaunay网格生成器 (版本号：1.0Belta)
 * 3D Isotropic Delaunay Mesh Generation (Version 1.0Belta)
 *
 * 陈建军 中国 浙江大学工程与科学计算研究中心
 * 版权所有	  2012年12月30日
 * Chen Jianjun  Center for Engineering & Scientific Computation,
 * Zhejiang University, P. R. China
 * Copyright reserved, 2012, 12, 30
 * 
 * 联系方式
 *   电话：+86-571-87951883
 *   传真：+86-571-87953167
 *   邮箱：chenjj@zju.edu.cn
 * For further information, please conctact
 *  Tel: +86-571-87951883
 *  Fax: +86-571-87953167
 * Mail: chenjj@zju.edu.cn
 *
 * ------------------------------------------------------------------------------
 * ------------------------------------------------------------------------------*/

#include <stdio.h>
#include <assert.h>
#include <functional>
#include <array>
#include <spdlog/spdlog.h> 
 #include "iso3d.h"
#include "iso3d_common.h"
#include "iso3d_error.h"
#include "iso3d_utility.h"
#include "locsmoothing.h"
#include "mshgen3d_def.h"

#ifdef _TIMING_PERFORMANCE
#include "spr.h"
#include "time_perform.h"
#endif

#define MAX_FILE_LEN 1024
#define _NGB 1

static const char *CASE_NAME = "embed_dtiso3d";
/**
 * @brief 生成三维各向同性四面体网格
 * 
 */
int meshGen3D_memo(
	/* ------------------------- input arguments --------------------------*/
	double		*pdSNX,	/* x coord. of boundary nodes */
	double		*pdSNY,	/* y coord. of boundary nodes */
	double		*pdSNZ,	/* z coord. of boundary nodes */
	int			nSN,		/* number of boundary nodes */
	int			*pnSFFm,	/* forming points of surface facets */
	int			*pnSFPt,	/* patches of surface facets */
	int         nSF,		/* number of facets */
	double		*pdBMNX,	/* x coord. of background mesh */
	double		*pdBMNY,	/* y coord. of background mesh */
	double		*pdBMNZ,	/* z coord. of background mesh */
	double		*pdBMNSpc,	/* space values of background mesh nodes */
	int			nBMN,		/* number of background mesh nodes */
	int			*pnBMEFm,	/* forming points of background mesh elements */
	int			*pnBMENg,	/* neighboring eles. of background mesh elements */
	int			nBME,		/* number of background mesh elements */
	double		*pdCX,		/* x coord. of center points of sources */
	double		*pdCY,		/* y coord. of center points of sources */
	double		*pdCZ,		/* z coord. of center points of sources */
	double		*pdIn,		/* inner radius of sources */
	double		*pdOu,		/* out radius of sources */
	double		*pdSp,		/* space values of sources */
	int			nPS,		/* number of point source */
	int			nLS,		/* number of line source */
	int			nTS,		/* number of triangle source */
	 /*-----------------------size funtion ----------------------------------*/
	std::function < double(std::array<double,3>) > func ,
				   /* ------------------------- output arguments -------------------------*/
				   double	  **ppdMNX,		/* x coord. of 3D mesh nodes */
				   double    **ppdMNY,		/* y coord. of 3D mesh nodes */
				   double    **ppdMNZ,		/* z coord. of 3D mesh nodes */
				   int        *pnMN,		/* number of 3D mesh nodes */
				   int       **ppnMEFm,	/* forming points of 3D mesh elements */
				   int       **ppnMENg,	/* neighboring eles. of 3D mesh elements */
				   int        *pnME,		/* number of 3D mesh elements */
				   int        **ppnPrt,		/* parents of surface facets */
					DTIso3DConfig config
				   )
{
	DTIso3D iso3d;

	iso3d.SetSizingFuntion(func);
    // DTIso3DConfig config;
    // g_bMeshingEnabled = true;  /* 网格生成功能是否被启动 -M */
    // g_bQuaImprEnabled = true;	/* 网格优化功能是否被启动 -Q */

#ifdef _ERROR_CHK
	int nErrCd;
#endif /* _ERROR_CHK */
	char strPLS[MAX_FILE_LEN], strBA3[MAX_FILE_LEN], strPL3[MAX_FILE_LEN], strQual[MAX_FILE_LEN];
#ifdef _NGB
	char strNGB[MAX_FILE_LEN];
#endif 

	/*
	 * 获取机器精度
	 */
	GEOM_FUNC::exactinit_threadSafe();
	/* 
	 * 设置SPR操作的参数 
	 */
	iso3d.sprimpl->spr_time_limit = 1.0; /* 时间限制，=0.0表示无限制 */

#if 0
	int iNod, iElem;
	int elems[1][4] = {{0, 1, 2, 3}};
	REAL qual, quals[1], qualgrads[6][3], volume, volumegrad[3], qualmeasure = 5;

	num_vertices = 4;
	vertices[0][0] = 0.021635;
	vertices[0][1] = -0.036032;
	vertices[0][2] = 0.003098;
	vertices[1][0] = 0.021751;
	vertices[1][1] = -0.036169;
	vertices[1][2] = 0.003115;
	vertices[2][0] = 0.021656;
	vertices[2][1] = -0.035976;
	vertices[2][2] = 0.003141;
	vertices[3][0] = 0.021518;
	vertices[3][1] = -0.035895;
	vertices[3][2] = 0.003080;
	quals[0] = tetquality(vertices[0], vertices[1], vertices[3], vertices[2], qualmeasure);

	::setMesh_LocSmoothing(elems, 1, quals);
	getoptinfo(0, 0, &qual, qualgrads, &volume, volumegrad, qualmeasure);
#endif

	SPRlogTime_Init();
#ifdef _TIMING_PERFORMANCE
	// iso3d.time_start[TOTAL_EXEC_TIME] =  SPRlogTime();
	iso3d.time_start[TOTAL_EXEC_TIME] = iso3d.time_start[READ_AND_INIT_TIME] = SPRlogTime();
#endif

	/*
	 * initialization
	 */

	if (iso3d.crtSources(pdCX, pdCY, pdCZ, pdIn, pdOu, pdSp, nPS, nLS, nTS) != 1)
	{
		spdlog::info("Error! fail to create sources\n");
		return 0;
	}

	if (iso3d.crtBkgMsh(pdBMNX, pdBMNY, pdBMNZ, pdBMNSpc, nBMN, 
		pnBMEFm, pnBMENg, nBME) != 1)
	{
		spdlog::info("Error! fail to create background mesh\n");
		return 0;
	}

	if (iso3d.crtSurface(pdSNX, pdSNY, pdSNZ, nSN, pnSFFm, pnSFPt, nSF) != 1)
	{
		spdlog::info("Error! fail to create surface\n");
		return 0;
	}

#ifdef _DEBUG
	const bool output = true;
#else
	const bool output = false;
#endif
	if (config.doMeshing)
	{
#ifdef _MULTITHREAD
		dtiso3d_meshing(iso3d, config.caseName, config.nMeshingThreads,
				/*embedded=*/true, /*output=*/output);
#else
		dtiso3d_meshing(iso3d, config.caseName,
				/*embedded=*/true);
#endif
	}

    if (config.doQualityImproving)
	{
		dtiso3d_improving(iso3d, config.caseName, config.nSmoothAttempts,
                /*embedded=*/true, /*output=*/output);
	}

#ifdef _TIMING_PERFORMANCE
	iso3d.time_start[WRIT_AND_FINA_TIME] = SPRlogTime();
#endif
	if (iso3d.outVolMsh(ppdMNX, ppdMNY, ppdMNZ, pnMN, ppnMEFm, ppnMENg, pnME, ppnPrt) != 1)
	{
		spdlog::info("Error. cannot output volume mesh\n");
		return 0;
	}
#ifdef _TIMING_PERFORMANCE
	iso3d.time_final[WRIT_AND_FINA_TIME] = SPRlogTime();
#endif

#ifdef _TIMING_PERFORMANCE
	iso3d.time_final[TOTAL_EXEC_TIME] = SPRlogTime();
#endif

#ifdef _ERROR_CHK
	nErrCd = iso3d.checkGlobalData(false);
 	if (nErrCd != 0)
	{
		printError(nErrCd);
		exit(1);
	}
#endif /* _ERROR_CHK */

	iso3d.printResult();
    printTimingProfile(config.doQualityImproving, iso3d);
	
	return 0;
}
