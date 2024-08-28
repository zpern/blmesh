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
#include <spdlog/spdlog.h> 
 #include "iso3d_error.h"
#include "iso3d_utility.h"
#include "iso3d.h"
#ifdef _TIMING_PERFORMANCE
#include "spr.h"
#include "time_perform.h"
#endif
#include "locsmoothing.h"
#include "iso3d_common.h"
#include "mshgen3d_def.h"
#ifdef _API_TEST_MODE
#include "tiger_volmtetrahedralizer.h"
#endif

#define _NGB 1

void printHead();
void printPrompt();

int main(int argc, char **argv)
{
	DTIso3D iso3d;
    DTIso3DConfig config;

    // g_bMeshingEnabled = true;  /* 网格生成功能是否被启动 -M */
    // g_bQuaImprEnabled = false;	/* 网格优化功能是否被启动 -Q */

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

	printHead();

    if (!config.parse(argc, argv))
    {
        printPrompt();
        // spdlog::info("Cannot initialize. Terminate\n");
        exit(1);
    }

    if (strlen(config.caseName) + 5 > MAX_FILE_LEN)
	{
        spdlog::info("too long file name %s", config.caseName);
		exit(1);
	}
	
	SPRlogTime_Init();
#ifdef _TIMING_PERFORMANCE
	iso3d.time_start[TOTAL_EXEC_TIME] =  SPRlogTime();
#endif

    if (config.doMeshing)
#ifdef _MULTITHREAD
        dtiso3d_meshing(iso3d, config.caseName, config.nMeshingThreads);
#else
        dtiso3d_meshing(iso3d, config.caseName);
#endif // _MULTITHREAD

	//目前不经过优化环节
	char savingName[MAX_FILE_LEN];
	config.getSavingName(savingName);
	dtiso3d_saving(iso3d, savingName);

    if (config.doQualityImproving)
	{
        if (!config.doMeshing)
		{
#ifdef _TIMING_PERFORMANCE
			iso3d.time_start[READ_AND_INIT_TIME] = SPRlogTime();
#endif
            if (dtiso3d_readVlmMesh(iso3d, config.caseName) != 1)
			{
				spdlog::info("Cannot read the volume mesh.\n");
				exit(1);
			}
#ifdef _TIMING_PERFORMANCE
			iso3d.time_final[READ_AND_INIT_TIME] = SPRlogTime();
#endif
			/* 建立哈希表，注意：目前只考虑二边流形问题；非二边流形问题是需要扩展 */
			iso3d.abstNodFirEle();
		}
        dtiso3d_improving(iso3d, config.caseName, config.nSmoothAttempts);
	}


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
    printTimingProfile(config.doQualityImproving,iso3d);
	
	return 0;
}
/*
 * 输出程序参数行输入格式要求
 * print the requirment for the command line
 */
void printPrompt()
{
	spdlog::info("Please check input arguments!\n");
	spdlog::info("Four options:\n");
	spdlog::info("\n\tOption 1: DTIso3D example_name");
	spdlog::info("\n\tOption 2: DTIso3D example_name -MQ\n");
	printf("\tDTIso3D reads a surface mesh (example_name.pls) &\n\ta mesh size map (optional, example_name.ba3)\n");
	printf("\tand then generates a volume mesh, and also improves \n\tthe volume mesh.\n");
	spdlog::info("\n\tOption 3: DTIso3D example_name -M\n");
	printf("\tDTIso3D reads a surface mesh (example_name.pls) &\n\ta mesh size map (optional, example_name.ba3)\n");
	printf("\tand then generates a volume mesh, but does not \n\timprove the volume mesh.\n");
	spdlog::info("\n\tOption 4: DTIso3D example_name -Q\n");
	spdlog::info("\tDTIso3D reads example_name.pl3 (volume mesh) and \n\tthen improves the volume mesh.\n\n");
}

/*
 * 输入程序题头信息
 * print the head info. of the program (PDMG)
 */
void printHead()
{
	spdlog::info(" ********************************************************************\n");
	spdlog::info("               3D Isotropic Delaunay Mesh Generation \n");
	spdlog::info("                           Version 1.0Belta\n\n");
	spdlog::info("         Center for Engineering & Scientific Computation\n");
	printf("                Zhejiang University, P. R. China\n");
	printf("	      CHEN Jianjun   Copyright reserved, 27/11/2014\n\n");
	printf("               For further information, please contact Dr. CHEN JJ\n");
	spdlog::info("                         Tel: +86-571-87951883\n");
	spdlog::info("                         Fax: +86-571-87953167\n");
	spdlog::info("                        Mail: chenjj@zju.edu.cn\n");
	spdlog::info(" ********************************************************************\n\n\n");
}
