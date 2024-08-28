#include <spdlog/spdlog.h> 
 #include "iso3d.h"
#include "iso3d_common.h"
#include "iso3d_error.h"
#include "iso3d_utility.h"
#include "mshgen3d_def.h"
#include "fileio.h"

#ifdef _TIMING_PERFORMANCE
#include "spr.h"
#include "time_perform.h"
#endif
#ifdef WIN32
//#include <windows.h> // TcharToChar
#endif

// bool g_bMeshingEnabled = true;  /* 网格生成功能是否被启动 -M */
// bool g_bQuaImprEnabled = true;	/* 网格优化功能是否被启动 -Q */
// int g_nThreads = 1; // Number of threads used to mesh.

static const char *EMBEDED_CASE_NAME = "embed_dtiso3d";

//#ifdef WIN32
//static void TcharToChar(const TCHAR * tchar, char * _char)
//{
//	int iLength;
//	//获取字节长度
//	iLength = WideCharToMultiByte(CP_ACP, 0, tchar, -1, nullptr, 0, nullptr, nullptr);
//	//将tchar值赋给_char
//	WideCharToMultiByte(CP_ACP, 0, tchar, -1, _char, iLength, nullptr, nullptr);
//}
//#endif

bool DTIso3DConfig::parse(int argc, char **argv)
{
	bool firstPath = true;
    for (int i = 1; i < argc; ++i)
    {
        // Start of an argument
        if (argv[i][0] == '-')
        {
            for (char *ptr = argv[i] + 1; *ptr; ++ptr)
            {
                switch (*ptr) {
                case 'M': doMeshing = true; break;
                case 'm': doMeshing = false; break;
                case 'Q': doQualityImproving = true; break;
                case 'q': doQualityImproving = false; break;
                case 'n':
                    if (ptr[1] != '\0' || i + 1 >= argc)
                        return false;
                    sscanf(argv[++i], "%d", &nMeshingThreads);
                    break;
                case 's':
                    if (ptr[1] != '\0' || i + 1 >= argc)
                        return false;
                    sscanf(argv[++i], "%d", &nSmoothAttempts);
                    break;
                default:
                    return false;
                }
            }
        }
		else {
			if (firstPath) {
				caseName = argv[i];
				firstPath = false;
			}
			else {
				saveName = argv[i];
				defaultSaving = false;
			}
		}
    }
    return true;
}

#ifdef _MULTITHREAD
int dtiso3d_meshing(DTIso3D& iso3d, const char chExampleName[],
                    int nThreads, bool embedded)
#else
int dtiso3d_meshing(DTIso3D& iso3d, const char chExampleName[], bool embedded)
#endif // _MULTITHREAD
{
#ifdef _ERROR_CHK
	int nErrCd;
#endif /* _ERROR_CHK */
	char strPLS[MAX_FILE_LEN], strBA3[MAX_FILE_LEN], strPL3[MAX_FILE_LEN];
	char strVTK[MAX_FILE_LEN];
#ifdef _NGB
	char strNGB[MAX_FILE_LEN];
#endif

	// FIXME: This license checking code is only valid for Windows.
#ifdef WIN
    TCHAR ExeFile[256];
    char licdir[256];

    memset(licdir, 0, sizeof(licdir));

    GetModuleFileName(nullptr, ExeFile, 256);
    TcharToChar(ExeFile, licdir);
//	spdlog::info("%s\n", licdir);
    (*strrchr(licdir, '\\')) = '\0';
//	spdlog::info("%s\n", licdir);
    strcat(licdir, "\\gbmesh.lic");
/*
	LiCheck liChk;
	//ret = liChk.CheckLic("gbmesh.lic");
	ret = liChk.CheckLic(licdir);
	if (ret == EXPIRED_LICENSE)
	{
		spdlog::info("License has expired.\n");
		exit(0);
	}
	else if (ret == N0_LICENSE_FILE)
	{
		spdlog::info("Cann't find license file.\n");
		exit(0);
	}
	else if (ret == INCORRECT_LICENSE_FILE)
	{
		spdlog::info("The license file is invalid.\n");
		exit(0);
	}
	else if (ret == INCORRECT_SYSTEM_TIME)
	{
		spdlog::info("The license file is invalid.\n");
		exit(0);
	}
	else if (ret == CORRECT_LICENSE_FILE)
	{
		//spdlog::info("The license file is valid.\n");
		;
	}*/
#endif
	if (embedded)
		chExampleName = EMBEDED_CASE_NAME;
	if (!embedded)
	{
		strcpy(strPLS, chExampleName);
		strcat(strPLS, ".pls");
		strcpy(strBA3, chExampleName);
		strcat(strBA3, ".ba3");
		strcpy(strPL3, chExampleName);
		strcat(strPL3, ".pl3");
		strcpy(strVTK, chExampleName);
		strcat(strVTK, ".vtk");
#ifdef _NGB
		strcpy(strNGB, chExampleName);
		strcat(strNGB, ".ngb");
#endif /* NGB */
		spdlog::info("Example Name: %s\n", chExampleName);

#ifdef _TIMING_PERFORMANCE
		iso3d.time_start[READ_AND_INIT_TIME] = SPRlogTime();
#endif
		//#TODO1 用宏区分BA3和VTK输入
		///*
		// * initialization
		// */
		//if (iso3d.readBA3(strBA3) <= 0)
		//{
		//	spdlog::info("fail to read BA3 file %s\n", strBA3);
		//}

#if 0
#ifdef DO_ZONE_DIVIDING
		std::shared_ptr<ZoneAttachInfo> zoneinfo = std::make_shared<ZoneAttachInfo>();
		zoneinfo->valid = false;
		iso3d.setZoneInfo(zoneinfo);
#endif
		if (iso3d.readPLS(strPLS) <= 0)
		{
			spdlog::info("fail to read PLS file %s\n", strPLS);
			exit(1);
		}
#else

		int bndPtNum = 0;
		double *bndPts = NULL;
		int bndFctNum = 0;
		int *bndFcts = NULL;
#ifdef DO_ZONE_DIVIDING
		std::shared_ptr<ZoneAttachInfo> zoneinfo = std::make_shared<ZoneAttachInfo>();
		int ret = readVtk(strVTK, &bndPtNum, &bndPts, &bndFctNum, &bndFcts, zoneinfo.get());

		iso3d.setZoneInfo(zoneinfo);
#else
		int ret = readVtk(strVTK, &bndPtNum, &bndPts, &bndFctNum, &bndFcts);
#endif // OUTPUT_ZONE_TAG
		if (ret < 0)
		{
			spdlog::info("fail to read VTK file %s\n", strVTK);
			exit(1);
		}
		iso3d.crtSurface(bndPtNum, bndPts, bndFctNum, bndFcts);
#endif

		//	if (iso3d.isSelfIntersect())
		//	{
		//		spdlog::info("Self-intersection found in the surface. \n");
		//		return 1;
		//	}

#ifdef _TIMING_PERFORMANCE
		iso3d.time_final[READ_AND_INIT_TIME] = SPRlogTime();
#endif
	}

#ifdef _TIMING_PERFORMANCE
	iso3d.time_start[BOUND_POINT_INST_TIME] = SPRlogTime();
#endif

#ifdef _DEBUG
    // {
    //     strcpy(strPLS, chExampleName);
    //     strcat(strPLS, ".pls");
    //     iso3d.writePLS(strPLS);
    // }
#endif
	/*
	 * insert boundary nodes
	 */
#ifdef _MULTITHREAD
	iso3d.bndPntInst(nThreads);
#else
	iso3d.bndPntInst();
#endif // _MULTITHREAD

#ifdef _DEBUG
    {
        strcpy(strPL3, chExampleName);
        strcat(strPL3, ".bndPntInst.pl3");
        iso3d.writePL3(strPL3);
    }
#endif

#ifdef _TIMING_PERFORMANCE
	iso3d.time_final[BOUND_POINT_INST_TIME] = SPRlogTime();
#ifdef _RECOVERY_FINAL
	iso3d.time_start[INNER_POINT_INST_TIME] = iso3d.time_final[BOUND_POINT_INST_TIME];
#endif
#endif

//	iso3d.readTetgenELE("mohne-a_cf-intmediate.ele", true);
//	iso3d.readTetgenELE("mohne-half-part-intmediate.ele", true);

#ifdef _ERROR_CHK
	// if (iso3d.isSelfIntersect())
	// {
	// 	spdlog::info("Self-intersection found in the surface. \n");
	// 	return 1;
	// }

    nErrCd = iso3d.checkGlobalData();
	if (nErrCd != 0)
	{
		printError(nErrCd);
		exit(1);
	}
#endif /* _ERROR_CHK*/

#ifdef _TIMING_PERFORMANCE
#ifdef _RECOVERY_FINAL
	iso3d.time_start[INNER_POINT_INST_TIME] = iso3d.time_final[BOUND_POINT_INST_TIME];
#endif
#endif

#ifndef _NO_FIELD_POINTS
#ifdef _RECOVERY_FINAL
	/*
	 * insert inner points
	 */
	iso3d.innerPntInst();
#endif /* _RECOVERY_FINAL */

#ifdef _TIMING_PERFORMANCE
#ifdef _RECOVERY_FINAL
	iso3d.time_final[INNER_POINT_INST_TIME] = SPRlogTime();
#endif /* _RECOVERY_FINAL */
#endif /* _TIMING_PERFORMANCE */
#endif /* #ifndef _NO_FIELD_POINTS */

#ifdef _RECOVER_WITHOUT_SPS

#ifdef _TIMING_PERFORMANCE
	iso3d.time_start[LSSP_RECV_EDGE_TIME] = SPRlogTime();
#endif /* _TIMING_PERFORMANCE */

	for (int i = 0; i < 0; i++)
	{
		iso3d.checkEdgeSP();
		iso3d.recvEdges_LessStnPnt();
	}

//	iso3d.writePL3("initmesh.pl3");
//	iso3d.writeTetgenELE("mohne-half-part-intmediate.1.ele");
//	exit(1);

#ifdef _ERROR_CHK
	nErrCd = iso3d.checkGlobalData();
	if (nErrCd != 0)
	{
		printError(nErrCd);
		exit(1);
	}
#endif /* _ERROR_CHK*/

#ifdef _TIMING_PERFORMANCE
	iso3d.time_final[LSSP_RECV_EDGE_TIME] = iso3d.time_start[NOSP_RECV_EDGE_TIME] = SPRlogTime();
#endif /* _TIMING_PERFORMANCE */

//	iso3d.checkEdgeSP();
//	exit(1);

#if 1

	iso3d.setBndProtectionFlag(true);
	iso3d.abstNodFirEle();	/* 获取首单元，建立哈希表 */
	//#zyj 先去掉
	for (int i = 0; i < 3; i++)
	{
		iso3d.recvBndEdge_NoStnPnt_V2();
		iso3d.recvBndFace_NoStnPnt_V2();
	}

	for (int i = 0; i < 0; i++)
	{
		iso3d.classifyBndEdges(nullptr);
		iso3d.recvBndEdge_NoStnPnt();
		iso3d.checkEdgeSP();
	}
//	iso3d.writeTetgenELE("cami1a-intmediate.1.ele");
//	exit(1);

#if 0
	iso3d.setBndProtectionFlag(true);
	iso3d.classifyBndEdges(nullptr);
	for (int i = 0; i < 0; i++)
	{
		iso3d.recvBndEdge_NoStnPnt();
		iso3d.checkEdgeSP();
	}
#endif

//	iso3d.writeTetgenELE("mohne_a-intmediate.1.ele");
#else
	iso3d.recvEdges_NoStnPnt();
#endif
//	exit(1);

#ifdef _ERROR_CHK
	nErrCd = iso3d.checkGlobalData();
	if (nErrCd != 0)
	{
		printError(nErrCd);
		exit(1);
	}
#endif /* _ERROR_CHK*/

#ifdef _TIMING_PERFORMANCE
	iso3d.time_final[NOSP_RECV_EDGE_TIME] = iso3d.time_start[NOSP_RECV_FACE_TIME] = SPRlogTime();
#endif /* _TIMING_PERFORMANCE */

#if 0
	iso3d.checkEdgeSP();

	iso3d.recvFaces_NoStnPnt();


#ifdef _ERROR_CHK
	nErrCd = iso3d.checkGlobalData();
	if (nErrCd != 0)
	{
		printError(nErrCd);
		exit(1);
	}
#endif /* _ERROR_CHK*/
#endif

#ifdef _TIMING_PERFORMANCE
	iso3d.time_final[NOSP_RECV_FACE_TIME] = SPRlogTime();
#endif /* _TIMING_PERFORMANCE */
#endif /* _RECOVER_WITHOUT_SPS */


#ifdef _TIMING_PERFORMANCE
	iso3d.time_start[CONF_RECV_EDGE_TIME] = SPRlogTime();
#endif /* _TIMING_PERFORMANCE */

//	iso3d.checkEdgeSP();
//	iso3d.writeTetgenELE("mohne-half-part-intmediate.ele");
	//exit(1);

	/*
	 * recover boundary edges
	 */
	iso3d.recvEdges();

//	if (iso3d.checkConfEdges() != 0)
//		exit(1);

#ifdef _ERROR_CHK
	nErrCd = iso3d.checkGlobalData();
	if (nErrCd != 0)
	{
		printError(nErrCd);
		exit(1);
	}
#endif /* _ERROR_CHK*/

#ifdef _TIMING_PERFORMANCE
	iso3d.time_final[CONF_RECV_EDGE_TIME] = iso3d.time_start[CONF_RECV_FACE_TIME] = SPRlogTime();
#endif

	/*
	 * recovery boudary faces
	 */
	iso3d.recvFaces();

 // if (iso3d.checkOuIn(true) != 0)
 // {
 //		exit(1);
 //	}

#ifdef _ERROR_CHK
	nErrCd = iso3d.checkGlobalData();
	if (nErrCd != 0)
	{
		printError(nErrCd);
		exit(1);
	}
#endif /* _ERROR_CHK*/

#ifdef _TIMING_PERFORMANCE
	iso3d.time_final[CONF_RECV_FACE_TIME] = iso3d.time_start[TYPE_DEL_ELEM_TIME] = SPRlogTime();
#endif

	exit(1);
	/*
	 * remove outer elements
	 */
	iso3d.typeEles();

#ifdef _ERROR_CHK
	nErrCd = iso3d.checkGlobalData(false);
	if (nErrCd != 0)
	{
		printError(nErrCd);
		exit(1);
	}
#endif /* _ERROR_CHK */

#ifdef _TIMING_PERFORMANCE
	iso3d.time_final[TYPE_DEL_ELEM_TIME] = SPRlogTime();
#endif

//	iso3d.printConformalRecvInfo();

#ifdef _CONSTRAINED_RECOVERY

#ifdef _TIMING_PERFORMANCE
	iso3d.time_start[CONS_RECV_EDGE_TIME] = SPRlogTime();
#endif
	iso3d.prepareRecvDesBnds();
	//iso3d.recvDesEdges();

#ifdef _ERROR_CHK
	nErrCd = iso3d.checkGlobalData(false);
	if (nErrCd != 0)
	{
		printError(nErrCd);
		exit(1);
	}
#endif /* _ERROR_CHK */

#ifdef _TIMING_PERFORMANCE
	iso3d.time_final[CONS_RECV_EDGE_TIME] = iso3d.time_start[CONS_RECV_FACE_TIME] = SPRlogTime();
#endif

	//iso3d.recvDesFacets();

#ifdef _ERROR_CHK
	nErrCd = iso3d.checkGlobalData(false);
	if (nErrCd != 0)
	{
		printError(nErrCd);
		exit(1);
	}
#endif /* _ERROR_CHK */

#ifdef _TIMING_PERFORMANCE
	iso3d.time_final[CONS_RECV_FACE_TIME] = iso3d.time_start[CONS_RECV_PNT_RMV_TIME] = SPRlogTime();
#endif

	for (int i = 0; i < 1; i++)
	{
		iso3d.removeInnerSteinerPoint_SPR();	/* 利用空腔重连算法消除Steiner点 */
		iso3d.removeInnerSteinerPoint_Flip();
	}

 //	iso3d.removeEdgeSteinerPoint_SPR();
 //	iso3d.removeFacSteinerPoint_SPR();
 //	iso3d.removeAllLeftNodes();
	iso3d.printSteinerPntInfo();

#ifdef _ERROR_CHK
	nErrCd = iso3d.checkGlobalData(false);
	if (nErrCd != 0)
	{
		printError(nErrCd);
		exit(1);
	}
#endif /* _ERROR_CHK */

#ifdef _TIMING_PERFORMANCE
	iso3d.time_final[CONS_RECV_PNT_RMV_TIME] = SPRlogTime();
#endif

#endif /* _CONSTRAINED_RECOVERY */

	iso3d.typeDelEles();

#ifdef _ERROR_CHK
	nErrCd = iso3d.checkGlobalData(false);
	if (nErrCd != 0)
	{
		printError(nErrCd);
		exit(1);
	}
#endif /* _ERROR_CHK */

#ifndef _NO_FIELD_POINTS
#ifndef _RECOVERY_FINAL
#ifdef _TIMING_PERFORMANCE
	iso3d.time_start[INNER_POINT_INST_TIME] = SPRlogTime();
#endif

	iso3d.innerPntInst();

#ifdef _ERROR_CHK
	nErrCd = iso3d.checkGlobalData(false);
	if (nErrCd != 0)
	{
		printError(nErrCd);
		exit(1);
	}
#endif /* _ERROR_CHK */

#ifdef _CONSTRAINED_RECOVERY
#ifdef _TIMING_PERFORMANCE
	iso3d.time_final[INNER_POINT_INST_TIME] = SPRlogTime();
#endif

	/* 在插入内部点后，可尝试再删除所有Steiner点 */
	iso3d.prepareRmvSteinerPntAftInstInnNode();
	iso3d.removeInnerSteinerPoint_SPR();	/* 利用空腔重连算法消除Steiner点 */
	iso3d.removeEdgeSteinerPoint_SPR(false);
	iso3d.removeFacSteinerPoint_SPR(false);
 //	iso3d.removeAllLeftNodes();
	iso3d.printSteinerPntInfo();


#ifdef _TIMING_PERFORMANCE
	iso3d.time_final[CONS_RECV_PNT_RMV_TIME] += SPRlogTime() - iso3d.time_final[INNER_POINT_INST_TIME];
#endif
#endif /* _CONSTRAINED_RECOVERY */

#endif /* #ifndef _RECOVERY_FINAL */
#endif/* _NO_FIELD_POINTS */


#ifdef _TIMING_PERFORMANCE
	iso3d.time_start[WRIT_AND_FINA_TIME] = SPRlogTime();
#endif


	iso3d.rmvNodsAndEles();

#ifdef _CONSTRAINED_RECOVERY
	//暂时方案，边界可能不一定是约束恢复的，modified by zyj
	iso3d.fillParents(true);
	//iso3d.fillParents(false);
#else
	iso3d.fillParents(true);
#endif

#ifdef DO_ZONE_DIVIDING
	iso3d.solveZone();
#endif // DO_ZONE_DIVIDING

#ifdef _TIMING_PERFORMANCE
	iso3d.time_final[WRIT_AND_FINA_TIME] = SPRlogTime();
#endif

	return 1;
}

int dtiso3d_saving(DTIso3D& iso3d, const char savePath[])
{
	char strVTK[MAX_FILE_LEN], extraVTK[MAX_FILE_LEN];
	strcpy(strVTK, savePath);
	strcat(strVTK, ".vtk");
	strcpy(extraVTK, savePath);
	strcat(extraVTK, "_visual.vtk");
#ifdef OUTPUT_ZONE_TAG
#ifdef OUTPUT_CENTRIC_POINT_INFO
	iso3d.writeVtkFile(strVTK, true, true);
#else
	iso3d.writeVtkFile(strVTK, true, false);
#endif // OUTPUT_CENTRIC_POINT_INFO
#else
#ifdef OUTPUT_CENTRIC_POINT_INFO
	iso3d.writeVtkFile(strVTK, false, true);
#else
	iso3d.writeVtkFile(strVTK, false, false);
#endif // OUTPUT_CENTRIC_POINT_INFO
#endif // OUTPUT_ZONE_TAG
	//输出可以可视化的vtk模型
	//iso3d.writeVisualVtkFile(extraVTK);
	return 0;
}

//#TODO2 移动improving部分的output模块到公共saving部分
int dtiso3d_improving(DTIso3D& iso3d, const char chExampleName[],
					  int nSmoothAttempts, bool embedded, bool output)
{
	if (embedded)
		chExampleName = EMBEDED_CASE_NAME;

#ifdef _ERROR_CHK
	int nErrCd;
#endif /* _ERROR_CHK */
	char strPL3[MAX_FILE_LEN], strQual[MAX_FILE_LEN];
#ifdef _NGB
	char strNGB[MAX_FILE_LEN];
#endif

#ifdef _TIMING_PERFORMANCE
	iso3d.time_start[QUAL_IMPRV_TIME] = SPRlogTime();
#endif
	/*added*/
	//iso3d.printResult();
	//printTimingProfile();
#if 0
	//优化
	iso3d.updatebdnpnts();
	iso3d.initoptimize();
//	iso3d.setsmoothpar(argv[1]);
#endif

	if (output)
	{
		strcpy(strQual, chExampleName);
		strcat(strQual, ".init.allangles.qual");
		iso3d.printElemQuality(strQual, QUALALLANGLES, 0.0, false);
	}

	iso3d.meshopt(nSmoothAttempts);

#ifdef _TIMING_PERFORMANCE
	iso3d.time_final[QUAL_IMPRV_TIME] = SPRlogTime();
#endif

	if (output)
	{
		strcpy(strQual, chExampleName);
		strcat(strQual, ".allangles.qual");
		iso3d.printElemQuality(strQual, QUALALLANGLES, 0.0, false);

		strcpy(strQual, chExampleName);
		strcat(strQual, ".sineangle.qual");
		iso3d.printElemQuality(strQual, QUALSINEANGLE, 30, false);

		strcpy(strQual, chExampleName);
		strcat(strQual, ".voledge.qual");
		iso3d.printElemQuality(strQual, QUALVLRMS3RATIO, 0.3, false);
	}
//
// 	iso3d.shelltransform(sin(ANGLE(30.0)));
// 	//iso3d.tranverseST(1);
// 	iso3d.smooth(8);
//
// 	iso3d.checkbadcells();
// 	iso3d.smooth(8);
//
// 	iso3d.shelltransform(sin(ANGLE(30.0)));
// 	//iso3d.tranverseST(1);
// 	iso3d.smooth(8);

// 	iso3d.checkbadcells();
// 	iso3d.smooth(8);

	if (output)
	{
		strcpy(strQual, chExampleName);
		strcat(strQual, ".bad_angles.pl3");
		iso3d.printBadCells(strQual, QUALMINSINE, sin(5.0*3.1415926/180.0));
	}

	iso3d.rmvNodsAndEles();
	iso3d.fillParents();

//	iso3d.evaluateQuality();
//	iso3d.printbadcells();

	if (output)
	{
		strcpy(strPL3, chExampleName);
		//	strcat(strPL3, ".pl3");
		strcat(strPL3, ".improved.pl3");
		iso3d.writePL3(strPL3);
	}

	return 1;
}

int dtiso3d_readVlmMesh(DTIso3D& iso3d, const char chExampleName[])
{
	char strPL3[MAX_FILE_LEN];

	strcpy(strPL3, chExampleName);
	strcat(strPL3, ".pl3");

	return iso3d.readPL3(strPL3);
}

/*
 * 输出必要的结果信息
 * print the resulting info. of the program (PDMG)
 */
void printTimingProfile(bool doQualityImproving, DTIso3D& iso3d)
{
#ifdef _TIMING_PERFORMANCE
	double totalTime = 0.0;
	double percentFactor = 0.0, computingRatio;
	int totalECCallings, succECCallings = 0, totalECCallingsSub = 0, succECCallingsSub = 0, i;

	for (i = TIME_START_IDX; i <= TYPE_DEL_ELEM_TIME; i++)
		iso3d.time_elaps[i] = iso3d.time_final[i] - iso3d.time_start[i];
	totalTime = iso3d.time_elaps[TOTAL_EXEC_TIME];

//	if (totalTime <= 0.0)
//		return;

	percentFactor = 100.0/totalTime;

	spdlog::info(" ********************************************************************\n");
	spdlog::info("                            Timing Profile\n");
	spdlog::info(" ********************************************************************\n");
	spdlog::info("%-24s%-12s%-12s\n", "Step Name", "Time(s)", "Percentage(%)");
	spdlog::info("%-24s%-12.2f%-12.2f\n", "Read & Init Input", iso3d.time_elaps[READ_AND_INIT_TIME], percentFactor*iso3d.time_elaps[READ_AND_INIT_TIME]);
	spdlog::info("%-24s%-12.2f%-12.2f\n", "Insert Boun. Pnts", iso3d.time_elaps[BOUND_POINT_INST_TIME], percentFactor*iso3d.time_elaps[BOUND_POINT_INST_TIME]);
    // spdlog::info("%-24s%-12.2f%-12.2f\n", "Insert Bound Pnts (Actual)", iso3d.time_elaps[BOUND_POINT_INST_TIME_ACTURL], percentFactor*iso3d.time_elaps[BOUND_POINT_INST_TIME_ACTURL]);
#ifndef _NO_FIELD_POINTS
#ifdef _RECOVERY_FINAL
	spdlog::info("%-24s%-12.2f%-12.2f\n", "Insert Innr. Pnts", iso3d.time_elaps[INNER_POINT_INST_TIME], percentFactor*iso3d.time_elaps[INNER_POINT_INST_TIME]);
#endif
#endif /* _NO_FIELD_POINTS */
#ifdef _RECOVER_WITHOUT_SPS
	spdlog::info("%-24s%-12.2f%-12.2f\n", "LsSP. Recv. Edges", iso3d.time_elaps[LSSP_RECV_EDGE_TIME], percentFactor*iso3d.time_elaps[LSSP_RECV_EDGE_TIME]);
	spdlog::info("%-24s%-12.2f%-12.2f\n", "NoSP. Recv. Edges", iso3d.time_elaps[NOSP_RECV_EDGE_TIME], percentFactor*iso3d.time_elaps[NOSP_RECV_EDGE_TIME]);
	spdlog::info("%-24s%-12.2f%-12.2f\n", "NoSP. Recv. Faces", iso3d.time_elaps[NOSP_RECV_FACE_TIME], percentFactor*iso3d.time_elaps[NOSP_RECV_FACE_TIME]);
#endif /* __RECOVER_WITHOUT_SPS */
	spdlog::info("%-24s%-12.2f%-12.2f\n", "Conf. Recv. Edges", iso3d.time_elaps[CONF_RECV_EDGE_TIME], percentFactor*iso3d.time_elaps[CONF_RECV_EDGE_TIME]);
	spdlog::info("%-24s%-12.2f%-12.2f\n", "Conf. Recv. Faces", iso3d.time_elaps[CONF_RECV_FACE_TIME], percentFactor*iso3d.time_elaps[CONF_RECV_FACE_TIME]);
	spdlog::info("%-24s%-12.2f%-12.2f\n", "Dele. Outer Elems.",iso3d.time_elaps[TYPE_DEL_ELEM_TIME], percentFactor*iso3d.time_elaps[TYPE_DEL_ELEM_TIME]);
#ifdef _CONSTRAINED_RECOVERY
	spdlog::info("%-24s%-12.2f%-12.2f\n", "Cons. Recv. Edges",iso3d.time_elaps[CONS_RECV_EDGE_TIME], percentFactor*   iso3d.time_elaps[CONS_RECV_EDGE_TIME]);
	spdlog::info("%-24s%-12.2f%-12.2f\n", "Cons. Recv. Faces",iso3d.time_elaps[CONS_RECV_FACE_TIME], percentFactor*   iso3d.time_elaps[CONS_RECV_FACE_TIME]);
	spdlog::info("%-24s%-12.2f%-12.2f\n", "Rmv. Steiner Pnts",iso3d.time_elaps[CONS_RECV_PNT_RMV_TIME], percentFactor*iso3d.time_elaps[CONS_RECV_PNT_RMV_TIME]);
#endif
#ifndef _NO_FIELD_POINTS
#ifndef _RECOVERY_FINAL
	spdlog::info("%-24s%-12.2f%-12.2f\n", "Insert Innr. Pnts", iso3d.time_elaps[INNER_POINT_INST_TIME], percentFactor*iso3d.time_elaps[INNER_POINT_INST_TIME]);
#endif
#endif /* _NO_FIELD_POINTS */
	spdlog::info("%-24s%-12.2f%-12.2f\n", "Finish and Output", iso3d.time_elaps[WRIT_AND_FINA_TIME], percentFactor*iso3d.time_elaps[WRIT_AND_FINA_TIME]);
	spdlog::info("%-24s%-12.2f%-12.2f\n", "Recover Edge Loop", iso3d.time_elaps[RECV_EDG_LOOP_TIME], percentFactor*iso3d.time_elaps[RECV_EDG_LOOP_TIME]);
	spdlog::info("%-24s%-12.2f%-12.2f\n", "Recover Face Loop", iso3d.time_elaps[RECV_FAC_LOOP_TIME], percentFactor*iso3d.time_elaps[RECV_FAC_LOOP_TIME]);
	spdlog::info("%-24s%-12.2f%-12.2f\n", "Remesh EDG. Shell", iso3d.time_elaps[REMESH_ASHELL_TIME], percentFactor*iso3d.time_elaps[REMESH_ASHELL_TIME]);
	spdlog::info("%-24s%-12.2f%-12.2f\n", "Remove Edge Recv.", iso3d.time_elaps[REMOV_AN_EDGE_TIME], percentFactor*iso3d.time_elaps[REMOV_AN_EDGE_TIME]);
	spdlog::info("%-24s%-12.2f%-12.2f\n", "Reflex Edge Check", iso3d.time_elaps[REFLEX_EDG_CH_TIME], percentFactor*iso3d.time_elaps[REFLEX_EDG_CH_TIME]);
	spdlog::info("%-24s%-12.2f%-12.2f\n", "Search ALL Shells", iso3d.time_elaps[FIND_ALLSHELL_TIME], percentFactor*iso3d.time_elaps[FIND_ALLSHELL_TIME]);
	spdlog::info("%-24s%-12.2f%-12.2f\n", "Dyn. Rm. Sh. Main", iso3d.time_elaps[DYN_RSHL_MAIN_TIME], percentFactor*iso3d.time_elaps[DYN_RSHL_MAIN_TIME]);
	spdlog::info("%-24s%-12.2f%-12.2f\n", "Chk. ED. Recovery", iso3d.time_elaps[CHK_EDGE_RECV_TIME], percentFactor*iso3d.time_elaps[CHK_EDGE_RECV_TIME]);
	spdlog::info("%-24s%-12.2f%-12.2f\n", "Chk. FC. Recovery", iso3d.time_elaps[CHK_FACE_RECV_TIME], percentFactor*iso3d.time_elaps[CHK_FACE_RECV_TIME]);
	spdlog::info("%-24s%-12.2f%-12.2f\n", "Qual. Improvement", iso3d.time_elaps[QUAL_IMPRV_TIME], percentFactor*iso3d.time_elaps[QUAL_IMPRV_TIME]);
	spdlog::info("%-24s%-12.2f%-12.2f\n", "Qual. Recurs. ST.", iso3d.time_elaps[QUAL_RECV_SH_TRANS_TIME], percentFactor*iso3d.time_elaps[QUAL_RECV_SH_TRANS_TIME]);
	spdlog::info("%-24s%-12.2f%-12.2f\n", "Qual. SPR Improv.", iso3d.time_elaps[QUAL_SPR_CALLING_TIME], percentFactor*iso3d.time_elaps[QUAL_SPR_CALLING_TIME]);
	spdlog::info("%-24s%-12.2f%-12.2f\n", "Opt. Mesh Smthing", iso3d.time_elaps[MESH_SMOOTHING_TIME], percentFactor*iso3d.time_elaps[MESH_SMOOTHING_TIME]);
	spdlog::info("%-24s%-12.2f%-12.2f\n", "Edge Contraction ", iso3d.time_elaps[EDGE_CONTRACT_TIME], percentFactor*iso3d.time_elaps[EDGE_CONTRACT_TIME]);
	spdlog::info("%-24s%-12.2f%-12.2f\n", "Edge Splitting   ", iso3d.time_elaps[EDGE_SPLITTING_TIME], percentFactor*iso3d.time_elaps[EDGE_SPLITTING_TIME]);
	spdlog::info("%-24s%-12.2f%-12.2f\n", "Qual. Evaluation ", iso3d.time_elaps[QUALITY_EVAL_TIME], percentFactor*iso3d.time_elaps[QUALITY_EVAL_TIME]);
	spdlog::info("%-24s%-12.2f%-12.2f\n", "Total Execu. Time", iso3d.time_elaps[TOTAL_EXEC_TIME], 100.0);
#endif /* _TIMING_PERFORMANCE */
#ifdef _VERBOSE
#ifdef _ROBUST_BW_KERNEL
	spdlog::info("%-24s%-12.2f%-12.2f\n", "Form Init. Cavity", iso3d.time_elaps[FORM_INIT_CAVITY_TIME], percentFactor*iso3d.time_elaps[FORM_INIT_CAVITY_TIME]);
	spdlog::info("%-24s%-12.2f%-12.2f\n", "Modify the Cavity", iso3d.time_elaps[MODIFY_CAVITY_TIME], percentFactor*iso3d.time_elaps[MODIFY_CAVITY_TIME]);
	spdlog::info("%-24s%-12.2f%-12.2f\n", "Update the Trian.", iso3d.time_elaps[UPDATE_TRIANGULATION_TIME], percentFactor*iso3d.time_elaps[UPDATE_TRIANGULATION_TIME]);
	spdlog::info("%-24s%-12.2f%-12.2f\n", "Create Inner Pnt.", iso3d.time_elaps[CREATE_INNER_POINT_TIME], percentFactor*iso3d.time_elaps[CREATE_INNER_POINT_TIME]);
#ifdef _INSIDE_TIMING_CALC
	spdlog::info("%-24s%-12.2f%-12.2f\n", "Calc. Elem. Para.", iso3d.time_elaps[CREATE_ELEM_PARA_TIME], percentFactor*iso3d.time_elaps[CREATE_ELEM_PARA_TIME]);
#endif /* _INSIDE_TIMING_CALC */
#endif /* _ROBUST_BW_KERNEL */
#endif /* _VERBOSE */
	spdlog::info("%-24s%-12.2f%-12.2f\n", "Call Poly. Recon.", iso3d.time_elaps[SPR_CALLING_TIME], percentFactor*iso3d.time_elaps[SPR_CALLING_TIME]);
#ifdef _VERBOSE
	printf("%-24s%-12.2f%-12.2f\t%fsec. per node\n", "Cln. af. Suc. Ad.", iso3d.time_elaps[ADD_SUCC_CLEAN_TIME],
		percentFactor*iso3d.time_elaps[ADD_SUCC_CLEAN_TIME],iso3d.time_elaps[ADD_SUCC_CLEAN_TIME]/(double)g_nInnInstPoint);
	printf("%-24s%-12.2f%-12.2f\t%fsec. per node\n", "Cln. af. FAI. Ad.", iso3d.time_elaps[ADD_FAIL_CLEAN_TIME],
		percentFactor*iso3d.time_elaps[ADD_FAIL_CLEAN_TIME],iso3d.time_elaps[ADD_FAIL_CLEAN_TIME]/(double)(g_nInnTryPoints - g_nInnInstPoint));
#endif /* _VERBOSE */

	printf("\nNo. of interior points that are attempted to be inserted and of which really inserted:: (%d,%d).\n", iso3d.g_nInnTryPoints, iso3d.g_nInnInstPoint);
	spdlog::info("Average cavity size for the really inserted points: {}.\n", iso3d.g_nInnInstPoint > 0 ? iso3d.g_nInnTotalCaviElem /(double)iso3d.g_nInnInstPoint : 0.0);
	printf("#Callings of the SPR routine and of which the successful ones: (%d,%d).\n", iso3d.g_nSPRCallings, iso3d.g_nSPRCallSucc);
	printf("Average Size of the reconnected polyhedra  & Average time cost per SPR calling. #face: %f, Time cost(s):%f.\n",
		(double)iso3d.g_nSPRPolySize/ iso3d.g_nSPRCallings, iso3d.time_elaps[SPR_CALLING_TIME]/ iso3d.g_nSPRCallings);
    if (doQualityImproving)
	{
		spdlog::info("#Callings of the smoothing routine in the smoothing passes: {}\n", iso3d.g_nGlobalSmoothingCallings);
		spdlog::info("#Successful Callings: {}\n", iso3d.g_nGlobalSmoothingCallSucc);
		spdlog::info("#Callings of the smoothing routine in the edge contraction passes: {}\n", iso3d.g_nEdContSmoothingCallings);
		spdlog::info("#Successful Callings: {}\n", iso3d.g_nEdContSmoothingCallSucc);
		spdlog::info("#Callings of the smoothing routine in the edge splitting passes: {}\n", iso3d.g_nEdSpltSmoothingCallings);
		spdlog::info("#Successful Callings: {}\n", iso3d.g_nEdSpltSmoothingCallSucc);

#if 0
		totalECCallings = succECCallings = 0;
		for (i = 0; i < MAX_EC_QUALITY_RATIO_DIVIDE; i++)
		{
			succECCallings += g_nECQualityRatio_SuccStatistics[i];
			totalECCallings += g_nECQualityRatio_TotalStatistics[i];
		}
		spdlog::info("Distribution of quality ratios for successful edge contraction callings.\n");
		printf("#Succ. Callings/#Total Callings:%d/%d = %f%%.\n", succECCallings, totalECCallings,
			totalECCallings > 0 ? 100.0*succECCallings/totalECCallings : 0.0);
		spdlog::info("Min = {}; Max = {}.\n", g_dECQualityRatio_MinRatio, g_dECQualityRatio_MaxRatio);
		spdlog::info("Distribution of the ratios.\n");
		totalECCallingsSub = 0;
		succECCallingsSub = 0;
		for (i = 0; i < MAX_EC_QUALITY_RATIO_DIVIDE; i++)
		{
			totalECCallingsSub += g_nECQualityRatio_TotalStatistics[i];
			succECCallingsSub += g_nECQualityRatio_SuccStatistics[i];
			printf(" [%f--%f]:%d/%d\t= %4.1f%%; Score VS Load: %4.1f%% : %4.1f%%.\n",
				g_dECQualityRatio_Divide[i], g_dECQualityRatio_Divide[i+1],
				g_nECQualityRatio_SuccStatistics[i], g_nECQualityRatio_TotalStatistics[i],
				g_nECQualityRatio_TotalStatistics[i] > 0 ?
				100.0*g_nECQualityRatio_SuccStatistics[i]/g_nECQualityRatio_TotalStatistics[i] :0.0,
				succECCallings > 0 ? 100.0*succECCallingsSub/succECCallings : 0.0,
				totalECCallings > 0 ? 100.0*totalECCallingsSub/totalECCallings : 0.0);
		}
#endif
	}
}

