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
#include "fileio.h"
#ifdef _TIMING_PERFORMANCE
#include "spr.h"
#include "time_perform.h"
#endif
#include "tiger_volmtetrahedralizer2.h"

//#define TESTPLS
//#define TESTVTK
//#define TESTMULTIDOMAIN
//#define TEST_MULTITHREAD
#define TEST_TIGER2

#if 0
int main(int argc, char **argv)
{
	const int MAX_FILE_LEN = 1024;
	char strPLS[MAX_FILE_LEN], strPL3[MAX_FILE_LEN];
	int bndPtNum = 0;
	double *bndPts = NULL;
	int bndFctNum = 0;
	int *bndFcts = NULL;
	int intCstNum[] = {0, 0, 0};
	int *intCsts = NULL;
	double (*sizingFunc)(double x, double y, double z) = NULL;
	int mshPtNum = 0;
	double *mshPts = NULL; 
	int volElmNum = 0;
	int *volElms = NULL;
	int errCode = 0;

	//printHead();

	if (argc < 2)
	{
		spdlog::info("Cannot initialize. Terminate\n");
		exit(1);
	}

	if (strlen(argv[1]) + 5 > MAX_FILE_LEN)
	{
		spdlog::info("too long file name %s", argv[1]);
		exit(1);
	}
	
//	SPRlogTime_Init();
//#ifdef _TIMING_PERFORMANCE
//	time_start[TOTAL_EXEC_TIME] =  SPRlogTime();
//#endif

	strcpy(strPLS, argv[1]);
	strcat(strPLS, ".pls");
	strcpy(strPL3, argv[1]);
	strcat(strPL3, ".pl3");

	readPLS(strPLS, &bndPtNum, &bndPts, &bndFctNum, &bndFcts);

	errCode = API_Tetrahedralize_Single_Domain(
		bndPtNum, bndPts, bndFctNum, bndFcts, 
		intCstNum, NULL,
		NULL,
		&mshPtNum, &mshPts,
		&volElmNum, &volElms);
	if (errCode != 0 && errCode != 1)
	{
		spdlog::info("Cannot tetrahedralize the domain.\n");
		goto FAIL;
	}

	writePL3(strPL3, mshPtNum, mshPts, volElmNum, volElms);

//#ifdef _TIMING_PERFORMANCE
//	time_final[TOTAL_EXEC_TIME] = SPRlogTime();
//#endif
//	spdlog::info("Elapsed time: {} seconds.\n", time_final[TOTAL_EXEC_TIME] - time_start[TOTAL_EXEC_TIME]); 
	
	goto SUCC;

FAIL:
SUCC:

	return errCode;
}
#endif

#ifdef TESTPLS
int main(int argc, char **argv)
{
	const int MAX_FILE_LEN = 1024;
	char strPLS[MAX_FILE_LEN], strPL3[MAX_FILE_LEN];
	char strPLS_2[MAX_FILE_LEN], strPL3_2[MAX_FILE_LEN];
	int bndPtNum = 0;
	double *bndPts = NULL;
	int bndFctNum = 0;
	int *bndFcts = NULL;
	int intCstNum[] = { 0, 0, 0 };
	int *intCsts = NULL;
	double(*sizingFunc)(double x, double y, double z) = NULL;
	int mshPtNum = 0;
	double *mshPts = NULL;
	int volElmNum = 0;
	int *volElms = NULL;
	int errCode = 0;

	//printHead();

	if (argc < 2)
	{
		spdlog::info("Cannot initialize. Terminate\n");
		exit(1);
	}

	if (strlen(argv[1]) + 5 > MAX_FILE_LEN)
	{
		spdlog::info("too long file name %s", argv[1]);
		exit(1);
	}

	//	SPRlogTime_Init();
	//#ifdef _TIMING_PERFORMANCE
	//	time_start[TOTAL_EXEC_TIME] =  SPRlogTime();
	//#endif

	strcpy(strPLS, argv[1]);
	strcat(strPLS, ".pls");
	strcpy(strPL3, argv[1]);
	strcat(strPL3, ".pl3");

	readPLS(strPLS, &bndPtNum, &bndPts, &bndFctNum, &bndFcts);

	//errCode = API_Tetrahedralize_Single_Domain(
	//	bndPtNum, bndPts, bndFctNum, bndFcts, 
	//	intCstNum, NULL,
	//	NULL,
	//	&mshPtNum, &mshPts,
	//	&volElmNum, &volElms);

	int obj, obj1, obj2;
	std::vector<int> vecObj;

	double vol;
	//errCode = API_Compute_Poly_Volume(bndPtNum, bndPts, bndFctNum, bndFcts, &vol);

	errCode = API_Tetrahedralize_Single_Domain(
		bndPtNum, bndPts, bndFctNum, bndFcts,
		intCstNum, NULL,
		NULL,
		&obj1);

	if (errCode != 0 && errCode != 1)
	{
		spdlog::info("Cannot tetrahedralize the domain.\n");
		return errCode;
	}
	vecObj.push_back(obj1);

#if 0
	strcpy(strPLS_2, "remesh_bnd.pls");
	strcpy(strPL3_2, "remesh_bnd.pl3");

	readPLS(strPLS_2, &bndPtNum, &bndPts, &bndFctNum, &bndFcts);

	errCode = API_Tetrahedralize_Single_Domain(
		bndPtNum, bndPts, bndFctNum, bndFcts,
		intCstNum, NULL,
		NULL,
		&obj2);
	if (errCode != 0 && errCode != 1)
	{
		spdlog::info("Cannot tetrahedralize the domain.\n");
		goto FAIL;
	}
	vecObj.push_back(obj2);
#endif

	double x, y, z;
	int a, b, c, d;
	for (int i = 0; i < vecObj.size(); i++)
	{
		obj = vecObj[i];
		mshPtNum = API_GetTetrahedraPntNum(obj);
		volElmNum = API_GetTetrahedraElemNum(obj);
		mshPts = (double *)malloc(mshPtNum * sizeof(double) * 3);
		volElms = (int *)malloc(volElmNum * sizeof(int) * 4);

		if (NULL == mshPts || NULL == volElms)
		{
			spdlog::info("Not enough memory.\n");
			errCode = -1;
			return errCode;
		}
		for (int j = 0; j < mshPtNum; j++)
		{
			API_GetTetrahedraPointCoord(obj, j + 1, mshPts[3 * j], mshPts[3 * j + 1], mshPts[3 * j + 2]);
		}

		for (int j = 0; j < volElmNum; j++)
		{
			API_GetTetrahedraElemPntIdx(obj, j + 1, volElms[4 * j], volElms[4 * j + 1], volElms[4 * j + 2], volElms[4 * j + 3]);
		}
		spdlog::info("mshPtNum = {}\n", mshPtNum);
		spdlog::info("volElmNum = {}\n", volElmNum);

		writePL3(strPL3, mshPtNum, mshPts, volElmNum, volElms);
		free(mshPts);
		free(volElms);
		mshPts = NULL;
		volElms = NULL;
	}

	int nMesh = API_TetrahedraObjNum();
	mshPtNum = API_GetTetrahedraPntNum(obj);
	volElmNum = API_GetTetrahedraElemNum(obj);
	// test 	
	API_DelTetrahedraObj(obj);
	int rtn = API_DelTetrahedraObj(obj);

	return errCode;
}
#endif

#ifdef TESTVTK
int main(int argc, char **argv)
{
	const int MAX_FILE_LEN = 1024;
	char strVtkIn[MAX_FILE_LEN], strVtkOut[MAX_FILE_LEN];
	int bndPtNum = 0;
	double *bndPts = NULL;
	int bndFctNum = 0;
	int *bndFcts = NULL;
	int intCstNum[] = { 0, 0, 0 };
	int *intCsts = NULL;
	double(*sizingFunc)(double x, double y, double z) = NULL;
	int mshPtNum = 0;
	double *mshPts = NULL;
	int volElmNum = 0;
	int *volElms = NULL;
	int errCode = 0;

	//printHead();

	if (argc < 2)
	{
		spdlog::info("Cannot initialize. Terminate\n");
		exit(1);
	}

	if (strlen(argv[1]) + 5 > MAX_FILE_LEN)
	{
		spdlog::info("too long file name %s", argv[1]);
		exit(1);
	}

	//	SPRlogTime_Init();
	//#ifdef _TIMING_PERFORMANCE
	//	time_start[TOTAL_EXEC_TIME] =  SPRlogTime();
	//#endif

	strcpy(strVtkIn, argv[1]);
	strcat(strVtkIn, ".vtk");
	strcpy(strVtkOut, argv[1]);
	strcat(strVtkOut, "_out.vtk");

	if (readVtk(strVtkIn, &bndPtNum, &bndPts, &bndFctNum, &bndFcts)) {
		spdlog::info("Read VTK failed! \n");
		return -1;
	}

	//errCode = API_Tetrahedralize_Single_Domain(
	//	bndPtNum, bndPts, bndFctNum, bndFcts, 
	//	intCstNum, NULL,
	//	NULL,
	//	&mshPtNum, &mshPts,
	//	&volElmNum, &volElms);

	int obj1, obj2;
	std::vector<int> vecObj;

	double vol;

	int inncst[3]{ 0,0,0 };
	errCode = API_Tetrahedralize_Single_Domain(bndPtNum, bndPts, bndFctNum, bndFcts, inncst, nullptr, nullptr, &obj1);
// 	errCode = API_Tetrahedralize_Multi_Domain (
// 		bndPtNum, bndPts, bndFctNum, bndFcts,
// 		&obj1);
	if (errCode != 0 && errCode != 1)
	{
		spdlog::info("Cannot tetrahedralize the domain.\n");
		return errCode;
	}
	vecObj.push_back(obj1);

	double x, y, z;
	int a, b, c, d;
	for (int i = 0; i < vecObj.size(); i++)
	{
		int obj = vecObj[i];
		mshPtNum = API_GetTetrahedraPntNum(obj);
		volElmNum = API_GetTetrahedraElemNum(obj);
		mshPts = (double *)malloc(mshPtNum * sizeof(double) * 3);
		volElms = (int *)malloc(volElmNum * sizeof(int) * 4);

		if (NULL == mshPts || NULL == volElms)
		{
			spdlog::info("Not enough memory.\n");
			errCode = -1;
			return errCode;
		}
		for (int j = 0; j < mshPtNum; j++)
		{
			API_GetTetrahedraPointCoord(obj, j + 1, mshPts[3 * j], mshPts[3 * j + 1], mshPts[3 * j + 2]);
		}

		for (int j = 0; j < volElmNum; j++)
		{
			API_GetTetrahedraElemPntIdx(obj, j + 1, volElms[4 * j], volElms[4 * j + 1], volElms[4 * j + 2], volElms[4 * j + 3]);
		}
		spdlog::info("mshPtNum = {}\n", mshPtNum);
		spdlog::info("volElmNum = {}\n", volElmNum);

		writeVtk(strVtkOut, mshPtNum, mshPts, volElmNum, volElms);
		free(mshPts);
		free(volElms);
		mshPts = NULL;
		volElms = NULL;

		int nMesh = API_TetrahedraObjNum();
		mshPtNum = API_GetTetrahedraPntNum(obj);
		volElmNum = API_GetTetrahedraElemNum(obj);
		// test 	
		API_DelTetrahedraObj(obj);
		int rtn = API_DelTetrahedraObj(obj);
	}

	return errCode;
}
#endif

#ifdef TESTMULTIDOMAIN
int main(int argc, char **argv)
{
	int errCode = 0;
	int obj;

	API_Tetrahedralize_Multi_Domain(2, nullptr, 2, nullptr, &obj);

	int mshPtNum = API_GetTetrahedraPntNum(obj);
	int volElmNum = API_GetTetrahedraElemNum(obj);
	double *mshPts = (double *)malloc(mshPtNum * sizeof(double) * 3);
	int *volElms = (int *)malloc(volElmNum * sizeof(int) * 4);

	if (NULL == mshPts || NULL == volElms)
	{
		spdlog::info("Not enough memory.\n");
		errCode = -1;
		return errCode;
	}
	for (int j = 0; j < mshPtNum; j++)
	{
		API_GetTetrahedraPointCoord(obj, j + 1, mshPts[3 * j], mshPts[3 * j + 1], mshPts[3 * j + 2]);
	}


	for (int j = 0; j < volElmNum; j++)
	{
		API_GetTetrahedraElemPntIdx(obj, j + 1, volElms[4 * j], volElms[4 * j + 1], volElms[4 * j + 2], volElms[4 * j + 3]);
	}
	spdlog::info("mshPtNum = {}\n", mshPtNum);
	spdlog::info("volElmNum = {}\n", volElmNum);

	writeVtk("./test.vtk", mshPtNum, mshPts, volElmNum, volElms);
	free(mshPts);
	free(volElms);
	mshPts = NULL;
	volElms = NULL;

	return errCode;
}
#endif

#ifdef  TEST_MULTITHREAD
int threadid = 0;
int finishedCnt = 0;
const int threadCnt = 10;
int g_argc;
char **g_argv;
std::mutex mtx;

int threadFunc() {
	mtx.lock();
	int tid = threadid++;
	mtx.unlock();
	const int MAX_FILE_LEN = 1024;
	char strVtkIn[MAX_FILE_LEN], strVtkOut[MAX_FILE_LEN];
	int bndPtNum = 0;
	double *bndPts = NULL;
	int bndFctNum = 0;
	int *bndFcts = NULL;
	int intCstNum[] = { 0, 0, 0 };
	int *intCsts = NULL;
	double(*sizingFunc)(double x, double y, double z) = NULL;
	int mshPtNum = 0;
	double *mshPts = NULL;
	int volElmNum = 0;
	int *volElms = NULL;
	int errCode = 0;

	if (g_argc < 2)
	{
		spdlog::info("Cannot initialize. Terminate\n");
		exit(1);
	}

	if (strlen(g_argv[1]) + 5 > MAX_FILE_LEN)
	{
		spdlog::info("too long file name %s", g_argv[1]);
		exit(1);
	}

	strcpy(strVtkIn, g_argv[1]);
	strcat(strVtkIn, ".vtk");
	strcpy(strVtkOut, g_argv[1]);
	std::string outstr = "_ out" + std::to_string(tid) + ".vtk";
	strcat(strVtkOut, outstr.c_str());

	if (readVtk(strVtkIn, &bndPtNum, &bndPts, &bndFctNum, &bndFcts)) {
		spdlog::info("Read VTK failed! \n");
		return -1;
	}

	int obj1, obj2;
	std::vector<int> vecObj;

	double vol;

	int inncst[3]{ 0,0,0 };
	errCode = API_Tetrahedralize_Single_Domain(bndPtNum, bndPts, bndFctNum, bndFcts, inncst, nullptr, nullptr, &obj1);

	if (errCode != 0 && errCode != 1)
	{
		spdlog::info("Cannot tetrahedralize the domain.\n");
		return errCode;
	}
	vecObj.push_back(obj1);

	double x, y, z;
	int a, b, c, d;
	for (int i = 0; i < vecObj.size(); i++)
	{
		int obj = vecObj[i];
		mshPtNum = API_GetTetrahedraPntNum(obj);
		volElmNum = API_GetTetrahedraElemNum(obj);
		mshPts = (double *)malloc(mshPtNum * sizeof(double) * 3);
		volElms = (int *)malloc(volElmNum * sizeof(int) * 4);

		if (NULL == mshPts || NULL == volElms)
		{
			spdlog::info("Not enough memory.\n");
			errCode = -1;
			return errCode;
		}
		for (int j = 0; j < mshPtNum; j++)
		{
			API_GetTetrahedraPointCoord(obj, j + 1, mshPts[3 * j], mshPts[3 * j + 1], mshPts[3 * j + 2]);
		}

		for (int j = 0; j < volElmNum; j++)
		{
			API_GetTetrahedraElemPntIdx(obj, j + 1, volElms[4 * j], volElms[4 * j + 1], volElms[4 * j + 2], volElms[4 * j + 3]);
		}
		spdlog::info("mshPtNum = {}\n", mshPtNum);
		spdlog::info("volElmNum = {}\n", volElmNum);

		writeVtk(strVtkOut, mshPtNum, mshPts, volElmNum, volElms);
		free(mshPts);
		free(volElms);
		mshPts = NULL;
		volElms = NULL;

		int nMesh = API_TetrahedraObjNum();
		mshPtNum = API_GetTetrahedraPntNum(obj);
		volElmNum = API_GetTetrahedraElemNum(obj);
		// test 	
		API_DelTetrahedraObj(obj);
		int rtn = API_DelTetrahedraObj(obj);
	}

	mtx.lock();
	finishedCnt++;
	spdlog::info("finish thread " << tid << " !");
	mtx.unlock();
	return 0;
}

int main(int argc, char **argv)
{
	g_argv = argv;
	g_argc = argc;

	std::vector<std::thread> vthread;

	for (int i = 0; i < threadCnt; i++) {
		vthread.emplace_back(threadFunc);
	}

	while (true) {
		mtx.lock();
		if(finishedCnt == threadCnt) break;
		mtx.unlock();
	}

	for (auto& t : vthread) {
		t.join();
	}

	return 0;
}


#endif //  TEST_MULTITHREAD

#ifdef TEST_TIGER2
int main(int argc, char** argv)
{
	char infilename[1024]; 
	char outfilename[1024];
	strcpy(infilename, argv[1]);
	strcpy(outfilename, argv[1]);
	strcat(infilename, ".vtk");
	strcat(outfilename, "_out.vtk");
	int inBndPntNum;
	double* inBndPnts;
	int inBndFctNum;
	int* inBndFcts;
	ZoneAttachInfo* inZoneAttachInfo = new ZoneAttachInfo;
	readVtk(infilename, &inBndPntNum, &inBndPnts, &inBndFctNum, &inBndFcts, inZoneAttachInfo);

	SurfMesh_t* surfMesh = new SurfMesh_t;
	API_CreateSurfaceMeshObj(inBndPntNum, inBndPnts, inBndFctNum, inBndFcts, surfMesh);

	delete[]inBndPnts;
	delete[]inBndFcts;

	if (inZoneAttachInfo->valid) {// domainTag tri2domain
		int* domainTag = new int[inZoneAttachInfo->tri2Domain.size()];
		for (int i = 0; i < inZoneAttachInfo->tri2Domain.size(); i++) {
			domainTag[i] = inZoneAttachInfo->tri2Domain[i];
		}
		// originSurfNums tri2FacesOffsets
		int* originSurfNums = new int[inBndFctNum];
		for (int i = 0; i < inBndFctNum; i++) {
			if (i == inBndFctNum - 1)
				originSurfNums[i] = inZoneAttachInfo->tri2Faces.size() - inZoneAttachInfo->tri2FaceOffsets[i];
			else
				originSurfNums[i] = inZoneAttachInfo->tri2FaceOffsets[i + 1] - inZoneAttachInfo->tri2FaceOffsets[i];
		}
		// originSurfs tri2Faces
		int* originSurfs = new int[inZoneAttachInfo->tri2Faces.size()];
		for (int i = 0; i < inZoneAttachInfo->tri2Faces.size(); i++)
			originSurfs[i] = inZoneAttachInfo->tri2Faces[i];
		// blockSurfNums block2FacesOffset
		int blockNum = inZoneAttachInfo->block2FaceOffsets.size();
		int* blockSurfNums = new int[blockNum];
		for (int i = 0; i < blockNum; i++) {
			if (i == blockNum - 1)
				blockSurfNums[i] = inZoneAttachInfo->block2Faces.size() - inZoneAttachInfo->block2FaceOffsets[i];
			else
				blockSurfNums[i] = inZoneAttachInfo->block2FaceOffsets[i + 1] - inZoneAttachInfo->block2FaceOffsets[i];
		}
		// blockSurfs block2Faces
		int* blockSurfs = new int[inZoneAttachInfo->block2Faces.size()];
		for (int i = 0; i < inZoneAttachInfo->block2Faces.size(); i++) {
			blockSurfs[i] = inZoneAttachInfo->block2Faces[i];
		}
		// pecBlocks
		int pecNum = inZoneAttachInfo->pecBlocks.size();
		int* pecBlocks = new int[pecNum];
		for (int i = 0; i < pecNum; i++) {
			pecBlocks[i] = inZoneAttachInfo->pecBlocks[i];
		}
		// sizeBlocks and sizes
		// blockSize
		int sizeNum = 0;
		std::vector<int> vecSizeBlocks;
		std::vector<double> vecSizes;
		for (int i = 0; i < inZoneAttachInfo->blockSize.size(); i++) {
			if (inZoneAttachInfo->blockSize[i] < 0)
				continue;
			else {
				sizeNum++;
				vecSizeBlocks.push_back(i);
				vecSizes.push_back(inZoneAttachInfo->blockSize[i]);
			}
		}
		int* sizeBlocks = new int[sizeNum];
		double* sizes = new double[sizeNum];
		for (int i = 0; i < sizeNum; i++) {
			sizeBlocks[i] = vecSizeBlocks[i];
			sizes[i] = vecSizes[i];
		}
		API_AddSurfaceOriginInformation(*surfMesh, domainTag,
			originSurfNums, originSurfs,
			blockNum, blockSurfNums, blockSurfs,
			pecNum, pecBlocks,
			sizeNum, sizeBlocks, sizes);
	}

	VolMesh_t* volMesh = new VolMesh_t;
	API_Tetrahedralize(*surfMesh, volMesh);

	int outPntNum = API_GetTetrahedraPntNum_High(*volMesh);
	//int outPntNum = API_GetTetrahedraPntNum(*volMesh, true);
	double* outPnts = new double[3 * outPntNum];
	for (int i = 0; i < outPntNum; i++){
		API_GetTetrahedraPointCoord(*volMesh, i, outPnts[3 * i], outPnts[3 * i + 1], outPnts[3 * i + 2]);
	}
	int outElmNum = API_GetTetrahedraElemNum(*volMesh);
	//int* outElm = new int[4 * outElmNum];
	int* outElm = new int[10 * outElmNum];
	for (int i = 0; i < outElmNum; i++) {
		//API_GetTetrahedraElemPntIdx(*volMesh, i, outElm[4 * i], outElm[4 * i + 1], outElm[4 * i + 2], outElm[4 * i + 3]);
		int tempElm[10];
		API_GetTetrahedraElemPntIdx_High(*volMesh, i, tempElm);
		for (int j = 0; j < 10; j++)
			outElm[10 * i + j] = tempElm[j];
	}

# if 1
	int outSurfTriNum = API_GetTetrahedraBndElemNum(*volMesh);
	int* outTriElems = new int[outSurfTriNum * 3];
	int* outTri2Domain = new int[outSurfTriNum];
	int* outTri2FaceNums = new int[outSurfTriNum];
	std::vector<int> vecTri2Face;
	for (int i = 0; i < outSurfTriNum; i++) {
		API_GetTetrahedraBndElemPts(*volMesh, i, outTriElems[3 * i], outTriElems[3 * i + 1], outTriElems[3 * i + 2]);
		outTri2Domain[i] = API_GetTetrahedraBndElemSubdomain(*volMesh, i);
		outTri2FaceNums[i] = API_GetTetrahedraBndElemSurfaceNum(*volMesh, i);
		for (int j = 0; j < outTri2FaceNums[i]; j++)
			vecTri2Face.push_back(API_GetTetrahedraBndElemSurface(*volMesh, i, j));
	}
	int* outTri2Face = new int[vecTri2Face.size()];
	for (int i = 0; i < vecTri2Face.size(); i++)
		outTri2Face[i] = vecTri2Face[i];
	int* outElem2Subdomain = new int[outElmNum];
	for (int i = 0; i < outElmNum; i++)
		outElem2Subdomain[i] = API_GetTetrahedraElemSubdomain(*volMesh, i);
	int outSubdomainNum = API_GetTetrahedraSubdomainNum(*volMesh);
	int* outSubdomainBlockNums = new int[outSubdomainNum];
	std::vector<int> vecSubdomain2Block;
	for (int i = 0; i < outSubdomainNum; i++) {
		outSubdomainBlockNums[i] = API_GetTetrahedraSubdomainBlockNum(*volMesh, i);
		for (int j = 0; j < outSubdomainBlockNums[i]; j++)
			vecSubdomain2Block.push_back(API_GetTetrahedraSubdomainBlock(*volMesh, i, j));
	}
	int* outSubdomain2Block = new int[vecSubdomain2Block.size()];
	for (int i = 0; i < vecSubdomain2Block.size(); i++)
		outSubdomain2Block[i] = vecSubdomain2Block[i];


	writeVtk_CenterPoint_new(outfilename, outPntNum, outPnts, outElmNum, outElm,
		outSurfTriNum, outTriElems, false,
		outTri2Face, outTri2FaceNums,
		outTri2Domain,
		outElem2Subdomain,
		outSubdomainNum, outSubdomain2Block, outSubdomainBlockNums);
# endif
	//writeVtk(outfilename, outPntNum, outPnts, outElmNum, outElm);
	//writeVtk_CenterPoint(outfilename, outPntNum, outPnts, outElmNum, outElm);
	return 0;
}
#endif