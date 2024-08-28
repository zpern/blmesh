/** ----------------------------------------------------------------------------
	* TIGER, a Trustable Intelligent GridER
	*
	* TIGER的四面体网格生成API集成接口定义 (版本号：1.0)
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

#include <spdlog/spdlog.h> 
 #include "tiger_volmtetrahedralizer.h"
#include "iso3d.h"
#include "tiger_volmtet.h"
#include "iso3d_multidomain.h"
#include <map>
#include <utility>
#include <memory>
#include <mutex>
#include "fileio.h"
#include "singleton_terminate.h"

namespace VOLM_TET
{
	std::map<int, TigerVolmTet> mapVolmTet;
	static int volmTetObj = 0;
	std::mutex mutex;
}

/** -------------------------------------------------------------------------------
	* Given a surface triangulation defining a single-connected domain and a user
	* defined sizing function for element spacing, try to discretize the domain into a
	* combination of tetrahedra
	*
	* The inputs include three parts:
	* The first part is the group of boundary facets with all their normals pointing to
	* the interior of the domain
	*
	*	The second part is the interior constraints, i.e., points, edges and/or facets (if any).
	*
	* The third part of the input is a user defined sizing function, which is expected to
	* output a sizing value for an arbitrary point inside the problem domain
	*
	* The outputs are a sequence of tetrahedra filling in the domain, with the input boundary
	* & interior constraints respected
	*
	* Input Paramerters:
	*	bndPtNum		# of boundary points defining the exterior boundary & interior constraints.
	*	bndPts			the coordinates of boundary points. dim = 3*bndPtNum
	*	bndFctNum		# of exterior boundary facets
	* bndFcts			the connections of exterior boundary facets. dim = 3*bndFctNum
	* intCstNum		# of interior constraints. dim=3.
	*							intCstNum[0] -- # of interior points
	*							intCstNum[1] -- # of interior edges
	*							intCstNum[2] -- # of interior facets
	*	intCsts			the connections of interior constraints
	*							dim = intCstNum[0] + 2 * intCstNum[1] + 3*intCstNum[2]
	* sizingFunc	the user specified sizing function (if any)
	* Output Paramerters:
	*	mshPtNums		# of mesh points
	*	mshPts			the coordinates of mesh points. dim = 3*mshPtNums
	*	volElmNum		# of volume elements
	*	volElms			the connections of volume elements. dim = 4*volElmNum
	* Returned value:
	*	0				Completely succeed
	*	1				Partially succeed
	*	otherwise		Fail. Error info. will be defined later
	* ------------------------------------------------------------------------------*/
int API_Tetrahedralize_Single_Domain(
	int bndPtNum, double bndPts[],
	int bndFctNum, int bndFcts[],
	int intCstNum[], int intCsts[],
	double (*sizingFunc)(double x, double y, double z),
	int* mshPtNum, double** mshPts,
	int* volElmNum, int** volElms
)
{
	DTIso3D iso3d;
	int errCode = 0;

	/*
	 * 获取机器精度
	 */
	GEOM_FUNC::exactinit_threadSafe();
	/*
	 * 设置SPR操作的参数
	 */
	iso3d.sprimpl->spr_time_limit = 1.0; /* 时间限制，=0.0表示无限制 */

	if (intCstNum[0] > 0 || intCstNum[1] > 0 && intCstNum[2] > 0)
	{
		spdlog::info("The present version does not support inputs with interior constraints.\n");
		errCode = 2;
		goto FAIL;
	}

	if (sizingFunc != nullptr)
	{
		spdlog::info("The present version does not support an input sizing function. \n");
		errCode = 2;
		goto FAIL;
	}

	/*
	 * initialization
	 */
	if (iso3d.crtSurface(bndPtNum, bndPts, bndFctNum, bndFcts) <= 0)
	{
		spdlog::info("fail to initialize the surface");
		errCode = 2;
		goto FAIL;
	}

	/* insert boundary nodes */
	iso3d.bndPntInst();

	/*
	 * recover boundaries without adding Steiner points
	 */
	iso3d.setBndProtectionFlag(true);
	iso3d.abstNodFirEle();	/* 获取首单元，建立哈希表 */
	iso3d.recvBndEdge_NoStnPnt_V2();
	iso3d.recvBndFace_NoStnPnt_V2();

	/*
	 * recover boundaries by adding Steiner points
	 */
	iso3d.recvEdges();
	iso3d.recvFaces();

	/*
	 * classify elements as INNER/OUTER ELEMENTS
	 */
	iso3d.typeEles();

	/*
	 * Constrained recovery
	 */
	iso3d.prepareRecvDesBnds();
	iso3d.recvDesEdges();
	iso3d.recvDesFacets();

	/* Attempt to remove Steiner points */
	for (int i = 0; i < 1; i++)
	{
		//iso3d.removeInnerSteinerPoint_SPR();	/* 利用空腔重连算法消除Steiner点 */
		iso3d.removeInnerSteinerPoint_Flip();
	}
	iso3d.printSteinerPntInfo();

	/* Classify those deleted elements */
	iso3d.typeDelEles();

	/* Insert field points */
	iso3d.innerPntInst();

	/* Attempt to remove Steiner points again */
	//iso3d.prepareRmvSteinerPntAftInstInnNode();
	//iso3d.removeInnerSteinerPoint_SPR();	/* 利用空腔重连算法消除Steiner点 */
	//iso3d.removeEdgeSteinerPoint_SPR(false);
	//iso3d.removeFacSteinerPoint_SPR(false);
	//	iso3d.removeAllLeftNodes();
	//iso3d.printSteinerPntInfo();

	/* Remove delete nodes & elements & add parent info. */
	iso3d.rmvNodsAndEles();
	iso3d.fillParents(false);

	/* ----------------------------------------------------------------------
	 * Now, we start to improve mesh quality by local operators
	 * --------------------------------------------------------------------*/
	iso3d.meshopt(10);

	/* Again, we remove delete nodes & elements & add parent info. */
	iso3d.rmvNodsAndEles();
	iso3d.fillParents();

	/* finally, we fill in the result */
	errCode = iso3d.outVolMsh(mshPtNum, mshPts, volElmNum, volElms);
	if (errCode != 0 && errCode != 1)
		goto FAIL;

	goto SUCC;
FAIL:
	// do nothing presently
SUCC:
	return errCode;
}
int API_Tetrahedralize_Single_Domain(
	int bndPtNum, double bndPts[],
	int bndFctNum, int bndFcts[],
	int intCstNum[], int intCsts[],
	double(*sizingFunc)(double x, double y, double z),
	int* objVolmTet
) {
	return API_Tetrahedralize_Single_Domain(bndPtNum, bndPts, bndFctNum, bndFcts, intCstNum, intCsts, 5, sizingFunc, objVolmTet);
}

int API_Tetrahedralize_Single_Domain(
	int bndPtNum, double bndPts[],
	int bndFctNum, int bndFcts[],
	int intCstNum[], int intCsts[],
	int optimizeCount,
	double(*sizingFunc)(double x, double y, double z),
	int* objVolmTet
)
{
#ifdef INFO
	spdlog::info("Call API_Tetrahedralize_Single_Domain\n");
#endif // INFO
	std::shared_ptr<DTIso3D> iso3d = std::shared_ptr<DTIso3D>(new DTIso3D());
	int errCode = 0;

	/*
	 * 获取机器精度
	 */
	GEOM_FUNC::exactinit_threadSafe();
	/*
	 * 设置SPR操作的参数
	 */
	iso3d->sprimpl->spr_time_limit = 1.0; /* 时间限制，=0.0表示无限制 */

	if (intCstNum[0] > 0 || intCstNum[1] > 0 && intCstNum[2] > 0)
	{
		spdlog::info("The present version does not support inputs with interior constraints.\n");
		errCode = 2;
		return errCode;
		//goto FAIL;
	}

#ifdef USE_SIZING_FUNC
	if (sizingFunc)
	{
		iso3d->setSizingFunc(sizingFunc);
		spdlog::info("Use indicated sizing function. \n");
	}
	else
	{
		spdlog::info("No sizing function input. Use default sizing control logic. \n");
	}
#else
	if (sizingFunc != NULL)

	{
		spdlog::info("The present version does not support an input sizing function. \n");
		errCode = 2;
		return errCode;
		//goto FAIL;
	}
#endif

	/*
	 * initialization
	 */
#ifdef INFO
	spdlog::info("Begin to initialize the surface\n");
#endif // INFO

	if (iso3d->crtSurface(bndPtNum, bndPts, bndFctNum, bndFcts) <= 0)
	{
		spdlog::info("fail to initialize the surface");
		errCode = 2;
		return errCode;
		//goto FAIL;
	}
	//	iso3d->writeVisualVtkFile("start.vtk");
		/* insert boundary nodes */
	iso3d->bndPntInst();

	if (checkterminate()) {
		iso3d->rmvNodsAndEles();
		return -1;
	}
	/*
	 * recover boundaries without adding Steiner points
	 */
	 // #3_不插点边界恢复
	iso3d->setBndProtectionFlag(true);
	iso3d->abstNodFirEle();	/* 获取首单元，建立哈希表 */
	iso3d->recvBndEdge_NoStnPnt_V2();
	iso3d->recvBndFace_NoStnPnt_V2();

	if (checkterminate()) {
		iso3d->rmvNodsAndEles();
		return -1;
	}

	/*
	 * recover boundaries by adding Steiner points
	 */
	 // #4_插点边界恢复
	iso3d->recvEdges();
	iso3d->recvFaces();
	//iso3d->writeVisualVtkFile("recover.vtk");
	/*
	 * classify elements as INNER/OUTER ELEMENTS
	 */

	if (checkterminate()) {
		iso3d->rmvNodsAndEles();
		return -1;
	}
	iso3d->typeEles();

	/*
	 * Constrained recovery
	 */
	 // #5_约束边界恢复
	iso3d->prepareRecvDesBnds();
	iso3d->recvDesEdges();
	iso3d->recvDesFacets();

	if (checkterminate()) {
		iso3d->rmvNodsAndEles();
		return -1;
	}
	//iso3d->writeVisualVtkFile("recover2.vtk");
	/* Attempt to remove Steiner points */
	// #6_移除边界点
	for (int i = 0; i < 1; i++)
	{
		//iso3d->removeInnerSteinerPoint_SPR();	/* 利用空腔重连算法消除Steiner点 */
		iso3d->removeInnerSteinerPoint_Flip();
	}
	iso3d->printSteinerPntInfo();
	if (checkterminate()) {
		iso3d->rmvNodsAndEles();
		return -1;
	}
	/* Classify those deleted elements */
	iso3d->typeDelEles();

	// #7_内部插点
	// 先恢复边界，划分出内部区域，再插内部点，可以防止在边界外部插入无用点
	/* Insert field points */
	iso3d->innerPntInst();
	if (checkterminate()) {
		iso3d->rmvNodsAndEles();
		return -1;
	}
	// #8_移除边界点
	/* Attempt to remove Steiner points again */
	//iso3d->prepareRmvSteinerPntAftInstInnNode();
	//iso3d->removeInnerSteinerPoint_SPR();	/* 利用空腔重连算法消除Steiner点 */
	//iso3d->removeEdgeSteinerPoint_SPR(false);
	//iso3d->removeFacSteinerPoint_SPR(false);
	//	iso3d.removeAllLeftNodes();
	iso3d->printSteinerPntInfo();
	if (checkterminate()) {
		iso3d->rmvNodsAndEles();
		return -1;
	}
	//iso3d->writeVisualVtkFile("remove_outer.vtk");
	/* Remove delete nodes & elements & add parent info. */
	iso3d->rmvNodsAndEles();
	iso3d->fillParents(false);

	/* ----------------------------------------------------------------------
	 * Now, we start to improve mesh quality by local operators
	 * --------------------------------------------------------------------*/
	 // #9_网格优化
	iso3d->meshopt(optimizeCount);
	if (checkterminate()) {
		iso3d->rmvNodsAndEles();
		return -1;
	}

	/* Again, we remove delete nodes & elements & add parent info. */
	iso3d->rmvNodsAndEles();
	iso3d->fillParents();

	VOLM_TET::mutex.lock();
	VOLM_TET::volmTetObj++;
	TigerVolmTet* ptet = &(VOLM_TET::mapVolmTet[VOLM_TET::volmTetObj]);
	VOLM_TET::mutex.unlock();

	/* finally, we fill in the result */
	errCode = iso3d->outVolMsh(ptet);
	if (errCode != 0 && errCode != 1)
	{
		API_DelTetrahedraObj(VOLM_TET::volmTetObj);
		return errCode;
		//goto FAIL;
	}
	//iso3d->writeVisualVtkFile("final.vtk");
	*objVolmTet = VOLM_TET::volmTetObj;

	goto SUCC;
FAIL:
	// do nothing presently
SUCC:
	return errCode;
}

int getTestObj(int* obj1, int* obj2)
{
	VOLM_TET::volmTetObj++;
	*obj1 = VOLM_TET::volmTetObj;
	TigerVolmTet& tet1 = VOLM_TET::mapVolmTet[*obj1];
	VOLM_TET::volmTetObj++;
	*obj2 = VOLM_TET::volmTetObj;
	TigerVolmTet& tet2 = VOLM_TET::mapVolmTet[*obj2];

	tet1.initPnt(5);
	tet2.initPnt(5);
	tet1.initTet(4);
	tet2.initTet(4);

	tet1.setCoord(1, 0, 0, 0);
	tet1.setCoord(2, 4, 0, 0);
	tet1.setCoord(3, 0, 4, 0);
	tet1.setCoord(4, 0, 0, 4);
	tet1.setCoord(5, 1, 1, 1);
	tet1.setTet(1, 0, 1, 3, 4);
	tet1.setTet(2, 0, 1, 2, 4);
	tet1.setTet(3, 0, 3, 2, 4);
	tet1.setTet(4, 1, 2, 3, 4);

	tet2.setCoord(1, 0, 0, 0);
	tet2.setCoord(2, -4, 0, 0);
	tet2.setCoord(3, 0, 4, 0);
	tet2.setCoord(4, 0, 0, 4);
	tet2.setCoord(5, -1, 1, 1);
	tet2.setTet(1, 0, 1, 3, 4);
	tet2.setTet(2, 0, 1, 2, 4);
	tet2.setTet(3, 0, 3, 2, 4);
	tet2.setTet(4, 1, 2, 3, 4);

	return 0;
}

#ifdef USE_TETGEN
int API_Tetrahedralize_Multi_Domain(int bndPtNum, double bndPts[], int bndFctNum, int bndFcts[], int* objVolmTet)
{
	int domainNum;
	int* domainBndPtNum, * domainBndFctNum;
	double* domainBndPts;
	int* domainBndFcts;
	int* ltowBndPtMap;
	splitMultiDomain(bndPtNum, bndPts, bndFctNum, bndFcts, &domainNum, &domainBndPtNum, &domainBndPts, &domainBndFctNum, &domainBndFcts, &ltowBndPtMap);

	int* objTets = new int[domainNum];
	int ptIndexOffset = 0;
	int fcIndexOffset = 0;
	int intCstNum[] = { 0, 0, 0 };
	for (int i = 0; i < domainNum; i++) {
		API_Tetrahedralize_Single_Domain(domainBndPtNum[i], domainBndPts + ptIndexOffset, domainBndFctNum[i], domainBndFcts + fcIndexOffset, intCstNum, nullptr, nullptr, &(objTets[i]));
		ptIndexOffset += domainBndPtNum[i] * 3;
		fcIndexOffset += domainBndFctNum[i] * 3;
	}

	API_Merge_Multi_Domain(domainNum, objTets, domainBndPtNum, ltowBndPtMap, objVolmTet);

	//释放合并前四面体
	for (int i = 0; i < domainNum; i++) {
		API_DelTetrahedraObj(objTets[i]);
	}
	// #TODO1 free memory
	return 0;
}
#endif
int API_Merge_Multi_Domain(int domainNum, int domainVolmTets[], int* bndPtNum, int* ltowBndPtMap, int* objVolmTet)
{
	VOLM_TET::mutex.lock();
	VOLM_TET::volmTetObj++;
	if (VOLM_TET::mapVolmTet.find(VOLM_TET::volmTetObj) != VOLM_TET::mapVolmTet.end()) {
		return -1;
	}
	*objVolmTet = VOLM_TET::volmTetObj;
	TigerVolmTet& objTet = VOLM_TET::mapVolmTet[VOLM_TET::volmTetObj];
	VOLM_TET::mutex.unlock();

	// get total boundary point number
	int totBndPtNum = -1;
	std::vector<int> domainOff;
	int curoff = 0;
	for (int i = 0; i < domainNum; i++) {
		domainOff.push_back(curoff);
		curoff += bndPtNum[i];
		for (int j = 0; j < bndPtNum[i]; j++) {
			int curBndIndex = ltowBndPtMap[j + domainOff[i]];
			if (curBndIndex > totBndPtNum) {
				totBndPtNum = curBndIndex;
			}
		}
	}
	totBndPtNum++;

	// get total tet number and point number
	int totPtNum = totBndPtNum, totTetNum = 0;
	for (int i = 0; i < domainNum; i++) {
		totPtNum += API_GetTetrahedraPntNum(domainVolmTets[i]) - bndPtNum[i];
		totTetNum += API_GetTetrahedraElemNum(domainVolmTets[i]);
	}

	objTet.initPnt(totPtNum);
	objTet.initTet(totTetNum);

	// set boundary point coord
	std::unique_ptr<bool[]> bndPtHash(new bool[totBndPtNum]);
	for (int i = 0; i < totBndPtNum; i++) {
		bndPtHash[i] = false;
	}
	for (int i = 0; i < domainNum; i++) {
		for (int j = 0; j < bndPtNum[i]; j++) {
			int wIdx = ltowBndPtMap[j + domainOff[i]];
			if (bndPtHash[wIdx]) continue;
			bndPtHash[wIdx] = true;
			double x, y, z;
			API_GetTetrahedraPointCoord(domainVolmTets[i], j + 1, x, y, z);
			objTet.setCoord(wIdx + 1, x, y, z);
		}
	}
	bndPtHash.reset();

	// set inner point coord and element index
	int curFreeIdx = totBndPtNum;
	int curFreeTet = 0;
	auto lToW = [&domainOff](int localIndex, int* curBndPtMap, int domainIndex, int beginFreeIdx, int curBndPtNum) {
		if (localIndex >= curBndPtNum) {
			return localIndex - curBndPtNum + beginFreeIdx;
		}
		else {
			return curBndPtMap[domainOff[domainIndex] + localIndex];
		}
		};
	for (int i = 0; i < domainNum; i++) {
		int curTotPtNum = API_GetTetrahedraPntNum(domainVolmTets[i]);
		int curBndPtNum = bndPtNum[i];
		int beginFreeIdx = curFreeIdx;
		// set inner point coord
		for (int j = curBndPtNum; j < curTotPtNum; j++) {
			double x, y, z;
			API_GetTetrahedraPointCoord(domainVolmTets[i], j + 1, x, y, z);
			objTet.setCoord((curFreeIdx++) + 1, x, y, z);
		}

		// set tet element indices
		int curTotTetNum = API_GetTetrahedraElemNum(domainVolmTets[i]);
		for (int j = 0; j < curTotTetNum; j++) {
			int p1, p2, p3, p4;
			API_GetTetrahedraElemPntIdx(domainVolmTets[i], j + 1, p1, p2, p3, p4);
			p1 = lToW(p1, ltowBndPtMap, i, beginFreeIdx, curBndPtNum);
			p2 = lToW(p2, ltowBndPtMap, i, beginFreeIdx, curBndPtNum);
			p3 = lToW(p3, ltowBndPtMap, i, beginFreeIdx, curBndPtNum);
			p4 = lToW(p4, ltowBndPtMap, i, beginFreeIdx, curBndPtNum);
			objTet.setTet((curFreeTet++) + 1, p1, p2, p3, p4);
		}
	}

	return 0;
}

int API_GetTetrahedraPntNum(int objVolmTet)
{
	VOLM_TET::mutex.lock();
	auto iter = VOLM_TET::mapVolmTet.find(objVolmTet);
	if (iter == VOLM_TET::mapVolmTet.end())
	{
		printf("Cann't find the Teterahedral Object, please check the Object Index!\n");
		return 0;
	}

	//TigerVolmTet *pVolmTet = iter->second;
	TigerVolmTet& VolmTet = iter->second;
	VOLM_TET::mutex.unlock();

	VolmTet.lock();
	int rtn = VolmTet.getPntNum();
	VolmTet.unlock();
	return rtn;
}

int API_GetTetrahedraElemNum(int objVolmTet)
{
	VOLM_TET::mutex.lock();
	auto iter = VOLM_TET::mapVolmTet.find(objVolmTet);
	if (iter == VOLM_TET::mapVolmTet.end())
	{
		printf("Cann't find the Teterahedral Object, please check the Object Index!\n");
		return 0;
	}

	TigerVolmTet& VolmTet = iter->second;
	VOLM_TET::mutex.unlock();

	VolmTet.lock();
	int rtn = VolmTet.getTetNum();
	VolmTet.unlock();
	return rtn;
}

/* 每次都需要查找，速度可能有影响 */
int API_GetTetrahedraPointCoord(int objVolmTet, int index, double& x, double& y, double& z)
{
	VOLM_TET::mutex.lock();
	auto iter = VOLM_TET::mapVolmTet.find(objVolmTet);
	if (iter == VOLM_TET::mapVolmTet.end())
	{
		printf("Cann't find the Teterahedral Object, please check the Object Index!\n");
		return 1;
	}

	TigerVolmTet& VolmTet = iter->second;
	VOLM_TET::mutex.unlock();

	VolmTet.lock();
	VolmTet.getCoord(index, x, y, z);
	VolmTet.unlock();

	return  0;
}

/* 每次都需要查找，速度可能有影响 */
int API_GetTetrahedraElemPntIdx(int objVolmTet, int index, int& a, int& b, int& c, int& d)
{
	VOLM_TET::mutex.lock();
	auto iter = VOLM_TET::mapVolmTet.find(objVolmTet);
	if (iter == VOLM_TET::mapVolmTet.end())
	{
		printf("Cann't find the Teterahedral Object, please check the Object Index!\n");
		return 1;
	}

	TigerVolmTet& VolmTet = iter->second;
	VOLM_TET::mutex.unlock();

	VolmTet.lock();
	VolmTet.getTet(index, a, b, c, d);
	VolmTet.unlock();

	return  0;
}
/*
* @note: 修改：删除内存泄漏 --yhf
*/
int API_DelTetrahedraObj(int objVolmTet)
{
	VOLM_TET::mutex.lock();
	if (VOLM_TET::mapVolmTet.find(objVolmTet) != VOLM_TET::mapVolmTet.end()) {
		VOLM_TET::mapVolmTet[objVolmTet].freeMemory();
		VOLM_TET::mapVolmTet.erase(objVolmTet);
		VOLM_TET::mutex.unlock();
		return 0;
	}
	VOLM_TET::mutex.unlock();
	return 1;
}

int API_TetrahedraObjNum()
{
	return VOLM_TET::mapVolmTet.size();
}

/* --------------------------------------------------------------------------------
 * Compute the volume of the polyhedron defined by a set of facets
 * -------------------------------------------------------------------------------*/
int API_Compute_Poly_Volume(
	int bndPtNum, double bndPts[],
	int bndFctNum, int bndFcts[],
	double* vol
)
{
	DTIso3D iso3d;
	int errCode = 0;

	/*
	 * 获取机器精度
	 */
	GEOM_FUNC::exactinit_threadSafe();
	/*
	 * 设置SPR操作的参数
	 */
	iso3d.sprimpl->spr_time_limit = 1.0; /* 时间限制，=0.0表示无限制 */

	if (!vol)
	{
		spdlog::info("nullptr POINTER INPUT.\n");
		errCode = 2;
		goto FAIL;
	}
	/*
	 * initialization
	 */
	if (iso3d.crtSurface(bndPtNum, bndPts, bndFctNum, bndFcts) <= 0)
	{
		spdlog::info("fail to initialize the surface");
		errCode = 2;
		goto FAIL;
	}

	/* insert boundary nodes */
	iso3d.bndPntInst();

	/*
	 * recover boundaries without adding Steiner points
	 */
	iso3d.setBndProtectionFlag(true);
	iso3d.abstNodFirEle();	/* 获取首单元，建立哈希表 */
	iso3d.recvBndEdge_NoStnPnt_V2();
	iso3d.recvBndFace_NoStnPnt_V2();

	/*
	 * recover boundaries by adding Steiner points
	 */
	iso3d.recvEdges();
	iso3d.recvFaces();

	/*
	 * classify elements as INNER/OUTER ELEMENTS
	 */
	iso3d.typeEles();

	iso3d.computeInteriorVolume(vol);

	goto SUCC;
FAIL:
	// do nothing presently
SUCC:
	return errCode;
}