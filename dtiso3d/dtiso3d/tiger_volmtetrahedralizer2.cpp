#include "tiger_volmtetrahedralizer2.h"
#include "iso3d.h"
#include "tiger_volmtet.h"
#include "iso3d_multidomain.h"
#include <spdlog/spdlog.h> 
 #include <map>
#include <utility>
#include <memory>
#include <mutex>
#include "fileio.h"
#include <vector>
#include <unordered_map>
#include "mmg/mmgs/libmmgs.h"
#include "mmg/mmg3d/libmmg3d.h"
namespace VOLM_TET
{
	std::map<int, TigerMesh> mapMesh;
	static int meshObj = 0;
	std::mutex mutex;
	inline TigerMesh* getMeshObj(int obj) {
		if (mapMesh.find(obj) == mapMesh.end())
			return nullptr;
		else
			return &(mapMesh[obj]);
	}
}

DECL_VOLTET int API_CreateSurfaceMeshObj(int bndPtNum, double bndPts[], int bndFctNum, int bndFcts[], SurfMesh_t* surfMeshObj)
{
	VOLM_TET::mutex.lock();
	*surfMeshObj = ++VOLM_TET::meshObj;
	auto& meshobj = VOLM_TET::mapMesh[*surfMeshObj];
	VOLM_TET::mutex.unlock();
	meshobj.setCoords(bndPtNum, bndPts);
	meshobj.initTri(bndFctNum);
	for (int fc = 0; fc < bndFctNum; fc++) {
		meshobj.setTri(fc, bndFcts + 3 * fc);
	}
	spdlog::info("Successful!\n");
	return 0;
}

DECL_VOLTET int API_CreateVolumeMeshObj(int ptNum, double pts[], int bndFctNum, int bndFcts[], int tetNum, int tets[], TigerMesh_t *volMeshObj){
	VOLM_TET::mutex.lock();
	*volMeshObj = ++VOLM_TET::meshObj;
	auto &meshobj = VOLM_TET::mapMesh[*volMeshObj];
	VOLM_TET::mutex.unlock();


	meshobj.initTri(bndFctNum);

	if (ptNum){
		meshobj.setCoords(ptNum, pts);		
	}

	if (tetNum){
		meshobj.setTets(tetNum, tets);		
	}

	for (int fc = 0; fc < bndFctNum; fc++){
		meshobj.setTri(fc, bndFcts + 3*fc);
	}

	spdlog::info("Successful!\n");	
	return 0;
}

DECL_VOLTET int API_AddSurfaceOriginInformation(SurfMesh_t surfMeshObj, int *domainTag, int *originSurfNums, int *originSurfs, int blockNum, int *blockSurfNums, int *blockSurfs, int pecNum, int *pecBlocks, int sizeNum, int *sizeBlocks, double* sizes)
{
	std::shared_ptr<ZoneAttachInfo> info = std::make_shared<ZoneAttachInfo>();
	int ntri = VOLM_TET::mapMesh[surfMeshObj].getSurTriNum();
	int npnt = VOLM_TET::mapMesh[surfMeshObj].getPntNum(false);
	// if (!ntri || !npnt) return -1;
	info->tri2Domain.insert(info->tri2Domain.begin(), domainTag, domainTag + ntri);
	spdlog::info("tri domains: {} \n", info->tri2Domain.size());
	
	int off = 0;
	for (int i = 0; i < ntri+1; i++) {
		info->tri2FaceOffsets.push_back(off);
		off += originSurfNums[i];
	}
	info->tri2Faces.insert(info->tri2Faces.begin(), originSurfs, originSurfs + info->tri2FaceOffsets.back());

	off = 0;
	for (int i=0;i<blockNum+1;i++)
	{
		info->block2FaceOffsets.push_back(off);
		off += blockSurfNums[i];
	}
	info->block2Faces.insert(info->block2Faces.begin(), blockSurfs, blockSurfs + info->block2FaceOffsets.back());
	info->pecBlocks.insert(info->pecBlocks.begin(), pecBlocks, pecBlocks + pecNum);
	
	info->blockSize = std::vector<double>(blockNum, -1);
	for (int i = 0; i < sizeNum; i++) {
		info->blockSize[sizeBlocks[i]] = sizes[i];
	}

	info->valid = true;

	VOLM_TET::mapMesh[surfMeshObj].setZoneInfo(info);
	return 0;
}

DECL_VOLTET int API_Tetrahedralize(SurfMesh_t surfMeshObj, VolMesh_t* volMeshObj)
{
	spdlog::info("Start Tetrahedralizing...\n");
	fflush(stdout);
	std::unique_ptr<DTIso3D> pIso3d = std::make_unique<DTIso3D>();
	DTIso3D& iso3d = *pIso3d;
	int errCode = 0;

	/*
	 * 获取机器精度
	 */
	GEOM_FUNC::exactinit_threadSafe();
	/*
	 * 设置SPR操作的参数
	 */
	iso3d.sprimpl->spr_time_limit = 1.0; /* 时间限制，=0.0表示无限制 */

	auto& surfMeshobj = VOLM_TET::mapMesh[surfMeshObj];
	/*
	 * initialization
	 */
	iso3d.setZoneInfo(surfMeshobj.getZoneInfo());
	if (iso3d.crtSurface(surfMeshobj.getPntNum(false), surfMeshobj.getCoords(), surfMeshobj.getSurTriNum(), surfMeshobj.getTris()) <= 0)
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
#ifdef _RECV_EDGE_BY_MID_SP
	iso3d.recvBndEdge_MidStnPnt();
#else // INSERT_MID_SP
	iso3d.recvBndEdge_NoStnPnt_V2();
#endif
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
	//TODO1 启用约束恢复
	//iso3d.recvDesEdges();
	//iso3d.recvDesFacets();

	/* Attempt to remove Steiner points */
	for (int i = 0; i < 1; i++)
	{
		iso3d.removeInnerSteinerPoint_SPR();	/* 利用空腔重连算法消除Steiner点 */
		iso3d.removeInnerSteinerPoint_Flip();
	}
	iso3d.printSteinerPntInfo();

	/* Classify those deleted elements */
	iso3d.typeDelEles();

	/* Insert field points */
	//#TODO0 要开启该功能之前，首先要实现插内部点时对于边界子面的保护（目前已实现对于原边界面的保护，未考虑conforming的情况）
	//iso3d.innerPntInst();


	/* Attempt to remove Steiner points again */
	iso3d.prepareRmvSteinerPntAftInstInnNode();
	iso3d.removeInnerSteinerPoint_SPR();	/* 利用空腔重连算法消除Steiner点 */
	iso3d.removeEdgeSteinerPoint_SPR(false);
	iso3d.removeFacSteinerPoint_SPR(false);
	//	iso3d.removeAllLeftNodes();
	iso3d.printSteinerPntInfo();

	/* Remove delete nodes & elements & add parent info. */
	iso3d.rmvNodsAndEles();
	iso3d.fillParents(true);

	/* ----------------------------------------------------------------------
	 * Now, we start to improve mesh quality by local operators
	 * --------------------------------------------------------------------*/
	//#UP20210611 启用网格优化
	//iso3d.meshopt(3);

	iso3d.solveZone();

	/* Again, we remove delete nodes & elements & add parent info. */
	//iso3d.rmvNodsAndEles();
	//iso3d.fillParents();

	/* finally, we fill in the result */
	VOLM_TET::mutex.lock();
	*volMeshObj = ++VOLM_TET::meshObj;
	auto& meshobj = VOLM_TET::mapMesh[*volMeshObj];
	VOLM_TET::mutex.unlock();

	errCode = iso3d.outVolMsh(&meshobj);
	if (errCode != 0 && errCode != 1)
	{
		API_DelTetrahedraObj(*volMeshObj);
		return errCode;
		//goto FAIL;
	}

	meshobj.testQuality();

	errCode = iso3d.outSurMsh(&meshobj);
	if (errCode != 0 && errCode != 1)
	{
		API_DelTetrahedraObj(*volMeshObj);
		return errCode;
		//goto FAIL;
	}

	errCode = iso3d.outVolMshZoneInfo(&meshobj);
	if (errCode != 0 && errCode != 1)
	{
		API_DelTetrahedraObj(*volMeshObj);
		return errCode;
		//goto FAIL;
	}
	//目前自动生成中点信息
	meshobj.initMidPoint();
	goto SUCC;
FAIL:
	// do nothing presently
SUCC:
	return errCode;

}

DECL_VOLTET int API_GetTetrahedraPntNum(VolMesh_t objVolmTet)
{
	auto obj = VOLM_TET::getMeshObj(objVolmTet);
	if (!obj) return -1;
	else return obj->getPntNum(false);
}

DECL_VOLTET int API_GetTetrahedraPntNum_High(TigerMesh_t objVolmTet)
{
	auto obj = VOLM_TET::getMeshObj(objVolmTet);
	if (!obj) return -1;
	else return obj->getPntNum(true);
}

DECL_VOLTET int API_GetTetrahedraElemNum(VolMesh_t objVolmTet)
{
	auto obj = VOLM_TET::getMeshObj(objVolmTet);
	if (!obj) return -1;
	else return obj->getTetNum();
}

DECL_VOLTET int API_GetTetrahedraPointCoord(VolMesh_t objVolmTet, int index, double &x, double &y, double &z)
{
	auto obj = VOLM_TET::getMeshObj(objVolmTet);
	if (!obj) return -1;
	obj->getCoord(index, x, y, z);
	return 0;
}

DECL_VOLTET int API_GetTetrahedraElemPntIdx(VolMesh_t objVolmTet, int index, int &a, int &b, int &c, int &d)
{
	auto obj = VOLM_TET::getMeshObj(objVolmTet);
	if (!obj) return -1;
	obj->getTet(index+1, a,b,c,d);
	return 0;
}

DECL_VOLTET int API_GetTetrahedraElemPntIdx_High(VolMesh_t objVolmTet, int index, int pts[10])
{
	auto obj = VOLM_TET::getMeshObj(objVolmTet);
	if (!obj) return -1;
	obj->getTet(index+1, pts[0], pts[1], pts[2], pts[3]);
	obj->getTetMP(index, pts + 4);
	return 0;
}

DECL_VOLTET int API_GetTetrahedraSubdomainNum(VolMesh_t objVolmTet)
{
	auto obj = VOLM_TET::getMeshObj(objVolmTet);
	if (!obj) return -1;
	return obj->getZoneOut()->domain2Blocks.size();
}

DECL_VOLTET int API_GetTetrahedraSubdomainBlockNum(VolMesh_t objVolmTet, int index)
{
	auto obj = VOLM_TET::getMeshObj(objVolmTet);
	if (!obj) return -1;
	return obj->getZoneOut()->domain2Blocks[index].size();
}

DECL_VOLTET int API_GetTetrahedraSubdomainBlock(VolMesh_t objVolmTet, int index, int blkIndex)
{
	auto obj = VOLM_TET::getMeshObj(objVolmTet);
	if (!obj) return -1;
	return obj->getZoneOut()->domain2Blocks[index][blkIndex];
}

DECL_VOLTET int API_GetTetrahedraElemSubdomain(VolMesh_t objVolmTet, int index)
{
	auto obj = VOLM_TET::getMeshObj(objVolmTet);
	if (!obj) return -1;
	return obj->getZoneOut()->elem2Domains[index];
}

DECL_VOLTET int API_GetTetrahedraBndElemNum(VolMesh_t objVolmTet)
{
	auto obj = VOLM_TET::getMeshObj(objVolmTet);
	if (!obj) return -1;
	return obj->getSurTriNum();
}

DECL_VOLTET int API_GetTetrahedraBndElemPts(VolMesh_t objVolmTet, int index, int& p1, int& p2, int& p3)
{
	auto obj = VOLM_TET::getMeshObj(objVolmTet);
	if (!obj) return -1;
	return obj->getTri(index, p1, p2, p3);
}

DECL_VOLTET int API_GetTetrahedraBndElemSubdomain(VolMesh_t objVolmTet, int index)
{
	auto obj = VOLM_TET::getMeshObj(objVolmTet);
	if (!obj) return -1;
	return obj->getZoneOut()->tri2Domain[index];
}

DECL_VOLTET int API_GetTetrahedraBndElemSurfaceNum(VolMesh_t objVolmTet, int index)
{
	auto obj = VOLM_TET::getMeshObj(objVolmTet);
	if (!obj) return -1;
	return obj->getZoneOut()->tri2FaceOffsets[index + 1] - obj->getZoneOut()->tri2FaceOffsets[index];
}

DECL_VOLTET int API_GetTetrahedraBndElemSurface(VolMesh_t objVolmTet, int index, int surfIndex)
{
	auto obj = VOLM_TET::getMeshObj(objVolmTet);
	if (!obj) return -1;
	return obj->getZoneOut()->tri2Face[obj->getZoneOut()->tri2FaceOffsets[index] + surfIndex];
}

// DECL_VOLTET int API_InsertPoint(VolMesh_t volMeshObj, double coord[3])
// {
// 	//Do Nothing
// 	return 0;
// }

//TODO:Add attach information
DECL_VOLTET int API_Adaptation(
	TigerMesh_t volMeshObjIn,
	TigerMesh_t *volMeshObjOut,
	double sizingValues[],
	int requiredTetNum,
	int requiredTets[],
	int requiredFctNum,
	int requiredFcts[]){

	MMG5_pMesh mmgMesh = nullptr;
	MMG5_pSol mmgSol = nullptr;
	MMG3D_Init_mesh(MMG5_ARG_start,
				 	MMG5_ARG_ppMesh, &mmgMesh, MMG5_ARG_ppMet, &mmgSol,
					 MMG5_ARG_end);

	auto &objIn = VOLM_TET::mapMesh[volMeshObjIn];
	int nv = objIn.getPntNum(false);
	double *coords = objIn.getCoords();
	
	int ntet = objIn.getTetNum();
	int *tets = objIn.getTets(); 
	
	int ntri = objIn.getSurTriNum();
	int *tris = objIn.getTris();

	std::shared_ptr<ZoneAttachInfo> zoneInfo = objIn.getZoneInfo();
	std::shared_ptr<ZoneOutInfo>  zoneOut  = objIn.getZoneOut();
	bool hasZone = false;
	std::vector<int> tri2Domain;
	std::vector<int> elem2Domain;
	std::vector<std::vector<int>> domain2Blocks;
	std::unordered_map<int, std::vector<int>> domain2Surface;
	if (!zoneOut->tri2Domain.empty()){
		hasZone = true;
		for (int i=0; i<zoneOut->tri2Domain.size(); i++){
			int domain = zoneOut->tri2Domain[i];
			tri2Domain.push_back(domain);
			if (!domain2Surface.count(domain)){
				for (int j = zoneOut->tri2FaceOffsets[i]; j<zoneOut->tri2FaceOffsets[i+1]; j++){
					domain2Surface[domain].push_back(zoneOut->tri2Face[j]);				
				}

			}
		}

		for(auto d: zoneOut->elem2Domains){
			elem2Domain.push_back(d);
		}
		domain2Blocks = zoneOut->domain2Blocks;
		
	}

	MMG3D_Set_meshSize(mmgMesh, nv, ntet, 0, ntri, 0, 0);
	for (int i=0; i<nv; i++){
		int base = 3*i;
		MMG3D_Set_vertex(mmgMesh, coords[3*i], coords[3*i+1], coords[3*i+2], 0, i+1);
	}
	
	for (int i=0; i<ntet; i++){
		int base = 4*i;
		int ref = 0;
		if(hasZone){
			ref = elem2Domain[i];
		}
		MMG3D_Set_tetrahedron(mmgMesh, tets[base]+1, tets[base+1]+1, tets[base+2]+1, tets[base+3]+1, ref, i+1);
	}

	for (int i=0; i<ntri; i++){
		int base = i*3;
		int ref = 0;
		if (hasZone){
			ref = tri2Domain[i];
		}
		MMG3D_Set_triangle(mmgMesh, tris[base]+1, tris[base+1]+1, tris[base+2]+1, ref, i+1);
	}

	MMG3D_Set_solSize(mmgMesh, mmgSol, MMG5_Vertex, nv, MMG5_Scalar);
	for (int i=0; i<nv; i++){
		MMG3D_Set_scalarSol(mmgSol, sizingValues[i], i+1);		
	}

	for(int i=0; i<requiredFctNum; i++){
		MMG3D_Set_requiredTriangle(mmgMesh, requiredFcts[i]+1);
	}

	for(int i=0; i<requiredTetNum; i++){
		MMG3D_Set_requiredTetrahedron(mmgMesh, requiredTets[i]+1);
	}

	MMG3D_mmg3dlib(mmgMesh, mmgSol);



	VOLM_TET::mutex.lock();
	*volMeshObjOut = ++VOLM_TET::meshObj;
	auto &meshobjOut = VOLM_TET::mapMesh[*volMeshObjOut];
	VOLM_TET::mutex.unlock();
	int nnv;
	int nntet;
	int nntri;
	int nnprism;
	int nnquad;
	int na;
	std::vector<int> newTri2Domain;
	std::vector<int> newElem2Domain;
	std::vector<int> newTri2Surface;
	std::vector<int> newTri2SurfaceOffset;
	newTri2SurfaceOffset.push_back(0);
	MMG3D_Get_meshSize(mmgMesh, &nnv, &nntet, &nnprism, &nntri, &nnquad, &na);

	meshobjOut.initPnt(nnv);
	meshobjOut.initTri(nntri);
	meshobjOut.initTet(nntet);
	for (int i=0; i<nnv; i++){
		double x, y, z;
		int ref;
		MMG3D_Get_vertex(mmgMesh, &x, &y, &z, &ref, nullptr, nullptr);
		meshobjOut.setCoord(i+1, x, y, z);
	}

	for (int i=0; i<nntri; i++){
		int n0, n1, n2, ref;
		MMG3D_Get_triangle(mmgMesh, &n0, &n1, &n2, &ref, nullptr);
		if (hasZone){
			newTri2Domain.push_back(ref);
		}
		meshobjOut.setTri(i, (n0-1), (n1-1), (n2-1));
	}

	for (int i=0; i<nntet; i++){
		int n0, n1, n2, n3, ref;
		MMG3D_Get_tetrahedron(mmgMesh, &n0, &n1, &n2, &n3, &ref, nullptr);
		if (hasZone){
			newElem2Domain.push_back(ref);
		}
		meshobjOut.setTet(i+1, n0-1, n1-1, n2-1, n3-1);
	}
	meshobjOut.initMidPoint();

	if (hasZone){
		meshobjOut.getZoneOut()->domain2Blocks = zoneOut->domain2Blocks;
		meshobjOut.getZoneOut()->tri2Domain = newTri2Domain;
		meshobjOut.getZoneOut()->elem2Domains = newElem2Domain;
		for (int i=0; i<newTri2Domain.size(); i++){
			int domain = newTri2Domain[i];
			newTri2SurfaceOffset.push_back(domain2Surface[domain].size() + newTri2SurfaceOffset[i]);
			for (auto f: domain2Surface[domain]){
				newTri2Surface.push_back(f);
			}
		}
		meshobjOut.getZoneOut()->tri2Face = newTri2Surface;
		meshobjOut.getZoneOut()->tri2FaceOffsets = newTri2SurfaceOffset;
	}
	return 0;
}



DECL_VOLTET int API_InsertPoints(TigerMesh_t volMeshObjIn, TigerMesh_t *volMeshObjOut, int addInsidePtNum, double addInsideCoords[], 
	int addBndPtNum, double addBndCoords[], int addBndFts[]){	
	spdlog::info("Start Points inserting...\n");
	fflush(stdout);	
	
	std::unique_ptr<DTIso3D> pIso3d = std::make_unique<DTIso3D>();
	DTIso3D& iso3d = *pIso3d;
	int errCode = 0;	

	/*
	 * 获取机器精度
	 */
	GEOM_FUNC::exactinit_threadSafe();
	/*
	 * 设置SPR操作的参数
	 */
	iso3d.sprimpl->spr_time_limit = 1.0; /* 时间限制，=0.0表示无限制 */

	auto &objIn = VOLM_TET::mapMesh[volMeshObjIn];
	// if (!objIn) return -1;



	std::vector<double> fakeCoords;//原始体网格的点+体内插入的点+面上插入的点
	double *originCoords = objIn.getCoords();
	int numOriginPnt = objIn.getPntNum(false);

	for (int i=0; i<numOriginPnt; i++){
		int base = 3*i;
		fakeCoords.emplace_back(originCoords[base]);
		fakeCoords.emplace_back(originCoords[base+1]);
		fakeCoords.emplace_back(originCoords[base+2]);
	}

	for(int i=0; i<addBndPtNum; i++){
		int base = 3*i;
		fakeCoords.emplace_back(addBndCoords[base]);
		fakeCoords.emplace_back(addBndCoords[base+1]);
		fakeCoords.emplace_back(addBndCoords[base+2]);
	}

	for (int i=0; i<addInsidePtNum; i++){
		int base = 3*i;
		fakeCoords.emplace_back(addInsideCoords[base]);
		fakeCoords.emplace_back(addInsideCoords[base+1]);
		fakeCoords.emplace_back(addInsideCoords[base+2]);
	}
	std::shared_ptr<ZoneAttachInfo> zoneInfo = objIn.getZoneInfo();
	std::unordered_map<int, std::vector<int>> domain2Surface;
	std::unordered_map<int, int> facetPnts2fakePnts;
	std::unordered_map<int, int> fakePnts2facetPnts;
	std::vector<std::vector<int>> triangleFacets;
	std::vector<int> triangleDomains;
	std::vector<std::vector<double>> facetPoints;
	int newIndex = 1;
	auto getNode=
	[&facetPnts2fakePnts, &fakePnts2facetPnts, &facetPoints, &fakeCoords, &newIndex]
	(int index){
		int rst = 0;
		if (fakePnts2facetPnts.count(index)){
			rst = fakePnts2facetPnts[index];
		}
		else{
			rst = newIndex++;
			fakePnts2facetPnts[index] = rst;
			facetPnts2fakePnts[rst] = index;
			int base = 3*index;
			facetPoints.emplace_back();
			facetPoints.back().push_back(fakeCoords[base]);
			facetPoints.back().push_back(fakeCoords[base+1]);
			facetPoints.back().push_back(fakeCoords[base+2]);		
		}
		return rst;
	};
	int numfct = objIn.getSurTriNum();
	int *fctTris = objIn.getTris();
	std::vector<bool> keepTris;

	for(int i=0; i<numfct; i++){
		int base = 3*i;
		triangleFacets.emplace_back();
		int domainTag = 0;
		if(zoneInfo->valid){
			domainTag = zoneInfo->tri2Domain[i];
		}
		triangleDomains.push_back(domainTag);
		if(zoneInfo->valid){
			if (!domain2Surface.count(domainTag)){

				for (int j=zoneInfo->tri2FaceOffsets[i]; j<zoneInfo->tri2FaceOffsets[i+1]; j++){
					domain2Surface[domainTag].push_back(zoneInfo->tri2Faces[j]);
				}			
			}				
		}
	
		
		triangleFacets.back().push_back(getNode(fctTris[base]));
		triangleFacets.back().push_back(getNode(fctTris[base+1]));
		triangleFacets.back().push_back(getNode(fctTris[base+2]));
		keepTris.push_back(1);
	}

	for(int i=0; i<addBndPtNum; i++){
		int triIndex = addBndFts[i];
		int addPntIndex  = numOriginPnt+i;
		keepTris[triIndex] = 0;
		int n0, n1, n2;
		n0 = triangleFacets[triIndex][0];
		n1 = triangleFacets[triIndex][1];
		n2 = triangleFacets[triIndex][2];
		int domainTag = triangleDomains[i];
		triangleDomains.push_back(domainTag);
		triangleFacets.emplace_back();
		triangleFacets.back().push_back(getNode(addPntIndex));
		triangleFacets.back().push_back(n0);
		triangleFacets.back().push_back(n1);

		triangleDomains.push_back(domainTag);
		triangleFacets.emplace_back();
		triangleFacets.back().push_back(getNode(addPntIndex));
		triangleFacets.back().push_back(n1);
		triangleFacets.back().push_back(n2);

		triangleDomains.push_back(domainTag);
		triangleFacets.emplace_back();
		triangleFacets.back().push_back(getNode(addPntIndex));
		triangleFacets.back().push_back(n2);
		triangleFacets.back().push_back(n0);
		
		keepTris.push_back(1);
		keepTris.push_back(1);
		keepTris.push_back(1);
	}
	int totalTris =0;
	for(auto i: keepTris){
		totalTris+=i;
	}
	spdlog::info("Start using mmgs!");
	MMG5_pMesh mmgMesh = nullptr;
	MMG5_pSol mmgSol = nullptr;
	MMGS_Init_mesh(MMG5_ARG_start,
					MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppMet,&mmgSol,
					MMG5_ARG_end);
	MMGS_Set_meshSize(mmgMesh, facetPoints.size(), triangleFacets.size(), 0);
	for(int i=0; i<facetPoints.size(); i++){
		MMGS_Set_vertex(mmgMesh, facetPoints[i][0], facetPoints[i][1], facetPoints[i][2], 0, i+1);
	}

	for(int i=0; i<keepTris.size(); i++){
		if (keepTris[i]){
			MMGS_Set_triangle(mmgMesh, triangleFacets[i][0], triangleFacets[i][1], triangleFacets[i][2], triangleDomains[i], i+1);
		}
	}

	MMGS_Set_iparameter(mmgMesh, mmgSol, MMGS_IPARAM_noinsert, 1);
	MMGS_Set_iparameter(mmgMesh, mmgSol, MMGS_IPARAM_nomove, 1);
	MMGS_mmgslib(mmgMesh, mmgSol);
	int np, nt, na;
	MMGS_Get_meshSize(mmgMesh, &np, &nt, &na);
	for(int i=0; i<np; i++){
		double x, y, z;
		int tmp0, tmp1, tmp2;
		MMGS_Get_vertex(mmgMesh, &x, &y, &z, &tmp0, &tmp1, &tmp2);
		int fakeIndex = facetPnts2fakePnts[i+1];
		int base = fakeIndex*3;
		fakeCoords[base] = x;
		fakeCoords[base+1] = y;
		fakeCoords[base+2] = z;
	}	
	std::vector<int> outputTris;
	std::vector<int> newDomainTags;
	for(int i=0; i<nt; i++){
		int n0, n1, n2, tag, tmp1;
		MMGS_Get_triangle(mmgMesh, &n0, &n1, &n2, &tag, &tmp1);
		newDomainTags.push_back(tag);
		outputTris.push_back(facetPnts2fakePnts[n0]);
		outputTris.push_back(facetPnts2fakePnts[n1]);
		outputTris.push_back(facetPnts2fakePnts[n2]);	
	}

	std::shared_ptr<ZoneAttachInfo> newZoneInfo = std::make_shared<ZoneAttachInfo>();
	if(zoneInfo->valid){
		newZoneInfo->block2FaceOffsets = objIn.getZoneInfo()->block2FaceOffsets;
		newZoneInfo->block2Faces = objIn.getZoneInfo()->block2Faces;
		newZoneInfo->blockSize = objIn.getZoneInfo()->blockSize;
		newZoneInfo->tri2FaceOffsets.push_back(0);
		for(int i=0; i<newDomainTags.size(); i++){
			int tag = newDomainTags[i];
			newZoneInfo->tri2Domain.push_back(tag);
			newZoneInfo->tri2FaceOffsets.push_back(domain2Surface[tag].size() + newZoneInfo->tri2FaceOffsets[i]);
			for(auto f: domain2Surface[tag]){
				newZoneInfo->tri2Faces.push_back(f);
			}
		}		
	}


	iso3d.setZoneInfo(newZoneInfo);
	if (iso3d.crtSurface(fakeCoords.size()/3, fakeCoords.data(), outputTris.size()/3, outputTris.data()) <= 0)
	{
		spdlog::info("fail to initialize the surface with additional points");
		errCode = 2;
		goto FAIL;
	};

	/* insert boundary nodes */
	iso3d.bndPntInst();

	/*
	 * recover boundaries without adding Steiner points
	 */
	iso3d.setBndProtectionFlag(true);
	iso3d.abstNodFirEle();	/* 获取首单元，建立哈希表 */
#ifdef _RECV_EDGE_BY_MID_SP
	iso3d.recvBndEdge_MidStnPnt();
#else // INSERT_MID_SP
	iso3d.recvBndEdge_NoStnPnt_V2(); 
#endif
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
	//TODO1 启用约束恢复
	//iso3d.recvDesEdges();
	//iso3d.recvDesFacets();

	/* Attempt to remove Steiner points */
	for (int i = 0; i < 1; i++)
	{
		iso3d.removeInnerSteinerPoint_SPR();	/* 利用空腔重连算法消除Steiner点 */
		iso3d.removeInnerSteinerPoint_Flip();
	}
	iso3d.printSteinerPntInfo();

	/* Classify those deleted elements */
	iso3d.typeDelEles();

	/* Remove delete nodes & elements & add parent info. */
	iso3d.rmvNodsAndEles();
	iso3d.fillParents(true);
	iso3d.solveZone();

	
	VOLM_TET::mutex.lock();
	*volMeshObjOut = ++VOLM_TET::meshObj;
	auto &meshobjOut = VOLM_TET::mapMesh[*volMeshObjOut];
	VOLM_TET::mutex.unlock();
	
	/* finally, we fill in the result */
	errCode = iso3d.outVolMsh(&meshobjOut);
	if (errCode != 0 && errCode != 1)
	{
		API_DelTetrahedraObj(volMeshObjIn);
		return errCode;
		//goto FAIL;
	}
	


	errCode = iso3d.outVolMshZoneInfo(&meshobjOut);
	if (errCode != 0 && errCode != 1)
	{
		API_DelTetrahedraObj(volMeshObjIn);
		return errCode;
		//goto FAIL;
	}
	//目前自动生成中点信息
	meshobjOut.initMidPoint();	
	
	goto SUCC;

	FAIL:
		// do nothing presently
	SUCC:
		return errCode;

}



DECL_VOLTET int API_DelTetrahedraObj(VolMesh_t objVolmTet)
{
	VOLM_TET::mutex.lock();
	if (VOLM_TET::mapMesh.find(objVolmTet) != VOLM_TET::mapMesh.end()) {
		VOLM_TET::mapMesh[objVolmTet].freeMemory();
		VOLM_TET::mapMesh.erase(objVolmTet);
		VOLM_TET::mutex.unlock();
		return 0;
	}
	VOLM_TET::mutex.unlock();
	return 1;
}

DECL_VOLTET int API_DelSurfMeshObj(SurfMesh_t objVolmTet)
{
	VOLM_TET::mutex.lock();
	if (VOLM_TET::mapMesh.find(objVolmTet) != VOLM_TET::mapMesh.end()) {
		VOLM_TET::mapMesh[objVolmTet].freeMemory();
		VOLM_TET::mapMesh.erase(objVolmTet);
		VOLM_TET::mutex.unlock();
		return 0;
	}
	VOLM_TET::mutex.unlock();
	return 0;
}

