#ifndef __tiger_volmtetrahedralizer2_h__
#define __tiger_volmtetrahedralizer2_h__

#ifdef WIN32
#ifdef libvolmt_EXPORTS 
#define DECL_VOLTET __declspec(dllexport)
#else
#define DECL_VOLTET //__declspec(dllimport)
#endif libvolmt_EXPORTS
#else
#define DECL_VOLTET //__declspec(dllimport)
#endif WIN32

typedef int SurfMesh_t;
typedef int VolMesh_t;
typedef int TigerMesh_t;

/**
 *  @brief	创建面网格对象
 *  @param[in]		bndPtNum		面网格点数量
 *  @param[in]		bndPts			面网格点坐标. dim = 3*bndPtNum
 *  @param[in]		bndFctNum		面网格单元数量
 *  @param[in]		bndFcts			面网格单元索引. dim = 3*bndFctNum
 *  @param[out]		surfMeshObj		输出面网格对象
 *  @return 
 *	0				Completely succeed
 *	otherwise		Fail. Error info. will be defined later
 */
DECL_VOLTET int API_CreateSurfaceMeshObj(
	int bndPtNum, double bndPts[],
	int bndFctNum, int bndFcts[],
	TigerMesh_t* surfMeshObj
);

/**
 * @brief  创建体网格对象，便于多次自适应循环中可以使用
 * 
 * @param[in]		ptNum           网格点数量
 * @param[in]		pts				网格点坐标. dim = 3*pts
 * @param[in]		bndFctNum		面网格单元数量
 * @param[in]		bndFcts			面网格单元索引. dim = 3*bndFctNum
 * @param[in]		tetNum			四面体网格单元数目
 * @param[in]		tets 			四面体网格单元索引. dim = 4*tetNum
 * @param[out]		volMeshObj 		输出体网格对象
 *	0				Completely succeed
 *	otherwise		Fail. Error info. will be defined later
 */
DECL_VOLTET int API_CreateVolumeMeshObj(
	int ptNum, double pts[],
	int bndFctNum, int bndFcts[], 
	int tetNum, int tets[],
	TigerMesh_t *volMeshObj
);





/**
 *  @brief		在面网格中添加额外的域信息
 *  @param[in]		surfMeshObj		面网格对象
 *  @param[in]		domainTag		面网格单元对应的子域，dim = bndFctNum
 *  @param[in]		originSurfNums	面网格单元对应的原曲面数量，dim = bndFctNum
 *  @param[in]		originSurfs		面网格单元对应的原曲面，dim = originSurfNums[0] + originSurfNums[1] + ...
 *  @param[in]		blockNum		体数量
 *  @param[in]		blockSurfNums	体对应的曲面数量，dim = blockNum
 *  @param[in]		blockSurfs		体对应的曲面，dim = blockSurfNums[0] + blockSurfNums[1] + ...
 *  @param[in]		pecNum			pec体数量
 *  @param[in]		pecBlocks		所有pec体, dim = pecNum
 *  @param[in]		sizeNum			尺寸标识数量
 *  @param[in]		sizeBlocks		需要设置尺寸的体，dim = sizeNum
 *  @param[in]		sizes			尺寸值，dim = sizeNum
 *  @return
 *	0				Completely succeed
 *	otherwise		Fail. Error info. will be defined later
 */
DECL_VOLTET int API_AddSurfaceOriginInformation(
	TigerMesh_t surfMeshObj,
	int *domainTag,
	int *originSurfNums, int *originSurfs,
	int blockNum, int *blockSurfNums, int *blockSurfs,
	int pecNum, int *pecBlocks,
	int sizeNum, int *sizeBlocks, double* sizes
);

/**
 *  @brief		在添加额外的尺寸场函数
 *  @param[in]		surfMeshObj		面网格对象
 *  @param[in]		sizingFunc		尺寸场函数
 *  @return
 *	0				Completely succeed
 *	otherwise		Fail. Error info. will be defined later
 */
DECL_VOLTET int API_AddSizingFunction(
	TigerMesh_t surfMeshObj,
	double(*sizingFunc)(double x, double y, double z)
);

/**
 *  @brief		四面体化面网格
 *  @param[in]		surfMeshObj		输入面网格
 *  @param[out]		volMeshObj		输出体网格
 *  @return
 *	0				Completely succeed
 *	otherwise		Fail. Error info. will be defined later
 */
//#TODO1 volMeshObj和surfMeshObj分开，把surfout收缩后放回到volMeshObj的surtri里
DECL_VOLTET int API_Tetrahedralize(
	TigerMesh_t surfMeshObj,
	TigerMesh_t* volMeshObj
);

/**
 *  @brief		获取四面体网格点数量
 *  @param[in]		objVolmTet		体网格对象
 *  @return
 *	>=0	网格点数量
 *	<0	出错
 */
DECL_VOLTET int API_GetTetrahedraPntNum(TigerMesh_t objVolmTet);

/**
 *  @brief		获取四面体网格点数量（包含高阶中点）
 *  @param[in]		objVolmTet		体网格对象
 *  @return
 *	>=0	网格点数量
 *	<0	出错
 */
DECL_VOLTET int API_GetTetrahedraPntNum_High(TigerMesh_t objVolmTet);

/**
 *  @brief		获取四面体网格单元数量
 *  @param[in]		objVolmTet		体网格对象
 *  @return		
 *	>=0	网格单元数量
 *	<0	出错
 */
DECL_VOLTET int API_GetTetrahedraElemNum(TigerMesh_t objVolmTet);

/**
 *  @brief		获取四面体网格点坐标
 *  @param[in]		objVolmTet		体网格对象
 *  @param[in]		index			点索引
 *  @param[out]		x				x坐标值
 *  @param[out]		y				y坐标值
 *  @param[out]		z				z坐标值
 *  @return
 *	0				Completely succeed
 *	otherwise		Fail. Error info. will be defined later
 */
DECL_VOLTET int API_GetTetrahedraPointCoord(TigerMesh_t objVolmTet, int index, double &x, double &y, double &z);

/**
 *  @brief		获取四面体单元的点索引
 *  @param[in]		objVolmTet		体网格对象
 *  @param[in]		index			单元索引
 *  @param[out]		a				点1索引
 *  @param[out]		b				点2索引
 *  @param[out]		c				点3索引
 *  @param[out]		d				点4索引
 *  @return
 *	0				Completely succeed
 *	otherwise		Fail. Error info. will be defined later
 */
DECL_VOLTET int API_GetTetrahedraElemPntIdx(TigerMesh_t objVolmTet, int index, int &a, int &b, int &c, int &d);

/**
 *  @brief		获取四面体单元的点索引（包含高阶中点）
 *  @param[in]		objVolmTet		体网格对象
 *  @param[in]		index			单元索引
 *  @param[out]		pts				点索引，所有高阶中点在常规点之后
 *  @return
 *	0				Completely succeed
 *	otherwise		Fail. Error info. will be defined later
 */
DECL_VOLTET int API_GetTetrahedraElemPntIdx_High(TigerMesh_t objVolmTet, int index, int pts[10]);

/**
 *  @brief		获取四面体子域数量
 *  @param[in]		objVolmTet		体网格对象
 *  @return		
 *	>=0	子域数量
 *	<0	出错
 */
DECL_VOLTET int API_GetTetrahedraSubdomainNum(TigerMesh_t objVolmTet);

/**
 *  @brief		获取四面体子域所属体数量
 *  @param[in]		objVolmTet		体网格对象
 *  @param[in]		index			子域索引
 *  @return		
 *	>=0	体数量
 *	<0	出错
 */
DECL_VOLTET int API_GetTetrahedraSubdomainBlockNum(TigerMesh_t objVolmTet, int index);

/**
 *  @brief		获取四面体子域的单个所属体编号
 *  @param[in]		objVolmTet		体网格对象
 *  @param[in]		index			子域索引
 *  @param[in]		blkIndex		体索引
 *  @return		
 *	>=0	原始几何体编号
 *	<0	出错
 */
DECL_VOLTET int API_GetTetrahedraSubdomainBlock(TigerMesh_t objVolmTet, int index, int blkIndex);

/**
 *  @brief		获取四面体单元对应子域
 *  @param[in]		objVolmTet		体网格对象
 *  @param[in]		index			单元索引
 *  @return		
 *	>=0	子域索引
 *	<0	出错
 */
DECL_VOLTET int API_GetTetrahedraElemSubdomain(TigerMesh_t objVolmTet, int index);

/**
 *  @brief		获取四面体网格的边界三角形单元数量
 *  @param[in]		objVolmTet		体网格对象
 *  @return
 *	>=0	边界三角形单元数量
 *	<0	出错
 */
DECL_VOLTET int API_GetTetrahedraBndElemNum(TigerMesh_t objVolmTet);

/**
 *  @brief		获取四面体网格的边界三角形单元点索引
 *  @param[in]		objVolmTet		体网格对象
 *  @param[in]		index			边界单元索引
 *  @param[in]		p1				点1索引
 *  @param[in]		p2				点2索引
 *  @param[in]		p3				点3索引
 *  @return
 *	0				Completely succeed
 *	otherwise		Fail. Error info. will be defined later
 */
DECL_VOLTET int API_GetTetrahedraBndElemPts(TigerMesh_t objVolmTet,int index, int& p1, int& p2, int& p3);

/**
 *  @brief		获取四面体网格的边界三角形单元对应子域编号
 *  @param[in]		objVolmTet		体网格对象
 *  @param[in]		index			边界单元索引
 *  @return
 *	>=0	子域编号
 *	<0	出错
 */
DECL_VOLTET int API_GetTetrahedraBndElemSubdomain(TigerMesh_t objVolmTet, int index);

/**
 *  @brief		获取四面体网格的边界三角形单元对应原始曲面数量
 *  @param[in]		objVolmTet		体网格对象
 *  @param[in]		index			边界单元索引
 *  @return
 *	>=0	单元对应原始曲面数量
 *	<0	出错
 */
DECL_VOLTET int API_GetTetrahedraBndElemSurfaceNum(TigerMesh_t objVolmTet, int index);

/**
 *  @brief		获取四面体网格的边界三角形单元的单个原始曲面的编号
 *  @param[in]		objVolmTet		体网格对象
 *  @param[in]		index			边界单元索引
 *  @param[in]		surfIndex		原始曲面索引
 *  @return
 *	>=0	单元对应指定原始曲面的编号
 *	<0	出错
 */
DECL_VOLTET int API_GetTetrahedraBndElemSurface(TigerMesh_t objVolmTet, int index, int surfIndex);

// /**
//  *  @brief		往四面体网格中插点
//  *  @param[in]		volMeshObj	体网格对象
//  *  @param[out]		coord		待插入点坐标
//  *  @return
//  *	0				Completely succeed
//  *	otherwise		Fail. Error info. will be defined later
//  */
// DECL_VOLTET int API_InsertPoint(
// 	TigerMesh_t volMeshObj,
// 	double coord[3]
// );

/**
 * @brief 在体网格中插入用户给定点
 * 
 * @param[in] volMeshObjIn 输入的体网格（用API_CreateVolumeMeshObj定义）
 * @param[out] volMeshObjOut 加点之后的体网格
 * @param[in] addInsidePtNum 插入内部点的数量
 * @param[in] addInsideCoords 插入内部点坐标. dim = 3*addInsidePtNum
 * @param[in] addBndPtNum 插入边界点的数量
 * @param[in] addBndCoords 插入边界点坐标 dim = 3*addBndPtNum
 * @param[in] addBndFts 插入边界点所在面片索引(对三角面片内部分裂，一个面片一个点，索引从0开始)
 * @return
 *  0				Completely succeed
 *	otherwise		Fail. Error info. will be defined later 
 */
DECL_VOLTET int API_InsertPoints(
	TigerMesh_t volMeshObjIn,
	TigerMesh_t *volMeshObjOut,
	int addInsidePtNum,
	double addInsideCoords[],
	int addBndPtNum,
	double addBndCoords[],
	int addBndFts[]
);


/**
 * @brief 体网格自适应
 * 
 * @param[in] volMeshObjIn 输入的体网格（包含边界面，用API_CreateVolumeMeshObj定义）
 * @param[out] volMeshObjOut 自适应之后的体网格
 * @param[in] sizingValues 输入体网格中每个节点上的尺寸值. dim = nPnt
 * @param[in] requiredTetNum 需要保留的四面体网格的数目
 * @param[in] requiredTets 需要保留的四面体网格 dim = requiredTetNum
 * @param[in] requiredFctNum 需要保留的边界面网格的数目
 * @param[in] requiredFcts 需要保留的边界面网格 dim = requiredFctNum
 * @return  
 *  0				Completely succeed
 *	otherwise		Fail. Error info. will be defined later 
 */

DECL_VOLTET int API_Adaptation(
	TigerMesh_t volMeshObjIn,
	TigerMesh_t *volMeshObjOut,
	double sizingValues[],
	int requiredTetNum,
	int requiredTets[],
	int requiredFctNum,
	int requiredFcts[]
);



/**
 *  @brief		释放四面体网格对象
 *  @param[in]		objVolmTet		体网格对象
 *  @return
 *	0				Completely succeed
 *	otherwise		Fail. Error info. will be defined later
 */
DECL_VOLTET int API_DelTetrahedraObj(TigerMesh_t objVolmTet);

/**
 *  @brief		释放面网格对象
 *  @param[in]		objSurfTet		面网格对象
 *  @return
 *	0				Completely succeed
 *	otherwise		Fail. Error info. will be defined later
 */
DECL_VOLTET int API_DelSurfMeshObj(TigerMesh_t objSurfTet);
#endif