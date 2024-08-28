//yhf 2021-03
#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <array>
#include <iomanip>
#include "vmesh.h"
#include "gtest/gtest.h"


TEST(API_Gen_Tetra_Mesh, sphere) {

	int obj_elli, obj_cube, obj_other;/* 源标识ID */
	/* 加密区域局部坐标系坐标轴 其中local_coordz方向为local_coordx和local_coordy的叉乘方向 */
	std::array<double, 3> local_coordx = { 2.0,0.0,0.0 };
	std::array<double, 3> local_coordy = { 0.0,2.0,0 };
	std::array<double, 3> local_coordz = { 0.0,0.0,2.0 };

	/* 测试椭圆的中心坐标（局部坐标系原点） */
	std::array<double, 3> local_coord_origin_Ellipsoid = { 0.4559882388381533,5,5 };
	/* 测试长方体的中心坐标（局部坐标系原点） */
	std::array<double, 3> local_coord_origin_Cuboid = { 4.762906708777568,4.5,2.5 };
	/* 测试要删除的obj中心坐标（局部坐标系原点） */
	std::array<double, 3> local_coord_origin_other = { 2.5, 2.4, 1.0 };
	double siz_Ellipsoid = 0.1;/* 椭圆内部尺寸 */
	double siz_Cuboid = 0.05;/* 长方体内部尺寸 */
	double volmGrowthRatio = 1.1;/* 尺寸拓展率 ，大于1，越接近1过渡越慢，建议1.3 */
	//API_Create_Cuboid_Sources_SF(
	//	local_coordx.data(), local_coordy.data(), local_coordz.data(), local_coord_origin_Ellipsoid.data(), siz_Ellipsoid, volmGrowthRatio,
	//	&obj_elli);/* 生成椭圆源 */
	//API_Create_Cuboid_Sources_SF(
	//	local_coordx.data(), local_coordy.data(), local_coordz.data(), local_coord_origin_Cuboid.data(), siz_Cuboid, volmGrowthRatio,
	//	&obj_cube);/* 生成长方体源 */

	//API_Create_Cuboid_Sources_SF(
	//	local_coordx.data(), local_coordy.data(), local_coordz.data(), local_coord_origin_other.data(), siz_Cuboid, volmGrowthRatio,
	//	&obj_other);/* 生成某一源 */

	//API_Delete_SF_Object(obj_other);/* 删除最后一个源 */

					 /* 测试类型 */
	std::function<double(std::array<double, 3>)> pfun = nullptr; /*  尺寸函数指针，用户也可以自定义函数，如上面的test函数 */


	/*  ------------------------- 输入参数 -------------------------- */
	double		*pdSNC = nullptr;			/*  曲面网格节点坐标，coord[3*i], coord[3*i+1], coord[3*i+2]为第i个节点x,y,z坐标  */
	int			nSN = 0;			/*  曲面网格边界节点数目  */
	int			*pnSFFm = nullptr;		/*  曲面网格单元节点编号  */
	int			*pnSFTp = nullptr;;		/*  曲面网格单元类型。当前仅支持三角形单元 目前暂时用不到，可以给空指针  */
	int			*pnSFPt = nullptr;;		/*  曲面网格单元所在几何面编号 注意，从1开始，不是0  */
	int			nSF = 0;			/*  曲面网格单元数目  */
	bool		bisostop = true;      /*  各向同性停止 */
	int nopt =10;/* 优化次数*/
	/*  ------------------------- 输出参数 ------------------------- */
	double	   *ppdMNC = nullptr;		/*  体网格节点坐标  */
	int         pnMN = 0;	/*  体网格节点数目  */
	int         *ppnMEFm = nullptr;		/*  体网格单元节点编号  */
	int         *ppnMETp = nullptr;		/*  体网格单元类型。当前支持单元类型：三棱柱、金字塔、四面体 参照typedef  */
	int        pnME = 0;			/*  体网格单元数目  */
	bool b_output_IO = false;
	/* read pls file */
	/* read pls file */

	std::string filename = "sphere/sphere.pls";
	std::ifstream fin(filename);
	std::string trash;
	fin >> nSF >> nSN >> trash >> trash >> trash >> trash;

	ASSERT_GT(nSF, 1)<<"No such file!"<<std::endl<<"[filename] "<<filename;
	ASSERT_GT(nSN, 3);

	pdSNC = new double[nSN * 3];
	for (int i = 0; i < nSN; i++) {
		fin >> trash >> pdSNC[3 * i] >> pdSNC[3 * i + 1] >> pdSNC[3 * i + 2];
	}
	pnSFFm = new int[3 * nSF];
	pnSFPt = new int[nSF];
	int maxF = 0;
	for (int i = 0; i < nSF; i++) {
		fin >> trash >> pnSFFm[3 * i + 1] >> pnSFFm[3 * i + 0] >> pnSFFm[3 * i + 2] >> pnSFPt[i];
		maxF = std::max(maxF, pnSFPt[i]);
	}

	ASSERT_NO_THROW(
		API_Gen_Tetra_Mesh(
			/*  ------------------------- 输入参数 -------------------------- */
			pdSNC,			/*  曲面网格节点坐标，coord[3*i], coord[3*i+1], coord[3*i+2]为第i个节点x,y,z坐标  */
			nSN,			/*  曲面网格边界节点数目  */
			pnSFFm,		/*  曲面网格单元节点编号  */
			pnSFTp,		/*  曲面网格单元类型。当前仅支持三角形单元  */
			pnSFPt,		/*  曲面网格单元所在几何面编号  */
			nSF,			/*  曲面网格单元数目  */
			nopt,
			pfun,
			/*  ------------------------- 输出参数 ------------------------- */
			&ppdMNC,		/*  体网格节点坐标  */
			&pnMN,			/*  体网格节点数目  */
			&ppnMEFm,		/*  体网格单元节点编号  */
			&ppnMETp,		/*  体网格单元类型。当前支持单元类型：三棱柱、金字塔、四面体  */
			&pnME,			/*  体网格单元数目  */
			b_output_IO,
			filename,
			volmGrowthRatio
		);
	);
	ASSERT_GT(pnME, nSF);
	ASSERT_GT(pnMN, nSN);
	ASSERT_NE(ppnMEFm, nullptr);
	ASSERT_NE(ppnMETp, nullptr);
	ASSERT_NE(ppdMNC,  nullptr);

	DELETE_ARRAY(pdSNC);
	DELETE_ARRAY(pnSFFm);
	DELETE_ARRAY(pnSFTp);
	DELETE_ARRAY(pnSFPt);
	DELETE_ARRAY(ppnMEFm);
	DELETE_ARRAY(ppnMETp);
	DELETE_ARRAY(ppdMNC);

}