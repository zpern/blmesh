//yhf 2021-03
#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <array>
#include <iomanip>
#include <string>
#include "vmesh.h"
#include "./writeblmesh.h"
#include "gtest/gtest.h"
#include "common_assert.h"

TEST(API_Gen_Vol_Mesh, intersection) {

	int obj_elli, obj_cube, obj_other;/* 源标识ID */
	std::function<double(std::array<double, 3>)> pfun = nullptr; /* 尺寸函数指针，用户也可以自定义函数，如上面的test函数 */


	/* ------------------------- 输入参数 -------------------------- */
	double		*pdSNC = nullptr;			/* 曲面网格节点坐标，coord[3*i], coord[3*i+1], coord[3*i+2]为第i个节点x,y,z坐标 */
	int			nSN = 0;			/* 曲面网格边界节点数目 */
	int			*pnSFFm = nullptr;		/* 曲面网格单元节点编号 */
	int			*pnSFTp = nullptr;;		/* 曲面网格单元类型。当前仅支持三角形单元 目前暂时用不到，可以给空指针 */
	int			*pnSFPt = nullptr;;		/* 曲面网格单元所在几何面编号 注意，从1开始，不是0 */
	int			nSF = 0;			/* 曲面网格单元数目 */
	int			*pnFT = nullptr;;			/* 几何面类型： 0为远场； 1 为物面； 2为对称面 pnFT[0]=1代表1号面为物面 */
	int			nLN = 100;            /* 边界层层数 */
	double		dLen = 0.01;			/* 边界层第一层厚度 */
	double		dRto = 1.25;			/* 边界层厚度增长因子 */
	bool		bisostop = true;      /* 各向同性停止 */
	int nopt = 10;/* 优化次数*/
	/* ------------------------- 输出参数 ------------------------- */
	double	   *ppdMNC = nullptr;		/* 体网格节点坐标 */
	int         pnMN = 0;	/* 体网格节点数目 */
	int         *ppnMEFm = nullptr;		/* 体网格单元节点编号 */
	int         *ppnMETp = nullptr;		/* 体网格单元类型。当前支持单元类型：三棱柱、金字塔、四面体 参照typedef */
	int        pnME = 0;			/* 体网格单元数目 */
			/* ------------------------- 边界信息 -------------------------*/
	int         num_boundary_face; /*边界面网格数量 */
	int         *boundary_mesh; /*边界面网格，注意每四个为一组，而且注意如果为三角形，最后一项为-1，id从0开始 */
	int         *boundary_face; /*边界面网格对应的面id，长度为num_boundary_face */
	/* ------------------------- 其他参数 -------------------------*/
	bool b_use_multiple_normals = false; /*是否启用多法向*/
	bool b_output_io_file = false; /*是否将api的输入和输出都输出到文件系统中（仅用于DEBUG）*/

		/* ------------------------- 控制参数 -------------------------*/
	double expan_beta = 1.2;      /* 尺寸过渡因子*/

	bool buse_multiple_normals = false;

	bool b_output_IO = true;
	int boundary_info;
	/* read pls file */
	/* read pls file */
	std::string filename = "VolumeMeshing.txt";
	//std::string filename = "intersection_vol/VolumeMeshing.txt";
	std::ifstream fin(filename);
	std::string trash;
	fin >>trash >> expan_beta >> trash >> nLN >> trash >> dLen >> trash >> dRto >> trash >> bisostop >> trash >> boundary_info;
	pnFT = new int[boundary_info];
	for (int i = 0; i < boundary_info; i++) {
		fin >> pnFT[i];
	}
	fin >> trash >> buse_multiple_normals;
	fin >> nSF >> nSN;
	ASSERT_GT(nSF, 1) << "No such file!" << std::endl << "[filename] " << filename;
	ASSERT_GT(nSN, 3);
	pdSNC = new double[nSN * 3];
	for (int i = 0; i < nSN; i++) {
		fin >> trash >> pdSNC[3 * i] >> pdSNC[3 * i + 1] >> pdSNC[3 * i + 2];
	}
	pnSFFm = new int[3 * nSF];
	pnSFPt = new int[nSF];
	for (int i = 0; i < nSF; i++) {
		fin >> trash >> pnSFFm[3 * i + 1] >> pnSFFm[3 * i + 0] >> pnSFFm[3 * i + 2] >> pnSFPt[i];
	}
	ASSERT_NO_THROW(
		TiGER::API_Gen_Vol_Mesh(
			/*  ------------------------- 输入参数 -------------------------- */
			pdSNC,			/*  曲面网格节点坐标，coord[3*i], coord[3*i+1], coord[3*i+2]为第i个节点x,y,z坐标  */
			nSN,			/*  曲面网格边界节点数目  */
			pnSFFm,		/*  曲面网格单元节点编号  */
			pnSFTp,		/*  曲面网格单元类型。当前仅支持三角形单元  */
			pnSFPt,		/*  曲面网格单元所在几何面编号  */
			nSF,			/*  曲面网格单元数目  */
			pnFT,           /* 面信息 */
			nLN,
			dLen,
			dRto, bisostop, nopt,
			pfun,
			/*  ------------------------- 输出参数 ------------------------- */
			&ppdMNC,		/*  体网格节点坐标  */
			&pnMN,			/*  体网格节点数目  */
			&ppnMEFm,		/*  体网格单元节点编号  */
			&ppnMETp,		/*  体网格单元类型。当前支持单元类型：三棱柱、金字塔、四面体  */
			&pnME,			/*  体网格单元数目  */
			&num_boundary_face,
			&boundary_mesh,
			&boundary_face,
			false,
			b_output_IO,
			filename,
			expan_beta
		);
	);

	writeBLVTK(filename+"_out.vtk",ppdMNC, pnMN, ppnMEFm, ppnMETp, pnME);
	writeBoundary(filename + "_out.bdry", num_boundary_face, boundary_mesh, boundary_face);
	


	DELETE_ARRAY(pdSNC);
	DELETE_ARRAY(pnSFFm);
	DELETE_ARRAY(pnSFTp);
	DELETE_ARRAY(pnSFPt);
	DELETE_ARRAY(ppnMEFm);
	DELETE_ARRAY(ppnMETp);
	DELETE_ARRAY(ppdMNC);

	//while (true) {}
}
