//yhf 2021-03
#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <array>
#include <iomanip>
#include "vmesh.h"
#include "gtest/gtest.h"

TEST(API_Gen_Tetra_Mesh, anistropic) {

	std::function<double(std::array<double, 3>)> pfun = nullptr; /*尺寸函数指针，用户也可以自定义函数，如上面的test函数 */


	/* ------------------------- 输入参数 -------------------------- */
	double		*pdSNC = nullptr;			/* 曲面网格节点坐标，coord[3*i], coord[3*i+1], coord[3*i+2]为第i个节点x,y,z坐标 */
	int			nSN = 0;			/* 曲面网格边界节点数目 */
	int			*pnSFFm = nullptr;		/* 曲面网格单元节点编号 */
	int			*pnSFTp = nullptr;;		/* 曲面网格单元类型。当前仅支持三角形单元 目前暂时用不到，可以给空指针 */
	int			*pnSFPt = nullptr;;		/* 曲面网格单元所在几何面编号 注意，从1开始，不是0 */
	int			nSF = 0;			/* 曲面网格单元数目 */
	int			*pnFT = nullptr;;			/* 几何面类型： 0为远场； 1 为物面； 2为对称面 pnFT[0]=1代表1号面为物面 */
	int			nLN = 200;            /* 边界层层数 */
	double		dLen = 0.01;			/* 边界层第一层厚度 */
	double		dRto = 1.2;			/* 边界层厚度增长因子 */
	bool		bisostop = true;      /* 各向同性停止 */
	int			nopt=0; /* 优化次数*/
	/* ------------------------- 输出参数 ------------------------- */
	double	   *ppdMNC = nullptr;		/* 体网格节点坐标 */
	int         pnMN = 0;	/* 体网格节点数目 */
	int         *ppnMEFm = nullptr;		/* 体网格单元节点编号 */
	int         *ppnMETp = nullptr;		/* 体网格单元类型。当前支持单元类型：三棱柱、金字塔、四面体 参照typedef */
	int        pnME = 0;			/* 体网格单元数目 */
		/* ------------------------- 边界层网格顶面层面网格参数 ------------------------- */
	double		*ppdSNC0 = nullptr;		/* 顶面网格节点坐标，coord[3*i], coord[3*i+1], coord[3*i+2]为第i个节点x,y,z坐标 */
	int			pnSN0;			/* 顶面网格边界节点数目 */
	int			pnSEO;			/* 顶面网格单元数目 */
	int			*ppnSFTpO = nullptr;		/* 顶面网格单元类型。当前仅支持三角形单元 */
	int			*ppnSFFmO = nullptr;	/* 顶面网格单元节点编号 */
	bool b_use_multiple_normals = false; /*是否启用多法向 */
	bool b_have_pyramid = false; /*是否有金字塔*/
	bool b_output_io_file = true; /*是否将api的输入和输出都输出到文件系统中（仅用于DEBUG）*/
	/*read pls file */
	/*read pls file */

	//std::string filename = "anistropic/anistropic.pls";
	std::string filename = "ans.o.pls";
	std::ifstream fin(filename);
	std::string trash;
	fin >> nSF >> nSN >> trash >> trash >> trash >> trash;
	ASSERT_GT(nSF, 1) << "No such file!" << std::endl << "[filename] " << filename;
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
	pnFT = new int[maxF];//其余均为物面
	for (int i = 0; i < maxF; i++)
		pnFT[i] = 1;
	for (int i = 0; i < 1; i++) {
		pnFT[i] = 1;
	}

	try {
		API_Gen_Boundary_Mesh(
			/* ------------------------- 输入参数 -------------------------- */
			pdSNC,			/* 曲面网格节点坐标，coord[3*i], coord[3*i+1], coord[3*i+2]为第i个节点x,y,z坐标 */
			nSN,			/* 曲面网格边界节点数目 */
			pnSFFm,		/* 曲面网格单元节点编号 */
			pnSFTp,		/* 曲面网格单元类型。当前仅支持三角形单元 */
			pnSFPt,		/* 曲面网格单元所在几何面编号 */
			nSF,			/* 曲面网格单元数目 */
			pnFT,			/* 几何面类型： 0为远场； 1 为物面； 2为对称面 */
			nLN,			/* 边界层层数 */
			dLen,			/* 边界层第一层厚度 */
			dRto,			/* 边界层厚度增长因子 */
			bisostop,      /* 各向同性停止 */
			/* ------------------------- 输出参数 ------------------------- */
			&ppdMNC,		/* 体网格节点坐标 */
			&pnMN,			/* 体网格节点数目 */
			&ppnMEFm,		/* 体网格单元节点编号 */
			&ppnMETp,		/* 体网格单元类型。当前支持单元类型：三棱柱、金字塔、四面体 */
			&pnME,			/* 体网格单元数目 */
			/* ------------------------- 边界层网格顶面层面网格参数 ------------------------- */
			&ppdSNC0,		/* 顶面网格节点坐标，coord[3*i], coord[3*i+1], coord[3*i+2]为第i个节点x,y,z坐标 */
			&pnSN0,			/* 顶面网格边界节点数目 */
			&pnSEO,			/* 顶面网格单元数目 */
			&ppnSFTpO,		/* 顶面网格单元类型。当前仅支持三角形单元 */
			&ppnSFFmO,		/* 顶面网格单元节点编号 */
			b_have_pyramid,
			b_use_multiple_normals, 
			b_output_io_file 
		);
	}
	catch (std::exception e) {
		std::cout << e.what();
	}
	//delete[]pnSFPt;
	//pnSFPt = new int[pnSEO];
	//double expan_beta = 1.2;
	//for (int i = 0; i < pnSEO; i++)
	//	pnSFPt[i] = 1;
	//for (int i = 0; i < 3 * pnSEO; i++)
	//	ppnSFFmO[i]++;
	//try {
	//	API_Gen_Tetra_Mesh(
	//		/* ------------------------- 输入参数 -------------------------- */
	//		ppdSNC0,			/* 曲面网格节点坐标，coord[3*i], coord[3*i+1], coord[3*i+2]为第i个节点x,y,z坐标 */
	//		pnSN0,			/* 曲面网格边界节点数目 */
	//		ppnSFFmO,		/* 曲面网格单元节点编号 */
	//		ppnSFTpO,		/* 曲面网格单元类型。当前仅支持三角形单元 */
	//		pnSFPt,		/* 曲面网格单元所在几何面编号,从1开始 */
	//		pnSEO,			/* 曲面网格单元数目 */
	//		nopt,
	//		nullptr, /*尺寸函数，输入为点坐标，输出为尺寸值 */
	//		/* ------------------------- 输出参数 ------------------------- */
	//		&ppdMNC,		/* 体网格节点坐标 */
	//		&pnMN,			/* 体网格节点数目 */
	//		&ppnMEFm,		/* 体网格单元节点编号 */
	//		&ppnMETp,		/* 体网格单元类型。当前支持单元类型：三棱柱、金字塔、四面体 */
	//		&pnME,			/* 体网格单元数目 */
	//		false,  /*是否将api的输入和输出都输出到文件系统中（仅用于DEBUG） */
	//		"virtualmesh",   /*几何文件名，缺省为virtualmesh */

	//		/* ------------------------- 控制参数 ------------------------- */
	//		expan_beta       /* 尺寸过渡因子 */
	//	);
	//}
	//catch (std::exception e) {
	//	std::cout << e.what();
	//}	
	//ASSERT_GT(pnME, nSF);
	//ASSERT_GT(pnMN, nSN);
	//ASSERT_GE(pnSN0, nSN);
	//ASSERT_GE(pnSEO, nSF);
	//ASSERT_NE(ppnMEFm, nullptr);
	//ASSERT_NE(ppnMETp, nullptr);
	//ASSERT_NE(ppdMNC, nullptr);
	//ASSERT_NE(ppdSNC0, nullptr);
	//ASSERT_NE(ppnSFTpO, nullptr);
	//ASSERT_NE(ppnSFFmO, nullptr);



}