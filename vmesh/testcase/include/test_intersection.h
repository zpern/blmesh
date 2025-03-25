//yhf 2021-03
#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <array>
#include <iomanip>
#include "vmesh.h"
#include "gtest/gtest.h"
#include "tiger_sizingfunction.h"
#include "tiger_sizingsmoother.h"
using namespace std;
TEST(API_Gen_Boundary_Mesh, intersection) {

	std::function<double(std::array<double, 3>)> pfun = nullptr; /*尺寸函数指针，用户也可以自定义函数，如上面的test函数 */

   /* ------------------------- 输入参数 -------------------------- */
	double		*pdSNC = nullptr;			/* 曲面网格节点坐标，coord[3*i], coord[3*i+1], coord[3*i+2]为第i个节点x,y,z坐标  */
	int			nSN = 0;			/* 曲面网格边界节点数目  */
	int			*pnSFFm = nullptr;		/* 曲面网格单元节点编号  */
	int			*pnSFTp = nullptr;;		/* 曲面网格单元类型。当前仅支持三角形单元 目前暂时用不到，可以给空指针  */
	int			*pnSFPt = nullptr;;		/* 曲面网格单元所在几何面编号 注意，从1开始，不是0  */
	int			nSF = 0;			/* 曲面网格单元数目  */
	std::map<int,int>			pnFT;			/* 几何面类型： 0为远场； 1 为物面； 2为对称面 pnFT[0]=1代表1号面为物面  */
	int			nLN = 20;            /* 边界层层数  */
	double		dLen = 0.01;			/* 边界层第一层厚度  */
	double		dRto = 1.3;			/* 边界层厚度增长因子  */
	double		max_skewnwass,  /* 各向异性停止**/
	bool		bisostop = true;      /* 各向同性停止 */
	int nopt = 0; /* 优化次数 */
	/* ------------------------- 输出参数 ------------------------- */
	double	   *ppdMNC1 = nullptr;		/* 体网格节点坐标  */
	int         pnMN1 = 0;	/* 体网格节点数目  */
	int         *ppnMEFm1 = nullptr;		/* 体网格单元节点编号  */
	int         *ppnMETp1 = nullptr;		/* 体网格单元类型。当前支持单元类型：三棱柱、金字塔、四面体 参照typedef  */
	int        pnME1 = 0;			/* 体网格单元数目  */

	double	   *ppdMNC2 = nullptr;		/* 体网格节点坐标  */
	int         pnMN2 = 0;	/* 体网格节点数目  */
	int         *ppnMEFm2 = nullptr;		/* 体网格单元节点编号  */
	int         *ppnMETp2 = nullptr;		/* 体网格单元类型。当前支持单元类型：三棱柱、金字塔、四面体 参照typedef  */
	int        pnME2 = 0;			/* 体网格单元数目  */

	double	   *ppdMNC = nullptr;		/* 体网格节点坐标  */
	int         pnMN = 0;	/* 体网格节点数目  */
	int         *ppnMEFm = nullptr;		/* 体网格单元节点编号  */
	int         *ppnMETp = nullptr;		/* 体网格单元类型。当前支持单元类型：三棱柱、金字塔、四面体 参照typedef  */
	int        pnME = 0;			/* 体网格单元数目  */
		/* ------------------------- 边界层网格顶面层面网格参数 ------------------------- */
	double		*ppdSNC0 = nullptr;		/* 顶面网格节点坐标，coord[3*i], coord[3*i+1], coord[3*i+2]为第i个节点x,y,z坐标  */
	int			pnSN0;			/* 顶面网格边界节点数目  */
	int			pnSEO;			/* 顶面网格单元数目  */
	int			*ppnSFTpO = nullptr;		/* 顶面网格单元类型。当前仅支持三角形单元  */
	int			*ppnSFFmO = nullptr;	/* 顶面网格单元节点编号  */
	bool buse_multiple_normals = false; /*是否启用多法向 */
	bool b_output_IO = false;

	/*Boundary Meshing*/

	std::string filename = "BoundaryMeshing.txt";
	//
	std::ifstream fin(filename);
	std::string trash;
	int boundary_info;
	fin >> trash >> nLN >> trash >> dLen >> trash >> dRto >> trash >> bisostop >> trash >> boundary_info;
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
	bool b_have_pyramid = true;

	pnSFFm = new int[3 * nSF];
	pnSFPt = new int[nSF];
	int maxF = 0;
	for (int i = 0; i < nSF; i++) {
		fin >> trash >> pnSFFm[3 * i + 1] >> pnSFFm[3 * i + 0] >> pnSFFm[3 * i + 2] >> pnSFPt[i];
		maxF = std::max(maxF, pnSFPt[i]);
	}
	int* l2g;             /*  顶面网格id到全局点id的映射  */
	double* ppnsizing; /* 顶面网格目标尺寸 */

	int num_boundary_face;  /////* 边界面网格数量 **/ 
	int* boundary_mesh;  /////* 边界面网格，注意每四个为一组，而且注意如果为三角形，最后一项为-1，id从0开始 **/ 
	int* boundary_face;  /////* 边界面网格对应的面id，长度为num_boundary_face **/ 


	//std::array<double, 12> per_face = {1,0,0,0,984808,0.173648,0,-0.173648,0.984808,0,0,0};
	std::array<double, 12> per_face = { 1,0,0,0,0.959493,0.281733,0,-0.281733,0.959493,0,0,0 };
	try {
		TiGER::API_Gen_Boundary_Mesh(
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
			&ppdMNC1,		/* 体网格节点坐标 */
			&pnMN1,			/* 体网格节点数目 */
			&ppnMEFm1,		/* 体网格单元节点编号 */
			&ppnMETp1,		/* 体网格单元类型。当前支持单元类型：三棱柱、金字塔、四面体 */
			&pnME1,			/* 体网格单元数目 */
			/* ------------------------- 边界层网格顶面层面网格参数 ------------------------- */
			&ppdSNC0,		/* 顶面网格节点坐标，coord[3*i], coord[3*i+1], coord[3*i+2]为第i个节点x,y,z坐标 */
			&pnSN0,			/* 顶面网格边界节点数目 */
			&pnSEO,			/* 顶面网格单元数目 */
			&ppnSFTpO,		/* 顶面网格单元类型。当前仅支持三角形单元 */
			&ppnSFFmO,		/* 顶面网格单元节点编号 */
			&l2g,			/*  顶面网格id到全局点id的映射  */
			&ppnsizing,		/* 顶面网格目标尺寸 */
			/*-------------------------边界信息 ---------------------*/
			&num_boundary_face,  /////* 边界面网格数量 **/ 
			&boundary_mesh, /////* 边界面网格，注意每四个为一组，而且注意如果为三角形，最后一项为-1，id从0开始 **/ 
			&boundary_face,  /////* 边界面网格对应的面id，长度为num_boundary_face **/ 
			b_have_pyramid, /*是否有金字塔*/
			buse_multiple_normals, /*是否启用多法向 */
			b_output_IO,
			"virtualmesh",
			per_face
		);
	}
	catch (std::exception e) {
		std::cout << e.what();
	}
	std::ofstream fout("top_layer.vtk");
	fout << "# vtk DataFile Version 2.0" << endl;
	fout << "boundary layer mesh" << endl;
	fout << "ASCII" << endl;
	fout << "DATASET UNSTRUCTURED_GRID" << endl;
	fout << "POINTS " << pnSN0 << " float" << endl;
	for (int i = 0; i < pnSN0; i++) {
		fout << (ppdSNC0)[3 * i] << " " << (ppdSNC0)[3 * i + 1] << " " << (ppdSNC0)[3 * i + 2] << endl;
	}

	fout << "CELLS " << pnSEO << " " << 4* pnSEO << endl;
	//cout << "pnSNO=" << pnSNO << endl;


	for (int i = 0; i < pnSEO; i++) {
		int type = ppnMETp1[i];
		fout << 3;
		for (int k = 0; k < 3; k++)
			fout << " " << (ppnSFFmO)[3*i+k];
		fout << endl;
	}

	fout << "CELL_TYPES " << pnSEO << endl;
	for (int i = 0; i < pnSEO; i++) {
		fout << 5 << endl;
	}



	fout.close();
	fout= std::ofstream("prismes.vtk");
	fout << "# vtk DataFile Version 2.0" << endl;
	fout << "boundary layer mesh" << endl;
	fout << "ASCII" << endl;
	fout << "DATASET UNSTRUCTURED_GRID" << endl;
	fout << "POINTS " << pnMN1 << " float" << endl;
	for (int i = 0; i < pnMN1; i++) {
		fout << (ppdMNC1)[3 * i] << " " << (ppdMNC1)[3 * i + 1] << " " << (ppdMNC1)[3 * i + 2] << endl;
	}


	int count2 = 0;
	for (int i = 0; i < pnME1; i++) {
		int type = ppnMETp1[i];
		switch (type)
		{
		case 11://tetra
			count2 += 4;
			break;
		case 13://prism
			count2 += 6;
			break;
		case 14://pyramid
			count2 += 5;
			break;
		default:
			break;
		}
	}


	fout << "CELLS " << pnME1 << " " << count2 + pnME1 << endl;
	//cout << "pnSNO=" << pnSNO << endl;


	int count = 0;
	for (int i = 0; i < pnME1; i++) {
		int type = ppnMETp1[i];
		switch (type)
		{
		case 11://tetra
			fout << 4;
			for (int k = 0; k < 4; k++)
				fout << " " << (ppnMEFm1)[count++];
			fout << endl;
			break;
		case 13://prism
			fout << 6;
			for (int k = 0; k < 6; k++)
				fout << " " << (ppnMEFm1)[count++];
			fout << endl;
			break;
		case 14://pyramid
			fout << 5;
			for (int k = 0; k < 5; k++)
				fout << " " << (ppnMEFm1)[count++];
			fout << endl;
			break;
		default:
			break;
		}
	}

	fout << "CELL_TYPES " << pnME1 << endl;
	for (int i = 0; i < pnME1; i++) {
		int type = ppnMETp1[i];
		switch (type)
		{
		case 11://tetra
			fout << 10 << endl;
			break;
		case 13://prism
			fout << 13 << endl;
			break;
		case 14://pyramid
			fout << 14 << endl;
			break;
		default:
			break;
		}
	}

	fout.close();










	return ;
	delete[]pnSFPt;
	pnSFPt = new int[pnSEO];
	double expan_beta = 1.2;
	for (int i = 0; i < pnSEO; i++)
		pnSFPt[i] = 1;
	for (int i = 0; i < 3 * pnSEO; i++)
		ppnSFFmO[i]++;
	try {
		TiGER::API_Gen_Tetra_Mesh(
			/* ------------------------- 输入参数 -------------------------- */
			ppdSNC0,			/* 曲面网格节点坐标，coord[3*i], coord[3*i+1], coord[3*i+2]为第i个节点x,y,z坐标 */
			pnSN0,			/* 曲面网格边界节点数目 */
			ppnSFFmO,		/* 曲面网格单元节点编号 */
			ppnSFTpO,		/* 曲面网格单元类型。当前仅支持三角形单元 */
			pnSFPt,		/* 曲面网格单元所在几何面编号,从1开始 */
			pnSEO,			/* 曲面网格单元数目 */
			0,

			nullptr, /*尺寸函数，输入为点坐标，输出为尺寸值 */
			/* ------------------------- 输出参数 ------------------------- */
			&ppdMNC2,		/* 体网格节点坐标 */
			&pnMN2,			/* 体网格节点数目 */
			&ppnMEFm2,		/* 体网格单元节点编号 */
			&ppnMETp2,		/* 体网格单元类型。当前支持单元类型：三棱柱、金字塔、四面体 */
			&pnME2,			/* 体网格单元数目 */
			false,  /*是否将api的输入和输出都输出到文件系统中（仅用于DEBUG） */
			"virtualmesh",   /*几何文件名，缺省为virtualmesh */

			/* ------------------------- 控制参数 ------------------------- */
			expan_beta       /* 尺寸过渡因子 */
		);
	}
	catch (std::exception e) {
		std::cout << e.what();
	}
	TiGER::API_Mesh_Merge(
		/* ------------------------- 输入参数 -------------------------- */
		ppdMNC1,		/* 体网格节点坐标 */
		pnMN1,			/* 体网格节点数目 */
		ppnMEFm1,		/* 体网格单元节点编号 */
		ppnMETp1,		/* 体网格单元类型。当前支持单元类型：三棱柱、金字塔、四面体 */
		pnME1,			/* 体网格单元数目 */

		ppdMNC2,		/* 体网格节点坐标 */
		pnMN2,			/* 体网格节点数目 */
		ppnMEFm2,		/* 体网格单元节点编号 */
		ppnMETp2,		/* 体网格单元类型。当前支持单元类型：三棱柱、金字塔、四面体 */
		pnME2,			/* 体网格单元数目 */
	/* ------------------------- 输出参数 ------------------------- */
		&ppdMNC,		/* 体网格节点坐标 */
		&pnMN,			/* 体网格节点数目 */
		&ppnMEFm,		/* 体网格单元节点编号 */
		&ppnMETp,		/* 体网格单元类型。当前支持单元类型：三棱柱、金字塔、四面体 */
		&pnME			/* 体网格单元数目 */
	);

	std::cout<< pnMN1+ pnMN2- pnMN << " points Merged , where the number of point in tetra api input is "<< pnSN0 << std::endl;
	fout= std::ofstream("MergedMesh.vtk");
	//	fout.setf(ios::scientific);

	fout << "# vtk DataFile Version 2.0" << endl;
	fout << "boundary layer mesh" << endl;
	fout << "ASCII" << endl;
	fout << "DATASET UNSTRUCTURED_GRID" << endl;
	fout << "POINTS " << pnMN << " float" << endl;
	for (int i = 0; i < pnMN; i++) {
		fout << (ppdMNC)[3 * i] << " " << (ppdMNC)[3 * i + 1] << " " << (ppdMNC)[3 * i + 2] << endl;
	}


	count2 = 0;
	for (int i = 0; i < pnME; i++) {
		int type = ppnMETp[i];
		switch (type)
		{
		case 11://tetra
			count2 += 4;
			break;
		case 13://prism
			count2 += 6;
			break;
		case 14://pyramid
			count2 += 5;
			break;
		default:
			break;
		}
	}


	fout << "CELLS " << pnME << " " << count2 + pnME << endl;
	//cout << "pnSNO=" << pnSNO << endl;


	count = 0;
	for (int i = 0; i < pnME; i++) {
		int type = ppnMETp[i];
		switch (type)
		{
		case 11://tetra
			fout << 4;
			for (int k = 0; k < 4; k++)
				fout << " " << (ppnMEFm)[count++];
			fout << endl;
			break;
		case 13://prism
			fout << 6;
			for (int k = 0; k < 6; k++)
				fout << " " << (ppnMEFm)[count++];
			fout << endl;
			break;
		case 14://pyramid
			fout << 5;
			for (int k = 0; k < 5; k++)
				fout << " " << (ppnMEFm)[count++];
			fout << endl;
			break;
		default:
			break;
		}
	}

	fout << "CELL_TYPES " << pnME << endl;
	for (int i = 0; i < pnME; i++) {
		int type = ppnMETp[i];
		switch (type)
		{
		case 11://tetra
			fout << 10 << endl;
			break;
		case 13://prism
			fout << 13 << endl;
			break;
		case 14://pyramid
			fout << 14 << endl;
			break;
		default:
			break;
		}
	}

	fout.close();

ASSERT_GT(pnME, nSF);
ASSERT_GT(pnMN, nSN);
ASSERT_GE(pnSN0, nSN);
ASSERT_GE(pnSEO, nSF);
ASSERT_NE(ppnMEFm, nullptr);
ASSERT_NE(ppnMETp, nullptr);
ASSERT_NE(ppdMNC, nullptr);
ASSERT_NE(ppdSNC0, nullptr);
ASSERT_NE(ppnSFTpO, nullptr);
ASSERT_NE(ppnSFFmO, nullptr);



}