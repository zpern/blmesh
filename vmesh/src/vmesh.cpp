#include "vmesh.h"
#include "Config.h"
#include "blpre.h"
#include "blmeshapi.h"
#include "../common/singleton_terminate.h"
#include "../postprocess/postprocess.h"
#include <spdlog/spdlog.h> 
 #include <vector>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <string>
#include <iostream>
#include <functional>
#include <iomanip>
#include <array>
#include<igl/remove_duplicate_vertices.h>
#ifdef WIN32
#include "windows.h"
#endif
#include<array>


#include <thread>
using namespace std;
namespace TiGER {
	/*type def
	enum EntityTopology
	{
		NODE = 5,
		LINE = 6,
		POLYGON =7,
		TRIANGLE =8,
		QUADRILATERAL =9,
		POLYHEDRON =10,
		TETRAHEDRON =11,
		HEXAHEDRON =12,
		PRISM =13,
		PYRAMID =14,
		SEPTAHEDRON =15,
		MIXED
	};
	**/
	void API_Terminate()
	{ 
		spdlog::info("try to stop the program！**************************************************************");
		SingletonTerminate::GetSingletonPtr()->terminate = true;
	}
	int API_Gen_Tetra_Mesh(

		/* ------------------------- 输入参数 --------------------------**/
		double		*pdSNC,			/* 曲面网格节点坐标，coord[3*i], coord[3*i+1], coord[3*i+2]为第i个节点x,y,z坐标 **/
		int			nSN,			/* 曲面网格边界节点数目 **/
		int			*pnSFFm,		/* 曲面网格单元节点编号 **/
		int			*pnSFTp,		/* 曲面网格单元类型。当前仅支持三角形单元 **/
		int			*pnSFPt,		/* 曲面网格单元所在几何面编号,从1开始 **/
		int			nSF,			/* 曲面网格单元数目 **/
		int			nopt,           /* 优化次数 **/
		std::function<double(std::array<double, 3>)> sizefunction, /* 尺寸函数，输入为点坐标，输出为尺寸值 **/
		/* ------------------------- 输出参数 ------------------------- **/ 
		double	   **ppdMNC,		/* 体网格节点坐标 **/ 
		int         *pnMN,			/* 体网格节点数目 **/ 
		int         **ppnMEFm,		/* 体网格单元节点编号 **/ 
		int         **ppnMETp,		/* 体网格单元类型。当前支持单元类型：四面体 **/ 
		int         *pnME,			/* 体网格单元数目 **/ 
		bool b_output_io_file,      /* 是否将api的输入和输出都输出到文件系统中（仅用于DEBUG）**/ 
		std::string filename,       /*  几何文件名，缺省为virtualmesh **/ 

		/* ------------------------- 控制参数 ------------------------- **/
		double expan_beta       /* 尺寸过渡因子**/

	) {


		VM surface;
		surface.ppdSNC0 = pdSNC;
		surface.nSN0 = nSN;
		surface.pnSNO = nSF;
		surface.ppnSFFmO = new int[nSF * 3]; ;
		for (int i = 0; i < 3 * surface.pnSNO; i++)
			surface.ppnSFFmO[i] = pnSFFm[i] - 1;
		surface.convert2NonManifold();



		tiger::Config cf;
		cf.SetDefaultConfig();
		cf.isotropic_stop = 1;
		double  dLen = 0.01;
		double dRto = 1.30;
		int nLN = 0;
		cf.dStepLen = dLen;
		cf.dStepLenRatio = dRto;
		cf.nLayerNum = nLN;
		cf.sGeoFileName = filename;


		vector<int> symm;
		vector<int> box;
		int numFace = 0;
		for (int i = 0; i < nSF; i++) {
			numFace = max(numFace, pnSFPt[i]);
		}
		cf.vecBoxFc = box;
		cf.vecSymmFc = symm;

		int npt = surface.nSN0;
		int nelm = surface.pnSNO;

		stringstream fout;

		fout << nelm << " " << npt << " " << 0 << " " << 0 << " " << 0 << " " << 0 << endl;
		std::vector<std::array<double, 3>> points;
		for (int i = 0; i < npt; i++) {
			points.push_back(std::array<double, 3>{surface.ppdSNC0[3 * i + 0], surface.ppdSNC0[3 * i + 1], surface.ppdSNC0[3 * i + 2]});
		}
		for (int i = 0; i < nelm; i++) {
			fout << i + 1 << " " << surface.ppnSFFmO[3 * i + 1] + 1 << " " << surface.ppnSFFmO[3 * i + 0] + 1 << " " << surface.ppnSFFmO[3 * i + 2] + 1 << " " << pnSFPt[i] << endl;
		}
		//fout.close();
		string input = fout.str();

		if (b_output_io_file) {
			ofstream ftest("TetraMeshing.txt");
			ftest.setf(ios::fixed, ios::floatfield);  // 设定为 fixed 模式，以小数点表示浮点数
			ftest.precision(17);
			ftest << "expan_beta: " << expan_beta << endl;
			ftest << nelm << " " << npt << " " << endl;
			std::vector<std::array<double, 3>> points;
			for (int i = 0; i < npt; i++) {
				ftest << i + 1 << " " << pdSNC[3 * i + 0] << " " << pdSNC[3 * i + 1] << " " << pdSNC[3 * i + 2] << endl;
			}
			for (int i = 0; i < nelm; i++) {
				ftest << i + 1 << " " << pnSFFm[3 * i + 1] << " " << pnSFFm[3 * i + 0] << " " << pnSFFm[3 * i + 2] << " " << pnSFPt[i] << endl;
			}

			cf.WriteConfigFile();
			ftest.close();
		}

		blpreConfig blconfig;
		blconfig.box = box;
		blconfig.symm = symm;
		blconfig.n = nLN;
		blconfig.Ro = dRto;
		blconfig.len = dLen;
		blconfig.use_multiple_normals = false;
		ControlVolume cv;
		auto bdyfile = PRE::blpre(input, blconfig, points,cv);
		VM v = blmesh(bdyfile, blconfig, true, false, true, false, b_output_io_file, sizefunction, expan_beta, nopt);
		*ppdMNC = v.ppdMNC;
		*pnMN = v.pnMN;
		*pnME = v.pnME;
		*ppnMEFm = v.ppnMEFm;
		*ppnMETp = v.ppnMETp;
		delete[]v.boundary_face;
		v.boundary_face = nullptr;
		delete[]v.boundary_mesh;
		v.boundary_mesh = nullptr;
		delete[]surface.ppnSFFmO;

		if (cv.f.size()) {
			int old_point_num = *pnMN;
			double* nppdMNC = new double[(*pnMN + cv.lower_point_num) * 3];
			for (int i = 0; i < *pnMN * 3; i++) {
				nppdMNC[i] = *ppdMNC[i];
			}
			for (int i = 0; i < cv.lower_point_num; i++) {
				for (int k = 0; k < 3; k++) {
					nppdMNC[*pnMN * 3 + 3 * i + k] = cv.v[i][k];
				}
			}
			*pnMN = *pnMN + cv.lower_point_num;
			//delete* ppdMNC[];
			ppdMNC = &nppdMNC;
			int* nppnMEFm = new int[6 * (*pnME) + 6 * cv.f.size()];
			int cnt = 0;
			for (int i = 0; i < *pnMN; i++) {
				int k = 0;
				if (*ppnMETp[i] == 11) {
					k = 4;
				}
				if (*ppnMETp[i] == 13) {
					k = 6;
				}
				if (*ppnMETp[i] == 14) {
					k = 5;
				}
				for (int j = 0; j < k; j++) {
					nppnMEFm[cnt] = *ppnMEFm[cnt];
					cnt++;
				}
			}

			for (int i = 0; i < cv.f.size(); i++) {
				int k = cv.f[i].size();
				for (int j = 0; j < k; j++) {
					if (cv.f[i][k] > cv.lower_point_num)
						nppnMEFm[cnt] = cv.f[i][k] - cv.lower_point_num;
					else
						nppnMEFm[cnt] = cv.f[i][k] + old_point_num;
					cnt++;
				}
			}
			int* nppnMETp = new int[*pnME + cv.f.size()];
			for (int i = 0; i < *pnME; i++) {
				nppnMETp[i] = *ppnMETp[i];
			}
			for (int i = 0; i < cv.f.size(); i++) {
				int k = cv.f[i].size();
				if (k == 4) { nppnMETp[*pnME + i] = 11; }
				if (k == 5) { nppnMETp[*pnME + i] = 14; }
				if (k == 6) { nppnMETp[*pnME + i] = 13; }
			}

			*pnME = *pnMN + cv.f.size();
			ppnMEFm = &nppnMEFm;
			ppnMETp = &nppnMETp;
		}

		return 0;

	}

	int API_Gen_Vol_Mesh(
		/* ------------------------- 输入参数 --------------------------**/
		double		*pdSNC,			/* 曲面网格节点坐标，coord[3*i], coord[3*i+1], coord[3*i+2]为第i个节点x,y,z坐标 **/
		int			nSN,			/* 曲面网格边界节点数目 **/
		int			*pnSFFm,		/* 曲面网格单元节点编号 **/
		int			*pnSFTp,		/* 曲面网格单元类型。当前仅支持三角形单元 **/
		int			*pnSFPt,		/* 曲面网格单元所在几何面编号,从1开始 **/
		int			nSF,			/* 曲面网格单元数目 **/
		std::map<int, int> pnFT,			/* 几何面类型： 0为远场； 1 为物面； 2为对称面 ,3为不长边界层的面,4为周期性面**/
		int			nLN,			/* 边界层层数 **/
		std::vector<int> layer_vec, /* 边界层层数数组 **/
		double		dLen,			/* 边界层第一层厚度 **/
		std::vector<double> length_vec, /* 边界层第一层厚度数组 **/
		double		dRto,			/* 边界层厚度增长因子 **/
		bool		bisostop,       /* 各向同性停止**/
		int			nopt,           /* 优化次数 **/
		std::function<double(std::array<double, 3>)> sizefunction, /*尺寸函数，输入为点坐标，输出为尺寸值**/
		/* ------------------------- 输出参数 -------------------------**/
		double	   **ppdMNC,		/* 体网格节点坐标 **/
		int         *pnMN,			/* 体网格节点数目 **/
		int         **ppnMEFm,		/* 体网格单元节点编号 **/
		int         **ppnMETp,		/* 体网格单元类型。当前支持单元类型：三棱柱、金字塔、四面体 **/
		int         *pnME,			/* 体网格单元数目 **/
		/* ------------------------- 边界信息 -------------------------**/
		int         *num_boundary_face, /* 边界面网格数量 **/
		int         **boundary_mesh, /* 边界面网格，注意每四个为一组，而且注意如果为三角形，最后一项为-1，id从0开始 **/
		int         **boundary_face, /* 边界面网格对应的面id，长度为num_boundary_face **/
		/* ------------------------- 其他参数 -------------------------**/
		bool b_use_multiple_normals, /* 是否启用多法向 **/
		bool b_output_io_file,  /* 是否将api的输入和输出都输出到文件系统中（仅用于DEBUG）**/
		std::string filename,   /* 几何文件名，缺省为virtualmesh **/
		/* ------------------------- 控制参数 -------------------------**/
		double expan_beta       /* 尺寸过渡因子**/
	) {
		VM v;
		if (checkterminate()) {
			*ppdMNC = v.ppdMNC;
			*pnMN = v.pnMN;
			*pnME = v.pnME;
			*ppnMEFm = v.ppnMEFm;
			*ppnMETp = v.ppnMETp;
			*boundary_face = v.boundary_face;
			*boundary_mesh = v.boundary_mesh;
			*num_boundary_face = v.num_boundary_face;
			return -1;
		}
		if (nLN == 0) {
			/** 退化为四面体程序 ***/
			API_Gen_Tetra_Mesh(

				/* ------------------------- 输入参数 --------------------------**/
				pdSNC,			/* 曲面网格节点坐标，coord[3*i], coord[3*i+1], coord[3*i+2]为第i个节点x,y,z坐标 **/
				nSN,			/* 曲面网格边界节点数目 **/
				pnSFFm,		/* 曲面网格单元节点编号 **/
				pnSFTp,		/* 曲面网格单元类型。当前仅支持三角形单元 **/
				pnSFPt,		/* 曲面网格单元所在几何面编号,从1开始 **/
				nSF,			/* 曲面网格单元数目 **/
				nopt,           /* 优化次数 **/
				sizefunction, /* 尺寸函数，输入为点坐标，输出为尺寸值 **/
				/* ------------------------- 输出参数 -------------------------**/
				ppdMNC,		/* 体网格节点坐标 **/
				pnMN,			/* 体网格节点数目 **/
				ppnMEFm,		/* 体网格单元节点编号 **/
				ppnMETp,		/* 体网格单元类型。当前支持单元类型：四面体 **/
				pnME,			/* 体网格单元数目 **/ 
				b_output_io_file,    /* 是否将api的输入和输出都输出到文件系统中（仅用于DEBUG）**/ 
				filename,      /* 几何文件名，缺省为virtualmesh **/ 

				/* ------------------------- 控制参数 -------------------------**/
				expan_beta       /* 尺寸过渡因子 **/

			);
			*num_boundary_face = 0;
			*boundary_mesh = nullptr;
			*boundary_face = nullptr;
			return 0;

		}
		tiger::Config cf;
		cf.SetDefaultConfig();
		cf.isotropic_stop = 1;
		cf.dStepLen = dLen;
		cf.dStepLenRatio = dRto;
		cf.nLayerNum = nLN;
		cf.sGeoFileName = filename;
		vector<vector<int>> bcs;
		bcs.resize(10);
		//vector<int> symm;
		//vector<int> box;
		//vector<int> match;
		//vector<int> per;
		for (auto i : pnFT) {
			bcs[i.second].push_back(i.first);
		}
		cf.vecBoxFc = bcs[0];
		cf.vecSymmFc = bcs[2];
		cf.vecMatchFc = bcs[4];
		cf.vecAdjacentFc = bcs[5];
		//cf.WriteConfigFile();

		int npt = nSN;
		int nelm = nSF;
		//ofstream fout("virtualmesh.pls");
		stringstream fout;

		fout << nelm;
		fout << " " << npt << " ";
		fout << 0 << " " << 0 << " " << 0 << " " << 0 << endl;
		std::vector<std::array<double, 3>> points;
		for (int i = 0; i < npt; i++) {
			points.push_back(std::array<double, 3>{pdSNC[3 * i + 0], pdSNC[3 * i + 1], pdSNC[3 * i + 2]});
		}
		for (int i = 0; i < nelm; i++) {
			fout << i + 1 << " " << pnSFFm[3 * i + 1] << " " << pnSFFm[3 * i + 0] << " " << pnSFFm[3 * i + 2] << " " << pnSFPt[i] << endl;
		}
		//fout.close();
		string input(fout.str());
		if (b_output_io_file) {
			ofstream ftest("VolumeMeshing.txt");
			ftest.setf(ios::fixed, ios::floatfield);  // 设定为 fixed 模式，以小数点表示浮点数
			ftest.precision(17);
			ftest << "expan_beta: " << expan_beta << endl;
			ftest << "nLN: " << nLN << endl;
			ftest << "dLen: " << dLen << endl;
			ftest << "dRto: " << dRto << endl;
			ftest << "bisostop: " << bisostop << endl;
			ftest << "boundary_info: "  << endl;


			for (auto i:pnFT) {
				ftest << i.first << " " << i.second << " " << std::endl;;
			}
			ftest << endl;
			ftest << "b_use_multiple_normals: " << b_use_multiple_normals << endl;
			ftest << nelm << " " << npt << " " << endl;
			std::vector<std::array<double, 3>> points;
			for (int i = 0; i < npt; i++) {
				ftest << i + 1 << " " << pdSNC[3 * i + 0] << " " << pdSNC[3 * i + 1] << " " << pdSNC[3 * i + 2] << endl;
			}
			for (int i = 0; i < nelm; i++) {
				ftest << i + 1 << " " << pnSFFm[3 * i + 1] << " " << pnSFFm[3 * i + 0] << " " << pnSFFm[3 * i + 2] << " " << pnSFPt[i] << endl;
			}

			ftest.close();

		}


		blpreConfig blconfig;
		blconfig.box = bcs[0];
		blconfig.symm = bcs[2];
		blconfig.n = nLN;
		blconfig.Ro = dRto;
		blconfig.len = dLen;
		blconfig.match = bcs[3];
		blconfig.per = bcs[4];
		blconfig.adjacent = bcs[5];
		blconfig.use_multiple_normals = b_use_multiple_normals;
		ControlVolume cv;
		auto bdyfile = PRE::blpre((input), blconfig, points,cv);

		if (checkterminate()) {
			*ppdMNC = v.ppdMNC;
			*pnMN = v.pnMN;
			*pnME = v.pnME;
			*ppnMEFm = v.ppnMEFm;
			*ppnMETp = v.ppnMETp;
			*boundary_face = v.boundary_face;
			*boundary_mesh = v.boundary_mesh;
			*num_boundary_face = v.num_boundary_face;
			return -1;
		}

		try {
			v = blmesh(bdyfile, blconfig, true, false, true, bisostop, b_output_io_file, sizefunction, expan_beta, nopt);
		}
		catch (std::exception e) {
			*ppdMNC = nullptr;
			*pnMN = 0;
			*pnME = 0;
			*ppnMEFm = nullptr;
			*ppnMETp = nullptr;
			throw e;
		}
		*ppdMNC = v.ppdMNC;
		*pnMN = v.pnMN;
		*pnME = v.pnME;
		*ppnMEFm = v.ppnMEFm;
		*ppnMETp = v.ppnMETp;
		*boundary_face = v.boundary_face;
		*boundary_mesh = v.boundary_mesh;
		*num_boundary_face = v.num_boundary_face;

		///merge cv
		///...
		if (cv.f.size()) {
			int old_point_num = *pnMN;
			double* nppdMNC = new double[(*pnMN + cv.lower_point_num) * 3];
			for (int i = 0; i < *pnMN * 3; i++) {
				nppdMNC[i] = *ppdMNC[i];
			}
			for (int i = 0; i < cv.lower_point_num; i++) {
				for (int k = 0; k < 3; k++) {
					nppdMNC[*pnMN * 3 + 3 * i + k] = cv.v[i][k];
				}
			}
			*pnMN = *pnMN + cv.lower_point_num;
			//delete* ppdMNC[];
			ppdMNC = &nppdMNC;
			int* nppnMEFm = new int[6 * (*pnME) + 6 * cv.f.size()];
			int cnt = 0;
			for (int i = 0; i < *pnMN; i++) {
				int k = 0;
				if (*ppnMETp[i] == 11) {
					k = 4;
				}
				if (*ppnMETp[i] == 13) {
					k = 6;
				}
				if (*ppnMETp[i] == 14) {
					k = 5;
				}
				for (int j = 0; j < k; j++) {
					nppnMEFm[cnt] = *ppnMEFm[cnt];
					cnt++;
				}
			}

			for (int i = 0; i < cv.f.size(); i++) {
				int k = cv.f[i].size();
				for (int j = 0; j < k; j++) {
					if (cv.f[i][k] > cv.lower_point_num)
						nppnMEFm[cnt] = cv.f[i][k] - cv.lower_point_num;
					else
						nppnMEFm[cnt] = cv.f[i][k] + old_point_num;
					cnt++;
				}
			}
			int* nppnMETp = new int[*pnME + cv.f.size()];
			for (int i = 0; i < *pnME; i++) {
				nppnMETp[i] = *ppnMETp[i];
			}
			for (int i = 0; i < cv.f.size(); i++) {
				int k = cv.f[i].size();
				if (k == 4) { nppnMETp[*pnME + i] = 11; }
				if (k == 5) { nppnMETp[*pnME + i] = 14; }
				if (k == 6) { nppnMETp[*pnME + i] = 13; }
			}

			*pnME = *pnMN + cv.f.size();
			ppnMEFm = &nppnMEFm;
			ppnMETp = &nppnMETp;
		}
		return 0;

	}
	int API_Gen_Boundary_Mesh(
		/* ------------------------- 输入参数 --------------------------**/
		double* pdSNC,			/* 曲面网格节点坐标，coord[3*i], coord[3*i+1], coord[3*i+2]为第i个节点x,y,z坐标 **/
		int			nSN,			/* 曲面网格边界节点数目 **/
		int* pnSFFm,		/* 曲面网格单元节点编号 **/
		int* pnSFTp,		/* 曲面网格单元类型。当前仅支持三角形单元 **/
		int* pnSFPt,		/* 曲面网格单元所在几何面编号 ,从1开始 **/
		int			nSF,			/* 曲面网格单元数目 **/
		std::map<int, int> pnFT,			/* 几何面类型： 0为远场； 1 为物面； 2为对称面 ,3为不长边界层的面,4为周期性面**/
		int			nLN,			/* 边界层层数 **/
		std::vector<int> layer_vec, /* 边界层层数数组 **/
		int			max_layer_diff, /*相邻网格边界层层数差*/
		double		dLen,			/* 边界层第一层厚度 **/
		std::vector<double> length_vec, /* 边界层第一层厚度数组 **/
		double		dRto,			/* 边界层厚度增长因子 **/
		double		max_skewnwass,  /* 各向异性停止**/
		double		bisostop,       /* 各向同性停止**/
		/* ------------------------- 输出参数 -------------------------**/
		double** ppdMNC,		/* 体网格节点坐标 **/
		int* pnMN,			/* 体网格节点数目 **/
		int** ppnMEFm,		/* 体网格单元节点编号 **/
		int** ppnMETp,		/* 体网格单元类型。当前支持单元类型：三棱柱，金字塔 **/
		int* pnME,			/* 体网格单元数目 **/
		/* ------------------------- 边界层网格顶面层面网格参数 -------------------------**/
		double** ppdSNC0,		/* 顶面网格节点坐标，coord[3*i], coord[3*i+1], coord[3*i+2]为第i个节点x,y,z坐标 **/
		int* pnSN0,			/* 顶面网格边界节点数目 **/
		int* pnSEO,			/* 顶面网格单元数目 **/
		int** ppnSFTpO,		/* 顶面网格单元类型。当前仅支持三角形单元 **/
		int** ppnSFFmO,			/* 体网格单元数目 **/
		int** l2g,             /*  顶面网格id到全局点id的映射  */
		double** ppnsizing, /* 顶面网格目标尺寸 */
		/////* ------------------------- 边界信息 ------------------------- **/ 
		int* num_boundary_face,  /////* 边界面网格数量 **/ 
		int** boundary_mesh,  /////* 边界面网格，注意每四个为一组，而且注意如果为三角形，最后一项为-1，id从0开始 **/ 
		int** boundary_face,  /////* 边界面网格对应的面id，长度为num_boundary_face **/ 
		/* ------------------------- 其他参数 -------------------------**/
		bool b_have_pyramid, /* 是否有金字塔 **/
		bool b_use_multiple_normals, /* 是否启用多法向 缺省为false **/
		bool b_output_io_file,  /* 是否将api的输入和输出都输出到文件系统中（仅用于DEBUG）缺省为false **/
		std::string filename,   /* 几何文件名，缺省为virtualmesh **/
		std::array<double, 12> per_matrix  /* 周期性面控制矩阵,前9位为旋转矩阵 m00，m01，m02 .... ，后三位为位移向量xyz **/
	) {
		tiger::Config cf;
		cf.SetDefaultConfig();
		cf.isotropic_stop = 1;
		cf.dStepLen = dLen;
		cf.dStepLenRatio = dRto;
		cf.nLayerNum = nLN;
		cf.nInitLayerNum_vec = layer_vec;
		cf.dStepLen_vec = length_vec;
		cf.sGeoFileName = filename;
		vector<vector<int>> bcs;
		bcs.resize(10);
		//vector<int> symm;
		//vector<int> box;
		//vector<int> match;
		//vector<int> per;
		for (auto i: pnFT) {
			bcs[i.second].push_back(i.first);
		}

		cf.vecBoxFc = bcs[0];

		
		cf.vecSymmFc = bcs[2];
		cf.vecMatchFc = bcs[4];
		cf.vecAdjacentFc = bcs[5];
		//cf.WriteConfigFile();

		int npt = nSN;
		int nelm = nSF;
		//ofstream fout("virtualmesh.pls");
		stringstream* fout = new stringstream();


		*fout << nelm << " " << npt << " " << 0 << " " << 0 << " " << 0 << " " << 0 << endl;

		std::vector<std::array<double, 3>> points;
		for (int i = 0; i < npt; i++) {
			points.push_back(std::array<double, 3>{pdSNC[3 * i + 0], pdSNC[3 * i + 1], pdSNC[3 * i + 2]});
		}
		for (int i = 0; i < nelm; i++) {
			*fout << i + 1 << " " << pnSFFm[3 * i + 1] << " " << pnSFFm[3 * i + 0] << " " << pnSFFm[3 * i + 2] << " " << pnSFPt[i] << endl;
		}
		//fout.close();
		string input(fout->str());

		if (b_output_io_file) {
			ofstream ftest("BoundaryMeshing.txt");
			ftest.setf(ios::fixed, ios::floatfield);  // 设定为 fixed 模式，以小数点表示浮点数
			ftest.precision(17);
			ftest << "nLN: " << nLN << endl;
			ftest << "dLen: " << dLen << endl;
			ftest << "dRto: " << dRto << endl;
			ftest << "max_skewnwass: " << max_skewnwass << endl;
			ftest << "bisostop: " << bisostop << endl;
			ftest << "boundary_info: " << endl;
			for (auto i : pnFT) {
				ftest << i.first << " " << i.second << " " << std::endl;;
			}
			ftest << "b_use_multiple_normals: " << b_use_multiple_normals << endl;
			ftest << nelm << " " << npt << " " << endl;
			std::vector<std::array<double, 3>> points;
			for (int i = 0; i < npt; i++) {
				ftest << i + 1 << " " << pdSNC[3 * i + 0] << " " << pdSNC[3 * i + 1] << " " << pdSNC[3 * i + 2] << endl;
			}
			for (int i = 0; i < nelm; i++) {
				ftest << i + 1 << " " << pnSFFm[3 * i + 1] << " " << pnSFFm[3 * i + 0] << " " << pnSFFm[3 * i + 2] << " " << pnSFPt[i] << endl;
			}

			ftest.close();
		}

		blpreConfig blconfig;
		blconfig.n = nLN;
		blconfig.Ro = dRto;
		blconfig.len = dLen;
		blconfig.layer_vec = layer_vec;
		blconfig.len_vec = length_vec;
        blconfig.box = bcs[0];
        blconfig.wall = bcs[1];
        blconfig.symm = bcs[2];
		blconfig.match = bcs[3];
		blconfig.per = bcs[4];
		blconfig.adjacent = bcs[5];
		blconfig.use_multiple_normals = b_use_multiple_normals;
		blconfig.max_equal_skewnwass = max_skewnwass;
		blconfig.max_layer_diff = max_layer_diff;
		ControlVolume cv;
		auto bdyfile = PRE::blpre(input, blconfig, points,cv);
		delete fout;



		VM v = blmesh(bdyfile, blconfig, false, true, b_have_pyramid, bisostop, b_output_io_file, nullptr, 1.2, 0, per_matrix);

		
		*ppdMNC = v.ppdMNC;
		*pnMN = v.pnMN;
		*pnME = v.pnME;
		*ppnMEFm = v.ppnMEFm;
		*ppnMETp = v.ppnMETp;
		*pnSEO = v.pnSNO;
		*ppnSFTpO = v.ppnSFTpO;
		*ppnSFFmO = v.ppnSFFmO;
		*ppdSNC0 = v.ppdSNC0;
		*pnSN0 = v.nSN0;
                *l2g = v.l2g;
                *ppnsizing = v.sizing;


		*boundary_face = v.boundary_face;
		*boundary_mesh = v.boundary_mesh;
		*num_boundary_face = v.num_boundary_face;

		if (cv.f.size()) {
			int old_point_num = *pnMN;
			double* nppdMNC = new double[(*pnMN + cv.lower_point_num+cv.add_point_num) * 3];
			std::map<std::array<double, 3>, int> coord_to_id;
			for (int i = 0; i < (*pnMN) ; i++) {
				coord_to_id[{(*ppdMNC)[3 * i + 0], (*ppdMNC)[3 * i + 1], (*ppdMNC)[3 * i + 2]}] = i;
			}
			for (int i = 0; i < (*pnMN) * 3; i++) {
				nppdMNC[i] = (*ppdMNC)[i];
			}
			for (int i = 0; i < cv.lower_point_num; i++) {
				for (int k = 0; k < 3; k++) {
					nppdMNC[(*pnMN) * 3 + 3 * i + k] = cv.v[i][k];
				}
			}
			for (int i = 0; i < cv.add_point_num; i++) {
				for (int k = 0; k < 3; k++) {
					nppdMNC[(*pnMN) * 3+3* cv.lower_point_num + 3 * i + k] = cv.v[cv.v.size()-cv.add_point_num+i][k];
				}
			}
			*pnMN = *pnMN + cv.lower_point_num + cv.add_point_num;
			//delete* ppdMNC[];
			*ppdMNC = nppdMNC;
			int* nppnMEFm = new int[6 * (*pnME) + 6 * cv.f.size()];
			int cnt = 0;
			for (int i = 0; i < *pnME; i++) {
				int k = 0;
				if ((*ppnMETp)[i] == 11) {
					k = 4;
				}
				if ((*ppnMETp)[i] == 13) {
					k = 6;
				}
				if ((*ppnMETp)[i] == 14) {
					k = 5;
				}
				for (int j = 0; j < k; j++) {
					nppnMEFm[cnt] = (*ppnMEFm)[cnt];
					cnt++;
				}
			}
			int* nppnMETp = new int[*pnME + cv.f.size()];
			for (int i = 0; i < *pnME; i++) {
				nppnMETp[i] = (*ppnMETp)[i];
			}
			for (int i = 0; i < cv.f.size(); i++) {
				int k = cv.f[i].size();
				for (int j = 0; j < k; j++) {
					if (cv.f[i][j] >= cv.lower_point_num) {
						if(cv.f[i][j]>= cv.v.size() -cv.add_point_num)
							nppnMEFm[cnt] = cv.f[i][j] -cv.v.size() +cv.add_point_num + old_point_num +cv.lower_point_num;
						else {
							//nppnMEFm[cnt] = cv.f[i][j] - cv.lower_point_num;
							nppnMEFm[cnt] = coord_to_id[cv.v[cv.f[i][j]]];
							//if (cv.v[cv.f[i][j]][0] != (*ppdMNC)[3 * cv.f[i][j] - cv.lower_point_num]) {
							//	std::cout << "debug here";
							//}
						}
					}
					else
						nppnMEFm[cnt] = cv.f[i][j] + old_point_num;
					cnt++;
				}

				if (k == 4) { nppnMETp[*pnME + i] = 11; }
				if (k == 5) { nppnMETp[*pnME + i] = 14; }
				if (k == 6) { nppnMETp[*pnME + i] = 13; }
			}

			*pnME = *pnME + cv.f.size();
			*ppnMEFm = nppnMEFm;
			*ppnMETp = nppnMETp;
		}
		return 0;

	}

	int API_Mesh_Optimize(double * pdMNC, int num_boundary_face, int* boundary_mesh, int nMN, int * pnMEFm, int * pnMETp, int nME)
	{
		spdlog::info("Golbal Point Optimization processing...");
		std::vector<bool> point_map(nMN, false);
		std::vector<BLVector> point_coordinate;
		for (int i = 0; i < nMN; i++) {
			point_coordinate.push_back(pdMNC + 3 * i);
		}
		std::vector<vector<int>> ele_map(nME);
		std::set<int> near_tetra;
		int count = 0;
		for (int i = 0; i < nME; i++) {
			int type = pnMETp[i];

			switch (type)
			{
			case 11://tetra
				for (int j = 0; j < 4; j++)
					ele_map[i].push_back(pnMEFm[count++]);
				break;
			case 13://prism
				for (int j = 0; j < 6; j++)
					ele_map[i].push_back(pnMEFm[count++]);
				break;
			case 14://pyramid
				for (int j = 0; j < 5; j++)
					ele_map[i].push_back(pnMEFm[count++]);
				break;
			default:
				break;
			}
		}
		std::set<int> point_id_s;
		for (int i = 0; i < num_boundary_face; i++) {
			for (int k = 0; k < 4; k++)
				point_id_s.insert(boundary_mesh[4 * i + k]);
		}
		//twice optimization
		while (MeshOptimize(point_coordinate, ele_map, point_id_s)) {}
		for (int i = 0; i < nMN; i++) {
			pdMNC[3 * i + 0] = point_coordinate[i].x;
			pdMNC[3 * i + 1] = point_coordinate[i].y;
			pdMNC[3 * i + 2] = point_coordinate[i].z;

		}

		spdlog::info("Golbal Point Optimization Done");
		if (checkterminate())
			return -1;
		return 0;
	}
	static bool is_point_equal(double* p1, double* p2) {
		return p1[0] == p2[0] && p1[1] == p2[1] && p1[2] == p2[2];
	}
	int API_Mesh_Merge(double *  pdMNC1, int  nMN1, int *  pnMEFm1, int *  pnMETp1, int   nME1, 
				       double *  pdMNC2, int  nMN2, int *  pnMEFm2, int *  pnMETp2, int   nME2,
					   double ** pdMNC3, int *nMN3, int ** pnMEFm3, int ** pnMETp3, int * nME3)
	{
		Eigen::MatrixXd V1;
		Eigen::MatrixXd V3;
		V1.resize(nMN1 + nMN2, 3);
		for (int i = 0; i < nMN1; i++) {
			for (int j = 0; j < 3; j++)
				V1(i, j) = pdMNC1[3 * i + j];
		}
		//	V2.resize(nMN2, 3);
		for (int i = nMN1; i < nMN1 + nMN2; i++) {
			for (int j = 0; j < 3; j++)
				V1(i, j) = pdMNC2[3 * (i-nMN1) + j];
		}
		Eigen::MatrixXi SVI;
		Eigen::MatrixXi SVJ;
		igl::remove_duplicate_vertices(V1, 0, V3, SVI, SVJ);
		
		*pdMNC3 = new double[V3.rows()*3];
		for (int i = 0; i < V3.rows(); i++) {
			for (int j = 0; j < 3; j++)
				(*pdMNC3)[3 * i + j] = V3(i, j);
		}
		*nMN3 = V3.rows();
		*nME3 = nME1 + nME2;
		*pnMETp3 = new int[*nME3];
		for (int i = 0; i < nME1; i++) {
			(*pnMETp3)[i] = pnMETp1[i];
		}
		for (int i = nME1; i < nME1+ nME2; i++) {
			(*pnMETp3)[i] = pnMETp2[i- nME1];
		}
	
		

		int count1 = 0;
		for (int i = 0; i < nME1; i++) {
			int type = pnMETp1[i];
			switch (type)
			{
			case 11://tetra
				count1 += 4;		
				break;
			case 13://prism
				count1 += 6;
				break;
			case 14://pyramid
				count1 += 5;
				break;
            case 12:
                count1 += 8;
                break;
			default:
				break;
			}
		}



		int count2 = 0;
		for (int i = 0; i < nME2; i++) {
			int type = pnMETp2[i];
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
            case 12:
                count2 += 8;
                break;
			default:
				break;
			}
		}

		*pnMEFm3 = new int[count1+count2];
		for (int i = 0; i < count1; i++) {
			(*pnMEFm3)[i] = SVJ(pnMEFm1[i]);
		}
		for (int i = count1; i < count1+count2; i++) {
			(*pnMEFm3)[i] = SVJ(pnMEFm2[i-count1]+nMN1);
		}
		return 0;
	}


}
