//  Copyright 2020 hfye@zju.edu.cn
/* 
 * main.cpp
 *
 *  Created on: Sep 23, 2015
 *      Author: Zhoufang Xiao
 **/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <omp.h>
#include <memory>

#include <spdlog/spdlog.h> 
 #include "./BLMesh.h"
#include "./BLFront3D.h"
#include "./BLMesh_define.h"
#include "./blmeshapi.h"
#include "./blpre.h"

#ifdef WIN32 
#include <windows.h>
#include <psapi.h>
#pragma comment(lib,"psapi.lib")
void showMemoryInfo(void)
{
	HANDLE handle = GetCurrentProcess();
	PROCESS_MEMORY_COUNTERS pmc;
	GetProcessMemoryInfo(handle, &pmc, sizeof(pmc));

	//cout << "内存使用:" << pmc.WorkingSetSize / 1024/1024 << "M/" << pmc.PeakWorkingSetSize / 1024 / 1024 << "M + " << pmc.PagefileUsage / 1024 / 1024 << "M/" << pmc.PeakPagefileUsage / 1024 / 1024 << "M" << endl;
}

BOOL WINAPI signal_handler(DWORD dwCtrlType) {
	if (dwCtrlType == CTRL_C_EVENT) {
		//cout << "Received signal " << dwCtrlType << ", stopping program..." << endl;
	//	throw std::runtime_error("Received signal stopping program...");
		exit(dwCtrlType);
	}
	return TRUE;
}


//#include <vld.h>
#endif


#define _CONFIG_FILE

//#define _2D_BLMESH
#define _GEN_BLMESH
#include "Octree.h"
ConfigArgc cf;
int main(int argc, char* argv[]) {
	double start_time = clock();
#ifndef USE_OPENMP
	//omp_set_num_threads(1);
#endif // USE_OPENMP

#ifdef WIN32 
	SetConsoleCtrlHandler(signal_handler, TRUE);
#endif
	cout << "Program started. Press Ctrl+C to stop." << endl;
	string bdryfile;
	string finconfig, bdry_file, volumn_mesh_file_vtk;
	parsecommand(argc, argv, cf);
	finconfig = std::string(cf.filename) + "_o";
	bdryfile = finconfig + ".pls";
	bdry_file = string(cf.filename) + ".bdry";
	volumn_mesh_file_vtk = string(cf.filename) + ".vtk";
	std::shared_ptr<BLMesh> blmesh(new BLMesh(cf.layer_num, 0, cf.isotropic_stop));
	if (argc == 3) {
		blmesh->max_depth_ = stoi(string(argv[1]));
		blmesh->max_obj_ = stoi(string(argv[2]));
	}
	ifstream inputfile(string(cf.filename + string(".pls")));
	stringstream ss;
	ss << inputfile.rdbuf();
	bool use_multi_normal = false;
	string smesh = ss.str();
	inputfile.close();
	blpreConfig blcf;
	blcf.box = cf.boxfc;
	blcf.len = cf.step_len;
	blcf.n = cf.layer_num;
	blcf.Ro = cf.ratio;
	blcf.symm = cf.symfc;
	blcf.use_multiple_normals = false;
	blcf.match = cf.matchfc;
	blcf.per = cf.perfc;
	auto ofile = PRE::blpre(smesh, blcf);
	use_multi_normal = false;
	/*删除内存**/
	ss.clear();
	inputfile.clear();
	smesh.clear();

	blmesh->InitBLMesh();
	PrintConfig(&cf);
	blmesh->cf = cf;
	blmesh->SetBoundary(ofile);

	std::get<0>(ofile).clear();

	double gen_mesh_start_time = clock();

	if (use_multi_normal) {
		blmesh->CalZeroNorm();
	}

	try {
		blmesh->GenerateBLMesh();
	}
	catch (std::exception error) {
		spdlog::info(error.what() );
		return -1;
	}
	double gen_bdymesh_end_time = clock();	
#ifdef WIN32
	showMemoryInfo();
#endif // WIN32

	std::cout.setf(ios::left, ios::adjustfield);
	std::cout.fill(' ');
	spdlog::info("|Generate Boundary Layer Mesh Cost|{}",  (gen_bdymesh_end_time - gen_mesh_start_time) / CLOCKS_PER_SEC);// '|' );
	spdlog::info("|Total Check Insertion Cost       |{}",  blmesh->check_prism_time / CLOCKS_PER_SEC);// );
	spdlog::info("|Total Normal Smoothing Cost      |{}", blmesh->smooth_time / CLOCKS_PER_SEC );//);
	std::cout.setf(ios::left, ios::adjustfield);
	std::cout.fill(' ');
	//spdlog::info(std::setiosflags(ios::left) << std::setfill(' ') );
	spdlog::info("             Mesh Information    " );
	spdlog::info("+================================================+");
	spdlog::info("|             Item                |number of cell|" );
	spdlog::info("-------------------------------------------------" );
	spdlog::info("|            Pyramid              |{}" , blmesh->m_nPyramid);// );
	spdlog::info("|            Tetrahedron          |{}" , blmesh->m_nTet );//);
	spdlog::info("|         Triangular prism        |{}" , blmesh->m_nPrism );//);
	spdlog::info("|            Triangle             |{}"  ,blmesh->m_nTri );//);
	spdlog::info("|          Quadrilateral          |{}"  ,blmesh->m_nQuad );//);

#ifdef WIN32
	showMemoryInfo();
#endif
	spdlog::info("+================================================+");
	


	/*auto nodes = blmesh->GetNarrowConstrainedPointidx();
	for (auto i : nodes) {
		cout << i << " ";
	}**/


	if (use_multi_normal)
		blmesh->RemoveOverlapNodeAndElement();

#ifdef _DEBUG
	//if(cf.surface_normal == 1)
	blmesh->SaveSurfGrid2("surfacemesh.vtk");
#endif

#ifdef _DEBUG
	//cf.smooth_attempt = 0;
#endif

#ifdef OUTPUT_SPERATELY
	blmesh->SaveBLMesh(const_cast<char*>(("_boundary_"+volumn_mesh_file_vtk).c_str()));
#endif
	//VM v;
	//blmesh->GenTopMesh(v);
	try {
		blmesh->GenOuterMesh(cf.smooth_attempt);
	}
	catch (string error) {
		spdlog::info(error );
		return 0;
	}


	
	blmesh->RemvNodElm();

	blmesh->GlobalOptmize();

	blmesh->FreeMemoryInFrontAndNode();

	double gen_mesh_end_time = clock();
#ifdef OUTPUT_SPERATELY
	blmesh->SaveTetMesh(const_cast<char*>(("_tet_" + volumn_mesh_file_vtk).c_str()));
#endif
	blmesh->SaveBLMesh(const_cast<char*>(volumn_mesh_file_vtk.c_str()));
	blmesh->SaveBdry(const_cast<char*>(bdry_file.c_str()));

	spdlog::info("Average number of layers{}",  blmesh->GetAvergeLayer() );
	double end_time = clock();



	std::cout.setf(ios::left, ios::adjustfield);
	std::cout.fill(' ');
	//spdlog::info(std::setiosflags(ios::left) << std::setfill(' ') );
	spdlog::info("             Mesh Information    " );
	spdlog::info("+================================================+" );
	spdlog::info("|             Item                |number of cell|" );
	spdlog::info("-------------------------------------------------" );
	spdlog::info("|            Pyramid              |{}", blmesh->m_nPyramid );
	spdlog::info("|            Tetrahedron          |{}", blmesh->m_nTet );
	spdlog::info("|         Triangular prism        |{}", blmesh->m_nPrism );
	spdlog::info("|            Triangle             |{}", blmesh->m_nTri );
	spdlog::info("|          Quadrilateral          |{}", blmesh->m_nQuad );
	spdlog::info("+================================================+" );





	std::cout.setf(ios::left, ios::adjustfield);
	std::cout.fill(' ');
	//spdlog::info(setiosflags(ios::left) << setfill(' ') );
	spdlog::info("              Timing Cost Table" );
	spdlog::info("+================================================+" );
	spdlog::info("|             Item                |Time(second)  |" );
	spdlog::info("-------------------------------------------------" );
	spdlog::info("|Generate Boundary Layer Mesh Cost|{}", (gen_bdymesh_end_time - gen_mesh_start_time) / CLOCKS_PER_SEC );
	spdlog::info("|Total Prism Check Insertion Cost |{}", blmesh->check_prism_time / CLOCKS_PER_SEC );
	spdlog::info("|Total Prism Check valid Cost     |{}", (blmesh->check_prism_quality_time - blmesh->check_prism_time) / CLOCKS_PER_SEC );
	spdlog::info("|Total Normal Smoothing Cost      |{}", blmesh->smooth_time / CLOCKS_PER_SEC );
	spdlog::info("|Total Pyramid generation Cost    |{}", blmesh->generate_pyramid_time / CLOCKS_PER_SEC );
	spdlog::info("|Other process In bdymesh Cost    |{}", (gen_bdymesh_end_time - gen_mesh_start_time - blmesh->smooth_time - blmesh->check_prism_quality_time - blmesh->generate_pyramid_time) / CLOCKS_PER_SEC );
	spdlog::info("|Generate Tetrahedral Mesh Cost   |{}", (gen_mesh_end_time - gen_bdymesh_end_time) / CLOCKS_PER_SEC );
	spdlog::info("|Generate Mesh Cost        (no IO)|{}", (gen_mesh_end_time - gen_mesh_start_time) / CLOCKS_PER_SEC );
	spdlog::info("|Total Cost         (including IO)|{}", (end_time - start_time) / CLOCKS_PER_SEC );
	spdlog::info("+================================================+" );



#ifdef MEMORY_DEBUG
	showMemoryInfo()// If you want to test memory leak --- use dead loop: "while (1) {}"
#endif
}
