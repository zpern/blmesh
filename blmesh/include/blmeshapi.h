#pragma once
#ifndef _BLMESH_API_H
#define _BLMESH_API_H
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#include <string>
#include <vector>
#include <functional>
#include <array>
#include <tuple>
#include <map>
#include <set>
#include "BLMesh_define.h"
#include "VirtualMesh.h"
#include "blpre.h"
//#define MEMORY_DEBUG //on - check memory leak  default off
//#define VIRTUAL_MESH
#define USE_MEMORY_POOL //on-use memory pool ,but still has bug
#define APIFUNC
//#define OUTPUT_SPERATELY
#define CHECK_SKEWNWSS
#define FRONT_ORDER 1 //1:hirbert sort 2:nature order 3:random
#define REFER_DS 1 //1:array 2:hashtable 3:balance tree
#define ONED_FILTER 1//1:enable 0:disable
#define CHECK_VOLUMN
#define CHECK_EVERY_VOLUMN
//#define USE_OPENMP
#define USE_DYNAMIC_ARRAY
//#undef DEBUG
//#define BRUTE_FORCE_CHECK //on - check valid but time-consuming only used in debug 
//#define OUTPUT_SPERATELY
//#define _DEBUG //on -  print debug infomation to cmd windows default - off
//#define USE_TETGEN
//#define USE_TRIANGLE
//#define ORTHO
//#ifdef NDEBUG
//#define DEBUG
//#endif
//#define AVERC
//#define CHANGE_STEP_BY_DISTANCE
//#define SMOOTH_FRONT_SIZE
#define SMOOTH_HOR_NORMAL
#define GEN_PYRAMID
#define GEN_TETRA
//#define GEN_OUTER_MESH



#define VERSION "2.0"

//multi-pyramid option
#define MAX_DIFF_RATIO 0.25 //define the max height of pyramid height/prism length
#define MIN_LAYER 5





#ifdef APIFUNC 
VM blmesh(std::tuple<std::string, double*, int*, int*, std::vector<double>> sfile,
	blpreConfig blconfig,
	bool genoutmesh,
	bool genboundarymesh,
	bool havepyramid,
	double isostop,
	bool outputIO,
	std::function<double(std::array<double, 3>)> sizefunction,
	double expan_beta,
	int number_opt,
	std::array<double, 12> per_matrix = std::array<double, 12>());
#endif

#endif
