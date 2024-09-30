//#include <vld.h>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <array>
#include <iomanip>
#include <spdlog/spdlog.h> 
 #include "common_assert.h"
#include "vmesh.h"

#include "gtest/gtest.h"
//#include "test_sphere_tetra.h"
//#include "test_sphere_boundary.h"
//#include "test_f6_box_tetra.h"
//#include "test_f6_box_boundary.h"
//#include "test_cubic_tetra.h"
//#include "test_cubic_boundary.h"
//#include "test_volume_density.h"
//#include "test_non_manifold.h"
//#include "test_anistropic.h"
//#include "test_f6_box_vol.h"

//#include "test_f6_box_symm.h"
//#include "test_intersection_tetra.h"
//#include "test_intersection_vol.h"
//#include "test_intersection.h"
//#include "test_two_cubic.h"
#include "test_intersection.h"
//#include "test_postprocess_vol.h"

#define TESTVOLMESH 0
#define TESTTETRAMESH 1
#define TESTBOUNDARYMESH 2
using namespace std;

// enum BLmeshEntityTopology
//{
//	NODE = 5,
//	LINE = 6,
//	POLYGON = 7,
//	TRIANGLE = 8,
//	QUADRILATERAL = 9,
//	POLYHEDRON = 10,
//	TETRAHEDRON = 11,
//	HEXAHEDRON = 12,
//	PRISM = 13,
//	PYRAMID = 14,
//	SEPTAHEDRON = 15,
//	MIXED
//};

/*
调用的两个API，在tiger_sizingfunction.h中定义
TigerAdasize_API int API_Create_Cuboid_Sources_SF(
double local_coordinate_x[3],double local_coordinate_y[3],double local_coordinate_z[3], double local_coordinate_origin[3], double soValues,
int *sfObjID);

TigerAdasize_API int API_Create_Ellipsoid_Sources_SF(
double local_coordinate_x[3], double local_coordinate_y[3], double local_coordinate_z[3],double local_coordinate_origin[3], double soValues,
int *sfObjID);
*/


int main(int argc, char* argv[]) {
	
	testing::InitGoogleTest(&argc, argv);
	
	int obj_elli, obj_cube, obj_other;/*源标识ID*/
	//::testing::UnitTest::GetInstance()->run;
	RUN_ALL_TESTS();

	return 0;
//
//
//	int sum = 0, npr = 0, npy = 0, nte = 0;
//	for (int i = 0; i < pnME; i++) {
//		int type = (ppnMETp)[i];
//		if (type == EntityTopology::PRISM) {
//			sum += 6; npr++;
//		}
//		if (type == EntityTopology::PYRAMID) {
//			sum += 5; npy++;
//		}
//		if (type == EntityTopology::TETRAHEDRON) {
//			sum += 4; nte++;
//		}
//		//cout << type;
//	}
//
//	ofstream fout(filename.substr(0, filename.size() - 4) + ".vtk");
//	fout.setf(ios::scientific);
//	fout << "# vtk DataFile Version 2.0" << endl;
//	fout << "boundary layer mesh" << endl;
//	fout << "ASCII" << endl;
//	fout << "DATASET UNSTRUCTURED_GRID" << endl;
//	fout << "POINTS " << pnMN << " float" << endl;
//	for (int i = 0; i < pnMN; i++) {
//		fout << (ppdMNC)[3 * i] << " " << (ppdMNC)[3 * i + 1] << " " << (ppdMNC)[3 * i + 2] << endl;
//	}
//	fout << "CELLS " << pnME << " " << npr * 7 + nte * 5 + npy * 6 << endl;
//
//	int t = 0;
//	for (int i = 0; i < pnME; i++) {
//		int type = (ppnMETp)[i];
//		if (type == EntityTopology::PRISM) {
//			fout << setprecision(10) << 6 << " " << (ppnMEFm)[t] << " " << (ppnMEFm)[t + 1] << " " << (ppnMEFm)[t + 2] << " " << (ppnMEFm)[t + 3] << " " << (ppnMEFm)[t + 4] << " " << (ppnMEFm)[t + 5] << endl;
//			t += 6;
//		}
//		else if (type == EntityTopology::PYRAMID) {
//			fout << 5 << " " << (ppnMEFm)[t] << " " << (ppnMEFm)[t + 1] << " " << (ppnMEFm)[t + 2] << " " << (ppnMEFm)[t + 3] << " " << (ppnMEFm)[t + 4] << endl;
//			t += 5;
//		}
//		else if (type == EntityTopology::TETRAHEDRON) {
//			fout << 4 << " " << (ppnMEFm)[t] << " " << (ppnMEFm)[t + 1] << " " << (ppnMEFm)[t + 2] << " " << (ppnMEFm)[t + 3] << endl;
//			t += 4;
//		}
//	}
//	fout << "CELL_TYPES " << npr + nte + npy << endl;
//	for (int i = 0; i < pnME; i++)
//	{
//		int type = (ppnMETp)[i];
//		if (type == EntityTopology::PRISM)
//			fout << 13 << endl;
//		else if (type == EntityTopology::PYRAMID)
//			fout << 14 << endl;
//		else if (type == EntityTopology::TETRAHEDRON)
//			fout << 10 << endl;
//	}
//
//	fout.close();
//	if (testType == TESTBOUNDARYMESH) {
//		ofstream fout("toplayer.pls");
//
//		fout.setf(ios::scientific);
//		//*
//		fout << pnSEO << " " << pnSN0 << " 0 0 0 0" << endl;
//		for (int i = 0; i < pnSN0; i++) {
//			fout << i + 1 << " " << (ppdSNC0)[3 * i] << " " << (ppdSNC0)[3 * i + 1] << " " << (ppdSNC0)[3 * i + 2] << endl;
//		}
//
//		for (int i = 0; i < pnSEO; i++) {
//			//cout << i;
//			fout << i + 1 << " " << (ppnSFFmO)[3 * i] + 1 << " " << (ppnSFFmO)[3 * i + 1] + 1 << " " << (ppnSFFmO)[3 * i + 2] + 1 << " 1" << endl;
//		}//*/
//		/*
//		fout << "# vtk DataFile Version 2.0" << endl;
//		fout << "boundary layer mesh" << endl;
//		fout << "ASCII" << endl;
//		fout << "DATASET UNSTRUCTURED_GRID" << endl;
//		fout << "POINTS " <<  pnSN0<< " float" << endl;
//		for (int i = 0; i < pnSN0; i++) {
//			fout << (ppdSNC0)[3 * i] << " " << (ppdSNC0)[3 * i + 1] << " " << (ppdSNC0)[3 * i + 2] << endl;
//		}
//		fout << "CELLS " << pnSEO << " " << pnSEO*4 << endl;
//		//cout << "pnSNO=" << pnSNO << endl;
//		for (int i = 0; i < pnSEO; i++) {
//			//cout << i;
//				fout << 3 << " " << (ppnSFFmO)[3*i] << " " << (ppnSFFmO)[3*i + 1] << " " << (ppnSFFmO)[3*i + 2] << " " << endl;
//		}
//		fout << "CELL_TYPES " << pnSEO << endl;
//		for (int i = 0; i < pnSEO; i++)
//		{
//				fout << 5 << endl;
//		}
//		//*/
//		fout.close();
//	}
//	return 0;
}
