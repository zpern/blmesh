#include "tiger_sizingfunction.h"
#include "inputdata.h"
#include "createsol.h"
#include <spdlog/spdlog.h> 
 #include <algorithm>
#include <chrono>
#include <iostream>
#include <fstream>
#include <ctime>
//#define GET_STL
#ifdef GET_STL
#include "getStlName.h"
string stlfileName;
void setStlFileName(const char* filename){
    stlfileName=filename;
}
#endif
using namespace std;
map<int,InputData*> sfObjMp;
int g_sfObjId=0;
//使用外部尺寸场
//#define USE_OUTTERSF
// 测量所有过程花费的时间
#define MEASURE_TIME
int API_Create_SurfBKG_SF(
	const char *vtkFile,
	double volmGrowthRatio,
	int *sfObjID
) {
	InputData* tmp = new InputData;
	std::string filename = string(vtkFile);

	using namespace AABBSER;
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	Eigen::VectorXd Psize;

	readfile(filename, V, F, Psize);

	//remove((tmpFilePath + "\\" + filename).c_str());
	tmp->meshForQuery = new bkMesh(V, F, Psize);
	spdlog::info("Sizing Function ok\n");
	++g_sfObjId;
	sfObjMp[g_sfObjId] = tmp;
	//cout << "sfObjMp size:" << sfObjMp.size() << endl;
	*sfObjID = g_sfObjId;
	return 0;
}
TigerAdasize_API int API_Create_Cuboid_Sources_SF(double local_coordinate_x[3], double local_coordinate_y[3], double local_coordinate_z[3], double local_coordinate_origin[3], double soValues, double volmGrowthRatio, int * sfObjID)
{
	InputData* tmp = new InputData;
	tmp->coord[0] = Vector(local_coordinate_x);
	tmp->coord[1] = Vector(local_coordinate_y);
	tmp->coord[2] = Vector(local_coordinate_z);
	tmp->source_t = InputData::CUBOID_SOLID_SOUCE;
	tmp->coord_zero = Vector(local_coordinate_origin);
	tmp->basic_source = soValues;
	tmp->volmGrowthRatio = volmGrowthRatio;

	++g_sfObjId;
	sfObjMp[g_sfObjId] = tmp;
	//cout << "sfObjMp size:" << sfObjMp.size() << endl;
	*sfObjID = g_sfObjId;
	return 0;

}
TigerAdasize_API int API_Create_Ellipsoid_Sources_SF(double local_coordinate_x[3], double local_coordinate_y[3], double local_coordinate_z[3], double local_coordinate_origin[3], double soValues, double volmGrowthRatio, int * sfObjID)
{
	InputData* tmp = new InputData;
	tmp->coord[0] = Vector(local_coordinate_x);
	tmp->coord[1] = Vector(local_coordinate_y);
	tmp->coord[2] = Vector(local_coordinate_z);
	tmp->source_t = InputData::ELLIOPSOID_SOLID_SOUCE;
	tmp->coord_zero = Vector(local_coordinate_origin);
	tmp->basic_source = soValues;
	tmp->volmGrowthRatio = volmGrowthRatio;

	++g_sfObjId;
	sfObjMp[g_sfObjId] = tmp;
	//cout << "sfObjMp size:" << sfObjMp.size() << endl;
	*sfObjID = g_sfObjId;
	return 0;
}
//*Returned value :
//*0				Completely succeed
//* 1				Partially succeed
//*	otherwise		Fail.Error info.will be defined later
int API_Query_SF_Objects(int *num, int **sfObjIDs) {
	*num = sfObjMp.size();
	int cnt = 0;
	*sfObjIDs = (int*)malloc(sizeof(int)*sfObjMp.size());
	for (auto it = sfObjMp.begin(); it != sfObjMp.end();it++) {
		*sfObjIDs[cnt++] = it->first;
	}
	return 0;
}
int API_Enable_SF_Object(int sfObjID) {
	if (sfObjMp.find(sfObjID) != sfObjMp.end()) {
		sfObjMp[sfObjID]->setEnable();
		return 0;
	}
	return 1;
}
int API_Disable_SF_Object(int sfObjID) {
	if (sfObjMp.find(sfObjID) != sfObjMp.end()) {
		sfObjMp[sfObjID]->setDisable();
		return 0;
	}
	return 1;
}
int API_Query_SF_Object_Validity(int sfObjID, int *valid) {
	*valid = sfObjMp.find(sfObjID) != sfObjMp.end();
	return 0;
}
int API_Delete_SF_Object(int sfObjID) {
	if (sfObjMp.find(sfObjID) != sfObjMp.end()) {
		cout << "del" << sfObjMp[sfObjID]<<endl;
		delete sfObjMp[sfObjID];
		sfObjMp.erase(sfObjID);
		return 0;
	}
	return 1;
}
int API_Query_SF_Object_Type(int sfObjID, int *type) {
	if (sfObjMp.find(sfObjID) != sfObjMp.end()) {
		*type = sfObjMp[sfObjID]->getType();
		return 0;
	}
	return 1;
}
double API_Sizing_Query( double x, double y, double z) {
	double ret=0;
	bool flg = true;
	
	//cout << sfObjMp.size()<<endl;
	
	for (auto obj: sfObjMp) {
		if (obj.second->isEnabled()) {
			if (obj.second->source_t == InputData::MESH_SF_SOUCE) {
				if (flg) {
					ret = obj.second->meshForQuery->getSizeSpatial(x, y, z);
					flg = false;
				}
				else ret = min(ret, obj.second->meshForQuery->getSizeSpatial(x, y, z));
			}
			else {
				double coord[3]={x,y,z};
				if (flg) {
					ret = obj.second->getSolidSouceSize(coord);
					flg = false;
				}
				else ret = min(ret, obj.second->getSolidSouceSize(coord));
			}
		}
	}
	//cout << "cal size query: " << "x:" << x << " y:" << y << " z:" << z <<" ret:"<<ret << endl;
	return ret;
}

double API_Sizing_Query(double p[3])
{
	return API_Sizing_Query(p[0], p[1], p[2]);
}
