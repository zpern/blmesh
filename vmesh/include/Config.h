#pragma once
#include <map>
#include <vector>
#include <string>
#include <cstring>
namespace tiger {


class Config
{
public:
    Config(void);
    ~Config(void);

    int SetDefaultConfig();

    int ReadConfigFile(const char* filename);

    int WriteConfigFile(const int iLv = 0,std::string filename="");

    int CopyConfig(const Config *cf);

private:
    int WriteConfigGeomPart(FILE *fp);

    int WriteConfigMeshSizePart(FILE *fp);

    int WriteConfigAdapPart(FILE *fp);

    int WriteConfigLayerMeshPart(FILE *fp);

    int WriteConfigMeshOptPart(FILE *fp);

    /* public data */
public:
    std::string sConfigName;
    int iWrtLvl; // 0: all menbers; 1: geometry members; 2: mesh members
    double dElemNum;    // total element number

    /* basic members */
    //std::string filename;  	// the project name

    /* geometry members */
    int iGeoType;
    std::string sKwGeoType;
    std::string sCmtGeoType;

    std::string sGeoFileName;	// the geometry file name
    std::string sKwGeoFileName;	// the key word
    std::string sCmtGeoFileName; // the comment

    double dSTolerance;				//set sTolerance, default value 1.0e-8
    std::string sKwSTolerance;
    std::string sCmtSTolerance;

    double dGTolerance;				//set gTolerance, default value 1.0e-8
    std::string sKwGTolerance;
    std::string sCmtGTolerance;

    bool   bDelSmallEdges;			//set del_small_edges, default value false
    std::string sKwDelSmallEdges;
    std::string sCmtDelSmallEdges;

    std::string sKwisotropic_stop;
    std::string sCmtisotropic_stop;
    int isotropic_stop;
    /* parameters for mesh size */
    double dGlblSize;				//the global mesh size
    std::string sKwGlblSize;
    std::string sCmtGlblSize;

    double dMinSize;				//the minimize mesh size
    std::string sKwMinSize;
    std::string sCmtMinSize;

    double dExpandRatio;
    std::string sKwExpandRatio;
    std::string sCmtExpandRatio;

    std::map<int, double> mapPtSize;
    std::string sKwPtSize;
    std::string sCmtPtSize;

    std::map<int, double> mapCvSize;
    std::string sKwCvSize;
    std::string sCmtCvSize;

    std::map<int, double> mapFcSize;
    std::string sKwFcSize;
    std::string sCmtFcSize;

    /* parameters for adaptive */
    int iAdaptive;
    std::string sKwAdaptive;
    std::string sCmtAdaptive;

    int  iProximityNum;
    std::string sKwProximityNum;
    std::string sCmtProximityNum;

    double dCurvatureAngle;
    std::string sKwCurvatureAngle;
    std::string sCmtCurvatureAngle;

    int iuNum;
    std::string sKwuNum;
    std::string sCmtuNum;

    int ivNum;
    std::string sKwvNum;
    std::string sCmtvNum;

    std::map<int, int> mapFcBku;
    std::string sKwFcBku;
    std::string sCmtFcBku;

    std::map<int, int> mapFcBkv;
    std::string sKwFcBkv;
    std::string sCmtFcBkv;

    /* layer mesh */
    int nLayerNum;
    std::string sKwLayerNum;
    std::string sCmtLayerNum;

	std::string sKwSymmFc;
	std::string sCmtSymmFc;
    std::vector<int> vecBoxFc;


	std::vector<int> vecSymmFc;
	std::vector<int> vecMatchFc;
    std::string sKwBoxFc;
    std::string sCmtBoxFc;

	std::vector<int> vecAdjacentFc;

    double dStepLen;
    std::string sKwStepLen;
    std::string sCmtStepLen;

    double dStepLenRatio;
    std::string sKwStepLenRatio;
    std::string sCmtStepLenRatio;

    double dStepLenInitRatio;
    //std::string sKwStepLenInitRatio;
    //std::string sCmtStepLenInitRatio;

    int nInitLayerNum;
    //std::string sKwInitLayerNum;
    //std::string sCmtInitLayerNum;

    /* mesh improvement */
    int nSmoothLoop;
    std::string sKwSmoothLoop;
    std::string sCmtSmoothLoop;


private:

};
}
