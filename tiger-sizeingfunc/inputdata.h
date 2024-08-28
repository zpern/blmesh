#ifndef INPUTDATA_H
#define INPUTDATA_H
#include <string>
#include <vector>
#include <cmath>
#include <map>
#include "createsol.h"
using namespace std;
struct vertex {
    float x,y,z;
    int idx;
};
using namespace AABBSER;
class InputData
{
public:
	InputData();
	~InputData() {
		delete meshForQuery;
	}
	double volmGrowthRatio;
	AABBSER::Vector coord[3];
	AABBSER::Vector coord_zero;
	double basic_source;
    InputData(string g,string s="",string t="");
	enum SOUCETYPE {
		MESH_SF_SOUCE=0,
		CUBOID_SOLID_SOUCE=1,
		ELLIOPSOID_SOLID_SOUCE=2

	}source_t;
    string gm3file;
    string Tran_ConfigFile;
    string stlplsfile;

	double getSolidSouceSize(double p[3]);
    int cvNums; vector<int> cvPtNums; vector<double> cvPts;
    int fcNums; vector<int> fcPtNums; vector<double> fcPts;
    vector<int> lpCvNums; vector<int> lpCvs;
    int bndPtNum; vector<double> bndPts;
    int bndFctNum; vector<int> bndFcts;vector<int> fctFaces;
    vector<int> spSettingNum;
    vector<double> spGlSettings; vector<double> spPtSettings;
    vector<double> spCvSettings; vector<double> spFcSettings;

    int ptNumsAll;
    int cvPtNumsAll;
    int lpCvNumsAll;
    int fcPtsNumsAll;

    //从文件中读取相应数据到API的输入中
    bool readGm3();
    bool readStlPls();
    bool readTranConfig();
    //从配置数据中读取配置信息
    bool readTF();
    //将API输入数据转化成adasize所需格式文件
    bool writeGm3(string fileName);
    bool writeStlPls(string fileName);
    bool writeTranConfig(string fileName);
    //从cf结构体中直接导出配置文件
    bool writeTranConfigFromCF(string fileName);

    bool writeTigerGm3(string fileName);
    bool writeTigerStlPls(string fileName);

    bool writeGm3String(string str);
	string getFileName(const string& s);
	string getPathName(const string& s);

    vector<vertex> gm3V;
    struct curve {
        int idx,ncp,type;
        vector<int> v;
    };
    vector<curve> gm3C;
    struct surface {
        int idx,type,nsp1,nsp2;
        vector<vector<int>> v;
    };
    vector<surface> gm3S;
    struct loop {
        int idx,sfidx,nlc;
        vector<int> lc;
    };
    vector<loop> gm3L;

    vector<vertex> TigerV;
    vector<curve> TigerC;
    vector<surface> TigerS;
    vector<loop> TigerL;

    struct element{
        int idx;
        int sfidx;
        int a,b,c;
    };
    vector<element> stlplsE;
    vector<vertex> stlplsV;

    struct ConfigArgc{
        char filenam[1024];		// the project name
        int step=-1;
        // hmin
        double min_space;
        // hmax
        double global_space;
        // theta
        double curvature_angle;
        // beta
        double expand_ratio;
        // 间隔 nsegPt,nsegCv,nsegFc
        double proximity_num;
        // 点的hmax
        std::map<int, double> ptsize;
        // 点的hmax
        std::map<int, double> lnsize;
        // 点的hmax
        std::map<int, double> fcsize;
        ConfigArgc() {
            step=-1;
            min_space=-1;
            global_space=-1;
            curvature_angle=-1;
            expand_ratio=-1;
            proximity_num=-1;
        }
    };
    ConfigArgc cf;


//    int* spSettingNum;
//    double* spGlSettings; double* spPtSettings;
//    double* spCvSettings; double* spFcSettings;
//    double volmGrowthRatio;
    struct TigerConfig{
        struct globalSettings{
            double hmax,hmin,beta,theta,nsegPt,nsegCv,nsegFc;
            //hmax:单元允许的空间最大值
            //hmin:单元允许的空间最小值
            //beta:相邻单元间增长比率
            //theta:用户参数 曲率
            //nsegPt:临近点间隔
            //nsegCv:临近曲线间隔
            //nsegFc:临近面间隔
        }GlSettings;
        struct pointSettings{
            //第一个参数为点id
            double id,hmax,hmin,theta,nsegPt,nsegCv,nsegFc;
        };
        vector<pointSettings> ptSettings;
        struct curveSettings{
            //第一个参数为曲线id
            double id,hmax,hmin,theta,nsegCv,nsegFc;
        };
        vector<curveSettings> cvSettings;
        struct faceSettings{
            //第一个参数为面id
            double id,hmax,hmin,beta,theta,nsegFc;
        };
        vector<faceSettings> fcSettings;
        double volumeGrowthRatio=0;
        int pntSetNums=0;
        int cvSetNums=0;
        int fcSetNums=0;
    }tf;
    string Gm3AndStlFileName;
    string TranConfigFileName;
	
	bkMesh* meshForQuery=nullptr;
	inline bool isEnabled() {
		return enabled;
	}
	inline void setEnable() {
		enabled = true;
	}
	inline void setDisable() {
		enabled = false;
	}
	enum SF_OBJ_TYPE
	{
		SF_OBJ_GLOB_VALUE = 0,
		SF_OBJ_CURVE_MESH,
		SF_OBJ_SURF_MESH,
		SF_OBJ_VOLM_MESH,
		SF_OBJ_SIMPLEX_SOURCES
	};
	inline void setType(SF_OBJ_TYPE t) {
		type = t;
	}
	inline int getType() {
		return type;
	}
private:

	bool enabled = true;
	SF_OBJ_TYPE type;
};
//全局变量，记录所有sfObj
extern map<int,InputData*> sfObjMp;

#endif // INPUTDATA_H
