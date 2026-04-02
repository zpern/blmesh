#ifndef _BLMesh_Define_
#define _BLMesh_Define_

#include "BLVector.h"
#include <Eigen/Dense>
#include <cstring>
#include <limits>
#include <map>
#include <memory>
#include <spdlog/spdlog.h>
#include <vector>
const int MAX_NCONN = 20;
const int MAX_POINT_IN_ELEMENT = 6;
#define MAX_FRONT_NODES      3
#define NUM_ADD_RATIO        0.25
#define MAX_NORMAL_COMPONENT 3
#define COORD_DIM2           2
#define COORD_DIM3           3
#define DIM2                 2
#define DIM3                 3
#define rad(a)               (PI * ((double)a / 180))
#ifndef RETURN_INSERTION
#define RETURN_INSERTION   2
#define RETURN_NONMANIFOLD 3
#endif
#ifndef PI
const double PI = 3.141592653589793;
#endif
#ifdef NULL
#undef NULL
#define NULL nullptr
#endif // !NULL

enum BLEntityTopology {
    NODE = 5,
    LINE = 6,
    POLYGON = 7,
    TRIANGLE = 8,
    QUADRILATERAL = 9,
    POLYHEDRON = 10,
    TETRAHEDRON = 11,
    HEXAHEDRON = 12,
    PRISM = 13,
    PYRAMID = 14,
    SEPTAHEDRON = 15,
    MIXED = 16
};

struct Elem {
    Elem()
        : topo(BLEntityTopology::TETRAHEDRON)
        , nconn(4)
    {}
    int conn[MAX_POINT_IN_ELEMENT];
    short nconn;

    int iLayer;
    short igom;

    void *pointer;
    BLEntityTopology topo;
};
// class AbstractElem {
// public:
//	int* conn;
//	//virtual int nConn() = 0;
// };
// class SurfaceElem:public AbstractElem {
// public:
//	short igom;
//	void * pointer;
// };
// class TriElem :public SurfaceElem {
// public:
//	TriElem() { conn = new int[3]; }
// };
// class QuadElem :public SurfaceElem {
// public:
//	QuadElem() { conn = new int[4]; }
// };
// class VolumnElem :public AbstractElem {};
// class TetraElem :public VolumnElem {
//	TetraElem() { conn = new int[4]; }
// };
// class PyramidElem :public VolumnElem {
//	PyramidElem() { conn = new int[5]; }
// };
// class PrismElem :public VolumnElem {
//	PrismElem() { conn = new int[6]; }
// };
// class Elem {
// public:
//	AbstractElem& operator [](int index) {
//		return element_.at(index);
//	}
//	void AddTriElem(const int* conn,int igeom) {
//		TriElem tri=
//	}
// protected:
//	std::vector<AbstractElem> element_;
// };
#include <ctime>
#define DEFAULTCFNAME "TRAN_CONFIG"
#define FILENAMESIZE  120

static char config[FILENAMESIZE];

#define ERRORCONFIG     -1
#define DEFAULTCONFIG   1
#define INPUTCONFIG     0
#define NODEFAULTCONFIG -2

#define N0_LICENSE_FILE        -1
#define INCORRECT_LICENSE_FILE 0
#define CORRECT_LICENSE_FILE   1
#define EXPIRED_LICENSE        2
#define INCORRECT_SYSTEM_TIME  3
#include <array>
struct SearchTriangle {
    SearchTriangle()
        : time_stamp(std::numeric_limits<unsigned int>::max()) {};
    SearchTriangle(const std::array<int, 3> &i)
        : time_stamp(std::numeric_limits<unsigned int>::max())
        , indexs(i) {};
    std::array<int, 3> indexs;
    unsigned int time_stamp;
    int operator[](int i) { return indexs[i]; }
    int operator[](short i) { return indexs[i]; }
};
struct MBLNode {
    MBLNode()
    {
        isymfc.reserve(2);
        bsysm = false;
        bdrypt = false;
        bfarfield = false;
        badjacent = false;
        bwall = false;
    }
    double coord[3];
    bool bsysm;              // symmetry plane
    bool bdrypt;             // boundary point
    bool bfarfield;          // farfield
    bool badjacent;
    bool bwall;              // wall
    std::vector<int> isymfc; // symmetry face index

    void *pointer;
    int reserved;
    // double distance;//for test only yhf20190817
};

enum SymmetryPlane {
    SYSMMETRY_NON = 0,
    SYSMMETRY_X = 1,
    SYSMMETRY_Y = 2,
    SYSMMETRY_Z = 3,
    SYSMMETRY_DEFINE = 4
};

enum NormalGenMethod {
    GEOMETRICAL = 0,
    PHYSICAL = 1,
    GEOMETRICAL_COPY = 2, // first layer points geometrical, others copy
    PHYSICAL_COPY = 3
};

enum BoundaryType {
    farfield = 0,
    wall = 1,
    symmetry = 2,
    match = 3,
    per = 4,
    adjacent = 5
};

enum PotentialType {
    BPOTENTIAL = 1,
    BFLUX = 2
};

// define error type
#define BLM_SUCCESS     0
#define OUT_OF_MEMORY   -1
#define OPEN_FILE_FAILD -2
#define BLM_FAILED      -3

struct ConfigArgc {
    ConfigArgc() {}
    char filename[1024]; // the project name
    int layer_num;
    std::vector<int> layer_num_vec;
    int surface_normal;
    double isotropic_stop;
    double step_len;
    std::vector<double> step_len_vec;
    double ratio;
    int layer_ratio;
    double ratio1;
    double ratio2;
    int smooth_attempt;

    double min_volumn_eps;
    std::vector<double> max_skewness;
    std::vector<double> max_orth;
    int max_layer_diff;
    double max_ratio_diff;
    bool iscompresslen;
    double clearance;
    /*do not generate bdylayermesh */
    std::vector<int> matchfc;

    /*far field face*/
    std::vector<int> boxfc;

    /*symmetry face*/
    std::vector<int> symfc;

    /*perodioc face*/
    std::vector<int> perfc;
    Eigen::Matrix3d rotate_matrix;
    Eigen::Matrix3d reverse_matrix;
    Eigen::Vector3d shift_vec;

    /*adjacent grid face*/
    std::vector<int> adjacentfc;
};
static int PrintConfig(const struct ConfigArgc *cf)
{
    char filename[256]; // the project name
    int layer_num;
    double step_len;
    double ratio;
    std::vector<int> boxfc;
    spdlog::info("**************input parameters******************\n");
    spdlog::info("file name:          {}", cf->filename);
    spdlog::info("total layer number: {}", cf->layer_num);
    spdlog::info("initial step length:{}", cf->step_len);
    spdlog::info("growth rate:        {}", cf->ratio2);
    spdlog::info("isotropic stop      {}", cf->isotropic_stop);
    spdlog::info("smooth attempt times{}", cf->smooth_attempt);

    return 0;
}
static void SetDefaultConfig(ConfigArgc &cf)
{
    strcpy(cf.filename, "virtualmesh");
    cf.layer_num = 10;
    cf.step_len = 0.1;
    cf.ratio = 1.25;
    cf.surface_normal = 0;
    cf.smooth_attempt = 10;
    cf.ratio1 = 1.2;
    cf.ratio2 = 1.2;
    cf.layer_ratio = 5;
    cf.max_skewness = std::vector<double>{0.98, 0.98};
    cf.max_orth = std::vector<double>{0.98, 0.98};
}

static int GetConfigName(int argc, char **argv)
{
    // spdlog::info("Note : Use the default config.\n");
    strcpy(config, DEFAULTCFNAME);
    return DEFAULTCONFIG;
}
static int parsecommand(int argc, char **argv, struct ConfigArgc &cf)
{
    char word[100];                 /* variable name from config file */
    int value = 0;                  /* variable value from config file */
    double fvalue = 0.0, fv1 = 0.0; /* variable for floats in config file */
    char str[100];
    char fstr[100];
    char temp[100];

    /* temp stuff for line reading */
    int nbytes = 1024;
    char line[1024];
    int numassigned;
    cf.isotropic_stop = 1;

    GetConfigName(argc, argv);
    FILE *fp = fopen(config, "r");
    if (!fp) {
        spdlog::info("Config File does not exist.\n");
        throw(std::string("Config File does not exist."));
    }
    SetDefaultConfig(cf);

    while (fgets(line, nbytes, fp) != NULL) {
        /* attempt to fetch a variable name and value from the config file */
        numassigned = sscanf(line, "%s %s", word, temp);
        fvalue = atof(fstr);
        /* check if this is a comment */
        if (word[0] == '#' || word[0] == '\n' || numassigned < 2) {
            continue;
        }

        if (!strcmp(word, "filename")) {
            sscanf(line, "%s %s", word, str);
            strcpy(cf.filename, str);
        }
        if (!strcmp(word, "layer_num")) {
            sscanf(line, "%s %d", word, &value);
            cf.layer_num = value;
        }
        if (!strcmp(word, "box_face")) {
            int iface;
            sscanf(line, "%s %d", word, &iface);
            cf.boxfc.push_back(iface);
        }
        if (!strcmp(word, "match_face")) {
            int iface;
            sscanf(line, "%s %d", word, &iface);
            cf.matchfc.push_back(iface);
        }
        if (!strcmp(word, "sym_face") || !strcmp(word, "symm_face")) {
            int iface;
            sscanf(line, "%s %d", word, &iface);
            cf.symfc.push_back(iface);
        }
        if (!strcmp(word, "match_face")) {
            int iface;
            sscanf(line, "%s %d", word, &iface);
            cf.matchfc.push_back(iface);
        }
        if (!strcmp(word, "step_len")) {
            sscanf(line, "%s %lf", word, &fvalue);
            cf.step_len = fvalue;
        }
        if (!strcmp(word, "ratio")) {
            sscanf(line, "%s %d %lf %lf", word, &value, &fvalue, &fv1);
            cf.ratio1 = fvalue;
            cf.ratio2 = fv1;
            cf.layer_ratio = value;
        }
        if (!strcmp(word, "smooth_attempt")) {
            sscanf(line, "%s %d", word, &value);
            cf.smooth_attempt = value;
        }
        if (!strcmp(word, "surface_normal")) {
            int iface;
            sscanf(line, "%s %d", word, &iface);
            cf.surface_normal = iface;
        }

        if (!strcmp(word, "isotropic_stop")) {
            double iface;
            sscanf(line, "%s %lf", word, &iface);
            cf.isotropic_stop = iface;
        }

        if (!strcmp(word, "per_face")) {
            int iface1, iface2;
            sscanf(line, "%s %d %d", word, &iface1, &iface2);
            cf.perfc.push_back(iface1);
            cf.perfc.push_back(iface2);
            for (int j = 0; j < 3; j++) {
                fgets(line, nbytes, fp);
                sscanf(line,
                       "%lf %lf %lf\n",
                       &cf.rotate_matrix(j, 0),
                       &cf.rotate_matrix(j, 1),
                       &cf.rotate_matrix(j, 2));
            }
            cf.reverse_matrix = cf.rotate_matrix.inverse();
            fgets(line, nbytes, fp);
            sscanf(line, "%lf %lf %lf\n", &cf.shift_vec(0), &cf.shift_vec(1), &cf.shift_vec(2));
        }
    }
    fclose(fp);
    return 0;
}

static bool TimeCheck()
{
    bool ret = true;
    int year = 2020, month = 10, day = 1;
    time_t rawtime;
    struct tm *d;
    ;
    time(&rawtime);
    d = localtime(&rawtime);
    // spdlog::info("year:%d month: %d day:%d\n", d->tm_year+1900, d->tm_mon+1, d->tm_mday);

    if (d->tm_year + 1900 < year) {
        ret = true;
    } else if (d->tm_mon + 1 < month) {
        ret = true;
    } else if ((d->tm_year + 1900 == year) && (d->tm_mon + 1 == month) && d->tm_mday <= day) {
        ret = true;
    } else {
        ret = false;
    }

    if (!ret) {
        // spdlog::info("Unknown error: please contact the software provider of
        // chenjj@zju.edu.cn!\n"); exit(0);
        return ret;
    }
    return ret;
}
extern ConfigArgc cf;

/*
 *@author yhf
 *@brife 求解三角形外心
 */
static BLVector solveCenterPointOfCircle(BLVector A, BLVector B, BLVector C)
{
    BLVector pd[3] = {A, B, C};
    double a1, b1, c1, d1;
    double a2, b2, c2, d2;
    double a3, b3, c3, d3;

    double x1 = pd[0].x, y1 = pd[0].y, z1 = pd[0].z;
    double x2 = pd[1].x, y2 = pd[1].y, z2 = pd[1].z;
    double x3 = pd[2].x, y3 = pd[2].y, z3 = pd[2].z;

    a1 = (y1 * z2 - y2 * z1 - y1 * z3 + y3 * z1 + y2 * z3 - y3 * z2);
    b1 = -(x1 * z2 - x2 * z1 - x1 * z3 + x3 * z1 + x2 * z3 - x3 * z2);
    c1 = (x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2);
    d1 = -(x1 * y2 * z3 - x1 * y3 * z2 - x2 * y1 * z3 + x2 * y3 * z1 + x3 * y1 * z2 - x3 * y2 * z1);

    a2 = 2 * (x2 - x1);
    b2 = 2 * (y2 - y1);
    c2 = 2 * (z2 - z1);
    d2 = x1 * x1 + y1 * y1 + z1 * z1 - x2 * x2 - y2 * y2 - z2 * z2;

    a3 = 2 * (x3 - x1);
    b3 = 2 * (y3 - y1);
    c3 = 2 * (z3 - z1);
    d3 = x1 * x1 + y1 * y1 + z1 * z1 - x3 * x3 - y3 * y3 - z3 * z3;
    BLVector centerpoint;
    centerpoint[0] =
        -(b1 * c2 * d3 - b1 * c3 * d2 - b2 * c1 * d3 + b2 * c3 * d1 + b3 * c1 * d2 - b3 * c2 * d1) /
        (a1 * b2 * c3 - a1 * b3 * c2 - a2 * b1 * c3 + a2 * b3 * c1 + a3 * b1 * c2 - a3 * b2 * c1);
    centerpoint[1] =
        (a1 * c2 * d3 - a1 * c3 * d2 - a2 * c1 * d3 + a2 * c3 * d1 + a3 * c1 * d2 - a3 * c2 * d1) /
        (a1 * b2 * c3 - a1 * b3 * c2 - a2 * b1 * c3 + a2 * b3 * c1 + a3 * b1 * c2 - a3 * b2 * c1);
    centerpoint[2] =
        -(a1 * b2 * d3 - a1 * b3 * d2 - a2 * b1 * d3 + a2 * b3 * d1 + a3 * b1 * d2 - a3 * b2 * d1) /
        (a1 * b2 * c3 - a1 * b3 * c2 - a2 * b1 * c3 + a2 * b3 * c1 + a3 * b1 * c2 - a3 * b2 * c1);

    // centerpoint[3] = sqrt(pow(pd[0].x - centerpoint[0], 2) + pow(pd[0].y - centerpoint[1], 2) +
    // pow(pd[0].z - centerpoint[2], 2));

    return centerpoint.normalized();
}

#endif
