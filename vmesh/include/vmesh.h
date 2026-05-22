#pragma once
#ifndef _API_EE_H
#define _API_EE_H
#include <array>
#include <functional>
#include <map>
#include <string>
#include <vector>

#ifdef Tiger_EXPORTS
#define Tiger_API
#define Tiger_API __declspec(dllexport)
#else
#define Tiger_API
#define Tiger_API __declspec(dllimport)
#endif

namespace TiGER {
enum EntityTopology {
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
    MIXED = 16,
};

/**
 * @brife: you could call this function by anthor thread for emergency stop
 ***/
void API_Terminate();

int API_Gen_Vol_Mesh(
    /* ------------------------- 输入参数 --------------------------**/
    double *pdSNC,                  /* 曲面网格节点坐标，coord[3*i], coord[3*i+1], coord[3*i+2]为第i个节点x,y,z坐标 **/
    int nSN,                        /* 曲面网格边界节点数目 **/
    int *pnSFFm,                    /* 曲面网格单元节点编号  ,�?开�?*/
    int *pnSFTp,                    /* 曲面网格单元类型。当前仅支持三角形单�?**/
    int *pnSFPt,                    /* 曲面网格单元所在几何面编号 ,�?开�?**/
    int nSF,                        /* 曲面网格单元数目 **/
    std::map<int, int> pnFT,        /* 几何面类型： 0为远场； 1 为物面； 2为对称面 ,3为不长边界层的面,4为周期性面**/
    int nLN,                        /* 边界层层�?**/
    std::vector<int> layer_vec,     /* 边界层层数数�?**/
    double dLen,                    /* 边界层第一层厚�?**/
    std::vector<double> length_vec, /* 边界层第一层厚度数�?**/
    double dRto,                    /* 边界层厚度增长因�?**/
    bool bisostop,                  /* 各向同性停�?*/
    int nopt,                       /* 优化次数 **/
    std::function<double(std::array<double, 3>)> sizefunction, /*尺寸函数，输入为点坐标，输出为尺寸�?*/
    /* ------------------------- 输出参数 -------------------------**/
    double **ppdMNC,        /* 体网格节点坐�?**/
    int *pnMN,              /* 体网格节点数�?**/
    int **ppnMEFm,          /* 体网格单元节点编�?**/
    int **ppnMETp,          /* 体网格单元类型。当前支持单元类型：三棱柱、金字塔、四面体 **/
    int *pnME,              /* 体网格单元数�?**/
                            /* ------------------------- 边界信息 -------------------------**/
    int *num_boundary_face, /*边界面网格数�?**/
    int **boundary_mesh, /*边界面网格，注意每四个为一组，而且注意如果为三角形，最后一项为-1，id�?开�?
                          **/
    int **boundary_face, /*边界面网格对应的面id，长度为num_boundary_face **/
    /* ------------------------- 其他参数 -------------------------**/
    bool b_use_multiple_normals = false,  /*是否启用多法�?*/
    bool b_output_io_file = false,        /*是否将api的输入和输出都输出到文件系统中（仅用于DEBUG�?*/
    std::string filename = "virtualmesh", /*几何文件名，缺省为virtualmesh**/
    /* ------------------------- 控制参数 -------------------------**/
    double expan_beta = 1.2 /* 尺寸过渡因子**/

);
int API_Gen_Tetra_Mesh(

    /* ------------------------- 输入参数 --------------------------**/
    double *pdSNC, /* 曲面网格节点坐标，coord[3*i], coord[3*i+1], coord[3*i+2]为第i个节点x,y,z坐标 **/
    int nSN,       /* 曲面网格节点数目 **/
    int *pnSFFm,   /* 曲面网格单元节点编号,�?开�?**/
    int *pnSFTp,   /* 曲面网格单元类型。当前仅支持三角形单�?**/
    int *pnSFPt,   /* 曲面网格单元所在几何面编号,�?开�?**/
    int nSF,       /* 曲面网格单元数目 **/
    int nopt,      /* 优化次数 **/

    std::function<double(std::array<double, 3>)> sizefunction, /* 尺寸函数，输入为点坐标，输出为尺寸�?**/
    /* ------------------------- 输出参数 -------------------------**/
    double **ppdMNC,                      /* 体网格节点坐�?**/
    int *pnMN,                            /* 体网格节点数�?**/
    int **ppnMEFm,                        /* 体网格单元节点编�?**/
    int **ppnMETp,                        /* 体网格单元类型。当前支持单元类型：四面�?**/
    int *pnME,                            /* 体网格单元数�?**/
    bool b_output_io_file = false,        /*是否将api的输入和输出都输出到文件系统中（仅用于DEBUG�?*/
    std::string filename = "virtualmesh", /*几何文件名，缺省为virtualmesh**/

    /* ------------------------- 控制参数 -------------------------**/
    double expan_beta = 1.2 /* 尺寸过渡因子**/
);

int API_Gen_Boundary_ALM_Mesh(
    /* ------------------------- 输入参数 --------------------------**/
    double *pdSNC,                    /* 曲面网格节点坐标，coord[3*i], coord[3*i+1],coord[3*i+2]为第i个节点x,y,z坐标 **/
    int nSN,                          /* 曲面网格边界节点数目 **/
    int *pnSFFm,                      /* 曲面网格单元节点编号 **/
    int *pnSFTp,                      /* 曲面网格单元类型。当前仅支持三角形单�?�?开�?**/
    int *pnSFPt,                      /* 曲面网格单元所在几何面编号 ,�?开�?**/
    int nSF,                          /* 曲面网格单元数目 **/
    std::map<int, int> pnFT,          /* 几何面类型： 0为远场； 1 为物面； 2为对称面 ,3为不长边界层的面,4为周期性面**/
    int nLN,                          /* 边界层层�?**/
    std::vector<int> layer_vec,       /* 边界层层数数�?**/
    int max_layer_diff,               /* 相邻网格边界层层数差 **/
    double max_ratio_diff,            /* 边界层相邻层数高度尺寸比 **/
    double dLen,                      /* 边界层第一层厚�?**/
    std::vector<double> length_vec,   /* 边界层第一层厚度数�?**/
    double dRto,                      /* 边界层厚度增长因�?**/
    std::vector<double> dRto_vec,     /* 边界层厚度增长因子数�?**/
    std::vector<double> max_skewness, /* 等面积偏斜度**/
    std::vector<double> max_orth,     /* 最大非正交�?*/
    double bisostop,                  /* 各向同性停�?*/
    double clearance,                 /* 附面层留空是否留固定距离 **/
    /* ------------------------- 输出参数 -------------------------**/
    double **ppdMNC, /* 体网格节点坐�?**/
    int *pnMN,       /* 体网格节点数�?**/
    int **ppnMEFm,   /* 体网格单元节点编�?**/
    int **ppnMETp,   /* 体网格单元类型。当前支持单元类型：三棱柱，金字�?**/
    int *pnME,       /* 体网格单元数�?**/
    /* ------------------------- 边界层网格顶面层面网格参�?-------------------------**/
    double **ppdSNC0,      /* 曲面网格节点坐标，coord[3*i], coord[3*i+1],coord[3*i+2]为第i个节点x,y,z坐标 **/
    int *nSN0,             /* 曲面网格边界节点数目 **/
    int *pnSEO,            /* 顶面网格单元数目 **/
    int **ppnSFTpO,        /* 顶面网格单元类型。当前仅支持三角形单�?**/
    int **ppnSFFmO,        /* 顶面网格单元节点编号 **/
    int **l2g,             /* local点id到global点id的映�?*/
    double **point_sizing, /* 顶面网格点尺�?*/
    /* ---------------------------- 边界信息 ---------------------------**/
    int *num_boundary_face, /*边界面网格数�?**/
    int **boundary_mesh,    /*边界面网格，注意每四个为一组，而且注意如果为三角形，最后一项为-1，id�?开�?**/
    int **boundary_face,    /*边界面网格对应的面id，长度为num_boundary_face **/
    /* ------------------------- 其他参数 -------------------------**/
    bool b_have_pyramid = true,           /*是否有金字塔**/
    bool b_use_multiple_normals = false,  /*是否启用多法�?*/
    bool b_output_io_file = false,        /*是否将api的输入和输出都输出到文件系统中（仅用于DEBUG�?*/
    std::string filename = "virtualmesh", /*几何文件名，缺省为virtualmesh**/
    std::array<double, 12> per_matrix =
        std::array<double, 12>()          /* 周期性面控制矩阵,�?位为旋转矩阵 m00，m01，m02 ....，后三位为位移向量xyz**/
);

int API_Gen_Boundary_FullLayer_Mesh(
    /* ------------------------- 输入参数 --------------------------**/
    int multiple_numlayer,       /* 多法向层�?**/
    double multiple_steplength,  /* 多法向步�?**/
    std::vector<double> len_vec, /* 多法向步长数�?**/
    bool preMultiple,            /*是否对多法向角点做预处理**/
    bool fast_intersection,      /*是否启用快速相交检�?*/
    double *pdSNC,               /* 曲面网格节点坐标，coord[3*i], coord[3*i+1], coord[3*i+2]为第i个节点x,y,z坐标 **/
    int nSN,                     /* 曲面网格边界节点数目 **/
    int *pnSFFm,                 /* 曲面网格单元节点编号 **/
    int *pnSFTp,                 /* 曲面网格单元类型。当前仅支持三角形单�?**/
    int *pnSFPt,                 /* 曲面网格单元所在几何面编号 ,�?开�?**/
    int nSF,                     /* 曲面网格单元数目 **/
    std::map<int, int> pnFT,     /* 几何面类型： 0为远场； 1 为物面； 2为对称面 ,3为不长边界层的面,4为周期性面**/

    /* ------------------------- 输出参数 -------------------------**/
    double **ppdMNC, /* 体网格节点坐�?**/
    int *pnMN,       /* 体网格节点数�?**/
    int **ppnMEFm,   /* 体网格单元节点编�?**/
    int **ppnMETp,   /* 体网格单元类型。当前支持单元类型：三棱柱，金字�?**/
    int *pnME,       /* 体网格单元数�?**/

    /* ------------------------- 边界层网格顶面层面网格参�?-------------------------**/
    double **ppdSNC0,   /* 顶面网格节点坐标，coord[3*i], coord[3*i+1], coord[3*i+2]为第i个节点x,y,z坐标 **/
    int *pnSN0,         /* 顶面网格边界节点数目 **/
    int *pnSEO,         /* 顶面网格单元数目 **/
    int **ppnSFTpO,     /* 顶面网格单元类型。当前仅支持三角形单�?**/
    int **ppnSFFmO,     /* 体网格单元数�?**/
    int **l2g,          /*  ��������id��ȫ�ֵ�id��ӳ�� */
    double **ppnsizing, /* 顶面网格点尺�?*/

    /////* ------------------------- 边界信息 ------------------------- **/
    int *num_boundary_face, /* 边界面网格数�?**/
    int **boundary_mesh,    /* 边界面网格，注意每四个为一组，而且注意如果为三角形，最后一项为-1，id�?开�?**/
    int **boundary_face,    /* 边界面网格对应的面id，长度为num_boundary_face **/

    /* ------------------------- 其他参数 -------------------------**/
    bool b_output_io_file, /* 是否将api的输入和输出都输出到文件系统中（仅用于DEBUG）缺省为false **/
    std::string filename   /* 几何文件名，缺省为virtualmesh **/
);

int API_Mesh_Optimize(
    /* ------------------------- 输入输出参数 -------------------------**/
    double *pdMNC,         /* 体网格节点坐�?**/
    int num_boundary_face, /*边界面网格数�?**/
    int *boundary_mesh,    /*边界面网格，注意每四个为一组，而且注意如果为三角形，最后一项为-1，id�?开�?
                            **/
    int nMN,               /* 体网格节点数�?**/
    int *pnMEFm,           /* 体网格单元节点编�?**/
    int *pnMETp,           /* 体网格单元类型。当前支持单元类型：三棱柱，金字�?**/
    int nME                /* 体网格单元数�?**/
);

int API_Mesh_Merge(
    /* ------------------------- 输入输出参数 -------------------------**/
    double *pdMNC1, /* 体网格节点坐�?**/
    int nMN1,       /* 体网格节点数�?**/
    int *pnMEFm1,   /* 体网格单元节点编�?**/
    int *pnMETp1,   /* 体网格单元类型。当前支持单元类型：四面�?**/
    int nME1,       /* 体网格单元数�?**/

    double *pdMNC2, /* 体网格节点坐�?**/
    int nMN2,       /* 体网格节点数�?**/
    int *pnMEFm2,   /* 体网格单元节点编�?**/
    int *pnMETp2,   /* 体网格单元类型。当前支持单元类型：三棱柱，金字�?**/
    int nME2,       /* 体网格单元数�?**/
    /* ------------------------- 输出参数 -------------------------**/
    double **pdMNC3, /* 体网格节点坐�?**/
    int *nMN3,       /* 体网格节点数�?**/
    int **pnMEFm3,   /* 体网格单元节点编�?**/
    int **pnMETp3,   /* 体网格单元类型。当前支持单元类型：三棱柱，金字�?**/
    int *nME3        /* 体网格单元数�?**/
);
} // namespace TiGER
#endif
