/*
 * main.cpp
 * Created on: Sep 23, 2015
 * Author: Zhoufang Xiao
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iomanip>
#include <spdlog/spdlog.h> 
#include "BLMesh.h"
#include "BLMesh_define.h"
#include "blmeshapi.h"
#include "../common/singleton_terminate.h"

#define _CONFIG_FILE

 // #define _2D_BLMESH
#define _GEN_BLMESH

ConfigArgc cf;

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
    std::array<double, 12> per_matrix) {

    VM ret;
    double start_time = clock();
    int argc = 1;
    char* argv[2];
#else
void main(int argc, char* argv[]) {
#endif

    char bdryfile[256], bdryfileo[256], blmfile[256], blmfile2[256], gridfile[256], fieldfile[256];

#if 1
#ifdef _CONFIG_FILE
    char finconfig[1024];

#ifdef APIFUNC
    cf.boxfc = blconfig.box;
    cf.step_len = blconfig.len;
    cf.layer_num = blconfig.n;
    cf.isotropic_stop = isostop;
    cf.ratio = blconfig.Ro;
    cf.layer_ratio = 0;
    cf.ratio1 = blconfig.Ro;
    cf.ratio2 = blconfig.Ro;
    cf.adjacentfc = blconfig.adjacent;
    cf.smooth_attempt = number_opt;
    strcpy(cf.filenam, "virtualmesh");
    cf.layer_ratio = 0;

#ifdef _DEBUG
    cf.smooth_attempt = 0;
#endif

    bool use_per = blconfig.per.size();
    if (use_per) {
        for (auto i : blconfig.per) cf.perfc.push_back(i);
        for (int i = 0; i < 9; i++)
            cf.rotate_matrix(i / 3, i % 3) = per_matrix[i];
        for (int i = 0; i < 3; i++)
            cf.shift_vec[i] = per_matrix[i + 9];
        cf.reverse_matrix = cf.rotate_matrix.inverse();
        genboundarymesh = true;
        genoutmesh = false;
    }
#else
    parsecommand(argc, argv, cf);
#endif

    PrintConfig(&cf);
    memset(finconfig, 0, sizeof(finconfig));
    strcpy(finconfig, cf.filenam);
    strcat(finconfig, "_o");
    argv[1] = finconfig;

#endif
#endif

    memset(bdryfile, 0, sizeof(bdryfile));
    strcpy(bdryfile, argv[1]);
    memset(bdryfileo, 0, sizeof(bdryfileo));
    strcpy(bdryfileo, cf.filenam);

#ifdef _2D_BLMESH
    strcat(bdryfile, ".fr2");
    strcat(bdryfileo, "_o.fr2");
#else
    strcat(bdryfile, ".pls");
    strcat(bdryfileo, ".bdry");
#endif

    memset(blmfile, 0, sizeof(blmfile));
    strcpy(blmfile, argv[1]);
    strcat(blmfile, ".vtk");

    memset(blmfile2, 0, sizeof(blmfile2));
    strncpy(blmfile2, argv[1], strlen(argv[1]) - 2);
    strcat(blmfile2, ".vtk");

    memset(gridfile, 0, sizeof(gridfile));
    strcpy(gridfile, argv[1]);

#ifdef _2D_BLMESH
    strcat(gridfile, ".pl2");
#else
    strcat(gridfile, ".pl3");
#endif

    memset(fieldfile, 0, sizeof(fieldfile));
    strcpy(fieldfile, argv[1]);
    strcat(fieldfile, "_pl.vtk");

#ifdef _CONFIG_FILE
#ifdef _2D_BLMESH
    BLMesh* blmesh = new BLMesh(cf.layer_num, 0);
#else
    std::shared_ptr<BLMesh> blmesh(new BLMesh(cf.layer_num, 0, cf.isotropic_stop));
    blmesh->cf = cf;
    blmesh->expan_ratio = expan_beta;
#endif
#else
#ifdef _2D_BLMESH
    BLMesh* blmesh = new BLMesh(16, 0);
#else
    BLMesh* blmesh = new BLMesh(10, 0);
#endif
#endif

    blmesh->SetSizingFunction(sizefunction);
    double gen_mesh_start_time = clock();

#ifdef _2D_BLMESH
    blmesh->InitBLMesh(BLMType::blm2d);
#else
    blmesh->InitBLMesh();
#endif

#ifdef APIFUNC
    blmesh->SetBoundary(sfile);
#else
    blmesh->ReadBoundary(bdryfile, 1);
#endif

    blmesh->CalZeroNorm();
    blmesh->generate_pyramid = havepyramid;
    blmesh->GenerateBLMesh();

    if (checkterminate()) {
        blmesh->FreeMemoryInFrontAndNode();
        return ret;
    }

    blmesh->RemoveOverlapNodeAndElement();

#ifdef APIFUNC
    if (genboundarymesh) {
        blmesh->GenTopMesh(ret, genoutmesh);
    }
#endif

    double gen_bdymesh_end_time = clock();

#ifdef APIFUNC
    if (genoutmesh) {
        spdlog::info("opt times={}", cf.smooth_attempt);
        blmesh->GenOuterMesh(cf.smooth_attempt);
        if (checkterminate()) {
            blmesh->FreeMemoryInFrontAndNode();
            return ret;
        }
    }
    else {
        blmesh->FreeMemoryInFrontAndNode();
    }

    auto o_to_n = blmesh->RemvNodElm();
    for (int i = 0; i < ret.nSN0; i++) {
        ret.l2g[i] = o_to_n[ret.l2g[i]];
    }

    if (outputIO) {
        blmesh->SaveBLMesh(blmfile2);
    }
    blmesh->SaveBdry(ret);
    blmesh->SaveBLMeshAndFreeMemory(ret);

    double gen_mesh_end_time = clock();
    spdlog::info("Average number of layers{}", blmesh->GetAvergeLayer());
    double end_time = clock();

    spdlog::info("              Mesh infomation    ");
    spdlog::info("+================================================+");
    spdlog::info("|             Item                |number of cell|");
    spdlog::info("-------------------------------------------------");
    spdlog::info("|            Pyramid              |{}", blmesh->m_nPyramid);
    spdlog::info("|            Tetrahedron          |{}", blmesh->m_nTet);
    spdlog::info("|         Triangular prism        |{}", blmesh->m_nPrism);
    spdlog::info("|            Triangle             |{}", blmesh->m_nTri);
    spdlog::info("|          Quadrilateral          |{}", blmesh->m_nQuad);
    spdlog::info("+================================================+");

    spdlog::info("              Timing statistics");
    spdlog::info("+================================================+");
    spdlog::info("|             Item                |Time(second)   |");
    spdlog::info("-------------------------------------------------");
    spdlog::info("|Generate Boundary Layer Mesh Cost|{}", (gen_bdymesh_end_time - gen_mesh_start_time) / CLOCKS_PER_SEC);
    spdlog::info("|Total Check Insertion Cost       |{}", blmesh->check_prism_time / CLOCKS_PER_SEC);
    spdlog::info("|Total Normal Smoothing Cost      |{}", blmesh->smooth_time / CLOCKS_PER_SEC);
    spdlog::info("|Other process In bdymesh Cost    |{}", (gen_bdymesh_end_time - gen_mesh_start_time - blmesh->smooth_time - blmesh->check_prism_time) / CLOCKS_PER_SEC);
    spdlog::info("|Generate Tetrahedral Mesh Cost   |{}", (gen_mesh_end_time - gen_bdymesh_end_time) / CLOCKS_PER_SEC);
    spdlog::info("|Generate Mesh Cost        (no IO)|{}", (gen_mesh_end_time - gen_mesh_start_time) / CLOCKS_PER_SEC);
    spdlog::info("|Total Cost        (including IO) |{}", (end_time - start_time) / CLOCKS_PER_SEC);
    spdlog::info("+================================================+");

    return ret;
#endif
}
