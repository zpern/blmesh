#ifndef ISO3D_COMMON_H
#define ISO3D_COMMON_H

#include "iso3d.h"

void printTimingProfile();

void printTimingProfile(bool doQualityImproving);
int dtiso3d_improving(DTIso3D& iso3d, const char chExampleName[],
                      int nSmoothAttempts, bool embedded = false, bool output = true);
#ifdef _MULTITHREAD
int dtiso3d_meshing(DTIso3D& iso3d, const char chExampleName[],
                    int nThreads = 1, bool embedded = false, bool output = true);
#else
int dtiso3d_meshing(DTIso3D& iso3d, const char chExampleName[],
                    bool embedded = false, bool output = true);
#endif // _MULTITHREAD
int dtiso3d_readVlmMesh(DTIso3D& iso3d, const char chExampleName[]);

const int MAX_FILE_LEN = 1024;

// extern bool g_bMeshingEnabled;  /* 网格生成功能是否被启动 -M */
// extern bool g_bQuaImprEnabled;	/* 网格优化功能是否被启动 -Q */
// extern int g_nThreads; // Number of threads used to mesh.

#endif // ISO3D_COMMON_H
