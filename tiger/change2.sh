#!/bin/bash

# 定位到包含 .cpp 文件的目录
cd /home/gitlab-runner/builds/aVZ5M1Fu/0/tiger/tiger1.5/src/blmesh/MNormal/src 

# 循环遍历目录下的所有 .cpp 文件
for file in *.cpp; do
    # 使用 sed 命令替换文件内容
    sed -i 's@#include "abstractoptimizer.h"@#include "AbstractOptimizer.h"@g' "$file"
    sed -i 's@#include "./virtualspheremesh.h"@#include "./VirtualSpheremesh.h"@g' "$file"
    sed -i 's@#include "./VirtualSpheremesh.h"@#include "./VirtualSphereMesh.h"@g' "$file"
    sed -i 's@#include "../include/pointoptimizer.h"@#include "../include/PointOptimizer.h"@g' "$file"
    sed -i 's@#include "../include/topologyoptimizer.h"@#include "../include/TopologyOptimizer.h"@g' "$file"
    sed -i 's@#include "../include/virtualspheremeshhasher.h"@#include "../include/VirtualSphereMeshHasher.h"@g' "$file"
    sed -i 's@#include "../include/meshevaluation.h"@#include "../include/MeshEvaluation.h"@g' "$file"
    sed -i 's@#include "BLVector.h"@#include "blvector.h"@g' "$file"
    sed -i 's@#include "mnormalmesh.h"@#include "MNormalMesh.h"@g' "$file"
    sed -i 's@#include "chamferbehavior.h"@#include "ChamferBehavior.h"@g' "$file"
done

cd /home/gitlab-runner/builds/aVZ5M1Fu/0/tiger/tiger1.5/src/blmesh/MNormal/include 
# 循环遍历目录下的所有 .cpp 文件
for file in *.h; do
    # 使用 sed 命令替换文件内容
    sed -i 's@#include "abstractoptimizer.h"@#include "AbstractOptimizer.h"@g' "$file"
    sed -i 's@#include "./virtualspheremesh.h"@#include "./VirtualSpheremesh.h"@g' "$file"
    sed -i 's@#include "./VirtualSpheremesh.h"@#include "./VirtualSphereMesh.h"@g' "$file"
    sed -i 's@#include "../include/pointoptimizer.h"@#include "../include/PointOptimizer.h"@g' "$file"
    sed -i 's@#include "../include/topologyoptimizer.h"@#include "../include/TopologyOptimizer.h"@g' "$file"
    sed -i 's@#include "../include/virtualspheremeshhasher.h"@#include "../include/VirtualSphereMeshHasher.h"@g' "$file"
    sed -i 's@#include "../include/meshevaluation.h"@#include "../include/MeshEvaluation.h"@g' "$file"
    sed -i 's@#include "BLVector.h"@#include "blvector.h"@g' "$file"
    sed -i 's@#include "mnormalmesh.h"@#include "MNormalMesh.h"@g' "$file"
    sed -i 's@#include "chamferbehavior.h"@#include "ChamferBehavior.h"@g' "$file"
done

echo "All replacements are done."

