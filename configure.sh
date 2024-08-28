#!/bin/bash
cur_dir=$(cd "$(dirname "$0")"; pwd)
# step 1 build geogram 
git submodule update --init

cd ./third/geogram
if [ ! -d "build" ];
then
bash -f configure.sh Linux64-gcc
fi
cd build/Linux64-gcc-Release
if [ ! -f "Makefile" ];
then echo "please install libx11-dev libxrandr-dev libxinerama-dev  libxcursor-dev libxi-dev"
rm build -rf
fi
make -j

if [ ! -d "../../../../tiger/lib" ];
then
mkdir ../../../../tiger/lib
fi
cp lib/libgeogram.a ../../../../tiger/lib/libgeogram_r.a
cd ../Linux64-gcc-Debug
make -j
cp lib/libgeogram.a ../../../../tiger/lib/libgeogram_d.a
# step 2 build blmesh
cd $cur_dir
if [ ! -d "build" ];
then mkdir build
     mkdir build/Debug
     mkdir build/Release
fi
cd build/Debug
cmake ../.. -DCMAKE_BUILD_TYPE=Debug
make -j
cd ../Release
cmake ../.. -DCMAKE_BUILD_TYPE=Release
make -j
cp vmesh/test /tmp/test
