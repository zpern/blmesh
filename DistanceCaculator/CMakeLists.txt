cmake_minimum_required (VERSION 2.6)
project(DistanceCaculator)
#best opt
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0 -g")
#
message("building DistanceCaculator.lib")
#
add_definitions(-D _LIB)

add_library(DistanceCaculator STATIC 
./include/kdtree.h
./src/kdtree.cpp
./include/dataio.h
./src/dataio.cpp
./src/scendata.cpp
./include/scendata.h
../blmesh/include/BLVector.h
../blmesh/src/BLVector.cpp
)



target_include_directories(DistanceCaculator public ./include ./src ./geom ../blmesh/include ../blmesh/src)



