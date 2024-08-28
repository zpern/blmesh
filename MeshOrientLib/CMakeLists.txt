cmake_minimum_required (VERSION 2.6)
project(meshOrientLib)
#best opt
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0 -g")
#
message("building meshOrientLib")
#
add_definitions(-D _LIB)
include_directories(./include ./src)
add_library(meshOrientLib STATIC 
include/mesh_orient.h 
src/mesh_orient.cpp 
include/trimesh.h 
src/trimesh.cpp 
include/vector.h 
src/vector.cpp)


