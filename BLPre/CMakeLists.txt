cmake_minimum_required (VERSION 3.5)
project(blpre)
#best opt
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0 -g")
#
message("building blpre")
#
add_definitions(-D _LIB)
include_directories(./include ./src  ../MNormal)
add_library(blpre STATIC 
include/blpre.h
src/main.cpp
include/MeshInfo.h
src/MeshInfo.cpp 

)