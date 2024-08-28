cmake_minimum_required (VERSION 2.6)

project(MeshChamfer)
#best opt
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -g")
#
set(SOURCEFILE 
include/blvector.h
include/sphericalCap.h
include/complexnode.h
include/geometryfunction.h
include/MNormalMesh.h
include/VirtualSphereMesh.h
include/VirtualSphereMeshGenerator.h
include/VirtualSphereMeshHasher.h
include/MeshSplitter.h
include/singleton.h
include/BoundaryTriangulation.h
include/Common.h
include/Splitter.h
include/VirtualSphereMeshStrategy.h
include/ChamferBehavior.h
include/AbstractOptimizer.h
include/TopologyOptimizer.h
include/PointOptimizer.h
include/MeshEvaluation.h
include/mergeoptimizer.h
include/combineoptimizer.h
include/mostnormaloptimizer.h

geom/geom_func.h

src/combineoptimizer.cpp
src/mergeoptimizer.cpp
src/meshevaluation.cpp
src/TopologyOptimizer.cpp
src/VirtualSphereMeshStrategy.cpp
src/BoundaryTriangulation.cpp
src/VirtualSphereMeshHasher.cpp
src/MeshSplitter.cpp
src/pointoptimizer.cpp
src/VirtualSphereMeshGenerator.cpp
src/VirtualSphereMesh.cpp
src/BLVector.cpp
src/sphericalCap.cpp
src/complexnode.cpp
src/geometryfunction.cpp
src/MNormalMesh.cpp
src/mostnormaloptimizer.cpp
)
link_directories(${PROJECT_SOURCE_DIR}/lib)
IF(APPLE)
MESSAGE(STATUS "TODO： mac doesnot support opemp yet")
ELSE()
FIND_PACKAGE (OpenMP REQUIRED)
ENDIF()


#set(Boost_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/../third/boost_1_74_0")
#find_package(Boost COMPONENTS system filesystem REQUIRED)
#if(Boost_FOUND)
#    message(Boost_INCLUDE_DIRS " ${Boost_INCLUDE_DIRS}")
#    message(Boost_LIBRARY_DIRS " ${Boost_LIBRARY_DIRS}")
#   include_directories(${Boost_INCLUDE_DIRS})
#    link_directories(${Boost_LIBRARY_DIRS})
#endif()

#add_subdirectory(geom)

message("building MeshChamfer")
#
link_directories(${PROJECT_SOURCE_DIR}/../tiger/lib)
include_directories(./include ./src )
include_directories(./third/spdlog/include)
include_directories(${CMAKE_SOURCE_DIR}/extern/spdlog/include)
add_definitions(-D _CHECK_INTERSECTION)
add_library(libMeshChamfer STATIC 
${SOURCEFILE}
)
set_target_properties(libMeshChamfer PROPERTIES OUTPUT_NAME MeshChamfer)
add_executable(MeshChamfer
${SOURCEFILE}
src/main.cpp
)
target_link_libraries(MeshChamfer PUBLIC libgeom)


