cmake_minimum_required (VERSION 2.6)
project(blmesh)
#best opt
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -g")
#
set(SOURCEFILE 
include/algebra.h
include/BLFront.h
include/BLFrontList.h
include/BLMesh.h
include/BLMesh_define.h
include/NormalSmoothStrategy.h
include/IterationNormalSmoothStrategy.h
include/SimpleNormalSmoothStrategy.h
include/BLNode.h
include/common.h
include/PrismStorehouse.h
../dtiso2d/genmesh2d.h
../dtiso3d/dtiso3d/tiger_volmtetrahedralizer.h
../BLPre/include/blpre.h
../DistanceCaculator/include/dataio.h
../verdict/verdict.h
include/geom_func.h
include/mesh_orient.h
include/Octree.h
include/OctreeAgent.h
include/blmeshapi.h
include/BLVector.h
include/DynamicArray.h
include/LightMemoryPool.h
include/VirtualMesh.h
include/distancecalculation.h
include/singlesymmetry.h
include/VoronoiNormalSmoothStrategy.h
include/sphericalvoronoigeometry.h
include/WangzhiNormalSmoothStrategy.h
include/HeightSmoothingStrategy.h
include/symmetry.h

src/binary_tree.hpp
src/symmetry.cpp
src/HeightSmoothingStrategy.cpp
src/WangzhiNormalSmoothStrategy.cpp
src/VoronoiNormalSmoothStrategy.cpp
src/distancecalculation.cpp
src/VirtualMesh.cpp
src/algebra.cpp
src/BLFront.cpp
src/BLMesh.cpp
src/BLFrontList.cpp
src/BLMesh.cpp
src/BLNode.cpp
src/LightMemoryPool.cpp
src/NormalSmoothStrategy.cpp
src/IterationNormalSmoothStrategy.cpp
src/SimpleNormalSmoothStrategy.cpp
src/Octree.cpp
src/OctreeAgent.cpp
src/BLVector.cpp
src/DynamicArray.cpp
src/PrismStorehouse.cpp
)
link_directories(${PROJECT_SOURCE_DIR}/lib)



#set(Boost_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/../third/boost_1_74_0")
#find_package(Boost COMPONENTS system filesystem REQUIRED)
#if(Boost_FOUND)
#    message(Boost_INCLUDE_DIRS " ${Boost_INCLUDE_DIRS}")
#    message(Boost_LIBRARY_DIRS " ${Boost_LIBRARY_DIRS}")
#   include_directories(${Boost_INCLUDE_DIRS})
#    link_directories(${Boost_LIBRARY_DIRS})
#endif()

FIND_PACKAGE (OpenMP)
if(OPENMP_FOUND)
    message("OpenMP found")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

message("building blmesh")
if(EXISTS ${CMAKE_SOURCE_DIR}/extern/libigl/CMakeLists.txt) 
include_directories(${CMAKE_SOURCE_DIR}/extern/libigl/include)
include_directories(${CMAKE_SOURCE_DIR}/extern/eigen)
include_directories(${CMAKE_SOURCE_DIR}/extern/spdlog/include)
endif()

if(EXISTS ${CMAKE_SOURCE_DIR}/extern/spdlog/CMakeLists.txt) 
include_directories(${CMAKE_SOURCE_DIR}/extern/spdlog/include)
include_directories(${CMAKE_SOURCE_DIR}/extern/eigen)
endif()

link_directories(${PROJECT_SOURCE_DIR}/../tiger/lib)
link_directories(${PROJECT_BINARY_DIR}/../lib)
if(EXISTS ${CMAKE_SOURCE_DIR}/extern)
    include_directories(${CMAKE_SOURCE_DIR}/extern/libigl/include)
    include_directories(${CMAKE_SOURCE_DIR}/extern/eigen)
    link_directories(${CMAKE_SOURCE_DIR}/lib)
else()
    include_directories(../third/libigl/include)
    include_directories(../third/eigen)
endif()

include_directories(../third/triangle)
include_directories(./include ./src ../dtiso2d/ ../dtiso3d/dtiso3d/ ../BLPre/include ../DistanceCaculator/include )
add_definitions(-D _CHECK_INTERSECTION)
add_library(libblmesh STATIC 
${SOURCEFILE}
src/blmeshapi.cpp
)
set_target_properties(libblmesh PROPERTIES OUTPUT_NAME blmesh)
add_executable(blmesh
${SOURCEFILE}
src/main.cpp
)
if(WIN32)
target_link_libraries(
blmesh

-Wl,--start-group

libdtiso3d 
libdtiso2d 
libgeom 
meshOrientLib 
libspr 
blpre 
libMeshChamfer
libverdict 
DistanceCaculator 
libpostprocess
tetgen
2dremesh
triangle

-Wl,--end-group
)
else()
target_link_libraries(blmesh PUBLIC libdtiso3d libdtiso2d libgeom meshOrientLib libspr blpre libMeshChamfer libverdict DistanceCaculator libpostprocess)
get_directory_property(hasParent PARENT_DIRECTORY)
if(hasParent)
target_link_directories(blmesh PUBLIC ${hasparent}/bin)
else()
endif()
target_link_libraries(blmesh PUBLIC tetgen)
target_link_libraries(blmesh PUBLIC 2dremesh)
target_link_libraries(blmesh PUBLIC triangle)
endif()
