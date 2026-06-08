cmake_minimum_required (VERSION 3.5)
project(vmesh)
message("building vmesh")
option(VMESH_BUILD_EXE "Build the local vmesh API test executable" ON)
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++14" COMPILER_SUPPORTS_CXX14)

include_directories(./include ./src ../blmesh/include ../BLPre/include/ ../MNormal/include/)
include_directories(../dtiso3d/third/tetgen)


if(EXISTS ${CMAKE_SOURCE_DIR}/extern)
    include_directories(${CMAKE_SOURCE_DIR}/extern/libigl/include)
    include_directories(${CMAKE_SOURCE_DIR}/extern/eigen)
else()
    include_directories(../third/libigl/include)
    include_directories(../third/eigen)
endif()


include_directories(${CMAKE_SOURCE_DIR}/extern/googletest)
include_directories(../third/googletest/googletest/include)

include_directories(../tiger-sizeingfunc)
include_directories(./testcase/include)

link_directories(${PROJECT_SOURCE_DIR}/../tiger/lib)
link_directories(${PROJECT_BINARY_DIR}/../lib)

add_library(vmesh STATIC
../blmesh/include/blmeshapi.h
../BLPre/include/blpre.h
include/Config.h
include/vmesh.h
../common/singleton.h
../common/singleton_terminate.h
src/vmesh.cpp
src/Config.cpp
)

if(VMESH_BUILD_EXE)
add_executable(vmesh_test
src/main.cpp
)

target_link_libraries(vmesh_test PRIVATE
vmesh
libblmesh
tetgen
triangle
2dremesh
blpre
libdtiso3d
libdtiso2d
libgeom
meshOrientLib
libspr
libMeshChamfer
DistanceCaculator
libpostprocess
OpenMP::OpenMP_CXX
)
set_target_properties(vmesh_test PROPERTIES OUTPUT_NAME vmesh_test)
endif()
# if(WIN32)
# add_executable(test 
# include/vmesh.h
# ../tiger-sizeingfunc/tiger_sizingfunction.h
# ../tiger-sizeingfunc/tiger_sizingsmoother.h
# ./testcase/include/writeblmesh.h
# ./testcase/include/common_assert.h
# ./testcase/include/test_volume_density.h
# ./testcase/include/test_sphere_tetra.h
# ./testcase/include/test_sphere_boundary.h
# ./testcase/include/test_cubic_tetra.h
# ./testcase/include/test_cubic_boundary.h
# ./testcase/include/test_f6_box_tetra.h
# ./testcase/include/test_f6_box_boundary.h
# ./testcase/include/test_non_manifold.h
# ./testcase/include/test_anistropic.h
# ./testcase/include/test_f6_box_vol.h
# ./testcase/include/test_f6_box_symm.h
# ./testcase/include/test_intersection.h
# ./testcase/include/test_intersection_tetra.h
# ./testcase/include/test_intersection_vol.h
# ./testcase/include/test_postprocess_vol.h
# ./testcase/include/test_two_cubic.h

# src/main.cpp
# )
# if(UNIX)
# FIND_PACKAGE (OpenMP)
# if(OPENMP_FOUND)
#     message("OpenMP found")
#     set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#     set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
# endif()
# else()
# message("openmp are not support")
# endif()

# set_target_properties(
# test
# PROPERTIES EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/vmesh/testcase
# )
# target_link_libraries(test PUBLIC vmesh) 
# target_link_libraries(test PUBLIC gtest)
# target_link_libraries(test PUBLIC libblmesh)
# #target_link_libraries(test PUBLIC create_sol)
# target_link_libraries(test PUBLIC tetgen)
# target_link_libraries(test PUBLIC triangle)
# target_link_libraries(test PUBLIC 2dremesh)
# target_link_libraries(test PUBLIC blpre libdtiso3d libdtiso2d libgeom meshOrientLib libspr libMeshChamfer DistanceCaculator libpostprocess)
# endif()


