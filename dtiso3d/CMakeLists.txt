#############################################################
#                                                           #
# Author: Jianjing Zheng (zhengjianjing@zju.edu.cn)         #
#                                                           #
#############################################################

cmake_minimum_required (VERSION 3.5) 
option(DTISO3D_WITH_MMG               "Use mmg"                      OFF)
MESSAGE(STATUS "###################################")
MESSAGE("CMake config files are written by Jianjing Zheng")
MESSAGE("Config bug report: zhengjianjing@zju.edu.cn")
IF (WIN32)  
    set(CMAKE_SYSTEM_NAME Windows)
ELSEIF (APPLE)  
    MESSAGE(STATUS "Now is Apple systens.")  
ELSEIF (UNIX)  
    set(CMAKE_SYSTEM_NAME Linux) 
ENDIF () 
FIND_PACKAGE (OpenMP)
if(OPENMP_FOUND)
    message("OpenMP found")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
MESSAGE(STATUS "operation system is ${CMAKE_SYSTEM_NAME}")  
  
IF (CMAKE_SYSTEM_NAME MATCHES "Linux")  
    MESSAGE(STATUS "current platform: Linux ")  
ELSEIF (CMAKE_SYSTEM_NAME MATCHES "Windows")  
    MESSAGE(STATUS "current platform: Windows")  
ELSEIF (CMAKE_SYSTEM_NAME MATCHES "FreeBSD")  
    MESSAGE(STATUS "current platform: FreeBSD")  
ELSE ()  
    MESSAGE(STATUS "other platform: ${CMAKE_SYSTEM_NAME}")  
ENDIF (CMAKE_SYSTEM_NAME MATCHES "Linux")  
  
MESSAGE(STATUS "###################################")  

# projet name
project (DTIso3D)

# if (CMAKE_BUILD_TYPE STREQUAL Debug OR CMAKE_BUILD_TYPE STREQUAL RelWithDebInfo)
if (CMAKE_BUILD_TYPE STREQUAL Debug)
    message("In Debug Mode")
    add_definitions(
        -D_DEBUG
        # -D_PDMG_TEST
        )
endif()
# if(CMAKE_BUILD_TYPE STREQUAL Debug)
#     message("In Debug Mode")
#     add_definitions(
#         -D_DEBUG
#         # -D_PDMG_TEST
#         )
#     set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g")
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
#     set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -g")
#     if (OPTIMIZED)
#         set(OPTIMIZATION_LEVEL "-O3")
#         message("Optimized ${OPTIMIZATION_LEVEL}.")
#         set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OPTIMIZATION_LEVEL}")
#         set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPTIMIZATION_LEVEL}")
#         set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OPTIMIZATION_LEVEL}")
#     endif()
#     if (GLIBCXX_DEBUG)
#         message("Debug glibc++.")
#         add_definitions(-D_GLIBCXX_DEBUG)
#     endif()
# else()
#     message("In Release Mode")
#     add_definitions(-DNDEBUG)
#     set(OPTIMIZATION_LEVEL "-O3")
#     message("Optimized ${OPTIMIZATION_LEVEL}.")
#     set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OPTIMIZATION_LEVEL}")
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPTIMIZATION_LEVEL}")
#     set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OPTIMIZATION_LEVEL}")
# endif()

# openmp
MESSAGE(STATUS "TODO： mac doesnot support opemp yet")
if(DTISO3D_WITH_MMG)
	link_libraries(mmgs mmg3d)
	add_subdirectory(${PROJECT_SOURCE_DIR}/third/mmg)
	include_directories(${CMAKE_BINARY_DIR}/third/mmg/include)
endif()
link_directories(${CMAKE_BINARY_DIR}/lib)


# add subdirectories
add_subdirectory(geom)
#add_subdirectory(GSmooth)
add_subdirectory(spr)
add_subdirectory(dtiso3d)

add_subdirectory(third/tetgen)
file(COPY  "${PROJECT_SOURCE_DIR}/dtiso3d/tiger_volmtetrahedralizer2.h"  
DESTINATION  "${CMAKE_BINARY_DIR}/include")

# 版本号
#set (DTIso3D_VERSION_MAJOR 1)
#set (DTIso3D_VERSION_MINOR 0)
#set (DTIso3D_VERSION_REVISION Belta)

# 配置一个头文件，通过它向源代码中传递一些CMake设置。
#configure_file (
#  "${PROJECT_SOURCE_DIR}/DTIso3D.h.in"
#  "${PROJECT_BINARY_DIR}/DTIso3D.h"
#  )
if(MSVC)
# MSVC 编译器使用 /STACK
set_target_properties(DTIso3D PROPERTIES LINK_FLAGS "/STACK:200000000")
elseif(MINGW)
# MinGW 使用 -Wl,--stack
set_target_properties(DTIso3D PROPERTIES LINK_FLAGS "-Wl,--stack,200000000")
else()
# 其他编译器
endif()