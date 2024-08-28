set(CMAKE_DEBUG_POSTFIX "_d")
set(CMAKE_RELEASE_POSTFIX "_r")
cmake_minimum_required(VERSION 2.7.1 FATAL_ERROR)


project(Tiger-SizingFunc LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


SET(CMAKE_DEBUG_POSTFIX "d" CACHE STRING "add a postfix, usually d on windows")
SET(CMAKE_RELEASE_POSTFIX "" CACHE STRING "add a postfix, usually empty on windows")
SET(CMAKE_RELWITHDEBINFO_POSTFIX "rd" CACHE STRING "add a postfix, usually empty on windows")
SET(CMAKE_MINSIZEREL_POSTFIX "s" CACHE STRING "add a postfix, usually empty on windows")

include_directories(../third/eigen)
include_directories(../third)
include_directories(../third/eigen/libigl-2.1.0/include)
set(EXTERN_DIR "${CMAKE_SOURCE_DIR}/extern")

message("source = ${EXTERN_DIR}/geogram/src/lib")
# ÅÐ¶Ï¸¸Ä¿Â¼ÊÇ·ñ´æÔÚ
if(EXISTS "${EXTERN_DIR}/geogram/src/lib")
    include_directories(${EXTERN_DIR}/geogram/src/lib)
else()
    include_directories(../third/geogram/src/lib)
endif()
#set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")



set (COM_LIBS_DIR

    )
if (CMAKE_BUILD_TYPE MATCHES "Debug")
    link_directories(
        ../third_party_libs/win64Debug
        ../third_party_libs/mmgs/win64Debug
        ../third_party_libs/ipsolve55/libs/win64Debug
        ../third_party_libs/ipoptDebug/lib
        ../lib
        ${COM_LIBS_DIR}
        )
elseif(CMAKE_BUILD_TYPE MATCHES "Release")
    link_directories(
        ${COM_LIBS_DIR}
        ${PROJECT_SOURCE_DIR}/geogram_1.7.4/install/lib/Release
        )
endif()

if(WIN32)
add_compile_options("/utf-8" "/wd4828" "/openmp")
endif()
add_library(TigerAdasize STATIC
    createsol.cpp
    inputdata.cpp
    tiger_sizingfunction.cpp
    vector_own.cpp
    createsol.h
    inputdata.h
    tiger_sizingfunction.h
    vector_own.h
    )
link_directories(${PROJECT_SOURCE_DIR}/../tiger/lib)
link_directories(${PROJECT_BINARY_DIR}/../lib)
link_directories(${PROJECT_SOURCE_DIR}/lib)
target_link_libraries(TigerAdasize PUBLIC debug geogram_d optimized geogram_r)
set_target_properties( TigerAdasize
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/../tiger/lib"
    LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/../tiger/lib"
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/../tiger/dll"
)
