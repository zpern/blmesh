cmake_minimum_required (VERSION 2.6)
project(create_sol)
message("building create_sol")
add_definitions(-D _LIB)


if(EXISTS ${CMAKE_SOURCE_DIR}/extern)
    include_directories(${CMAKE_SOURCE_DIR}/extern/libigl/include)
    include_directories(${CMAKE_SOURCE_DIR}/extern/eigen)
    include_directories(${CMAKE_SOURCE_DIR}/geogram/src/lib)
else()
    include_directories(../third/libigl/include)
    include_directories(../third/eigen)
    include_directories(../third/geogram/src/lib)
endif()

include_directories(../third)
set(EXTERN_DIR "${CMAKE_SOURCE_DIR}/extern")

message("source = ${EXTERN_DIR}/geogram/src/lib")
if(EXISTS "${EXTERN_DIR}/geogram/src/lib")
    include_directories(${EXTERN_DIR}/geogram/src/lib)
else()
    include_directories(../third/geogram/src/lib)
endif()
include_directories(../third/eigen/libigl-2.1.0/include)

include_directories(./testcase/include)

add_library(create_sol STATIC 
bkmeshimpl.h
bkmeshimpl.cpp
createsol.h
createsol.cpp
vector_own.h
vector_own.cpp
)
set_target_properties( create_sol
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/../tiger/lib"
    LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/../tiger/lib"
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/../tiger/dll"
)
