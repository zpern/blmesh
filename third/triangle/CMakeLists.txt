cmake_minimum_required (VERSION 3.5)
project(triangle)
#best opt
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -g")
#
message("building triangle")



add_library(triangle STATIC 
triangle.h
triangle.cpp
)
set_target_properties( triangle
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/../../tiger/lib"
    LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/../../tiger/lib"
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/../../tiger/dll"
)