
cmake_minimum_required (VERSION 3.5)
project(tetgen)
#best opt
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -g")
#
message("building tetgen")
add_definitions(-D _LIB)



add_library(tetgen STATIC 
./tetgen.cpp
./tetgen.h
./predicates.cxx
)
set_target_properties( tetgen
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/../../../tiger/lib"
    LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/../../../tiger/lib"
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/../../../tiger/dll"
)