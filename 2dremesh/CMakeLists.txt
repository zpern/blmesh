cmake_minimum_required (VERSION 3.5)
project(2dremesh)
#best opt
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0 -g")
#
message("building 2dremesh.lib")
#
add_definitions(-D _LIB)




file(GLOB INC_FILES  ./*.h)
file(GLOB SRC_FILES  ./*.cpp)

set(ALL_FILES 
${INC_FILES}
${SRC_FILES} 
)
add_library(2dremesh STATIC ${ALL_FILES})

get_directory_property(hasParent PARENT_DIRECTORY)
if(hasParent)
target_include_directories(2dremesh PUBLIC ${hasParent}/third_party/eigen)
else()
target_include_directories(2dremesh PUBLIC ../third/eigen)
endif()
