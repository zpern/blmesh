project(postprocess)
set(POST_SRC
	postprocess.cpp
	postprocess.h
)

#AUX_SOURCE_DIRECTORY(. DIR_SRCS)

add_definitions(-DVERDICT_EXPORTS)
add_library(libpostprocess  STATIC ${POST_SRC})
set_target_properties(libpostprocess PROPERTIES OUTPUT_NAME postprocess)
# 设置库的输出路径
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)

# 设置可执行程序的输出路径
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
