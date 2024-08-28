# 设置编译的源文件
set(LIBSPR_SRC 
		celltopology.cpp celltopology.h
		gridedgehash.cpp gridedgehash.h
		instentlist.cpp instentlist.h
		IntIntMap.cpp IntIntMap.h
		IntStack.cpp IntStack.h
		SPR.cpp spr.h    )

#AUX_SOURCE_DIRECTORY(. DIR_SRCS)

# 定义宏
add_definitions(-D_DOUBLE -D_OPT_INTERSEC_ALGO)

# 设置头文件路径
include_directories (${PROJECT_SOURCE_DIR}/geom)

# 设置链接的库路径
link_directories(${PROJECT_BINARY_DIR}/lib)

# 设置生成库
if(WIN32)
    add_library(libspr STATIC ${LIBSPR_SRC})
else()
if(SHARED)
    add_library(libspr SHARED ${LIBSPR_SRC})
else()
    add_library(libspr STATIC ${LIBSPR_SRC})
endif()
endif()

# 设置链接的库
target_link_libraries(libspr libgeom)

# 设置库的输出路径
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)

# 设置库的输出名字
set_target_properties(libspr PROPERTIES OUTPUT_NAME "spr")

install(TARGETS libspr
	ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    )
