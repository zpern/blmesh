# 设置编译的源文件
set(COMMON_SRC
	blockedarray_instance.cpp
	earremoval.cpp
	iso3d_common.cpp
	iso3d.cpp
	iso3d_error.cpp
	iso3d_utility.cpp
	iso3d_multidomain.cpp
	locsmoothing.cpp
	matrix.cpp
	optdihangle.cpp
	skirtpolylist.cpp
	space.cpp
	tetramesh.cpp
	vector.cpp
	solvezone.cpp
	singleton_terminate.h
	tiger_volmtet.cpp
	tiger_volmtetrahedralizer.cpp #zyj 先去掉原来的接口  #TOYHF 需要改回原来的接口
#	tiger_volmtetrahedralizer2.cpp
	)
#IF (CMAKE_SYSTEM_NAME MATCHES "Windows")
IF (0)
	set(COMMON_SRC
            ${COMMON_SRC}
            ComputInfo.cpp
            HardDriveSerialNumer.cpp
			LiCheck.cpp
			md5.cpp
			)
ENDIF ()
set(SHARED false)
set(DTISO3D_SRC 
	${COMMON_SRC}
	Iso3D_Main.cpp)

set(libdtiso3d_SRC 
	${COMMON_SRC}
	iso3d_libmain.cpp)

set(libvolmt_SRC
	${COMMON_SRC}
#	tiger_volmtetrahedralizer.cpp #zyj 先去掉原来的接口
	tiger_volmtetrahedralizer.h
	)

set(volmt_test_SRC	
	 volmt_test.cpp
	 tiger_volmtetrahedralizer.h
	 )	
	
# 定义宏
add_definitions(
	-D_EXT_BOUND_CLASSIFIED
	-D_BR_MAIN_PROCEDURE_SPR
	-D_CONSTRAINED_RECOVERY
	-D_DOUBLE
	# -D_ERROR_CHK
	# -D_LOAD_BALANCING
	# -D_MULTITHREAD
	-D_NEW_HEAP_DATASTRUCTURE
	-D_RECOVER_WITHOUT_SPS
	-D_REPOSITION_STEINER_PNT
	-D_ROBUST_BW_KERNEL
	-D_SMOOTH_IN_QUEUE
	-D_TIMING_PERFORMANCE
	-D_WRITE_CPP_DOUBLE
	-D_RECV_EDGE_BY_MID_SP
	#_D_KEEP_ALL_INNER #保留所有内部网格（用于非流体算例）
   # -D_SPR_IMP  # 采用SPR优化网格
)

IF (CMAKE_SYSTEM_NAME MATCHES "Windows")
add_definitions(
	-DUNICODE
)
ENDIF ()

# 设置stack commit size.    (stack reserve size? How)
IF (WIN32)  
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /STACK:200000000")
ENDIF()

# 设置头文件路径
include_directories (
	${PROJECT_SOURCE_DIR}/spr
	${PROJECT_SOURCE_DIR}/geom
	${PROJECT_SOURCE_DIR}/meshmetric)

# 设置链接的库路径
link_directories(${CMAKE_BINARY_DIR}/lib)


########################################################
# DTIso3d
########################################################
# 设置生成可执行程序
add_executable(dtiso3d ${DTISO3D_SRC})

# 设置链接的库
target_link_libraries(dtiso3d libspr)
target_link_libraries(dtiso3d PUBLIC tetgen)
########################################################
# libdtiso3d
########################################################
if(SHARED)
    add_library(libdtiso3d SHARED ${libdtiso3d_SRC})
else()
    add_library(libdtiso3d STATIC ${libdtiso3d_SRC})
endif()

target_link_libraries(libdtiso3d libspr)


set_target_properties(libdtiso3d PROPERTIES OUTPUT_NAME dtiso3d)


########################################################
# volmt_test
########################################################
# 设置生成可执行程序
#add_executable(volmt_test ${volmt_test_SRC})

# 设置链接的库
#target_link_libraries(volmt_test libvolmt)


########################################################
# libvolmt
########################################################
if(SHARED)
    add_library(libvolmt SHARED ${libvolmt_SRC})
else()
    add_library(libvolmt STATIC ${libvolmt_SRC})
endif()

target_link_libraries(libvolmt libspr)
if(DTISO3D_WITH_MMG)
 target_link_libraries(libvolmt PUBLIC debug mmgs_d optimized mmgs_r)
 target_link_libraries(libvolmt PUBLIC debug mmg3d_d optimized mmg3d_r)
endif()

set_target_properties(libvolmt PROPERTIES OUTPUT_NAME volmt)

# 设置库的输出路径
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)

# 设置可执行程序的输出路径
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

install(TARGETS dtiso3d libdtiso3d libvolmt
	RUNTIME DESTINATION bin
	ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    )
	

