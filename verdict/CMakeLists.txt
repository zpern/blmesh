project(verdict)
set(LIBVerdict_SRC 
	V_EdgeMetric.cpp
	V_GaussIntegration.cpp
	V_HexMetric.cpp	
	V_KnifeMetric.cpp
	V_PrismMetric.cpp
	V_PyramidMetric.cpp
	V_QuadMetric.cpp
	V_TetMetric.cpp
	V_TriMetric.cpp
	VerdictVector.cpp
V_GaussIntegration.hpp
v_vector.h
verdict.h
Verdict.htm
verdict_defines.hpp
VerdictVector.hpp
)

#AUX_SOURCE_DIRECTORY(. DIR_SRCS)

add_definitions(-DVERDICT_EXPORTS)
add_library(libverdict  STATIC ${LIBVerdict_SRC})
set_target_properties(libverdict PROPERTIES OUTPUT_NAME verdict)
# 设置库的输出路径
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)

# 设置可执行程序的输出路径
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)


# ÉèÖÃÉú³É¶¯Ì¬¿â  


# ÉèÖÃŸ²Ì¬¿âµÄÊä³öÂ·Ÿ¶
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
