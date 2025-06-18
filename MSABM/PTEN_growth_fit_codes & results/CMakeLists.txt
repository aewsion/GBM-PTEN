### CMakeLists for p_ten_growth MSABM
project(p_ten_growth)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
aux_source_directory(src SRC_LIST)
include_directories(inc)
add_executable(p_ten_growth ${SRC_LIST})