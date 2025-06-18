### CMakeLists for MSABM p_ten_combine
cmake_minimum_required(VERSION 3.20)
project(p_ten_combine)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
aux_source_directory(src SRC_LIST)
include_directories(inc)
add_executable(p_ten_combine ${SRC_LIST})