cmake_minimum_required(VERSION 3.10)
project(EX1)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

find_package(SDPA REQUIRED)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)
include_directories( ${EIGEN3_INCLUDE_DIR} ) 

include_directories(${SDPA_INCLUDE_DIRS})

add_executable(gret-sdp-demo gret-sdp-demo.cpp)
target_link_libraries(gret-sdp-demo ${SDPA_LIBRARIES} Eigen3::Eigen)