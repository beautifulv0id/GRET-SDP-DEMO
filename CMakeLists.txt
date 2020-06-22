cmake_minimum_required(VERSION 3.10)
project(EX1)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

set(OpenGR_PROJECT_DIR "../OpenGR_Fork")

list(APPEND CMAKE_PREFIX_PATH "${OpenGR_PROJECT_DIR}/build/install/lib/cmake")

find_package(OpenGR REQUIRED)
include_directories(${OpenGR_INCLUDE_DIR})
link_directories(${OpenGR_LIB_DIR})

find_package(SDPA REQUIRED)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)
include_directories( ${EIGEN3_INCLUDE_DIR} ) 


find_package(Boost REQUIRED)
include_directories( ${Boost_INCLUDE_DIR} )

include_directories(${SDPA_INCLUDE_DIRS})

add_executable(gret-sdp-demo gret-sdp-demo.cpp)
target_link_libraries(gret-sdp-demo ${SDPA_LIBRARIES} ${Boost_LIBRARIES} Eigen3::Eigen ${OpenGR_LIBRARIES} gr::utils gr::algo) 