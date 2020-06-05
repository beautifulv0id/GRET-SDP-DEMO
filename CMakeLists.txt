cmake_minimum_required(VERSION 3.10)
project(EX1)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

find_package(SDPA REQUIRED)

include_directories(${SDPA_INCLUDE_DIRS})

add_executable(example1 examples/example1.cpp)
target_link_libraries(example1 ${SDPA_LIBRARIES})

add_executable(example5 examples/example5.cpp)
target_link_libraries(example5 ${SDPA_LIBRARIES})

add_executable(sdpa_exe examples/sdpa_exe.cpp)
target_link_libraries(sdpa_exe ${SDPA_LIBRARIES})