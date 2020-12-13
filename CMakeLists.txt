cmake_minimum_required(VERSION 3.10)
project(plane_finding)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_BUILD_TYPE "debug")

find_package(Pangolin REQUIRED)
include_directories(${Pangolin_INCLUDE_DIRS})
find_package(PCL 1.8 REQUIRED )
include_directories( ${PCL_INCLUDE_DIRS} )
add_definitions( ${PCL_DEFINITIONS} )
link_directories( ${PCL_LIBRARY_DIRS} )

add_executable(plane_finding main.cpp)
target_link_libraries(plane_finding ${Pangolin_LIBRARIES} ${PCL_LIBRARIES})