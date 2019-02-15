cmake_minimum_required(VERSION 3.5)
project(karman_vortex)

add_definitions(-std=c++14)

set(CXX_FLAGS "-Wall -Wextra -O3")
set(CMAKE_CXX_FLAGS "${CXX_FLAGS}")

add_executable(karman src/Karman.cpp src/main.cpp)
