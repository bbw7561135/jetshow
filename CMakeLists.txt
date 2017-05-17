cmake_minimum_required(VERSION 3.6)
project(jetshow)

include_directories("/usr/include/eigen3", "include")

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES src/main.cpp src/bfields.cpp src/nfields.cpp src/vfields.cpp src/geometry.cpp src/jet.cpp include/geometry.h include/jet.h include/vfields.h include/nfields.h)
add_executable(jetshow ${SOURCE_FILES})