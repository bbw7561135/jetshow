cmake_minimum_required(VERSION 3.6)
project(jetshow)

include_directories("/usr/include/eigen3", "include")

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -fopenmp -fext-numeric-literals")
set(SOURCE_FILES src/main.cpp src/BField.cpp src/Cone.cpp include/Cone.h src/Geometry.cpp include/Geometry.h src/Ray.cpp include/Ray.h src/Intersection.cpp include/Intersection.h include/Cell.h src/Cell.cpp include/linspace.h src/Pixel.cpp include/Pixel.h src/Image.cpp include/Image.h src/ImagePlane.cpp include/ImagePlane.h src/NField.cpp include/NField.h src/VField.cpp include/VField.h src/Jet.cpp include/Jet.h include/utils.h src/Cylinder.cpp include/Cylinder.h src/utils.cpp src/System.cpp include/System.h include/Observation.h src/Observation.cpp src/Cells.cpp include/Cells.h)
add_executable(jetshow ${SOURCE_FILES})

set(TEST_DIR tests)
set(TEST_FILES ${TEST_DIR}/simple_test.cpp include/Cone.h include/Intersection.h
        include/Ray.h src/Cone.cpp src/Intersection.cpp src/Ray.cpp
        include/Geometry.h src/Geometry.cpp include/utils.h)
add_executable(jetshow_test ${TEST_FILES})


# Require MPI for this project:
find_package(MPI REQUIRED)
set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
include_directories(${MPI_INCLUDE_PATH})
target_link_libraries(jetshow ${MPI_LIBRARIES})