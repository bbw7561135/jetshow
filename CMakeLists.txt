cmake_minimum_required(VERSION 3.6)
project(jetshow)

include_directories("/usr/include/eigen3", "include")

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES src/main.cpp src/Bfield.cpp src/Cone.cpp include/Cone.h src/Geometry.cpp include/Geometry.h src/Ray.cpp include/Ray.h src/Intersection.cpp include/Intersection.h include/Cell.h src/Cell.cpp include/linspace.h src/Pixel.cpp include/Pixel.h src/Image.cpp include/Image.h src/ImagePlane.cpp include/ImagePlane.h src/NField.cpp include/NField.h src/VField.cpp include/VField.h src/Jet.cpp include/Jet.h include/utils.h)
add_executable(jetshow ${SOURCE_FILES})