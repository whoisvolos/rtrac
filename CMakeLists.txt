cmake_minimum_required(VERSION 3.0.2)
project(MPI_MC)

set(CMAKE_INCLUDE_PATH "")

find_package(MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})

set(SOURCE_FILES
    hello3.cpp
    funcs.hpp
    time.hpp
    integrate.hpp
    integrate.cpp
    math/operations.h
    math/operations_ext.h
    math/plane.h
    math/plane.cpp
    math/triangle.h
    math/triangle.cpp
    math/types.h)

include_directories(${MPI_INCLUDE_PATH})

add_executable(hello ${SOURCE_FILES} math/plane.h math/operations_ext.h)
target_link_libraries(hello ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
  set_target_properties(hello PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(hello PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_BUILD_TYPE Release)
