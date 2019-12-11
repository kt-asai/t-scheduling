cmake_minimum_required(VERSION 3.9)

project(t-scheduling CXX)

set(TARGETS t-scheduling)

find_package(Boost)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    link_directories("${Boost_LIBRARY_DIRS}")
endif()

#add_subdirectory(src/tpar)

add_executable(${TARGETS} src/main.cpp
        src/io/circuit_reader.cpp
        src/util/util.cpp
        src/tpar/partition.cpp
        src/tpar/matroid.hpp
        src/circuit.cpp
        src/character.cpp
        src/synthesis.cpp
        src/circuit_builder.cpp
        src/gaussian_decomposer.cpp)


set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -O2")

add_definitions(${CMAKE_CXX_FLAGS} "-std=c++14 -Wall")