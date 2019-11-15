cmake_minimum_required(VERSION 3.9)

project(t-scheduling CXX)

set(TARGETS t-scheduling)

add_executable(${TARGETS} src/main.cpp src/circuit_reader.cpp src/circuit.cpp src/character.cpp)

set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_FLAGS "")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -O2")

add_definitions(${CMAKE_CXX_FLAGS} "-std=c++14 -Wall")