cmake_minimum_required(VERSION 3.5)

project(t-scheduling CXX)

set(TARGETS t-scheduling)

find_package(Z3)
if(Z3_FOUND)
    message(STATUS "Z3_FOUND: ${Z3_FOUND}")
    message(STATUS "Found Z3 ${Z3_VERSION_STRING}")
    message(STATUS "Z3_DIR: ${Z3_DIR}")

    include_directories(${Z3_CXX_INCLUDE_DIRS})
    link_libraries(${Z3_LIBRARIES})
endif()

find_package(Boost)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    link_directories("${Boost_LIBRARY_DIRS}")
endif()

add_executable(${TARGETS} src/main.cpp
        src/io/circuit_reader.cpp
        src/util/util.cpp
        src/tpar/partition.cpp
        src/tpar/matroid.hpp
        src/circuit/circuit.cpp
        src/character/character.cpp
        src/decomposer/gaussian_decomposer.cpp
        src/decomposer/parallel_decomposer.cpp
        src/synthesis/tpar_synthesis.cpp
        src/synthesis/tskd_synthesis.cpp
        src/synthesis/simple_circuit_builder.cpp
        src/synthesis/greedy_circuit_builder.cpp
        src/matrix/matrix_reconstructor.cpp
        src/parallel/parallelization_oracle.cpp
        src/layout/layout.cpp
        src/simulator/simulator.cpp)

set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -O2")

add_definitions(${CMAKE_CXX_FLAGS} "-std=c++14 -Wall")