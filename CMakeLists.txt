cmake_minimum_required(VERSION 3.2)

project(t-scheduling CXX)

add_executable(main main.cpp)

set(CMAKE_CXX_FLAGS "")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -O2")

add_definitions(${CMAKE_CXX_FLAGS} "-std=c++14 -Wall")