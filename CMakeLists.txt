cmake_minimum_required(VERSION 3.6)
project(geometry)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES main.cpp geometry.cpp templ_stack.h)
add_executable(geometry ${SOURCE_FILES})