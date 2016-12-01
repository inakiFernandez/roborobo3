cmake_minimum_required(VERSION 2.4)

PROJECT( auxlearnfunction )

file ( GLOB_RECURSE all_include auxlearnfunction.h )
file ( GLOB all_sources auxlearnfunction.cpp )
set(CMAKE_CXX_FLAGS "-Wall -Wextra -std=c++11 -Wno-switch") 
include_directories(. ../include/contrib/ ../include/core/Utilities/)

add_library(auxlearnfunction STATIC ${all_sources})

#ADD_EXECUTABLE( auxlearnfunction ${ex_sources} )
