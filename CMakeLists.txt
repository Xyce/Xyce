cmake_minimum_required(VERSION 3.10 FATAL_ERROR) 
project(Xyce VERSION 6.11 )   	
#cmake documentation says that C and CXX are enabled if no lang is specified.  
# however, check_include_file_cxx isn't found unless 'LANGUAGES CXX'
# is ALONE! 
# wut? Did I find a bug in cmake >3.12? 

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
include(cmake/config.cmake)

add_subdirectory(src) 
