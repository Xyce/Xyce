# SuiteSparse web page: http://suitesparse.com/
#
# The only part of SuiteSparse used by Xyce is AMD.
# When pointed to a SuiteSparse repository, this file is designed to build and
# install *only* the AMD library into a specified location.
#
# The recommended usage is to first obtain SuiteSparse. Either download a
# release tarball or make a shallow clone of the repository from GitHub:
#    git clone --depth 1 https://github.com/DrTimothyAldenDavis/SuiteSparse.git
#
# Next, use the following invocation from a build directory, making the install
# directory the *same one* as the Trilinos installation target:
#    cmake \
#    -D SuiteSparsePath=path/to/SuiteSparse \
#    -D CMAKE_INSTALL_PREFIX=/trilinos/install/path \
#    path/to/Xyce/cmake/trilinos/AMD
#
# Build/install AMD in the usual way; e.g., from the AMD build directory, run:
#    cmake --build . -t install
#
# Then be sure to use the following flags with the Trilnos CMake invocation:
#    -D CMAKE_INSTALL_PREFIX=/trilinos/install/path \
#    -D AMD_LIBRARY_DIRS=/trilinos/install/path/lib \
#    -D AMD_INCLUDE_DIRS=/trilinos/install/path/include \

cmake_minimum_required(VERSION 3.13 FATAL_ERROR)

project(AMD LANGUAGES C)

if(NOT DEFINED SuiteSparsePath)
     message(FATAL_ERROR "Must supply the path to SuiteSparse using:
     -D SuiteSparsePath=path/to/SuiteSparse")
endif()

set(SuiteSparsePath "" CACHE FILEPATH "Path to the SuiteSparse repository")
set(CMAKE_BUILD_TYPE Release)

add_library(amd
  STATIC
    ${SuiteSparsePath}/AMD/Source/amd_aat.c
    ${SuiteSparsePath}/AMD/Source/amd_1.c
    ${SuiteSparsePath}/AMD/Source/amd_2.c
    ${SuiteSparsePath}/AMD/Source/amd_dump.c
    ${SuiteSparsePath}/AMD/Source/amd_postorder.c
    ${SuiteSparsePath}/AMD/Source/amd_defaults.c
    ${SuiteSparsePath}/AMD/Source/amd_post_tree.c
    ${SuiteSparsePath}/AMD/Source/amd_order.c
    ${SuiteSparsePath}/AMD/Source/amd_control.c
    ${SuiteSparsePath}/AMD/Source/amd_info.c
    ${SuiteSparsePath}/AMD/Source/amd_valid.c
    ${SuiteSparsePath}/AMD/Source/amd_preprocess.c
    ${SuiteSparsePath}/SuiteSparse_config/SuiteSparse_config.c
  )

target_include_directories(amd
  PUBLIC
     ${SuiteSparsePath}/AMD/Include/
     ${SuiteSparsePath}/SuiteSparse_config/
  )

set_property(TARGET amd PROPERTY POSITION_INDEPENDENT_CODE ON)

install(TARGETS amd DESTINATION lib)
install(FILES
     ${SuiteSparsePath}/AMD/Include/amd.h
     ${SuiteSparsePath}/SuiteSparse_config/SuiteSparse_config.h
     DESTINATION include
  )


# NOTE 1:
#
# This could probably be combined with the Xyce/cmake/Patch_suitesparse.cmake
# file at some point, maybe as part of a modified superbuild?
#
#
# NOTE 2:
#
# This file does not build amd_global.c, which just contains a comment block
# stating the functions it provided are now in
# SuiteSparse_config/SuiteSparse_config.c.
#
#
# NOTE 3:
#
# The following fragments would more faithfully reproduce what the Makefiles
# do.  However, the AMD package has no ifdef's with DINT or DLONG, so it
# doesn't seem to make sense to do all the extra compiling with those flags.
# (SuiteSparse_config doesn't have them either.)  However, if AMD is ever
# changed to preferentially use ints vs longs, then something like the
# following will need to be done.
#
#    add_library(amd_integer
#      STATIC
#        amd_aat.c
#        amd_1.c
#        amd_2.c
#        amd_dump.c
#        amd_postorder.c
#        amd_defaults.c
#        amd_post_tree.c
#        amd_order.c
#        amd_control.c
#        amd_info.c
#        amd_valid.c
#        amd_preprocess.c
#      )
#    target_compile_definitions(amd_integer PUBLIC DINT)
#
#
#    add_library(amd_long
#      STATIC
#        amd_aat.c
#        amd_1.c
#        amd_2.c
#        amd_dump.c
#        amd_postorder.c
#        amd_defaults.c
#        amd_post_tree.c
#        amd_order.c
#        amd_control.c
#        amd_info.c
#        amd_valid.c
#        amd_preprocess.c
#      )
#    target_compile_definitions(amd_long PUBLIC DLONG)
#
#    target_link_libraries(amd suitesparseconfig amd_integer amd_long)

