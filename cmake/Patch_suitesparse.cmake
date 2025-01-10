cmake_minimum_required(VERSION 3.10)

project(suitesparse C CXX)

add_library(amd
  ${suitesparse_SOURCE_DIR}/AMD/Source/amd_1.c
  ${suitesparse_SOURCE_DIR}/AMD/Source/amd_2.c
  ${suitesparse_SOURCE_DIR}/AMD/Source/amd_aat.c
  ${suitesparse_SOURCE_DIR}/AMD/Source/amd_control.c
  ${suitesparse_SOURCE_DIR}/AMD/Source/amd_defaults.c
  ${suitesparse_SOURCE_DIR}/AMD/Source/amd_dump.c
  ${suitesparse_SOURCE_DIR}/AMD/Source/amd_global.c
  ${suitesparse_SOURCE_DIR}/AMD/Source/amd_info.c
  ${suitesparse_SOURCE_DIR}/AMD/Source/amd_order.c
  ${suitesparse_SOURCE_DIR}/AMD/Source/amd_post_tree.c
  ${suitesparse_SOURCE_DIR}/AMD/Source/amd_postorder.c
  ${suitesparse_SOURCE_DIR}/AMD/Source/amd_preprocess.c
  ${suitesparse_SOURCE_DIR}/AMD/Source/amd_valid.c
  ${suitesparse_SOURCE_DIR}/SuiteSparse_config/SuiteSparse_config.c)

target_include_directories(amd
  PUBLIC
    "${suitesparse_SOURCE_DIR}/AMD/Include/"
    "${suitesparse_SOURCE_DIR}/SuiteSparse_config"
    )

install(TARGETS amd DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
install(FILES ${suitesparse_SOURCE_DIR}/AMD/Include/amd.h DESTINATION ${CMAKE_INSTALL_PREFIX}/include)
install(FILES ${suitesparse_SOURCE_DIR}/SuiteSparse_config/SuiteSparse_config.h DESTINATION ${CMAKE_INSTALL_PREFIX}/include)
