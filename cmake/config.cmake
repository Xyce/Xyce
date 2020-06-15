
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules")
include(CheckIncludeFileCXX)
include(CheckCXXSymbolExists)
include(cmake/feature_modes.cmake)
include(cmake/tps.cmake)
include(cmake/check_include_file_cxx.cmake)
include(cmake/CPackConfig.cmake)
