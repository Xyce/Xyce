#
# a CMake script to compile verilogADMS 
# files for use in Xyce as a plugin library
#
# Requirements:
#  1. Xyce must be installed on the machine where this is run.
#  2. Xyce must have been built with plugin support enabled
#
#  If Xyce was built and installed properly, then
#  these required files will be in place
#
#    <Path to Xyce Install>/bin/admsXml
#    <Path to Xyce Install>/share/xyceAnalogFunction_nosac.xml
#    <Path to Xyce Install>/share/xxyceBasicTemplates_nosac.xml
#    <Path to Xyce Install>/xyceHeaderFile_nosac.xml
#    <Path to Xyce Install>/xyceImplementationFile_nosac.xml
#    <Path to Xyce Install>/xyceVersion_nosac.xml
#
# How to Run this Script
# 
# cmake -DVERILOG_FILES="file1.va;file2.va;..."  
#
# or i
#
# cmake -DCMAKE_PREFIX_PATH=<Path to Xyce Install> 
#
# Optional variables
#

cmake_minimum_required(VERSION 3.20)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(xyce_plugin)

# Xyce is required 
find_package(Xyce)

# Determine library name.
if(DEFINED PLUGIN_NAME)
  set(library_name "${PLUGIN_NAME}")
else()
  set(library_name "xyceplugin")
endif()

message(STATUS "Verilog files = ${VERILOG_FILES}")
message(STATUS "Xyce_DIR = ${Xyce_DIR}")
file(REAL_PATH "${Xyce_DIR}/../../../bin/" Xyce_BinDir )
file(REAL_PATH "${Xyce_DIR}/../../../share/" Xyce_XmlDir )
message(STATUS "Xyce_BinDir = ${Xyce_BinDir}")
message(STATUS "Xyce_XmlDir = ${Xyce_XmlDir}")
find_file( AdmsXMLCommand NAMES "admsXml" "admsXml.exe" PATHS ${Xyce_BinDir} REQUIRED )
message(STATUS "admsXml = ${AdmsXMLCommand}")

set(ModuleNames)
set(SourceFiles)
set(SourceHeaders)
# Loop over the input verilog files.
# to get source file names
foreach( vafile ${VERILOG_FILES} )
  # need to read the module name from each va input file
  file(READ ${vafile} vaFileText )
  string(REGEX MATCH "module[ \t]+[^ \t\r\n(]+" moduleAndName ${vaFileText} )
  string(REGEX REPLACE "module[ \t]+" "" module ${moduleAndName} )
  #string(REPLACE ".va" "" module "${vafile}" ) 
  string(CONCAT rawSource "N_DEV_ADMS" "${module}" ".C")
  string(CONCAT rawHeader "N_DEV_ADMS" "${module}" ".h")
  list(APPEND ModuleNames ${module})
  list(APPEND SourceFiles ${rawSource})
  list(APPEND SourceHeaders ${rawHeader})
  message( STATUS "working on ${rawSource} and ${rawHeader}")
  add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${rawSource} ${CMAKE_CURRENT_BINARY_DIR}/${rawHeader}
    COMMAND ${AdmsXMLCommand} 
            -D__XYCE__ -x 
            -e ${Xyce_XmlDir}/adms.implicit.xml 
            -e ${Xyce_XmlDir}/xyceVersion_nosac.xml  
            -e ${Xyce_XmlDir}/xyceBasicTemplates_nosac.xml 
            -e ${Xyce_XmlDir}/xyceAnalogFunction_nosac.xml 
            -e ${Xyce_XmlDir}/xyceHeaderFile_nosac.xml 
            -e ${Xyce_XmlDir}/xyceImplementationFile_nosac.xml ${vafile} 
    DEPENDS ${vaFile}
  )
endforeach()

message( STATUS  "ModuleNames ${ModuleNames}" )
message( STATUS  "SourceFiles ${SourceFiles}" )
message( STATUS  "SourceHeaders ${SourceHeaders}" )

# make the bootstrap code for the shared library
string( CONCAT BootstrapFileName "N_DEV_Module" "${library_name}" ".C" )
file( WRITE ${BootstrapFileName} "#include<Xyce_config.h>\n" )
foreach( header ${SourceHeaders} )
  file( APPEND ${BootstrapFileName} "#include<" "${header}" ">\n" )
endforeach()


file( APPEND ${BootstrapFileName}  "struct Bootstrap\n" )
file( APPEND ${BootstrapFileName}  "{\n" )
file( APPEND ${BootstrapFileName}  "  Bootstrap()\n" )
file( APPEND ${BootstrapFileName}  "  {\n" )
foreach( module ${ModuleNames})
  file( APPEND ${BootstrapFileName}  "    Xyce::Device::ADMS" "${module}" "::registerDevice();\n" )
endforeach()

file( APPEND ${BootstrapFileName}  "  }\n" )
file( APPEND ${BootstrapFileName}  "};\n" )
file( APPEND ${BootstrapFileName}  "Bootstrap s_bootstrap;\n" )
# done making the bootstrap file.

add_library("${library_name}" SHARED ${SourceFiles} ${BootstrapFileName} )
target_include_directories( "${library_name}" PUBLIC ${CMAKE_CURRENT_BINARY_DIR} )
target_link_libraries( "${library_name}" Xyce::XyceLib)

#target_include_directories("${library_name}")
#target_link_libraries("${library_name}" xyce)



