
# -- locate sources ---------------------------------------------------------- 

## DEBUG:  the order matters to final link! do not rearrange
# add standard top level source directories to list
set ( SOURCE_DIRS  
      AnalysisPKG 
      CircuitPKG
      DakotaLinkPKG      
      IOInterfacePKG 
      MultiTimePDEPKG 
      DeviceModelPKG 
      ErrorHandlingPKG 
      LoaderServicesPKG 
      NonlinearSolverPKG 
      LinearAlgebraServicesPKG 
      ParallelDistPKG 
      TimeIntegrationPKG 
      TopoManagerPKG 
      UtilityPKG )

# collect header and library information recursively
foreach( SRCDIR ${SOURCE_DIRS} )
  add_subdirectory( ${SRCDIR} )
endforeach( SRCDIR ${SOURCE_DIRS} )


# prune duplicate entries
list ( REMOVE_DUPLICATES XYCE_INCLUDE_DIRECTORIES )
list ( REMOVE_DUPLICATES XYCE_PACKAGE_LIBRARIES )


# add header directories to -Ilist 
include_directories ( ${XYCE_INCLUDE_DIRECTORIES} )


# explicity mark generated files
if ( XYCE_GENERATED_SOURCES )
  set_source_files_properties ( ${XYCE_GENERATED_SOURCES} 
                                PROPERTIES GENERATED TRUE )
endif ( XYCE_GENERATED_SOURCES )




# -- build package libraries (libMPDE.so, libMPDE.a, etc.) -------------------
# -- override optimiztion options for specific files (for debugging)
                                   
foreach ( PACKAGE_LIBRARY ${XYCE_PACKAGE_LIBRARIES} )


 ##############################################################################
 ## This is commented out but was used for a very long time because
 ## we used to build Xyce on a tiny virtual machine and the Intel compiler
 ## would crash trying to optimize ADMS devices.
 ## We now run on a beefier physical machine and should no longer need to
 ## downgrade the optimization like this:
 #if ( ${PACKAGE_LIBRARY} STREQUAL "ADMS" AND 
 #     COMPILER_VENDOR STREQUAL "Intel" )
 #  message ( STATUS "Changing optimization to -O1 for ADMS sources..." )
 #  set_source_files_properties ( ${${PACKAGE_LIBRARY}_SOURCES} 
 #                                 PROPERTIES COMPILE_FLAGS "-O1" )
 #endif ( ${PACKAGE_LIBRARY} STREQUAL "ADMS" AND 
 #        COMPILER_VENDOR STREQUAL "Intel" )
 ##############################################################################


 ## This is commented out but it can be useful in debugging problems where the
 ## compliler over optimizes something.  The code below lets one change the
 ## optimization level of individual files.
 ##  if ( ${PACKAGE_LIBRARY} STREQUAL "DeviceModel" AND 
 ##       COMPILER_VENDOR STREQUAL "Intel" )
 ##    
 ##    foreach ( DEV_SRC_FILE ${${PACKAGE_LIBRARY}_SOURCES} )
 ##      if ( ${DEV_SRC_FILE} MATCHES N_DEV_MutIndNonLin )
 ##       message ( WARNING "Changing optimization to -O1 for Nonlinear mutual inductor sources..." )
 ##        set_source_files_properties ( ${DEV_SRC_FILE} 
 ##                                   PROPERTIES COMPILE_FLAGS "-O1" )
 ##      endif ( ${DEV_SRC_FILE} MATCHES N_DEV_MutIndNonLin )
 ##    endforeach ( DEV_SRC_FILE ${${PACKAGE_LIBRARY}_SOURCES} )
 ##    
 ##  endif ( ${PACKAGE_LIBRARY} STREQUAL "DeviceModel" AND 
 ##          COMPILER_VENDOR STREQUAL "Intel" )
  ############################################################################


  #### DEBUG:  use this block for monolithic libxyce #########################
  set ( PACKAGE_SOURCES ${PACKAGE_SOURCES} ${${PACKAGE_LIBRARY}_SOURCES} )
  ############################################################################

  #### DEBUG:  use this block for subpackage libraries #######################
  # if ( Xyce_ENABLE_STATIC ) 
  #   add_library ( ${PACKAGE_LIBRARY}_static STATIC ${${PACKAGE_LIBRARY}_SOURCES} ) 
  #   set_target_properties ( ${PACKAGE_LIBRARY}_static PROPERTIES 
  #                           OUTPUT_NAME ${PACKAGE_LIBRARY} 
  #                           CLEAN_DIRECT_OUTPUT ON )
  # endif ( Xyce_ENABLE_STATIC ) 
  #
  # if ( Xyce_ENABLE_SHARED ) 
  #   add_library ( ${PACKAGE_LIBRARY}_shared SHARED ${${PACKAGE_LIBRARY}_SOURCES} ) 
  #   set_target_properties ( ${PACKAGE_LIBRARY}_shared PROPERTIES 
  #                           OUTPUT_NAME ${PACKAGE_LIBRARY} 
  #                           CLEAN_DIRECT_OUTPUT ON )
  #   DEBUG:  add subpackage target dependencies
  # endif ( Xyce_ENABLE_SHARED ) 
  #
  # DEBUG:  add code to enable installation of subpackage libs
  ############################################################################

endforeach ( PACKAGE_LIBRARY ${XYCE_PACKAGE_LIBRARIES} )




# -- build targets -----------------------------------------------------------

# create main library (libxyce.a or libxyce.so)
if ( Xyce_ENABLE_STATIC ) 

  #### DEBUG:  use this block for monolithic libxyce.a #######################
  add_library( lib_xyce_static STATIC ${Circuit_SOURCES} ${PACKAGE_SOURCES} ) 
  target_link_libraries ( lib_xyce_static 
                          ${ARCH_LIBRARIES} 
                          ${EXTERNAL_LIBRARIES} )
  ############################################################################

  ##### DEBUG:  use this block for subpackage libraries ######################
  # add_library( lib_xyce_static STATIC ${Circuit_SOURCES} ) 
  # target_link_libraries ( lib_xyce_static 
  #                         ${ARCH_LIBRARIES} 
  #                         ${XYCE_PACKAGE_LIBRARIES} 
  #                         ${EXTERNAL_LIBRARIES} )
  ############################################################################

  # rename library 
  set_target_properties ( lib_xyce_static PROPERTIES 
                          OUTPUT_NAME xyce          
                          CLEAN_DIRECT_OUTPUT ON )

  # bind top level subpackage dependencies
  if ( XYCE_TARGET_DEPENDENCIES )
    add_dependencies ( lib_xyce_static ${XYCE_TARGET_DEPENDENCIES} )
  endif ( XYCE_TARGET_DEPENDENCIES )

  # installation support
  install ( TARGETS lib_xyce_static 
            ARCHIVE DESTINATION lib )

endif ( Xyce_ENABLE_STATIC ) 


if ( Xyce_ENABLE_SHARED ) 

  #### DEBUG:  use this block for monolithic libxyce.so ######################
  add_library( lib_xyce_shared SHARED ${Circuit_SOURCES} ${PACKAGE_SOURCES} ) 
  target_link_libraries ( lib_xyce_shared 
                          ${ARCH_LIBRARIES} 
                          ${EXTERNAL_LIBRARIES} )
  
  ############################################################################

  #### DEBUG:  use this block for subpackage libraries #######################
  # add_library( lib_xyce_shared SHARED ${Circuit_SOURCES} ) 
  # target_link_libraries ( lib_xyce_shared 
  #                         ${ARCH_LIBRARIES} 
  #                         ${XYCE_PACKAGE_LIBRARIES} 
  #                         ${EXTERNAL_LIBRARIES} )
  ############################################################################

  # rename library 
  set_target_properties ( lib_xyce_shared PROPERTIES 
                          OUTPUT_NAME xyce 
                          CLEAN_DIRECT_OUTPUT ON )
                          
  # set RPATH to be used on MacOSX builds.  DYLD_LIBRARY_PATH AND LD_LIBRARY_PATH
  # are no longer supported on MacOSX.
  set_target_properties( lib_xyce_shared PROPERTIES MACOSX_RPATH True)
  
  # bind top level subpackage dependencies
  if ( XYCE_TARGET_DEPENDENCIES )
    add_dependencies ( lib_xyce_shared ${XYCE_TARGET_DEPENDENCIES} )
  endif ( XYCE_TARGET_DEPENDENCIES )

  # installation support
  install ( TARGETS lib_xyce_shared 
            LIBRARY DESTINATION lib 
            RUNTIME DESTINATION bin )

endif ( Xyce_ENABLE_SHARED ) 


# create executable (Xyce or Xyce.exe)
if ( NOT Xyce_DISABLE_BINARY )

  if ( Xyce_STATIC_BINARY )

    add_executable( Xyce_STATIC ${Xyce_SOURCE_DIR}/src/Xyce.C )
    target_link_libraries( Xyce_STATIC lib_xyce_static )

    # rename binary and enable static linking
    set_target_properties ( Xyce_STATIC PROPERTIES 
                            OUTPUT_NAME Xyce 
                            LINK_FLAGS "-static" 
                            LINK_SEARCH_END_STATIC ON )

    # installation support
    install ( TARGETS Xyce_STATIC 
              RUNTIME DESTINATION bin )

  else ( Xyce_STATIC_BINARY )

    add_executable( Xyce_SHARED ${Xyce_SOURCE_DIR}/src/Xyce.C )

    # prefer shared libs over static
    if ( Xyce_ENABLE_SHARED )
      target_link_libraries( Xyce_SHARED lib_xyce_shared )
    else ( Xyce_ENABLE_SHARED )
      target_link_libraries( Xyce_SHARED lib_xyce_static )
    endif ( Xyce_ENABLE_SHARED )

    # rename binary 
    set_target_properties ( Xyce_SHARED PROPERTIES OUTPUT_NAME Xyce ) 

    # installation support
    install ( TARGETS Xyce_SHARED 
              RUNTIME DESTINATION bin )

  endif ( Xyce_STATIC_BINARY )

endif ( NOT Xyce_DISABLE_BINARY )




# -- installation support ----------------------------------------------------

# header files
install ( FILES ${Xyce_SOURCE_DIR}/src/CircuitPKG/N_CIR_Xyce.h 
                ${Xyce_SOURCE_DIR}/src/IOInterfacePKG/N_IO_CmdParse.h 
                ${Xyce_SOURCE_DIR}/src/UtilityPKG/N_UTL_Xyce.h 
                ${Xyce_BINARY_DIR}/src/Xyce_config.h 
          DESTINATION include )


# documentation
install ( FILES ${Xyce_SOURCE_DIR}/doc/Reference_Guide/Xyce_RG.pdf
          DESTINATION doc 
          RENAME Reference_Guide.pdf 
          OPTIONAL )

install ( FILES ${Xyce_SOURCE_DIR}/doc/Users_Guide/Xyce_UG.pdf
          DESTINATION doc 
          RENAME Users_Guide.pdf 
          OPTIONAL )
         
  install ( FILES ${Xyce_SOURCE_DIR}/doc/Release_Notes_${Xyce_RELEASE_NOTES_VERSION}/Release_Notes_${Xyce_RELEASE_NOTES_VERSION}.pdf
          DESTINATION doc 
          OPTIONAL )

if ( Xyce_RAD_MODELS )
  install ( FILES ${Xyce_SOURCE_DIR}/doc/SandiaModelsDocs/Reference_Guide/Xyce_Sandia_Models.pdf
            DESTINATION doc 
            RENAME Reference_Guide_Sandia_Models.pdf
            OPTIONAL )
endif ( Xyce_RAD_MODELS ) 

# DEBUG build documentation instead of just including here?
# DEBUG: forcing regression path to documentation (remove later)
install ( FILES ${Xyce_BINARY_DIR}/../doc/Reference_Guide.pdf
                ${Xyce_BINARY_DIR}/../doc/Users_Guide.pdf 
		${Xyce_BINARY_DIR}/../doc/Release_Notes_${Xyce_RELEASE_NOTES_VERSION}.pdf
          DESTINATION doc 
          OPTIONAL )
if ( Xyce_RAD_MODELS )
  install ( FILES ${Xyce_BINARY_DIR}/../doc/Reference_Guide_Sandia_Models.pdf
            DESTINATION doc 
            OPTIONAL )
endif ( Xyce_RAD_MODELS )

#This is supposed to make Windows builds include appropriate system dlls
# IT DOES NOT WORK FOR XYCE, because on the machine where we build Xyce
# there is no full install of Visual Studio, and this thing only looks in
# the place where VS keeps the redistributable libraries.  VS Express doesn't
# have that directory, so we must look elsewhere.
include(InstallRequiredSystemLibraries)

# copyright, readme, license, etc.
install ( FILES ${Xyce_SOURCE_DIR}/distribution/README.TXT
            DESTINATION doc 
            OPTIONAL )


#Use a build-appropriate license file:
if ( Xyce_RAD_MODELS )
  set ( CPACK_RESOURCE_FILE_LICENSE "${Xyce_SOURCE_DIR}/distribution/CPack.ECILicense.txt" )
else ( Xyce_RAD_MODELS )
  if ( Xyce_NONFREE_MODELS )
    set ( CPACK_RESOURCE_FILE_LICENSE "${Xyce_SOURCE_DIR}/distribution/CPack.NonFreeLicense.txt" )
  else ( Xyce_NONFREE_MODELS )
    set ( CPACK_RESOURCE_FILE_LICENSE "${Xyce_SOURCE_DIR}/distribution/CPack.OSLicense.txt" )
  endif ( Xyce_NONFREE_MODELS )
endif ( Xyce_RAD_MODELS )  

set ( CPACK_RESOURCE_FILE_README "${Xyce_SOURCE_DIR}/distribution/CPack.Description.txt" )



# set packaging variables
set ( Xyce_INSTALL_NAME "Xyce ${Xyce_VERSION_STRING_LONG}" )
set ( CPACK_PACKAGE_NAME "${Xyce_INSTALL_NAME}" )  
set ( CPACK_PACKAGE_DESCRIPTION_SUMMARY "Xyce Parallel Electronic Simulator" )
set ( CPACK_PACKAGE_VENDOR "Sandia National Laboratories" )
set ( CPACK_PACKAGE_VERSION_MAJOR "${Xyce_VERSION_MAJOR}" )
set ( CPACK_PACKAGE_VERSION_MINOR "${Xyce_VERSION_MINOR}" )
if ( Xyce_VERSION_PATCH)
  set ( CPACK_PACKAGE_VERSION_PATCH "${Xyce_VERSION_PATCH}" )
else ( Xyce_VERSION_PATCH)
  set ( CPACK_PACKAGE_VERSION_PATCH "0" )
endif ( Xyce_VERSION_PATCH)

set ( CPACK_PACKAGE_INSTALL_DIRECTORY "${Xyce_INSTALL_NAME}" )

# generator specific settings


if ( CMAKE_HOST_UNIX )

if ( CMAKE_HOST_APPLE )

  # OSX bundle directives
  set ( CPACK_GENERATOR "PackageMaker" )

else ( CMAKE_HOST_APPLE )

  # rpm directives
  set ( CPACK_GENERATOR "RPM" )
endif ( CMAKE_HOST_APPLE )

  ### unset ( CPACK_PACKAGE_FILE_NAME )
  ### unset ( CPACK_PACKAGE_FILE_NAME CACHE )


  set ( CPACK_RPM_PACKAGE_DESCRIPTION "For more information, visit http://xyce.sandia.gov ." )
  set ( CPACK_RPM_PACKAGE_GROUP "Applications/Engineering" )
  set ( CPACK_RPM_PACKAGE_LICENSE "GPLv3" )

  if ( Xyce_VERSION_EXTRA )
    set ( CPACK_RPM_PACKAGE_NAME "Xyce-${Xyce_VERSION_EXTRA}" )
  else ( Xyce_VERSION_EXTRA )
    set ( CPACK_RPM_PACKAGE_NAME "Xyce" )
  endif ( Xyce_VERSION_EXTRA )


### FORCE MKL ###
#get_property ( IS64BIT GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS )


# NOTE: the use of "IS64BIT" above to determine if we should install
# 64-bit libraries is totally wrong.  FIND_LIBRARY_USE_LIB64_PATHS is almost
# always true, and means "search 64-bit library paths to resolve libraries
# if needed".  It does not mean "the compiler is actually using 64-bit 
# libraries".
# 
# We only use 32-bit on ONE platform anymore, and on that platform we always
# run a standard script that sets up the environment for Intel compiler
# builds.  So we just add an "export INTEL_COMPILER_32BIT=TRUE" to that
# script, and punt the autodetection.
#
# Ick.
#

if ( "$ENV{INTEL_COMPILER_32BIT}" STREQUAL "TRUE" )
#    message( WARNING "INTEL_COMPILER_32BIT is $ENV{INTEL_COMPILER_32BIT}.  Installing 32-bit libraries")
    set( MKLLIBPATH "$ENV{MKLROOT}/lib/32" )
    set( INTELLIBPATH "$ENV{INTEL_PATH}/cc/current/lib/ia32" )       
else ( "$ENV{INTEL_COMPILER_32BIT}" STREQUAL "TRUE" )
#   message( WARNING "INTEL_COMPILER_32BIT is $ENV{INTEL_COMPILER_32BIT}.  Installing 64-bit libraries")
    set( MKLLIBPATH "$ENV{MKLROOT}/lib/em64t" )
    if (CMAKE_HOST_APPLE)
        set( INTELLIBPATH "$ENV{INTEL_PATH}/$ENV{INTEL_COM_VER}/lib" )   
        set( MKLLIBPATH "$ENV{MKLROOT}/lib" )
    else (CMAKE_HOST_APPLE)

        set( INTELLIBPATH "$ENV{INTEL_PATH}/cce/current/lib/intel64" )   
    endif (CMAKE_HOST_APPLE)
endif ( "$ENV{INTEL_COMPILER_32BIT}" STREQUAL "TRUE" )

if (CMAKE_HOST_APPLE)
  set ( CPACK_PACKAGING_INSTALL_PREFIX "/" )
else (CMAKE_HOST_APPLE)
  set ( CPACK_PACKAGING_INSTALL_PREFIX "/usr/local/${CPACK_RPM_PACKAGE_NAME}-${Xyce_VERSION_STRING_SHORT}" )
endif (CMAKE_HOST_APPLE)

set ( CPACK_PACKAGE_DEFAULT_LOCATION "/usr/local/${CPACK_RPM_PACKAGE_NAME}-${Xyce_VERSION_STRING_SHORT}" ) 
if ( COMPILER_VENDOR STREQUAL "Intel" )
  if (CMAKE_HOST_APPLE)
   install ( DIRECTORY ${MKLLIBPATH}/ ${INTELLIBPATH}/  ${INTELLIBPATH}/intel64/
            DESTINATION lib 
            FILES_MATCHING PATTERN "*.dylib" )
  else (CMAKE_HOST_APPLE)
    install ( DIRECTORY ${MKLLIBPATH}/ ${INTELLIBPATH}/ 
            DESTINATION lib 
            FILES_MATCHING PATTERN "*.so" )
  endif (CMAKE_HOST_APPLE)
endif ( COMPILER_VENDOR STREQUAL "Intel" )
### ##### ### ###

  set ( CPACK_RPM_PACKAGE_VERSION "${Xyce_VERSION_SHORT}" )

  set ( CPACK_STRIP_FILES "bin/Xyce" )
  set ( CPACK_SOURCE_STRIP_FILES "" )

endif ( CMAKE_HOST_UNIX )


if ( CMAKE_HOST_WIN32 )

  # NSIS directives
  set ( CPACK_GENERATOR "NSIS" )
  set ( CPACK_NSIS_DISPLAY_NAME "${Xyce_INSTALL_NAME}" )
  set ( CPACK_NSIS_HELP_LINK "http:\\\\\\\\xyce.sandia.gov\\\\contact_us.html" )
  set ( CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\xyce.sandia.gov")
  set ( CPACK_NSIS_CONTACT "xyce-support@sandia.gov" )
  set ( CPACK_NSIS_MODIFY_PATH OFF )

  # registry settings
  set ( CPACK_PACKAGE_INSTALL_REGISTRY_KEY "${Xyce_INSTALL_NAME}" )

  # Add the required intel library.  This is permitted for redistribution.
  set ( INTELLIBPATH "$ENV{INTEL_DEV_REDIST}redist/intel64/compiler" )
  string ( REPLACE "\\" "/" INTELLIBPATH ${INTELLIBPATH} )
  if ( EXISTS "${INTELLIBPATH}/libmmd.dll" AND EXISTS "${INTELLIBPATH}/svml_dispmd.dll" )
    # For native Windows builds, check for these two dll's.  If they don't exist,
    # assume the package is being build on the legacy Windows build system and use that
    # dll and path instead.
    install ( FILES ${INTELLIBPATH}/libmmd.dll
	    DESTINATION bin
	    COMPONENT Runtime)
    install ( FILES ${INTELLIBPATH}/svml_dispmd.dll
	    DESTINATION bin
	    COMPONENT Runtime)
  else ( EXISTS "${INTELLIBPATH}/libmmd.dll" AND EXISTS "${INTELLIBPATH}/svml_dispmd.dll" )
    set ( INTELLIBPATH "$ENV{ICPP_COMPILER12_CYGWIN}/redist/ia32/compiler" )
    install ( FILES ${INTELLIBPATH}/libmmd.dll
	    DESTINATION bin
	    COMPONENT Runtime)
  endif ( EXISTS "${INTELLIBPATH}/libmmd.dll" AND EXISTS "${INTELLIBPATH}/svml_dispmd.dll" )

  #HACK:  This *should* be picked up by InstallRequiredSystemLibraries,
  # but is not because we have the "express" edition of MSVS --- that
  # does not come with the redistributable packet.  According to
  # posts from Microsoft on their msdn forum, a work-around for express
  # users is to take the frile from the windows/system32 or sysWOW64
  # directories (for x86 and 64-bit, respectively)

  install (FILES "C:/Windows/system32/msvcp100.dll" "C:/Windows/system32/msvcr100.dll"
	DESTINATION bin
	COMPONENT Runtime)

  # Start Menu entries
  set ( CPACK_NSIS_MENU_LINKS "doc" "Documentation" 
                              "http://xyce.sandia.gov" "${CPACK_PACKAGE_DESCRIPTION_SUMMARY}" )
  set ( CPACK_NSIS_CREATE_ICONS_EXTRA
        " CreateShortCut \\\"$SMPROGRAMS\\\\$STARTMENU_FOLDER\\\\Xyce ${Xyce_VERSION_STRING_LONG} Command Prompt.lnk\\\" \\\"%comspec%\\\" \\\"/k path=%PATH%;$INSTDIR\\\\bin\\\" 
          CreateShortCut  \\\"$DESKTOP\\\\Xyce ${Xyce_VERSION_STRING_LONG} Command Prompt.lnk\\\" \\\"%comspec%\\\" \\\"/k path=%PATH%;$INSTDIR\\\\bin\\\" " ) 

  set ( CPACK_NSIS_DELETE_ICONS_EXTRA
        " Delete \\\"$SMPROGRAMS\\\\$MUI_TEMP\\\\Xyce ${Xyce_VERSION_STRING_LONG} Command Prompt.lnk\\\"
          Delete \\\"$DESKTOP\\\\Xyce ${Xyce_VERSION_STRING_LONG} Command Prompt.lnk\\\" "
)

endif ( CMAKE_HOST_WIN32 )




# generate CPack configuration
include ( CPack )



# -- testing support ---------------------------------------------------------

# debug... move to ctest driven approach
## if ( DEFINED Xyce_SRCTEST )
  add_subdirectory ( test )
## endif ( DEFINED Xyce_SRCTEST )

