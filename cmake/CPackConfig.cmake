## Installation and CPack
# copyright, readme, license, etc.
install ( FILES ${Xyce_SOURCE_DIR}/distribution/README.TXT
            DESTINATION doc
            OPTIONAL COMPONENT Documentation )

# install optionally specified libraries (mainly to support Intel MKL issues)
install ( FILES ${Xyce_INSTALL_EXTRA_LIBS} DESTINATION lib COMPONENT Core )

#Use a build-appropriate license file:
if ( Xyce_RAD_MODELS )
  set ( CPACK_RESOURCE_FILE_LICENSE "${Xyce_SOURCE_DIR}/distribution/CPack.ECILicense.txt" )
  install ( FILES "${Xyce_SOURCE_DIR}/distribution/CPack.ECILicense.txt"
            DESTINATION doc
            OPTIONAL COMPONENT Documentation )
else()
  if ( Xyce_NONFREE_MODELS )
    set ( CPACK_RESOURCE_FILE_LICENSE "${Xyce_SOURCE_DIR}/distribution/CPack.NonFreeLicense.txt" )
    install ( FILES "${Xyce_SOURCE_DIR}/distribution/CPack.NonFreeLicense.txt"
            DESTINATION doc
            OPTIONAL COMPONENT Documentation )
  else ()
    set ( CPACK_RESOURCE_FILE_LICENSE "${Xyce_SOURCE_DIR}/distribution/CPack.OSLicense.txt" )
    install ( FILES "${Xyce_SOURCE_DIR}/distribution/CPack.OSLicense.txt"
            DESTINATION doc
            OPTIONAL COMPONENT Documentation )
  endif ()
endif ()

set ( CPACK_RESOURCE_FILE_README "${Xyce_SOURCE_DIR}/distribution/CPack.Description.txt" )

# set packaging variables
set( Xyce_BaseName "Xyce")
if (Xyce_RAD_MODELS)
  set( Xyce_BaseName "${Xyce_BaseName}Rad")
elseif (Xyce_NONFREE_MODELS)
  set( Xyce_BaseName "${Xyce_BaseName}")
else ()
  set( Xyce_BaseName "${Xyce_BaseName}OpenSource")
endif ()
if( Xyce_PARALLEL_MPI )
  set( Xyce_BaseName "${Xyce_BaseName}OMPI")
endif ()


set ( Xyce_INSTALL_NAME "${Xyce_BaseName}_${Xyce_VERSION_STRING_SHORT}" )
if( NOT DEFINED CPACK_PACKAGE_NAME)
  set ( CPACK_PACKAGE_NAME "${Xyce_INSTALL_NAME}" )
endif()
set ( CPACK_PACKAGE_DESCRIPTION_SUMMARY "Xyce Parallel Electronic Simulator" )
set ( CPACK_PACKAGE_VENDOR "Sandia National Laboratories" )
set ( CPACK_PACKAGE_VERSION_MAJOR "${Xyce_VERSION_MAJOR}" )
set ( CPACK_PACKAGE_VERSION_MINOR "${Xyce_VERSION_MINOR}" )
if ( Xyce_VERSION_PATCH)
  set ( CPACK_PACKAGE_VERSION_PATCH "${Xyce_VERSION_PATCH}" )
else()
  set ( CPACK_PACKAGE_VERSION_PATCH "0" )
endif()

# not all generators support this option.  In fact it may
# be only the Windows NSIS one that does.  So use with caution
set ( CPACK_PACKAGE_INSTALL_DIRECTORY "${Xyce_INSTALL_NAME}" )

# generator specific settings

if ( CMAKE_HOST_UNIX )
  # set default install location on unix to /usr/local
  set (CPACK_PACKAGING_INSTALL_PREFIX "/usr/local/${Xyce_INSTALL_NAME}")

  if ( CMAKE_HOST_APPLE )
    # OSX bundle directives
    if( NOT DEFINED CPACK_GENERATOR)
      set ( CPACK_GENERATOR "productbuild" )
    endif()
  else()
    # rpm directives if nothing specified
    SET( CPACK_GENERATOR "${GEN_TYPE}" )
  endif()

  set ( CPACK_RPM_PACKAGE_DESCRIPTION "For more information, visit http://xyce.sandia.gov ." )
  SET (CPACK_PACKAGE_DESCRIPTION_FILE "${Xyce_SOURCE_DIR}/distribution/CPack.Description.txt")
  set ( CPACK_RPM_PACKAGE_LICENSE "GPLv3" )
  SET( CPACK_RPM_PACKAGE_RELOCATABLE "true")
  SET( CPACK_DEBIAN_FILE_NAME "Xyce-${Xyce_VERSION_STRING_LONG}.deb" )
  SET( CPACK_DEBIAN_PACKAGE_ARCHITECHTURE "i386" )
  SET( CPACK_DEBIAN_PACKAGE_MAINTAINER "Sandia National Laboratories" )

  if ( Xyce_VERSION_EXTRA )
    set ( CPACK_RPM_PACKAGE_NAME "Xyce-${Xyce_VERSION_EXTRA}" )
    SET( CPACK_DEBIAN_PACKAGE_NAME "Xyce-${Xyce_VERSION_EXTRA}")
    set(CPACK_ARCHIVE_FILE_NAME "Xyce-${Xyce_VERSION_EXTRA}")
  else()
    set ( CPACK_RPM_PACKAGE_NAME "Xyce" )
    SET( CPACK_DEBIAN_PACKAGE_NAME "Xyce")
    set(CPACK_ARCHIVE_FILE_NAME "Xyce")
  endif()
  
endif()


if ( CMAKE_HOST_WIN32 )

  # NSIS directives
  set ( CPACK_GENERATOR "NSIS" )
  set ( CPACK_NSIS_DISPLAY_NAME "${Xyce_INSTALL_NAME}" )
  set ( CPACK_NSIS_HELP_LINK "http:\\\\\\\\xyce.sandia.gov\\\\contact_us.html" )
  set ( CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\xyce.sandia.gov")
  set ( CPACK_NSIS_CONTACT "xyce-support@sandia.gov" )
  set ( CPACK_NSIS_MODIFY_PATH OFF )

  #
  # All of the following extra install() calls
  # are probably not needed with the use of 
  # install(TARGETS Xyce ... RUNTIME_DEPENDENCIES ) in Xyce/src/CMakeLists.txt
  # However removal of these lines will need to be tested under Windows first
  # All lines from here to just before the "registry settings" line
  #
  
  # there is probably a better way to do this with
  # file(GET_RUNTIME_DEPENDENCIES...)
  # but that requires cmake 3.17.3 or higher.  For now
  # rely on the old method but we should transition to the
  # new way when more recent versions of cmake are available.

  file(TO_CMAKE_PATH "$ENV{ONEAPI_ROOT}" ONEAPI_ROOT)

  # Add the required intel library.  This is permitted for redistribution.
  find_file(SVML_PATH "svml_dispmd.dll" PATHS "${ONEAPI_ROOT}/bin" "${ONEAPI_ROOT}/compiler/latest/windows/redist/intel64_win/compiler" NO_DEFAULT_PATH)
  find_file(MMD_PATH "libmmd.dll" PATHS "${ONEAPI_ROOT}/bin" "${ONEAPI_ROOT}/compiler/latest/windows/redist/intel64_win/compiler" NO_DEFAULT_PATH)
  find_file(IOMP5_PATH "libiomp5md.dll" PATHS "${ONEAPI_ROOT}/bin" "${ONEAPI_ROOT}/compiler/latest/windows/redist/intel64_win/compiler" NO_DEFAULT_PATH)

  if(IS_SYMLINK ${SVML_PATH})
    file(READ_SYMLINK ${SVML_PATH} SVML_PATH)
    file(REAL_PATH ${SVML_PATH} SVML_PATH BASE_DIRECTORY "${ONEAPI_ROOT}/bin")
    set(SVML_PATH ${SVML_PATH} CACHE STRING "Converted from symbolic link to real path in ${CMAKE_CURRENT_LIST_FILE}" FORCE)
  endif()

  if(IS_SYMLINK ${MMD_PATH})
    file(READ_SYMLINK ${MMD_PATH} MMD_PATH)
    file(REAL_PATH ${MMD_PATH} MMD_PATH BASE_DIRECTORY "${ONEAPI_ROOT}/bin")
    set(MMD_PATH ${MMD_PATH} CACHE STRING "Converted from symbolic link to real path in ${CMAKE_CURRENT_LIST_FILE}" FORCE)
  endif()

  if(IS_SYMLINK ${IOMP5_PATH})
    file(READ_SYMLINK ${IOMP5_PATH} IOMP5_PATH)
    file(REAL_PATH ${IOMP5_PATH} IOMP5_PATH BASE_DIRECTORY "${ONEAPI_ROOT}/bin")
    set(IOMP5_PATH ${IOMP5_PATH} CACHE STRING "Converted from symbolic link to real path in ${CMAKE_CURRENT_LIST_FILE}" FORCE)
  endif()

  install ( FILES ${SVML_PATH} DESTINATION bin COMPONENT Runtime)  
  install ( FILES ${MMD_PATH} DESTINATION bin COMPONENT Runtime)
  install ( FILES ${IOMP5_PATH} DESTINATION bin COMPONENT Runtime)

  # For native Windows builds we also need assorted MS Visual Studio DLLs
  # This thing apparently takes care of it:
  include(InstallRequiredSystemLibraries)

  # registry settings
  set ( CPACK_PACKAGE_INSTALL_REGISTRY_KEY "${Xyce_INSTALL_NAME}" )

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

endif()
