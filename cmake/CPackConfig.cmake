## Installation and CPack
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

  # rpm directives if nothing specified
    SET( CPACK_GENERATOR "${GEN_TYPE}" )

  endif ( CMAKE_HOST_APPLE )

  set ( CPACK_RPM_PACKAGE_DESCRIPTION "For more information, visit http://xyce.sandia.gov ." )
  set ( CPACK_RPM_PACKAGE_LICENSE "GPLv3" )
  SET( CPACK_DEBIAN_FILE_NAME "Xyce-${Xyce_VERSION_STRING_LONG}.deb" )
  SET( CPACK_DEBIAN_PACKAGE_ARCHITECHTURE "i386" )
  SET( CPACK_DEBIAN_PACKAGE_MAINTAINER "Sandia National Laboratories" )

  if ( Xyce_VERSION_EXTRA )

    set ( CPACK_RPM_PACKAGE_NAME "Xyce-${Xyce_VERSION_EXTRA}" )
    SET( CPACK_DEBIAN_PACKAGE_NAME "Xyce-${Xyce_VERSION_EXTRA}")
    set(CPACK_ARCHIVE_FILE_NAME "Xyce-${Xyce_VERSION_EXTRA}")
  else ( Xyce_VERSION_EXTRA )

    set ( CPACK_RPM_PACKAGE_NAME "Xyce" )
    SET( CPACK_DEBIAN_PACKAGE_NAME "Xyce")
    set(CPACK_ARCHIVE_FILE_NAME "Xyce")

  endif ( Xyce_VERSION_EXTRA )

endif ( CMAKE_HOST_UNIX )


if ( CMAKE_HOST_WIN32 )

  # NSIS directives
  set ( CPACK_GENERATOR "NSIS" )
  set ( CPACK_NSIS_DISPLAY_NAME "${Xyce_INSTALL_NAME}" )
  set ( CPACK_NSIS_HELP_LINK "http:\\\\\\\\xyce.sandia.gov\\\\contact_us.html" )
  set ( CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\xyce.sandia.gov")
  set ( CPACK_NSIS_CONTACT "xyce-support@sandia.gov" )
  set ( CPACK_NSIS_MODIFY_PATH OFF )

  # there is probably a better way to do this with
  # file(GET_RUNTIME_DEPENDENCIES...)
  # but that requires cmake 3.17.3 or higher.  For now
  # rely on the old method but we should transition to the
  # new way when more recent versions of cmake are available.

  # Add the required intel library.  This is permitted for redistribution.
  find_file(SVNL_PATH
    "svml_dispmd.dll"
    PATHS
      "$ENV{INTEL_DEV_REDIST}redist/intel64/compiler"
      "$ENV{INTEL_DEV_REDIST}redist/intel64_win/compiler" )
  find_file(MMD_PATH
      "libmmd.dll"
      PATHS
        "$ENV{INTEL_DEV_REDIST}redist/intel64/compiler"
        "$ENV{INTEL_DEV_REDIST}redist/intel64_win/compiler" )
  if ( MMD_PATH AND SVNL_PATH)
    # For native Windows builds, check for these two dll's.  If they don't exist,
    # assume the package is being build on the legacy Windows build system and use that
    # dll and path instead.
    install ( FILES ${MMD_PATH}
	    DESTINATION bin
	    COMPONENT Runtime)
    install ( FILES ${SVNL_PATH}
	    DESTINATION bin
	    COMPONENT Runtime)
  else ()
    set ( INTELLIBPATH "$ENV{ICPP_COMPILER12_CYGWIN}/redist/ia32/compiler" )
    install ( FILES ${INTELLIBPATH}/libmmd.dll
	    DESTINATION bin
	    COMPONENT Runtime)
  endif()


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

endif ( CMAKE_HOST_WIN32 )
