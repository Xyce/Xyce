## With CMake 3.18, find_package(Python3) supports embedding
#find_package (Python3 COMPONENTS Development.Embed QUIET)
option(Xyce_PYMI "Build Python Model Interpreter (Xyce-PyMi)?" OFF)
if (Xyce_PYMI)
  find_package( Python COMPONENTS Development Interpreter QUIET )
  if ( EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/pybind11 )
    set(PYBIND11_INSTALL ON)
    add_subdirectory( pybind11 EXCLUDE_FROM_ALL )
    list(PREPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_BINARY_DIR}/pybind11")
    find_package(pybind11_dummy 2.6 NAMES pybind11 CONFIG QUIET)
    if (NOT pybind11_dummy_VERSION)
      message(FATAL_ERROR "Xyce_PYMI=ON and pybind11 found in source tree, but pybind11 version <2.6")
    endif()
  else()
    find_package( pybind11 2.6 CONFIG QUIET)
    if (NOT pybind11_FOUND)
      message(FATAL_ERROR "Xyce_PYMI=ON but pybind11 NOT found")
    endif()
  endif()
  
  if (pybind11_FOUND)
    if((NOT Python_FOUND) OR (Python_VERSION_MAJOR VERSION_LESS 3))
      message(FATAL_ERROR "Xyce_PYMI=ON and pybind11 found, but Python version used <3")
    else()
      get_target_property(python_target_type Python::Python TYPE)
      if (python_target_type STREQUAL STATIC_LIBRARY)
        message(FATAL_ERROR "Static Python library found but static Python libraries are not supported")
      else()
        configure_file( Xyce-PyMi.C.in ${CMAKE_CURRENT_BINARY_DIR}/Xyce-PyMi.C @ONLY)
        add_executable( Xyce-PyMi ${CMAKE_CURRENT_BINARY_DIR}/Xyce-PyMi.C )
        target_link_libraries(Xyce-PyMi PRIVATE XyceLib pybind11::embed)
        set_target_properties(Xyce-PyMi PROPERTIES
                              CXX_VISIBILITY_PRESET hidden
                              VISIBLITY_INLINES_HIDDEN ON)
        # See <https://github.com/pybind/pybind11/blob/master/docs/faq.rst#how-can-i-create-smaller-binaries>
        # for an explanation of why the "VISIBILITY" variables, above, are set.
        if (TARGET pybind11::thin_lto)
            target_link_libraries( Xyce-PyMi PRIVATE pybind11::thin_lto )
        elseif (TARGET pybind11::lto)
            target_link_libraries( Xyce-PyMi PRIVATE pybind11::lto )
        endif()
        # The above are, "An alternative to INTERPROCEDURAL_OPTIMIZATION for
        # adding link-time optimization." See
        # <https://github.com/pybind/pybind11/blob/master/docs/compiling.rst#pybind11_add_module>
        # for an explanation of why these are used.
  
        # If the location of the Python library contains other libraries such as libblas or liblapack,
        # then the RPATH provided will be given precedence over ldconfig at build/run time unless we also 
        # add the library locations for TPLs used by Trilinos to the RPATH 
        option(Xyce_PREPEND_TRILINOS_TPL_PATHS_TO_RPATH "Prepend Trilions TPLs' paths to Xyce-PyMi RPATH?" OFF)
        if (Xyce_PREPEND_TRILINOS_TPL_PATHS_TO_RPATH)
          get_target_property(PY_BUILD_RPATH Xyce-PyMi BUILD_RPATH)
          if (NOT PY_BUILD_RPATH)
            set(PY_BUILD_RPATH "")
          endif()
          foreach(tpl_loc ${Trilinos_TPL_LIBRARIES})
            get_filename_component(tpl_dirname ${tpl_loc} DIRECTORY)
            if (NOT(${tpl_dirname} IN_LIST PY_BUILD_RPATH))
              list(PREPEND PY_BUILD_RPATH ${tpl_dirname})
            endif()
          endforeach()
          set_target_properties(Xyce-PyMi PROPERTIES BUILD_RPATH "${PY_BUILD_RPATH}")
        endif()
  
        # ensures binary is in same folder as Xyce at build time and install time
        set_target_properties(Xyce-PyMi PROPERTIES 
            INSTALL_RPATH_USE_LINK_PATH ON
            RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/src"
        )
        install(TARGETS Xyce-PyMi DESTINATION bin)
        file(GLOB_RECURSE PYMI_PY_FILES ${PROJECT_SOURCE_DIR}/utils/PythonModelInterface/*.py)
        install(FILES ${PYMI_PY_FILES} DESTINATION include/pymi)
  
        message(STATUS "pybind11 found, Python>=3.4 found")
        message(STATUS "Python_EXECUTABLE: ${Python_EXECUTABLE}")
        message(STATUS "Python_LIBRARIES: ${Python_LIBRARIES}")
      endif()
    endif()
  endif()
endif(Xyce_PYMI)

add_executable( testGenCoup EXCLUDE_FROM_ALL testGenCoup.C )
target_link_libraries( testGenCoup XyceLib )

add_executable( TestGH19 EXCLUDE_FROM_ALL TestGH19.C )
target_link_libraries( TestGH19 XyceLib )

#if( Xyce_ENABLED_SHARED )
#  target_link_libraries( testGenCoup lib_xyce_shared )
#else( Xyce_ENABLED_SHARED )
#  target_link_libraries( testGenCoup lib_xyce_static )
#endif( Xyce_ENABLED_SHARED )

