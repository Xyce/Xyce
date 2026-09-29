## With CMake 3.18, find_package(Python3) supports embedding
#find_package (Python3 COMPONENTS Development.Embed QUIET)
option(Xyce_PYMI "Build Python Model Interpreter (Xyce-PyMi)?" OFF)
if (Xyce_PYMI)
  find_package(Python REQUIRED COMPONENTS Interpreter Development Development.Embed)
  if (Python_VERSION VERSION_LESS "3.10")
      MESSAGE(FATAL_ERROR "
      Xyce_PYMI:BOOL=ON but Python version <3.10 is not supported.
      ")
  endif()
  if (Python_VERSION VERSION_GREATER_EQUAL "3.12")
    find_package(Python REQUIRED COMPONENTS Development.SABIModule)
  endif()
  MESSAGE(STATUS "Python_EXECUTABLE: ${Python_EXECUTABLE}")

  # locate nanobind
  execute_process(
      COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
      OUTPUT_STRIP_TRAILING_WHITESPACE
      OUTPUT_VARIABLE nanobind_ROOT
  )
  find_package(nanobind CONFIG REQUIRED)
  nanobind_build_library(nanobind)
  
  get_target_property(python_target_type Python::Python TYPE)
  if (python_target_type STREQUAL STATIC_LIBRARY)
    message(FATAL_ERROR "Static Python library found but static Python libraries are not supported")
  else()
    configure_file( Xyce-PyMi.C.in ${CMAKE_CURRENT_BINARY_DIR}/Xyce-PyMi.C @ONLY)
    add_executable( Xyce-PyMi ${CMAKE_CURRENT_BINARY_DIR}/Xyce-PyMi.C )
    if (Python_VERSION VERSION_GREATER_EQUAL "3.12")
      target_compile_definitions(Xyce-PyMi PRIVATE NB_STABLE_ABI)
    endif()
    target_link_libraries(Xyce-PyMi PRIVATE XyceLib nanobind Python::Python)
    target_include_directories(Xyce-PyMi PRIVATE ${nanobind_ROOT}/../include)
    set_target_properties(Xyce-PyMi PROPERTIES
                          CXX_VISIBILITY_PRESET hidden
                          VISIBILITY_INLINES_HIDDEN ON)
    # See <https://github.com/pybind/pybind11/blob/master/docs/faq.rst#how-can-i-create-smaller-binaries>
    # for an explanation of why the "VISIBILITY" variables, above, are set.
    nanobind_lto(Xyce-PyMi)
    # See <https://nanobind.readthedocs.io/en/latest/api_cmake.html> for details on LTO
  
    # If the location of the Python library contains other libraries such as libblas or liblapack,
    # then the RPATH provided will be given precedence over ldconfig at build/run time unless we also 
    # add the library locations for TPLs used by Trilinos to the RPATH 
    option(Xyce_PREPEND_TRILINOS_TPL_PATHS_TO_RPATH "Prepend Trilinos TPLs' paths to Xyce-PyMi RPATH?" OFF)
    if (Xyce_PREPEND_TRILINOS_TPL_PATHS_TO_RPATH)
      get_target_property(PY_BUILD_RPATH Xyce-PyMi BUILD_RPATH)
      if (NOT PY_BUILD_RPATH)
        set(PY_BUILD_RPATH "")
      endif()
      foreach(tpl_loc ${Trilinos_TPL_LIBRARIES})
        get_filename_component(tpl_dirname ${tpl_loc} DIRECTORY)
        if (tpl_dirname)
          if (NOT(${tpl_dirname} IN_LIST PY_BUILD_RPATH))
            list(PREPEND PY_BUILD_RPATH ${tpl_dirname})
          endif()
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
    install(TARGETS nanobind DESTINATION lib)
    file(GLOB_RECURSE PYMI_PY_FILES ${PROJECT_SOURCE_DIR}/utils/PythonModelInterface/*.py)
    install(FILES ${PYMI_PY_FILES} DESTINATION include/pymi)
  
    message(STATUS "nanobind found, Python found")
    message(STATUS "Python_EXECUTABLE: ${Python_EXECUTABLE}")
    message(STATUS "Python_LIBRARIES: ${Python_LIBRARIES}")
  endif()
endif(Xyce_PYMI)

add_executable( testGenCoup testGenCoup.C )
target_link_libraries( testGenCoup XyceLib )

add_executable( TestGH19 TestGH19.C )
target_link_libraries( TestGH19 XyceLib )

add_executable( TestIssue29 TestIssue29.C )
target_link_libraries( TestIssue29 XyceLib )

