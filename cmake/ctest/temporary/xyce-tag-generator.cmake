# reads a "|" separated list of xyce build options and the associated
# tag(s) as referenced in tests and in the run_xyce_regression utility
file(STRINGS ${CTEST_SCRIPT_DIRECTORY}/build-opt-to-tag.txt rawOptTag)

if(NOT rawOptTag)
  message(FATAL_ERROR "ERROR: Unable to read file \"build-opt-to-tag.txt\"")
endif()

# parse the content of the file into a cmake LIST
foreach(tagLine IN LISTS rawOptTag)
  if(VERBOSITY GREATER 4)
    message("DEBUG[xtg]: original line: \"${tagLine}\"")
  endif()

  # remove any whitespace
  string(REGEX REPLACE "[ \t\r\n]" "" tagLine "${tagLine}")

  if(VERBOSITY GREATER 3)
    message("DEBUG[xtg]: stripped line: \"${tagLine}\"")
  endif()

  string(REGEX MATCH "^#" commentBool "${tagLine}")

  # finally, if the line isn't a comment process it
  if(NOT commentBool)
    list(APPEND buildOptToTag ${tagLine})
  endif()

endforeach()

if(VERBOSITY GREATER 2)
  message("DEBUG[xtg]: buildOptToTag = ${buildOptToTag}")
endif()

# now convert the LIST into a tag line for run_xyce_regression
foreach(tagEnt IN LISTS buildOptToTag)
  if(${Xyce_ADMS_MODELS})
    message("DEBUG[xtg]: YES")
  else()
    message("DEBUG[xtg]: NO")
  endif()
endforeach()
