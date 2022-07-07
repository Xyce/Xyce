# Generate timestamp.h
string(TIMESTAMP CONFIGURE_TIME %Y%m%d%H%M)
if(Git_FOUND AND EXISTS "${GIVEN_SOURCE_DIR}/.git")
     execute_process(
          COMMAND ${GIT_EXECUTABLE} -C ${GIVEN_SOURCE_DIR} describe --dirty
          OUTPUT_VARIABLE GITSHA
          ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
else()
     set(GITSHA "UNKNOWN")
endif()
if (NOT ${GITSHA} STREQUAL "")
     set(GITSTRING "-(${GITSHA})")
endif ()
file(WRITE timestamp.h "#define XYCEBUILDTIMESTAMP \"${CONFIGURE_TIME}${GITSTRING}\"")
