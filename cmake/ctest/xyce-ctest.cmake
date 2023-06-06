# hard-code most things for now, later convert to arguments passed in
# via the ctest command line

set(CTEST_PROJECT_NAME "Xyce")
set(CTEST_DROP_METHOD "http")
set(CTEST_DROP_SITE "joseki-srn.sandia.gov/CDash")
set(CTEST_DROP_LOCATION "/submit.php?project=Xyce")

set(CTEST_SOURCE_DIRECTORY "$ENV{WORKSPACE}/source/Xyce")
set(CTEST_BINARY_DIRECTORY "$ENV{WORKSPACE}/build")

find_program(HNAME NAMES hostname)
execute_process(COMMAND "${HNAME}"
  OUTPUT_VARIABLE myhost
  OUTPUT_STRIP_TRAILING_WHITESPACE)

set(CTEST_SITE "${myhost}")

set(CMAKE_BUILD_TYPE "Release")

set(CTEST_BUILD_NAME "glh_Intel64_RHEL7_Serial-cmake-cde-gcc")

set(CTEST_BUILD_FLAGS "-j16")

set(DONOTCLONE TRUE)
set(DONOTUPDATE TRUE)

set(CTEST_CMAKE_GENERATOR "Unix Makefiles")

set(MODEL "Experimental")

ctest_start(${MODEL} GROUP ${MODEL})
ctest_configure()
ctest_build()
ctest_submit(RETRY_COUNT 10 RETRY_DELAY 30)
