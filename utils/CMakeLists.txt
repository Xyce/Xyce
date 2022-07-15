
# add in subdirectories to be searched

add_subdirectory(XyceCInterface)
if (Xyce_SIMULINK)
  add_subdirectory(utils)
endif()
