
#This is being set here, because the check for iota in "check_include_file_cxx.cmake" is not working
set (HAVE_IOTA true)


if (EXISTS "${PROJECT_SOURCE_DIR}/src/DeviceModelPKG/ADMS" AND ( (NOT DEFINED Xyce_ADMS_MODELS) OR Xyce_ADMS_MODELS) )
     set (Xyce_ADMS_MODELS TRUE CACHE BOOL "Include the ADMS directory, if it exists")
     message(STATUS "Including the src/DeviceModelPKG/ADMS model directory")
elseif (NOT DEFINED Xyce_ADMS_MODELS)
     set (Xyce_ADMS_MODELS FALSE CACHE BOOL "Include the ADMS directory, if it exists")
elseif (NOT Xyce_ADMS_MODELS)
     message(STATUS "NOT including the src/DeviceModelPKG/ADMS model directory")
     set (Xyce_ADMS_MODELS FALSE CACHE BOOL "Include the ADMS directory, if it exists")
else ()
     set (Xyce_ADMS_MODELS FALSE CACHE BOOL "Include the ADMS directory, if it exists")
     message(WARNING "The src/DeviceModelPKG/ADMS directory is not in place - "
          "changing Xyce_ADMS_MODELS to FALSE")
endif ()

if (EXISTS "${PROJECT_SOURCE_DIR}/src/DeviceModelPKG/NeuronModels" AND ( (NOT DEFINED Xyce_NEURON_MODELS) OR Xyce_NEURON_MODELS) )
     set (Xyce_NEURON_MODELS TRUE CACHE BOOL "Include the Neuron directory, if it exists")
     message(STATUS "Including the src/DeviceModelPKG/NeuronModels model directory")
elseif (NOT DEFINED Xyce_NEURON_MODELS)
     set (Xyce_NEURON_MODELS FALSE CACHE BOOL "Include the Neuron directory, if it exists")
elseif (NOT Xyce_NEURON_MODELS)
     set (Xyce_NEURON_MODELS FALSE CACHE BOOL "Include the Neuron directory, if it exists")
     message(STATUS "NOT including the src/DeviceModelPKG/NeuronModels model directory")
else ()
     set (Xyce_NEURON_MODELS FALSE CACHE BOOL "Include the Neuron directory, if it exists")
     message(WARNING "The src/DeviceModelPKG/NeuronModels directory is not in place - "
          "changing Xyce_NEURON_MODELS to FALSE.")
endif ()

if (EXISTS "${PROJECT_SOURCE_DIR}/src/DeviceModelPKG/Xyce_NonFree" AND ( (NOT DEFINED Xyce_NONFREE_MODELS) OR Xyce_NONFREE_MODELS) )
     set (Xyce_NONFREE_MODELS TRUE CACHE BOOL "Include the NonFree directory, if it exists")
     message(STATUS "Including the src/DeviceModelPKG/Xyce_NonFree model directory")
elseif (NOT DEFINED Xyce_NONFREE_MODELS)
     set (Xyce_NONFREE_MODELS FALSE CACHE BOOL "Include the NonFree directory, if it exists")
elseif (NOT Xyce_NONFREE_MODELS)
     set (Xyce_NONFREE_MODELS FALSE CACHE BOOL "Include the NonFree directory, if it exists")
     message(STATUS "NOT including the src/DeviceModelPKG/Xyce_NonFree model directory")
else ()
     set (Xyce_NONFREE_MODELS FALSE CACHE BOOL "Include the NonFree directory, if it exists")
     message(WARNING "The src/DeviceModelPKG/NonFree directory is not in place - "
          "changing Xyce_NONFREE_MODELS to FALSE.")
endif ()

if (EXISTS "${PROJECT_SOURCE_DIR}/src/DeviceModelPKG/SandiaModels" AND ( (NOT DEFINED Xyce_RAD_MODELS) OR Xyce_RAD_MODELS) )
     set (Xyce_RAD_MODELS TRUE CACHE BOOL "Include the SandiaModels directory, if it exists")
     message(STATUS "Including the src/DeviceModelPKG/SandiaModels model directory")
elseif (NOT DEFINED Xyce_RAD_MODELS)
     set (Xyce_RAD_MODELS FALSE CACHE BOOL "Include the SandiaModels directory, if it exists")
elseif (NOT Xyce_RAD_MODELS)
     message(STATUS "NOT including the src/DeviceModelPKG/SandiaModels model directory")
     set (Xyce_RAD_MODELS FALSE CACHE BOOL "Include the SandiaModels directory, if it exists")
else ()
     set (Xyce_RAD_MODELS FALSE CACHE BOOL "Include the SandiaModels directory, if it exists")
     message(WARNING "The src/DeviceModelPKG/SandiaModels directory is not in place - "
          "changing Xyce_RAD_MODELS to FALSE.")
endif ()

#### Only on the systems that have bison and flex for now #######
# Also this is both a C++ and CMAKE variable set. There are both ifs 
# in the CmakeLists.txt in src/DeviceModelPKG/Core/CMakeLists.txt but 
# also #ifndefs inside the source code on the variable. 
#add_definitions( -DXyce_REACTION_PARSER) # This is a #cmakedefine in Xyce_config.h.cmake
set(Xyce_REACTION_PARSER true CACHE BOOL "The reaction parser depends on Flex and Bison. Turn this off if you don't need
this feature. This is a super sepcialized feature that is for chemical reactions." ) 
add_definitions( -DYYTEXT_POINTER) 
#
# If everything is in place, define `Xyce_REACTION_PARSER`
#* If everything is _not_ in place, give a warning (and make sure
#  `Xyce_REACTION_PARSER` is disabled)
#* Allow the developer to turn it off
#
# If Xyce_REACTION_PARSER is both a CMAKE variable and a C++ variable 
# the code in `src/DeviceModelPKG/Core` will
# through Bison and Flex. The `CMakeLists.txt` has an if on the CMAKE 
# variable Xyce_REACTION_PARSER which executes the Bison and Flex Find_package 
# macros names BISON_TARGET and FLEX_TARGET respectively. They're detailed in the
# CMAKE documentation for each find_package. 

