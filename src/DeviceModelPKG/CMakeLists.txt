# -- process remaining CMakeLists.txt files ----------------------------------

add_subdirectory (Core)
add_subdirectory (OpenModels)
add_subdirectory (TCADModels)
add_subdirectory (IBISModels)
add_subdirectory (EXTSC)

# optional source subdirectories
if (Xyce_ADMS_MODELS)
     add_subdirectory (ADMS)
endif ()

if (Xyce_NEURON_MODELS)
     add_subdirectory (NeuronModels)
endif ()

if (Xyce_NONFREE_MODELS)
     add_subdirectory (Xyce_NonFree)
endif ()

if (Xyce_RAD_MODELS)
     add_subdirectory (SandiaModels)
endif ()

# Currently, the EXTSC devices are required. To make them optional,
#    src/DeviceModelPKG/Core/N_DEV_RegisterDevices.C
# must be modified with the Xyce_EXTDEV ifdef.
# if ( Xyce_EXTDEV )
#   add_subdirectory( EXTSC )
# endif ( Xyce_EXTDEV )

#Catch the scope and put it into the top level
set (XYCE_SOURCES ${XYCE_SOURCES} PARENT_SCOPE)
set (XYCE_INCLUDE_DIRECTORIES ${XYCE_INCLUDE_DIRECTORIES} PARENT_SCOPE) 

