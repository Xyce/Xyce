# -- process remaining CMakeLists.txt files ----------------------------------

add_subdirectory (Core)
add_subdirectory (OpenModels)
add_subdirectory (TCADModels)
add_subdirectory (IBISModels)
add_subdirectory (EXTSC)

# optional source subdirectories
if (Xyce_ADMS_MODELS AND DEFINED Xyce_ADMS_MODELS_DIR)
     add_subdirectory (${Xyce_ADMS_MODELS_DIR} ${CMAKE_CURRENT_BINARY_DIR}/ADMS)
endif ()

if (Xyce_NEURON_MODELS AND DEFINED Xyce_NEURON_MODELS_DIR)
     add_subdirectory (${Xyce_NEURON_MODELS_DIR} ${CMAKE_CURRENT_BINARY_DIR}/NeuronModels)
endif ()

if (Xyce_NONFREE_MODELS AND DEFINED Xyce_NONFREE_MODELS_DIR)
     add_subdirectory (${Xyce_NONFREE_MODELS_DIR} ${CMAKE_CURRENT_BINARY_DIR}/Xyce_NonFree)
endif ()

if (Xyce_RAD_MODELS AND DEFINED Xyce_RAD_MODELS_DIR)
     add_subdirectory (${Xyce_RAD_MODELS_DIR} ${CMAKE_CURRENT_BINARY_DIR}/SandiaModels)
endif ()

# Currently, the EXTSC devices are required. To make them optional,
#    src/DeviceModelPKG/Core/N_DEV_RegisterDevices.C
# must be modified with the Xyce_EXTDEV ifdef.
# if ( Xyce_EXTDEV )
#   add_subdirectory( EXTSC )
# endif ( Xyce_EXTDEV )

