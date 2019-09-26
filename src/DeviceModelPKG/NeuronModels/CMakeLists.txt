
# -- process files for this library ------------------------------------------

# class source list
set ( NeuronModels_SOURCES
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_RegisterNeuronDevices.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_Neuron.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_Neuron2.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_Neuron3.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_Neuron4.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_Neuron5.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_Neuron6.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_Neuron7.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_Neuron8.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_Neuron9.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_NeuronPop1.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_Synapse.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_Synapse2.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_Synapse3.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_Synapse4.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_MembraneModel.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_MembranePassive.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_MembraneHH.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_MembraneCS.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_MembraneUserDefined.C
      CACHE INTERNAL "X_S" )

set (XYCE_SOURCES ${XYCE_SOURCES} ${NeuronModels_SOURCES} PARENT_SCOPE )
# add header files to project
set ( XYCE_INCLUDE_DIRECTORIES ${XYCE_INCLUDE_DIRECTORIES} 
      ${CMAKE_CURRENT_SOURCE_DIR}
      PARENT_SCOPE )




# -- build requested libraries -----------------------------------------------
set ( XYCE_PACKAGE_LIBRARIES ${XYCE_PACKAGE_LIBRARIES} 
      NeuronModels 
      CACHE INTERNAL "X_L" )
