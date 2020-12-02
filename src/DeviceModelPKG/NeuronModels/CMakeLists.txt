
# -- process files for this library ------------------------------------------

# class source list
target_sources(XyceLib PRIVATE
      N_DEV_RegisterNeuronDevices.C
      N_DEV_Neuron.C
      N_DEV_Neuron2.C
      N_DEV_Neuron3.C
      N_DEV_Neuron4.C
      N_DEV_Neuron5.C
      N_DEV_Neuron6.C
      N_DEV_Neuron7.C
      N_DEV_Neuron8.C
      N_DEV_Neuron9.C
      N_DEV_NeuronPop1.C
      N_DEV_Synapse.C
      N_DEV_Synapse2.C
      N_DEV_Synapse3.C
      N_DEV_Synapse4.C
      N_DEV_MembraneModel.C
      N_DEV_MembranePassive.C
      N_DEV_MembraneHH.C
      N_DEV_MembraneCS.C
      N_DEV_MembraneUserDefined.C)
# add header files to project
target_include_directories(XyceLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})




# -- build requested libraries -----------------------------------------------
set ( XYCE_PACKAGE_LIBRARIES ${XYCE_PACKAGE_LIBRARIES}
      NeuronModels
      CACHE INTERNAL "X_L" )
