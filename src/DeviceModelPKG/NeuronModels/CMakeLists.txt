#
# CMakeLists.txt for NeuronModels
#
# Author: Aaron Gibson (asgibso@sandia.gov)
#
# Builds the NeuronModels for Xyce.


add_to_xyce_library( NAME NeuronModels
        FILES
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
)

#install( FILES ${OpenModel_Headers}
#        DESTINATION include
#)

