#
# CMakeLists.txt for DeviceModelPKG
#
# Author: Aaron Gibson (asgibso@sandia.gov)
#
# This file determines which models should be compiled directly into Xyce.
#
# We first add in the directories that should always be built, then define
# some macros for managing additional models and their options.
#

# Note that we will print the diagnostic information about the devices that
# we are including in this build so that this information is available by
# running cmake again.
message( "Xyce Device Models:" )
message( "-- Core" )
message( "-- OpenModels" )
message( "-- TCADModels" )

# Include the models that we should always try and build.
add_xyce_subdirectory( "${CMAKE_CURRENT_SOURCE_DIR}/Core" )
add_xyce_subdirectory( "${CMAKE_CURRENT_SOURCE_DIR}/OpenModels" )
add_xyce_subdirectory( "${CMAKE_CURRENT_SOURCE_DIR}/TCADModels" )

# Now define macros for our optional models.
#
# add_optional_model
#
# This macro will add an option to this build to signal whether to include
# the given MODEL_NAME (located at MODEL_PATH) in the Xyce libraries and
# executable. MODEL_DESC stores the option description for use with ccmake
# or cmake-gui.
#
macro( add_optional_model MODEL_NAME MODEL_PATH MODEL_DESC )
	option( "Xyce_${MODEL_NAME}" "${MODEL_DESC}" ON )
	if( Xyce_${MODEL_NAME} )
		add_xyce_subdirectory( "${MODEL_PATH}" )
		
		# Add our diagnostic message indicating this model type is enabled.
		option_message( Xyce_${MODEL_NAME} "${MODEL_NAME}" )
	endif()
endmacro()

macro( add_model_from_other_repo MODEL_NAME MODEL_FILENAME MODEL_DESC )
	# Search for the path to the model in some default locations.
	set( SEARCH_PATHS 
		"${CMAKE_CURRENT_SOURCE_DIR}/${MODEL_FILENAME}"
		"${CMAKE_CURRENT_SOURCE_DIR}/Xyce_${MODEL_FILENAME}"
		"${CMAKE_SOURCE_DIR}/../${MODEL_FILENAME}"
		"${CMAKE_SOURCE_DIR}/../Xyce_${MODEL_FILENAME}"
	)
	set( SHOULD_ENABLE_MODEL OFF )
	set( MODEL_PATH "" )
	foreach( TMPPATH ${SEARCH_PATHS} )
		if( EXISTS ${TMPPATH} )
			set( SHOULD_ENABLE_MODEL ON )
			set( MODEL_PATH ${TMPPATH} )
			break()
		endif()
	endforeach()
	
	# Regardless of whether we found the path or not, we should add it
	# as an option so the developer knows it exists.
	set( Xyce_${MODEL_NAME} ${SHOULD_ENABLE_MODEL} CACHE BOOL "${MODEL_DESC}" )
	set( Xyce_${MODEL_NAME}_PATH ${MODEL_PATH} CACHE PATH "Path to ${MODEL_NAME}" )
	
	if( Xyce_${MODEL_NAME} )
		if( NOT EXISTS "${Xyce_${MODEL_NAME}_PATH}/CMakeLists.txt" )
			message( SEND_ERROR
				"Xyce_${MODEL_NAME} enabled, but Xyce_${MODEL_NAME}_PATH does not point to directory containing CMakeLists.txt!" )
		endif()
		add_xyce_subdirectory( ${Xyce_${MODEL_NAME}_PATH} )
	endif()

	# Add a diagnostic message.
	option_message( Xyce_${MODEL_NAME} "${MODEL_NAME}" )
	if( Xyce_${MODEL_NAME} )
		message( "-- * Path to Model: ${Xyce_${MODEL_NAME}_PATH}" )
	endif()
endmacro()

# Add some optional models that are a part of the Xyce repository.
add_optional_model( ADMS_MODELS "${CMAKE_CURRENT_SOURCE_DIR}/ADMS" "Enable ADMS Models" )
add_optional_model( EXTDEV "${CMAKE_CURRENT_SOURCE_DIR}/EXTSC" "Enable External Devices" )
add_optional_model( NEURON_MODELS "${CMAKE_CURRENT_SOURCE_DIR}/NeuronModels" "Enable Xyce Neuron Models" )

# Add some optional models that are not a part of the Xyce repository directly,
# but are stored elsewhere. 
add_model_from_other_repo( NONFREE_MODELS NonFree 
	"Path to Non-free models. If an invalid path, non-free models are disabled." )

add_model_from_other_repo( RAD_MODELS SandiaModels 
	"Path to SandiaModels. If an invalid path, SandiaModels will be disabled." )

