# -- process remaining CMakeLists.txt files ----------------------------------

add_subdirectory( Core )
add_subdirectory( OpenModels )
add_subdirectory( TCADModels )
add_subdirectory( IBISModels )

# optional source subdirectories
if ( Xyce_EXTDEV )
  add_subdirectory( EXTSC )
endif ( Xyce_EXTDEV )

if ( Xyce_RAD_MODELS )
  add_subdirectory( SandiaModels )
endif ( Xyce_RAD_MODELS )

if ( Xyce_NONFREE_MODELS )
  add_subdirectory( Xyce_NonFree )
endif ( Xyce_NONFREE_MODELS )

if ( Xyce_ADMS_MODELS )
  add_subdirectory( ADMS )
endif ( Xyce_ADMS_MODELS )

if ( Xyce_NEURON_MODELS )
  add_subdirectory( NeuronModels )
endif ( Xyce_NEURON_MODELS )

