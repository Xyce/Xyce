#
# -- process files for this library ------------------------------------------
#
set( DakotaLink_SOURCES
	N_DAK_DakotaController.C
	N_DAK_DakotaInterface.C
)
add_to_xyce_library( NAME DakotaLink
	FILES ${DakotaLink_SOURCES}
)
