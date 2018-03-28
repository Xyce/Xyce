#
# CMakeLists.txt for CircuitPKG
#

# Set the header files (to be installed) here.
add_to_xyce_library( NAME Circuit
	FILES
		N_CIR_Xyce.C
		N_CIR_SecondLevelSimulator.C
		N_CIR_MixedSignalSimulator.C
		N_CIR_Xygra.C
		N_CIR_GenCouplingSimulator.C
)

# Include the headers in the install.
set( CircuitPKG_Headers
	N_CIR_GenCouplingSimulator.h
	N_CIR_Xyce.h
	N_CIR_Xygra.h
)
install( FILES ${CircuitPKG_Headers}
	DESTINATION include
)
