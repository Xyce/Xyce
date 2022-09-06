
# -- process files for this library ------------------------------------------

# class source list
target_sources(XyceLib PRIVATE
      N_DEV_RegisterOpenDevices.C
      N_DEV_ArtificialParameters.C
      N_DEV_ACC.C
      N_DEV_BJT.C
      N_DEV_Battery.C
      N_DEV_Capacitor.C
      N_DEV_DAC.C
      N_DEV_ADC.C
      N_DEV_Bsrc.C
      N_DEV_Diode.C
      N_DEV_GeneralExternal.C
      N_DEV_ISRC.C
      N_DEV_Inductor.C
      N_DEV_JFET.C
      N_DEV_LTRA.C
      N_DEV_MemristorPEM.C
      N_DEV_MemristorTEAM.C
      N_DEV_MemristorYakopcic.C
      N_DEV_MESFET.C
      N_DEV_MOSFET1.C
      N_DEV_MOSFET2.C
      N_DEV_MOSFET3.C
      N_DEV_MOSFET6.C
      N_DEV_MOSFET_B3.C
      N_DEV_MOSFET_B4.C
      N_DEV_MOSFET_B4p61.C
      N_DEV_MOSFET_B4p70.C
      N_DEV_MOSFET_B3SOI.C
      N_DEV_MutIndLin.C
      N_DEV_MutIndNonLin.C
      N_DEV_MutIndNonLin2.C
      N_DEV_OpAmp.C
      N_DEV_Digital.C
      N_DEV_VDMOS.C
      N_DEV_Resistor.C
      N_DEV_Resistor3.C
      N_DEV_ThermalResistor.C
      N_DEV_SW.C
      N_DEV_TRA.C
      N_DEV_VCCS.C
      N_DEV_Vcvs.C
      N_DEV_Vsrc.C
      N_DEV_ROM.C
      N_DEV_TransLine.C
      N_DEV_PowerGrid.C
      N_DEV_PowerGridBranch.C
      N_DEV_PowerGridBusShunt.C
      N_DEV_PowerGridGenBus.C
      N_DEV_PowerGridTransformer.C
      N_DEV_AntiWindupLimiter.C
      N_DEV_YLin.C
      N_DEV_Delay.C)

# add header files to project
target_include_directories(XyceLib PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${REACTION_PARSER_INC}> )

# -- build requested libraries -----------------------------------------------
set ( XYCE_PACKAGE_LIBRARIES ${XYCE_PACKAGE_LIBRARIES}
      OpenModels
      CACHE INTERNAL "X_L" )
