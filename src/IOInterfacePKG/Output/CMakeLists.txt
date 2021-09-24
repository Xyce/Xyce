
# -- process files for this library ------------------------------------------

# class source list
target_sources(XyceLib PRIVATE
     N_IO_MeasureManager.C
     N_IO_MeasureBase.C
     N_IO_MeasureExtrema.C
     N_IO_MeasureStats.C
     N_IO_MeasureTranStats.C
     N_IO_MeasureRiseFallDelay.C
     N_IO_MeasureAverage.C
     N_IO_MeasureDuty.C
     N_IO_MeasureErrorFunctions.C
     N_IO_MeasureFrequency.C
     N_IO_MeasureMax.C
     N_IO_MeasureMin.C
     N_IO_MeasureOffTime.C
     N_IO_MeasureOnTime.C
     N_IO_MeasurePeakToPeak.C
     N_IO_MeasureRMS.C
     N_IO_MeasureWhenAT.C
     N_IO_MeasureFindWhen.C
     N_IO_MeasureEquationEvaluation.C
     N_IO_MeasureDerivativeEvaluation.C
     N_IO_MeasureIntegralEvaluation.C
     N_IO_MeasureError.C
     N_IO_MeasureFourier.C
     N_IO_MeasureFFT.C
     N_IO_Remeasure.C )

# add header files to project
target_include_directories(XyceLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})




# -- build requested libraries -----------------------------------------------
set ( XYCE_PACKAGE_LIBRARIES ${XYCE_PACKAGE_LIBRARIES}
      OutputInterface
      CACHE INTERNAL "X_L" )
