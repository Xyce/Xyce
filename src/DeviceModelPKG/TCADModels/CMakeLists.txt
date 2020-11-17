
# -- process files for this library ------------------------------------------

# class source list
target_sources(XyceStatic PRIVATE
      N_DEV_RegisterTCADDevices.C
      N_DEV_2DPDE_DAE.C
      N_DEV_2DPDE_Output.C
      N_DEV_2DPDEInstance.C
      N_DEV_2DPDEModel.C
      N_DEV_2DPDEParam.C
      N_DEV_2DPDESetup.C
      N_DEV_BernouliSupport.C
      N_DEV_DevicePDEInstance.C
      N_DEV_DiodePDE_DAE.C
      N_DEV_DiodePDEInstance.C
      N_DEV_DiodePDEModel.C
      N_DEV_DopeInfo.C
      N_DEV_MaterialLayer.C
      N_DEV_MaterialSupport.C
      N_DEV_PDEMeshContainer.C
      N_DEV_PDE_2DMesh.C
      N_DEV_PDE_Electrode.C)

# add header files to project
target_include_directories(XyceStatic PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR})

# -- build requested libraries -----------------------------------------------
set ( XYCE_PACKAGE_LIBRARIES ${XYCE_PACKAGE_LIBRARIES}
      TCADModels
      CACHE INTERNAL "X_L" )
