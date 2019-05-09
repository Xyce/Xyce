
# -- process files for this library ------------------------------------------

# class source list
set ( TCADModels_SOURCES  
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_RegisterTCADDevices.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_2DPDE_DAE.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_2DPDE_Output.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_2DPDEInstance.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_2DPDEModel.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_2DPDEParam.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_2DPDESetup.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_BernouliSupport.C  
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_DevicePDEInstance.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_DiodePDE_DAE.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_DiodePDEInstance.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_DiodePDEModel.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_DopeInfo.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_MaterialLayer.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_MaterialSupport.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_PDEMeshContainer.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_PDE_2DMesh.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_PDE_Electrode.C
      CACHE INTERNAL "X_S" )

# add header files to project
set ( XYCE_INCLUDE_DIRECTORIES ${XYCE_INCLUDE_DIRECTORIES} 
      ${CMAKE_CURRENT_SOURCE_DIR}
      CACHE INTERNAL "X_I" )

# -- build requested libraries -----------------------------------------------
set ( XYCE_PACKAGE_LIBRARIES ${XYCE_PACKAGE_LIBRARIES} 
      TCADModels
      CACHE INTERNAL "X_L" )
