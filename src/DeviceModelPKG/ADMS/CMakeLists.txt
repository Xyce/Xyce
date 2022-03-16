
# -- process files for this library ------------------------------------------

# class source list
target_sources(XyceLib PRIVATE
      N_DEV_RegisterADMSDevices.C
      N_DEV_ADMSvbic13.C
      N_DEV_ADMSvbic13_4t.C
      N_DEV_ADMSPSP103VA.C
      N_DEV_ADMSPSP103TVA.C
      N_DEV_ADMSPSP102VA.C
      N_DEV_ADMSJUNCAP200.C
      N_DEV_ADMSDIODE_CMC.C
      N_DEV_ADMSHBT_X.C
      N_DEV_ADMSbsimcmg.C
      N_DEV_ADMSbsimcmg_108.C
      N_DEV_ADMSbsimcmg_110.C
      N_DEV_ADMSbsim6.C
      N_DEV_ADMSbsimsoi.C
      N_DEV_ADMSbsimsoi450.C
      N_DEV_ADMSbjt504va.C
      N_DEV_ADMSbjt504tva.C
      N_DEV_ADMSmvs_2_0_0_etsoi.C
      N_DEV_ADMSmvs_2_0_0_hemt.C
      N_DEV_ADMSmvsg_cmc.C
      N_DEV_ADMShicumL2va.C
      N_DEV_ADMShic0_full.C
      N_DEV_ADMSl_utsoi.C)

# add header files to project
target_include_directories(XyceLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# -- build requested libraries -----------------------------------------------
set ( XYCE_PACKAGE_LIBRARIES ${XYCE_PACKAGE_LIBRARIES}
      ADMS
      CACHE INTERNAL "X_L" )
