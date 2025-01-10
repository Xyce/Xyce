#/bin/sh

# Be sure to trun this from the project root!

(cd utils/ADMS/examples/fbh_hbt-2.1 &&  make all-source) &
(cd utils/ADMS/examples/vbic_r1.3_prerelease &&  make all-source) &
(cd utils/ADMS/examples/bsimcmg_107.0.0/code &&  make all-source) &
(cd utils/ADMS/examples/BSIMCMG108.0.0_20140822/code &&  make all-source) &
(cd utils/ADMS/examples/BSIMCMG110.0.0_20160101/code &&  make all-source)&
(cd utils/ADMS/examples/BSIM-CMG_111.2.1_06062022/code &&  make all-source)&
(cd utils/ADMS/examples/BSIM6.1.1/code &&  make all-source)&
(cd utils/ADMS/examples/BSIM-SOI_4/bsimsoi4.6.1 &&  make all-source)&
(cd utils/ADMS/examples/BSIM-SOI_4/bsimsoi4.5.0 &&  make all-source)&
(cd utils/ADMS/examples/DIODE_CMC_2/diode_cmc_2.0.0 &&  make all-source)&
(cd utils/ADMS/examples/psp102 &&  make all-source)&
(cd utils/ADMS/examples/psp103 &&  make all-source)&
(cd utils/ADMS/examples/mextram_504.12.1 &&  make all-source)&
(cd utils/ADMS/examples/mvs_2.0.0 &&  make all-source)&
(cd utils/ADMS/examples/mvsg_cmc_v1.1.0/vacode &&  make all-source)&
(cd utils/ADMS/examples/hicum &&  make all-source)&
(cd utils/ADMS/examples/hicum_l0 &&  make all-source)&
(cd utils/ADMS/examples/L_UTSOI/L_UTSOI_102 && make all-source)&
(cd utils/ADMS/examples/EKV/ekv-2.6 && make all-source)&
(cd src/DeviceModelPKG/Xyce_NonFree/Verilog/ekv301_02 && ./make_ekv_usable.sh)&

wait

diff -u -w -Bb  src/DeviceModelPKG/ADMS/N_DEV_ADMSHBT_X.h utils/ADMS/examples/fbh_hbt-2.1/N_DEV_ADMSHBT_X.h
diff -u -w -Bb  src/DeviceModelPKG/ADMS/N_DEV_ADMSHBT_X.C utils/ADMS/examples/fbh_hbt-2.1/N_DEV_ADMSHBT_X.C
diff -u -w -Bb  src/DeviceModelPKG/ADMS/N_DEV_ADMSvbic13.C utils/ADMS/examples/vbic_r1.3_prerelease/N_DEV_ADMSvbic13.C
diff -u -w -Bb  src/DeviceModelPKG/ADMS/N_DEV_ADMSvbic13.h utils/ADMS/examples/vbic_r1.3_prerelease/N_DEV_ADMSvbic13.h
diff -u -w -Bb  src/DeviceModelPKG/ADMS/N_DEV_ADMSvbic13_4t.C utils/ADMS/examples/vbic_r1.3_prerelease/N_DEV_ADMSvbic13_4t.C
diff -u -w -Bb  src/DeviceModelPKG/ADMS/N_DEV_ADMSvbic13_4t.h utils/ADMS/examples/vbic_r1.3_prerelease/N_DEV_ADMSvbic13_4t.h
diff -u -w -Bb  src/DeviceModelPKG/ADMS/N_DEV_ADMSbsimcmg.h utils/ADMS/examples/bsimcmg_107.0.0/code/N_DEV_ADMSbsimcmg.h
diff -u -w -Bb  src/DeviceModelPKG/ADMS/N_DEV_ADMSbsimcmg.C utils/ADMS/examples/bsimcmg_107.0.0/code/N_DEV_ADMSbsimcmg.C
diff -u -w -Bb  src/DeviceModelPKG/ADMS/N_DEV_ADMSbsimcmg_108.h utils/ADMS/examples/BSIMCMG108.0.0_20140822/code/N_DEV_ADMSbsimcmg_108.h
diff -u -w -Bb  src/DeviceModelPKG/ADMS/N_DEV_ADMSbsimcmg_108.C utils/ADMS/examples/BSIMCMG108.0.0_20140822/code/N_DEV_ADMSbsimcmg_108.C
diff -u -w -Bb  src/DeviceModelPKG/ADMS/N_DEV_ADMSbsimcmg_110.h utils/ADMS/examples/BSIMCMG110.0.0_20160101/code/N_DEV_ADMSbsimcmg_110.h
diff -u -w -Bb  src/DeviceModelPKG/ADMS/N_DEV_ADMSbsimcmg_110.C utils/ADMS/examples/BSIMCMG110.0.0_20160101/code/N_DEV_ADMSbsimcmg_110.C
diff -u -w -Bb  src/DeviceModelPKG/ADMS/N_DEV_ADMSbsimcmg_111.h utils/ADMS/examples/BSIM-CMG_111.2.1_06062022/code/N_DEV_ADMSbsimcmg_111.h
diff -u -w -Bb  src/DeviceModelPKG/ADMS/N_DEV_ADMSbsimcmg_111.C utils/ADMS/examples/BSIM-CMG_111.2.1_06062022/code/N_DEV_ADMSbsimcmg_111.C
diff -u -w -Bb  src/DeviceModelPKG/ADMS/N_DEV_ADMSbsim6.h utils/ADMS/examples/BSIM6.1.1/code/N_DEV_ADMSbsim6.h
diff -u -w -Bb  src/DeviceModelPKG/ADMS/N_DEV_ADMSbsim6.C utils/ADMS/examples/BSIM6.1.1/code/N_DEV_ADMSbsim6.C
diff -u -w -Bb  src/DeviceModelPKG/ADMS/N_DEV_ADMSbsimsoi.h utils/ADMS/examples/BSIM-SOI_4/bsimsoi4.6.1/N_DEV_ADMSbsimsoi.h
diff -u -w -Bb  src/DeviceModelPKG/ADMS/N_DEV_ADMSbsimsoi.C utils/ADMS/examples/BSIM-SOI_4/bsimsoi4.6.1/N_DEV_ADMSbsimsoi.C
diff -u -w -Bb  src/DeviceModelPKG/ADMS/N_DEV_ADMSbsimsoi450.h utils/ADMS/examples/BSIM-SOI_4/bsimsoi4.5.0/N_DEV_ADMSbsimsoi450.h
diff -u -w -Bb  src/DeviceModelPKG/ADMS/N_DEV_ADMSbsimsoi450.C utils/ADMS/examples/BSIM-SOI_4/bsimsoi4.5.0/N_DEV_ADMSbsimsoi450.C
diff -u -w -Bb  src/DeviceModelPKG/ADMS/N_DEV_ADMSDIODE_CMC.h utils/ADMS/examples/DIODE_CMC_2/diode_cmc_2.0.0/N_DEV_ADMSDIODE_CMC.h
diff -u -w -Bb  src/DeviceModelPKG/ADMS/N_DEV_ADMSDIODE_CMC.C utils/ADMS/examples/DIODE_CMC_2/diode_cmc_2.0.0/N_DEV_ADMSDIODE_CMC.C
diff -u -w -Bb  src/DeviceModelPKG/ADMS/N_DEV_ADMSPSP102VA.h utils/ADMS/examples/psp102/N_DEV_ADMSPSP102VA.h
diff -u -w -Bb  src/DeviceModelPKG/ADMS/N_DEV_ADMSPSP102VA.C utils/ADMS/examples/psp102/N_DEV_ADMSPSP102VA.C
diff -u -w -Bb  src/DeviceModelPKG/ADMS/N_DEV_ADMSPSP103VA.h utils/ADMS/examples/psp103/N_DEV_ADMSPSP103VA.h
diff -u -w -Bb  src/DeviceModelPKG/ADMS/N_DEV_ADMSPSP103VA.C utils/ADMS/examples/psp103/N_DEV_ADMSPSP103VA.C
diff -u -w -Bb  src/DeviceModelPKG/ADMS/N_DEV_ADMSPSP103TVA.h utils/ADMS/examples/psp103/N_DEV_ADMSPSP103TVA.h
diff -u -w -Bb  src/DeviceModelPKG/ADMS/N_DEV_ADMSPSP103TVA.C utils/ADMS/examples/psp103/N_DEV_ADMSPSP103TVA.C
diff -u -w -Bb  src/DeviceModelPKG/ADMS/N_DEV_ADMSJUNCAP200.h utils/ADMS/examples/psp103/N_DEV_ADMSJUNCAP200.h
diff -u -w -Bb  src/DeviceModelPKG/ADMS/N_DEV_ADMSJUNCAP200.C utils/ADMS/examples/psp103/N_DEV_ADMSJUNCAP200.C
diff -u -w -Bb  src/DeviceModelPKG/ADMS/N_DEV_ADMSbjt504va.h utils/ADMS/examples/mextram_504.12.1/N_DEV_ADMSbjt504va.h
diff -u -w -Bb  src/DeviceModelPKG/ADMS/N_DEV_ADMSbjt504va.C utils/ADMS/examples/mextram_504.12.1/N_DEV_ADMSbjt504va.C
diff -u -w -Bb  src/DeviceModelPKG/ADMS/N_DEV_ADMSbjt504tva.h utils/ADMS/examples/mextram_504.12.1/N_DEV_ADMSbjt504tva.h
diff -u -w -Bb  src/DeviceModelPKG/ADMS/N_DEV_ADMSbjt504tva.C utils/ADMS/examples/mextram_504.12.1/N_DEV_ADMSbjt504tva.C
diff -u -w -Bb  src/DeviceModelPKG/ADMS/N_DEV_ADMSmvs_2_0_0_etsoi.h utils/ADMS/examples/mvs_2.0.0/N_DEV_ADMSmvs_2_0_0_etsoi.h
diff -u -w -Bb  src/DeviceModelPKG/ADMS/N_DEV_ADMSmvs_2_0_0_etsoi.C utils/ADMS/examples/mvs_2.0.0/N_DEV_ADMSmvs_2_0_0_etsoi.C
diff -u -w -Bb  src/DeviceModelPKG/ADMS/N_DEV_ADMSmvs_2_0_0_hemt.h utils/ADMS/examples/mvs_2.0.0/N_DEV_ADMSmvs_2_0_0_hemt.h
diff -u -w -Bb  src/DeviceModelPKG/ADMS/N_DEV_ADMSmvs_2_0_0_hemt.C utils/ADMS/examples/mvs_2.0.0/N_DEV_ADMSmvs_2_0_0_hemt.C
diff -u -w -Bb  src/DeviceModelPKG/ADMS/N_DEV_ADMSmvsg_cmc.h utils/ADMS/examples/mvsg_cmc_v1.1.0/vacode/N_DEV_ADMSmvsg_cmc.h
diff -u -w -Bb  src/DeviceModelPKG/ADMS/N_DEV_ADMSmvsg_cmc.C utils/ADMS/examples/mvsg_cmc_v1.1.0/vacode/N_DEV_ADMSmvsg_cmc.C
diff -u -w -Bb  src/DeviceModelPKG/ADMS/N_DEV_ADMShicumL2va.h utils/ADMS/examples/hicum/N_DEV_ADMShicumL2va.h
diff -u -w -Bb  src/DeviceModelPKG/ADMS/N_DEV_ADMShicumL2va.C utils/ADMS/examples/hicum/N_DEV_ADMShicumL2va.C
diff -u -w -Bb  src/DeviceModelPKG/ADMS/N_DEV_ADMShic0_full.h utils/ADMS/examples/hicum_l0/N_DEV_ADMShic0_full.h
diff -u -w -Bb  src/DeviceModelPKG/ADMS/N_DEV_ADMShic0_full.C utils/ADMS/examples/hicum_l0/N_DEV_ADMShic0_full.C
diff -u -w -Bb  src/DeviceModelPKG/ADMS/N_DEV_ADMSl_utsoi.h utils/ADMS/examples/L_UTSOI/L_UTSOI_102/N_DEV_ADMSl_utsoi.h
diff -u -w -Bb  src/DeviceModelPKG/ADMS/N_DEV_ADMSl_utsoi.C utils/ADMS/examples/L_UTSOI/L_UTSOI_102/N_DEV_ADMSl_utsoi.C
diff -u -w -Bb  src/DeviceModelPKG/ADMS/N_DEV_ADMSekv_va.h utils/ADMS/examples/EKV/ekv-2.6/N_DEV_ADMSekv_va.h
diff -u -w -Bb  src/DeviceModelPKG/ADMS/N_DEV_ADMSekv_va.C utils/ADMS/examples/EKV/ekv-2.6/N_DEV_ADMSekv_va.C
diff -u -w -Bb  src/DeviceModelPKG/Xyce_NonFree/N_DEV_ADMSekv3.h src/DeviceModelPKG/Xyce_NonFree/Verilog/ekv301_02/N_DEV_ADMSekv3.h
diff -u -w -Bb  src/DeviceModelPKG/Xyce_NonFree/N_DEV_ADMSekv3.C src/DeviceModelPKG/Xyce_NonFree/Verilog/ekv301_02/N_DEV_ADMSekv3.C
