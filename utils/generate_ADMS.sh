#/bin/sh

# Be sure to trun this from the project root!

(cd utils/ADMS/examples/fbh_hbt-2.1 &&  make all-source) &
(cd utils/ADMS/examples/vbic_r1.3_prerelease &&  make all-source) &
(cd utils/ADMS/examples/bsimcmg_107.0.0/code &&  make all-source) &
# Do not regenerate --- we use a hand-optimized version now
# (cd utils/ADMS/examples/BSIMCMG110.0.0_20160101/code &&  make all-source)&
(cd utils/ADMS/examples/BSIM6.1.1/code &&  make all-source)&
(cd utils/ADMS/examples/psp102 &&  make all-source)&
(cd utils/ADMS/examples/psp103 &&  make all-source)&
(cd utils/ADMS/examples/mextram_504.12.1 &&  make all-source)&
(cd utils/ADMS/examples/mvs_2.0.0 &&  make all-source)&
(cd utils/ADMS/examples/mvsg_cmc_v1.1.0/vacode &&  make all-source)&
(cd utils/ADMS/examples/hicum &&  make all-source)&
(cd utils/ADMS/examples/hicum_l0 &&  make all-source)&
(cd src/DeviceModelPKG/Xyce_NonFree/Verilog/ekv301_02 && ./make_ekv_usable.sh)&
(cd src/DeviceModelPKG/Xyce_NonFree/Verilog/ekv2.6 && ./make_ekv2.6_usable.sh)&

wait

diff -u -Bb utils/ADMS/examples/fbh_hbt-2.1/N_DEV_ADMSHBT_X.h src/DeviceModelPKG/ADMS/N_DEV_ADMSHBT_X.h
diff -u -Bb utils/ADMS/examples/fbh_hbt-2.1/N_DEV_ADMSHBT_X.C src/DeviceModelPKG/ADMS/N_DEV_ADMSHBT_X.C
#diff -u -Bb utils/ADMS/examples/vbic_r1.3_prerelease/N_DEV_ADMSvbic13.C src/DeviceModelPKG/ADMS/N_DEV_ADMSvbic13.C
#diff -u -Bb utils/ADMS/examples/vbic_r1.3_prerelease/N_DEV_ADMSvbic13.h src/DeviceModelPKG/ADMS/N_DEV_ADMSvbic13.h
diff -u -Bb utils/ADMS/examples/vbic_r1.3_prerelease/N_DEV_ADMSvbic13_4t.C src/DeviceModelPKG/ADMS/N_DEV_ADMSvbic13_4t.C
diff -u -Bb utils/ADMS/examples/vbic_r1.3_prerelease/N_DEV_ADMSvbic13_4t.h src/DeviceModelPKG/ADMS/N_DEV_ADMSvbic13_4t.h
diff -u -Bb utils/ADMS/examples/bsimcmg_107.0.0/code/N_DEV_ADMSbsimcmg.h src/DeviceModelPKG/ADMS/N_DEV_ADMSbsimcmg.h
diff -u -Bb utils/ADMS/examples/bsimcmg_107.0.0/code/N_DEV_ADMSbsimcmg.C src/DeviceModelPKG/ADMS/N_DEV_ADMSbsimcmg.C
#diff -u -Bb utils/ADMS/examples/BSIMCMG110.0.0_20160101/code/N_DEV_ADMSbsimcmg_110.h src/DeviceModelPKG/ADMS/N_DEV_ADMSbsimcmg_110.h
#diff -u -Bb utils/ADMS/examples/BSIMCMG110.0.0_20160101/code/N_DEV_ADMSbsimcmg_110.C src/DeviceModelPKG/ADMS/N_DEV_ADMSbsimcmg_110.C
diff -u -Bb utils/ADMS/examples/BSIM6.1.1/code/N_DEV_ADMSbsim6.h src/DeviceModelPKG/ADMS/N_DEV_ADMSbsim6.h
diff -u -Bb utils/ADMS/examples/BSIM6.1.1/code/N_DEV_ADMSbsim6.C src/DeviceModelPKG/ADMS/N_DEV_ADMSbsim6.C
diff -u -Bb utils/ADMS/examples/psp102/N_DEV_ADMSPSP102VA.h src/DeviceModelPKG/ADMS/N_DEV_ADMSPSP102VA.h
diff -u -Bb utils/ADMS/examples/psp102/N_DEV_ADMSPSP102VA.C src/DeviceModelPKG/ADMS/N_DEV_ADMSPSP102VA.C
diff -u -Bb utils/ADMS/examples/psp103/N_DEV_ADMSPSP103VA.h src/DeviceModelPKG/ADMS/N_DEV_ADMSPSP103VA.h
diff -u -Bb utils/ADMS/examples/psp103/N_DEV_ADMSPSP103VA.C src/DeviceModelPKG/ADMS/N_DEV_ADMSPSP103VA.C
diff -u -Bb utils/ADMS/examples/psp103/N_DEV_ADMSJUNCAP200.h src/DeviceModelPKG/ADMS/N_DEV_ADMSJUNCAP200.h
diff -u -Bb utils/ADMS/examples/psp103/N_DEV_ADMSJUNCAP200.C src/DeviceModelPKG/ADMS/N_DEV_ADMSJUNCAP200.C
diff -u -Bb utils/ADMS/examples/mextram_504.12.1/N_DEV_ADMSbjt504va.h src/DeviceModelPKG/ADMS/N_DEV_ADMSbjt504va.h
diff -u -Bb utils/ADMS/examples/mextram_504.12.1/N_DEV_ADMSbjt504va.C src/DeviceModelPKG/ADMS/N_DEV_ADMSbjt504va.C
diff -u -Bb utils/ADMS/examples/mextram_504.12.1/N_DEV_ADMSbjt504tva.h src/DeviceModelPKG/ADMS/N_DEV_ADMSbjt504tva.h
diff -u -Bb utils/ADMS/examples/mextram_504.12.1/N_DEV_ADMSbjt504tva.C src/DeviceModelPKG/ADMS/N_DEV_ADMSbjt504tva.C
diff -u -Bb utils/ADMS/examples/mvs_2.0.0/N_DEV_ADMSmvs_2_0_0_etsoi.h src/DeviceModelPKG/ADMS/N_DEV_ADMSmvs_2_0_0_etsoi.h
diff -u -Bb utils/ADMS/examples/mvs_2.0.0/N_DEV_ADMSmvs_2_0_0_etsoi.C src/DeviceModelPKG/ADMS/N_DEV_ADMSmvs_2_0_0_etsoi.C
diff -u -Bb utils/ADMS/examples/mvs_2.0.0/N_DEV_ADMSmvs_2_0_0_hemt.h src/DeviceModelPKG/ADMS/N_DEV_ADMSmvs_2_0_0_hemt.h
diff -u -Bb utils/ADMS/examples/mvs_2.0.0/N_DEV_ADMSmvs_2_0_0_hemt.C src/DeviceModelPKG/ADMS/N_DEV_ADMSmvs_2_0_0_hemt.C
diff -u -Bb utils/ADMS/examples/mvsg_cmc_v1.1.0/vacode/N_DEV_ADMSmvsg_cmc.h src/DeviceModelPKG/ADMS/N_DEV_ADMSmvsg_cmc.h
diff -u -Bb utils/ADMS/examples/mvsg_cmc_v1.1.0/vacode/N_DEV_ADMSmvsg_cmc.C src/DeviceModelPKG/ADMS/N_DEV_ADMSmvsg_cmc.C
diff -u -Bb utils/ADMS/examples/hicum/N_DEV_ADMShicumL2va.h src/DeviceModelPKG/ADMS/N_DEV_ADMShicumL2va.h
diff -u -Bb utils/ADMS/examples/hicum/N_DEV_ADMShicumL2va.C src/DeviceModelPKG/ADMS/N_DEV_ADMShicumL2va.C
diff -u -Bb utils/ADMS/examples/hicum_l0/N_DEV_ADMShic0_full.h src/DeviceModelPKG/ADMS/N_DEV_ADMShic0_full.h
diff -u -Bb utils/ADMS/examples/hicum_l0/N_DEV_ADMShic0_full.C src/DeviceModelPKG/ADMS/N_DEV_ADMShic0_full.C
diff -u -Bb src/DeviceModelPKG/Xyce_NonFree/Verilog/ekv301_02/N_DEV_ADMSekv3.h src/DeviceModelPKG/Xyce_NonFree/N_DEV_ADMSekv3.h
diff -u -Bb src/DeviceModelPKG/Xyce_NonFree/Verilog/ekv301_02/N_DEV_ADMSekv3.C src/DeviceModelPKG/Xyce_NonFree/N_DEV_ADMSekv3.C
diff -u -Bb  src/DeviceModelPKG/Xyce_NonFree/Verilog/ekv2.6/N_DEV_ADMSekv_va.h src/DeviceModelPKG/Xyce_NonFree/N_DEV_ADMSekv_va.h
diff -u -Bb  src/DeviceModelPKG/Xyce_NonFree/Verilog/ekv2.6/N_DEV_ADMSekv_va.C src/DeviceModelPKG/Xyce_NonFree/N_DEV_ADMSekv_va.C
