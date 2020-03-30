#!/bin/sh

admsXml -e ../../../../xyceVersion.xml -e ../../../../xyceBasicTemplates.xml -e ../../../../xyceImplementationFile.xml -e ../../../../xyceHeaderFile.xml vbic_3T_et_cf_Xyce_typed.vla
emacs N_DEV_ADMSvbic.C --batch --eval="(require 'cc-mode)" --eval="(c-set-offset 'substatement-open 0)" --eval="(c-set-offset 'arglist-intro 3)" --eval="(c-set-offset 'innamespace -2)" --eval="(setq-default indent-tabs-mode nil)" --eval='(indent-region (point-min) (point-max) nil)' -f save-buffer
emacs N_DEV_ADMSvbic.h --batch --eval="(require 'cc-mode)" --eval="(c-set-offset 'substatement-open 0)" --eval="(c-set-offset 'arglist-intro 3)" --eval="(c-set-offset 'innamespace -2)" --eval="(setq-default indent-tabs-mode nil)" --eval='(indent-region (point-min) (point-max) nil)' -f save-buffer
patch --ignore-whitespace < make_typed_vbic_usable.diff

# cp N_DEV_ADMSvbic.h ../../../../../../src/DeviceModelPKG/ADMS/include/
# cp N_DEV_ADMSvbic.C ../../../../../../src/DeviceModelPKG/ADMS/src/
