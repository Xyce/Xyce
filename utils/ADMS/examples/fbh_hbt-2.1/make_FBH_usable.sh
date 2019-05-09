#!/bin/sh

admsXml -e ../../xyceVersion.xml -e ../../xyceBasicTemplates.xml -e ../../xyceImplementationFile.xml -e ../../xyceHeaderFile.xml fbhhbt-2.1_nonoise_limited_inductors_typed.va
emacs N_DEV_ADMSHBT_X.C --batch --eval="(require 'cc-mode)" --eval="(c-set-offset 'substatement-open 0)" --eval="(c-set-offset 'arglist-intro 3)" --eval="(c-set-offset 'innamespace -2)" --eval="(setq-default indent-tabs-mode nil)" --eval='(indent-region (point-min) (point-max) nil)' -f save-buffer
emacs N_DEV_ADMSHBT_X.h --batch --eval="(require 'cc-mode)" --eval="(c-set-offset 'substatement-open 0)" --eval="(c-set-offset 'arglist-intro 3)" --eval="(c-set-offset 'innamespace -2)" --eval="(setq-default indent-tabs-mode nil)" --eval='(indent-region (point-min) (point-max) nil)' -f save-buffer
patch < limited_typed_source.patch
patch < limited_typed_header.patch

# cp N_DEV_ADMSHBT_X.h ../../../../src/DeviceModelPKG/ADMS/include/
# cp N_DEV_ADMSHBT_X.C ../../../../src/DeviceModelPKG/ADMS/src/
