#!/bin/sh

admsXml -e ../../xyceVersion.xml -e ../../xyceBasicTemplates.xml -e ../../xyceImplementationFile.xml -e ../../xyceHeaderFile.xml bjt504.va
emacs N_DEV_ADMSbjt504va.C --batch --eval="(require 'cc-mode)" --eval="(c-set-offset 'substatement-open 0)" --eval="(c-set-offset 'arglist-intro 3)" --eval="(c-set-offset 'innamespace -2)" --eval="(setq-default indent-tabs-mode nil)" --eval='(indent-region (point-min) (point-max) nil)' -f save-buffer
emacs N_DEV_ADMSbjt504va.h --batch --eval="(require 'cc-mode)" --eval="(c-set-offset 'substatement-open 0)" --eval="(c-set-offset 'arglist-intro 3)" --eval="(c-set-offset 'innamespace -2)" --eval="(setq-default indent-tabs-mode nil)" --eval='(indent-region (point-min) (point-max) nil)' -f save-buffer
