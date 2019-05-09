ADMSDIR = ../..

MEXTRAMDEPS = evaluate.inc \
        frontdef.inc \
        initialize.inc \
        opinfo.inc \
        opvars.inc \
        parameters.inc \
        tscaling.inc \
        variables.inc


XyceADMSFILES = $(ADMSDIR)/xyceVersion.xml \
        $(ADMSDIR)/xyceBasicTemplates.xml \
        $(ADMSDIR)/xyceHeaderFile.xml \
        $(ADMSDIR)/xyceImplementationFile.xml

CLEANFILES = N_DEV_ADMSbjt504tva.C N_DEV_ADMSbjt504tva.h  N_DEV_ADMSbjt504va.C N_DEV_ADMSbjt504va.h 

all-source: N_DEV_ADMSbjt504tva.C N_DEV_ADMSbjt504va.C 

clean:
	rm -f $(CLEANFILES)

N_DEV_ADMSbjt504va.h N_DEV_ADMSbjt504va.C: bjt504.va $(MEXTRAMDEPS) $(XyceADMSFILES)
	admsXml -e ../../xyceVersion.xml -e ../../xyceBasicTemplates.xml -e ../../xyceImplementationFile.xml -e ../../xyceHeaderFile.xml bjt504.va
	emacs N_DEV_ADMSbjt504va.C --batch --eval="(require 'cc-mode)" --eval="(c-set-offset 'substatement-open 0)" --eval="(c-set-offset 'arglist-intro 3)" --eval="(c-set-offset 'innamespace -2)" --eval="(setq-default indent-tabs-mode nil)" --eval='(indent-region (point-min) (point-max) nil)' -f save-buffer
	emacs N_DEV_ADMSbjt504va.h --batch --eval="(require 'cc-mode)" --eval="(c-set-offset 'substatement-open 0)" --eval="(c-set-offset 'arglist-intro 3)" --eval="(c-set-offset 'innamespace -2)" --eval="(setq-default indent-tabs-mode nil)" --eval='(indent-region (point-min) (point-max) nil)' -f save-buffer

N_DEV_ADMSbjt504tva.h N_DEV_ADMSbjt504tva.C: bjt504t.va $(MEXTRAMDEPS) $(XyceADMSFILES)
	admsXml -e ../../xyceVersion.xml -e ../../xyceBasicTemplates.xml -e ../../xyceImplementationFile.xml -e ../../xyceHeaderFile.xml bjt504t.va
	emacs N_DEV_ADMSbjt504tva.C --batch --eval="(require 'cc-mode)" --eval="(c-set-offset 'substatement-open 0)" --eval="(c-set-offset 'arglist-intro 3)" --eval="(c-set-offset 'innamespace -2)" --eval="(setq-default indent-tabs-mode nil)" --eval='(indent-region (point-min) (point-max) nil)' -f save-buffer
	emacs N_DEV_ADMSbjt504tva.h --batch --eval="(require 'cc-mode)" --eval="(c-set-offset 'substatement-open 0)" --eval="(c-set-offset 'arglist-intro 3)" --eval="(c-set-offset 'innamespace -2)" --eval="(setq-default indent-tabs-mode nil)" --eval='(indent-region (point-min) (point-max) nil)' -f save-buffer

