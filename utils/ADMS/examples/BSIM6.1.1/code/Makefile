ADMSDIR = ../../..

XyceADMSFILES = $(ADMSDIR)/adms.implicit.xml \
	$(ADMSDIR)/xyceVersion.xml \
	$(ADMSDIR)/xyceBasicTemplates.xml \
	$(ADMSDIR)/xyceHeaderFile.xml \
	$(ADMSDIR)/xyceImplementationFile.xml

CLEANFILES = N_DEV_ADMSbsim6.h N_DEV_ADMSbsim6.C disciplines.vams constants.vams

all-source: N_DEV_ADMSbsim6.C

clean:
	rm -f $(CLEANFILES)

N_DEV_ADMSbsim6.h N_DEV_ADMSbsim6.C: BSIM6.1.1.va $(XyceADMSFILES)
	admsXml -D__XYCE__=1 -x -e $(ADMSDIR)/adms.implicit.xml -e $(ADMSDIR)/xyceVersion.xml -e $(ADMSDIR)/xyceBasicTemplates.xml -e $(ADMSDIR)/xyceImplementationFile.xml -e $(ADMSDIR)/xyceHeaderFile.xml BSIM6.1.1.va
	emacs N_DEV_ADMSbsim6.C --batch --eval="(require 'cc-mode)" --eval="(c-set-offset 'substatement-open 0)" --eval="(c-set-offset 'arglist-intro 3)" --eval="(c-set-offset 'innamespace -2)" --eval="(setq-default indent-tabs-mode nil)" --eval='(indent-region (point-min) (point-max) nil)' -f save-buffer
	emacs N_DEV_ADMSbsim6.h --batch --eval="(require 'cc-mode)" --eval="(c-set-offset 'substatement-open 0)" --eval="(c-set-offset 'arglist-intro 3)" --eval="(c-set-offset 'innamespace -2)" --eval="(setq-default indent-tabs-mode nil)" --eval='(indent-region (point-min) (point-max) nil)' -f save-buffer
