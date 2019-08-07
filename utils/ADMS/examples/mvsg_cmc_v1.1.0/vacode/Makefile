ADMSDIR = ../../..

XyceADMSFILES = $(ADMSDIR)/xyceVersion.xml \
	$(ADMSDIR)/xyceBasicTemplates.xml \
	$(ADMSDIR)/xyceHeaderFile.xml \
	$(ADMSDIR)/xyceImplementationFile.xml

CLEANFILES = N_DEV_ADMSmvsg_cmc.h N_DEV_ADMSmvsg_cmc.C disciplines.vams constants.vams

all-source: N_DEV_ADMSmvsg_cmc.C

clean:
	rm -f $(CLEANFILES)

N_DEV_ADMSmvsg_cmc.h N_DEV_ADMSmvsg_cmc.C: mvsg_cmc_1.1.0.va $(XyceADMSFILES)
	admsXml -D__MVSG_HV_ORIG__=1 -e $(ADMSDIR)/xyceVersion.xml -e $(ADMSDIR)/xyceBasicTemplates.xml -e $(ADMSDIR)/xyceImplementationFile.xml -e $(ADMSDIR)/xyceHeaderFile.xml mvsg_cmc_1.1.0.va
	emacs N_DEV_ADMSmvsg_cmc.C --batch --eval="(require 'cc-mode)" --eval="(c-set-offset 'substatement-open 0)" --eval="(c-set-offset 'arglist-intro 3)" --eval="(c-set-offset 'innamespace -2)" --eval="(setq-default indent-tabs-mode nil)" --eval='(indent-region (point-min) (point-max) nil)' -f save-buffer
	emacs N_DEV_ADMSmvsg_cmc.h --batch --eval="(require 'cc-mode)" --eval="(c-set-offset 'substatement-open 0)" --eval="(c-set-offset 'arglist-intro 3)" --eval="(c-set-offset 'innamespace -2)" --eval="(setq-default indent-tabs-mode nil)" --eval='(indent-region (point-min) (point-max) nil)' -f save-buffer
