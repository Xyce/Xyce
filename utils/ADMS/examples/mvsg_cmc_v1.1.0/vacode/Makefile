ADMSDIR = ../../..

XyceADMSFILES = $(ADMSDIR)/adms.implicit.xml \
	$(ADMSDIR)/xyceVersion_nosac.xml \
	$(ADMSDIR)/xyceBasicTemplates_nosac.xml \
	$(ADMSDIR)/xyceAnalogFunction_nosac.xml \
	$(ADMSDIR)/xyceHeaderFile_nosac.xml \
	$(ADMSDIR)/xyceImplementationFile_nosac.xml \
	$(ADMSDIR)/xyceOutVarsDoc.xml

CLEANFILES = N_DEV_ADMSmvsg_cmc.h N_DEV_ADMSmvsg_cmc.C disciplines.vams constants.vams *.tex

all-source: N_DEV_ADMSmvsg_cmc.C

clean:
	rm -f $(CLEANFILES)

N_DEV_ADMSmvsg_cmc.h N_DEV_ADMSmvsg_cmc.C: mvsg_cmc_1.1.0.va $(XyceADMSFILES)
	admsXml -x -e $(ADMSDIR)/adms.implicit.xml -e $(ADMSDIR)/xyceVersion_nosac.xml -e $(ADMSDIR)/xyceBasicTemplates_nosac.xml -e $(ADMSDIR)/xyceAnalogFunction_nosac.xml -e $(ADMSDIR)/xyceImplementationFile_nosac.xml -e $(ADMSDIR)/xyceHeaderFile_nosac.xml -e $(ADMSDIR)/xyceOutVarsDoc.xml mvsg_cmc_1.1.0.va
	mv M_2002_OutputVars.tex M_2002_OutputVars.tex_borken
	sed -e 's|\([tla]\)_|\1\\_|g' < M_2002_OutputVars.tex_borken > M_2002_OutputVars.tex
	rm M_2002_OutputVars.tex_borken
#	emacs N_DEV_ADMSmvsg_cmc.C --batch --eval="(require 'cc-mode)" --eval="(c-set-offset 'substatement-open 0)" --eval="(c-set-offset 'arglist-intro 3)" --eval="(c-set-offset 'innamespace -2)" --eval="(setq-default indent-tabs-mode nil)" --eval='(indent-region (point-min) (point-max) nil)' -f save-buffer
#	emacs N_DEV_ADMSmvsg_cmc.h --batch --eval="(require 'cc-mode)" --eval="(c-set-offset 'substatement-open 0)" --eval="(c-set-offset 'arglist-intro 3)" --eval="(c-set-offset 'innamespace -2)" --eval="(setq-default indent-tabs-mode nil)" --eval='(indent-region (point-min) (point-max) nil)' -f save-buffer
