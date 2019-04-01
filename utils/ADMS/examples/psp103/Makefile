ADMSDIR = ../..

PSPDEPS = Common103_macrodefs.include \
	JUNCAP200_InitModel.include \
	JUNCAP200_macrodefs.include \
	JUNCAP200_parlist.include \
	JUNCAP200_varlist1.include \
	JUNCAP200_varlist2.include \
	PSP103_ChargesNQS.include \
	PSP103_InitNQS.include \
	PSP103_SPCalculation.include \
	PSP103_binning.include \
	PSP103_binpars.include \
	PSP103_macrodefs.include \
	PSP103_module.include \
	PSP103_nqs_macrodefs.include 

XyceADMSFILES = $(ADMSDIR)/xyceVersion.xml \
	$(ADMSDIR)/xyceBasicTemplates_nosac.xml \
	$(ADMSDIR)/xyceHeaderFile.xml \
	$(ADMSDIR)/xyceImplementationFile_nosac.xml

CLEANFILES = N_DEV_ADMSJUNCAP200.C N_DEV_ADMSPSP103VA.C N_DEV_ADMSPSP103TVA.C N_DEV_ADMSJUNCAP200.h N_DEV_ADMSPSP103VA.h N_DEV_ADMSPSP103TVA.h discipline.h

all-source: N_DEV_ADMSJUNCAP200.C N_DEV_ADMSPSP103VA.C N_DEV_ADMSPSP103TVA.C

clean:
	rm -f $(CLEANFILES)

N_DEV_ADMSPSP103VA.h N_DEV_ADMSPSP103VA.C: psp103.va $(PSPDEPS) $(XyceADMSFILES)
	admsXml -D_XYCE_HACK=1 -e $(ADMSDIR)/xyceVersion.xml -e $(ADMSDIR)/xyceBasicTemplates_nosac.xml -e $(ADMSDIR)/xyceImplementationFile_nosac.xml -e $(ADMSDIR)/xyceHeaderFile.xml psp103.va
	emacs N_DEV_ADMSPSP103VA.C --batch --eval="(require 'cc-mode)" --eval="(c-set-offset 'substatement-open 0)" --eval="(c-set-offset 'arglist-intro 3)" --eval="(c-set-offset 'innamespace -2)" --eval="(setq-default indent-tabs-mode nil)" --eval='(indent-region (point-min) (point-max) nil)' -f save-buffer
	emacs N_DEV_ADMSPSP103VA.h --batch --eval="(require 'cc-mode)" --eval="(c-set-offset 'substatement-open 0)" --eval="(c-set-offset 'arglist-intro 3)" --eval="(c-set-offset 'innamespace -2)" --eval="(setq-default indent-tabs-mode nil)" --eval='(indent-region (point-min) (point-max) nil)' -f save-buffer

N_DEV_ADMSJUNCAP200.h N_DEV_ADMSJUNCAP200.C: juncap200.va $(PSPDEPS)  $(XyceADMSFILES)
	admsXml -e $(ADMSDIR)/xyceVersion.xml -e $(ADMSDIR)/xyceBasicTemplates_nosac.xml -e $(ADMSDIR)/xyceImplementationFile_nosac.xml -e $(ADMSDIR)/xyceHeaderFile.xml juncap200.va
	emacs N_DEV_ADMSJUNCAP200.C --batch --eval="(require 'cc-mode)" --eval="(c-set-offset 'substatement-open 0)" --eval="(c-set-offset 'arglist-intro 3)" --eval="(c-set-offset 'innamespace -2)" --eval="(setq-default indent-tabs-mode nil)" --eval='(indent-region (point-min) (point-max) nil)' -f save-buffer
	emacs N_DEV_ADMSJUNCAP200.h --batch --eval="(require 'cc-mode)" --eval="(c-set-offset 'substatement-open 0)" --eval="(c-set-offset 'arglist-intro 3)" --eval="(c-set-offset 'innamespace -2)" --eval="(setq-default indent-tabs-mode nil)" --eval='(indent-region (point-min) (point-max) nil)' -f save-buffer

N_DEV_ADMSPSP103TVA.h N_DEV_ADMSPSP103TVA.C: psp103t.va $(PSPDEPS) $(XyceADMSFILES)
	admsXml -D_XYCE_HACK=1 -e $(ADMSDIR)/xyceVersion.xml -e $(ADMSDIR)/xyceBasicTemplates.xml -e $(ADMSDIR)/xyceImplementationFile.xml -e $(ADMSDIR)/xyceHeaderFile.xml psp103t.va
	emacs N_DEV_ADMSPSP103TVA.C --batch --eval="(require 'cc-mode)" --eval="(c-set-offset 'substatement-open 0)" --eval="(c-set-offset 'arglist-intro 3)" --eval="(c-set-offset 'innamespace -2)" --eval="(setq-default indent-tabs-mode nil)" --eval='(indent-region (point-min) (point-max) nil)' -f save-buffer
	emacs N_DEV_ADMSPSP103TVA.h --batch --eval="(require 'cc-mode)" --eval="(c-set-offset 'substatement-open 0)" --eval="(c-set-offset 'arglist-intro 3)" --eval="(c-set-offset 'innamespace -2)" --eval="(setq-default indent-tabs-mode nil)" --eval='(indent-region (point-min) (point-max) nil)' -f save-buffer
