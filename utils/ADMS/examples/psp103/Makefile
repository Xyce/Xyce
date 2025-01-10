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

XyceADMSFILES = $(ADMSDIR)/adms.implicit.xml \
	$(ADMSDIR)/xyceVersion_nosac.xml \
	$(ADMSDIR)/xyceBasicTemplates_nosac.xml \
	$(ADMSDIR)/xyceHeaderFile_nosac.xml \
	$(ADMSDIR)/xyceImplementationFile_nosac.xml \
	$(ADMSDIR)/xyceOutVarsDoc.xml

CLEANFILES = N_DEV_ADMSJUNCAP200.C N_DEV_ADMSPSP103VA.C N_DEV_ADMSPSP103TVA.C N_DEV_ADMSJUNCAP200.h N_DEV_ADMSPSP103VA.h N_DEV_ADMSPSP103TVA.h discipline.h *.tex

all-source: N_DEV_ADMSJUNCAP200.C N_DEV_ADMSPSP103VA.C N_DEV_ADMSPSP103TVA.C

clean:
	rm -f $(CLEANFILES)

N_DEV_ADMSPSP103VA.h N_DEV_ADMSPSP103VA.C: psp103.va $(PSPDEPS) $(XyceADMSFILES)
	admsXml -D_XYCE_HACK=1 -x -e $(ADMSDIR)/adms.implicit.xml -e $(ADMSDIR)/xyceVersion_nosac.xml -e $(ADMSDIR)/xyceBasicTemplates_nosac.xml -e $(ADMSDIR)/xyceImplementationFile_nosac.xml -e $(ADMSDIR)/xyceHeaderFile_nosac.xml -e $(ADMSDIR)/xyceOutVarsDoc.xml psp103.va
	mv M_103_OutputVars.tex M_103_OutputVars.tex_borken
	sed -e 's|lp_|lp\\_|g' < M_103_OutputVars.tex_borken > M_103_OutputVars.tex
	rm M_103_OutputVars.tex_borken
#	emacs N_DEV_ADMSPSP103VA.C --batch --eval="(require 'cc-mode)" --eval="(c-set-offset 'substatement-open 0)" --eval="(c-set-offset 'arglist-intro 3)" --eval="(c-set-offset 'innamespace -2)" --eval="(setq-default indent-tabs-mode nil)" --eval='(indent-region (point-min) (point-max) nil)' -f save-buffer
#	emacs N_DEV_ADMSPSP103VA.h --batch --eval="(require 'cc-mode)" --eval="(c-set-offset 'substatement-open 0)" --eval="(c-set-offset 'arglist-intro 3)" --eval="(c-set-offset 'innamespace -2)" --eval="(setq-default indent-tabs-mode nil)" --eval='(indent-region (point-min) (point-max) nil)' -f save-buffer

N_DEV_ADMSJUNCAP200.h N_DEV_ADMSJUNCAP200.C: juncap200.va $(PSPDEPS)  $(XyceADMSFILES)
	admsXml -x -e $(ADMSDIR)/adms.implicit.xml -e $(ADMSDIR)/xyceVersion_nosac.xml -e $(ADMSDIR)/xyceBasicTemplates_nosac.xml -e $(ADMSDIR)/xyceImplementationFile_nosac.xml -e $(ADMSDIR)/xyceHeaderFile_nosac.xml -e $(ADMSDIR)/xyceOutVarsDoc.xml juncap200.va
	mv D_200_OutputVars.tex D_200_OutputVars.tex_borken
	sed -e 's|lp_|lp\\_|g' < D_200_OutputVars.tex_borken > D_200_OutputVars.tex
	rm D_200_OutputVars.tex_borken
#	emacs N_DEV_ADMSJUNCAP200.C --batch --eval="(require 'cc-mode)" --eval="(c-set-offset 'substatement-open 0)" --eval="(c-set-offset 'arglist-intro 3)" --eval="(c-set-offset 'innamespace -2)" --eval="(setq-default indent-tabs-mode nil)" --eval='(indent-region (point-min) (point-max) nil)' -f save-buffer
#	emacs N_DEV_ADMSJUNCAP200.h --batch --eval="(require 'cc-mode)" --eval="(c-set-offset 'substatement-open 0)" --eval="(c-set-offset 'arglist-intro 3)" --eval="(c-set-offset 'innamespace -2)" --eval="(setq-default indent-tabs-mode nil)" --eval='(indent-region (point-min) (point-max) nil)' -f save-buffer

N_DEV_ADMSPSP103TVA.h N_DEV_ADMSPSP103TVA.C: psp103t.va $(PSPDEPS) $(XyceADMSFILES)
	admsXml -x -e $(ADMSDIR)/adms.implicit.xml -D_XYCE_HACK=1 -e $(ADMSDIR)/xyceVersion_nosac.xml -e $(ADMSDIR)/xyceBasicTemplates_nosac.xml -e $(ADMSDIR)/xyceImplementationFile_nosac.xml -e $(ADMSDIR)/xyceHeaderFile_nosac.xml -e $(ADMSDIR)/xyceOutVarsDoc.xml psp103t.va
	mv M_1031_OutputVars.tex M_1031_OutputVars.tex_borken
	sed -e 's|lp_|lp\\_|g' < M_1031_OutputVars.tex_borken > M_1031_OutputVars.tex
	rm M_1031_OutputVars.tex_borken
#	emacs N_DEV_ADMSPSP103TVA.C --batch --eval="(require 'cc-mode)" --eval="(c-set-offset 'substatement-open 0)" --eval="(c-set-offset 'arglist-intro 3)" --eval="(c-set-offset 'innamespace -2)" --eval="(setq-default indent-tabs-mode nil)" --eval='(indent-region (point-min) (point-max) nil)' -f save-buffer
#	emacs N_DEV_ADMSPSP103TVA.h --batch --eval="(require 'cc-mode)" --eval="(c-set-offset 'substatement-open 0)" --eval="(c-set-offset 'arglist-intro 3)" --eval="(c-set-offset 'innamespace -2)" --eval="(setq-default indent-tabs-mode nil)" --eval='(indent-region (point-min) (point-max) nil)' -f save-buffer
