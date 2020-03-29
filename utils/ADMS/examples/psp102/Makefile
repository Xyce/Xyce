ADMSDIR = ../..

PSPDEPS = JUNCAP200_InitModel.include \
	JUNCAP200_macrodefs.include \
	JUNCAP200_parlist.include \
	JUNCAP200_varlist1.include \
	JUNCAP200_varlist2.include \
	PSP102_ChargesNQS.include \
	PSP102_InitNQS.include \
	PSP102_binning.include \
	PSP102_binpars.include \
	PSP102_macrodefs.include \
	PSP102_module.include \
	PSP102_nqs_macrodefs.include \
	Common102_macrodefs.include \

XyceADMSFILES = $(ADMSDIR)/adms.implicit.xml \
	$(ADMSDIR)/xyceVersion_nosac.xml \
	$(ADMSDIR)/xyceBasicTemplates_nosac.xml \
	$(ADMSDIR)/xyceHeaderFile_nosac.xml \
	$(ADMSDIR)/xyceImplementationFile_nosac.xml \
	$(ADMSDIR)/xyceOutVarsDoc.xml

CLEANFILES = N_DEV_ADMSPSP102VA.C N_DEV_ADMSPSP102VA.h N_DEV_ADMSJUNCAP200.C N_DEV_ADMSJUNCAP200.h discipline.h *.tex

all-source: N_DEV_ADMSPSP102VA.C

clean:
	rm -f $(CLEANFILES)

N_DEV_ADMSPSP102VA.h N_DEV_ADMSPSP102VA.C: psp102.va $(PSPDEPS) $(XyceADMSFILES)
	admsXml -D_XYCE_HACK=1 -x -e $(ADMSDIR)/adms.implicit.xml -e $(ADMSDIR)/xyceVersion_nosac.xml -e $(ADMSDIR)/xyceBasicTemplates_nosac.xml -e $(ADMSDIR)/xyceImplementationFile_nosac.xml -e $(ADMSDIR)/xyceHeaderFile_nosac.xml 	-e $(ADMSDIR)/xyceOutVarsDoc.xml psp102.va
	mv M_102_OutputVars.tex M_102_OutputVars.tex_borken
	sed -e 's|lp_|lp\\_|g' < M_102_OutputVars.tex_borken > M_102_OutputVars.tex
	rm M_102_OutputVars.tex_borken
#	emacs N_DEV_ADMSPSP102VA.C --batch --eval="(require 'cc-mode)" --eval="(c-set-offset 'substatement-open 0)" --eval="(c-set-offset 'arglist-intro 3)" --eval="(c-set-offset 'innamespace -2)" --eval="(setq-default indent-tabs-mode nil)" --eval='(indent-region (point-min) (point-max) nil)' -f save-buffer
#	emacs N_DEV_ADMSPSP102VA.h --batch --eval="(require 'cc-mode)" --eval="(c-set-offset 'substatement-open 0)" --eval="(c-set-offset 'arglist-intro 3)" --eval="(c-set-offset 'innamespace -2)" --eval="(setq-default indent-tabs-mode nil)" --eval='(indent-region (point-min) (point-max) nil)' -f save-buffer
