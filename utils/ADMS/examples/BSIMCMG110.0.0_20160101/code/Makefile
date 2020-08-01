ADMSDIR = ../../..

XyceADMSFILES = $(ADMSDIR)/adms.implicit.xml \
	$(ADMSDIR)/xyceVersion_nosac.xml \
	$(ADMSDIR)/xyceBasicTemplates_nosac.xml \
	$(ADMSDIR)/xyceAnalogFunction_nosac.xml \
	$(ADMSDIR)/xyceHeaderFile_nosac.xml \
	$(ADMSDIR)/xyceImplementationFile_nosac.xml \
	$(ADMSDIR)/xyceOutVarsDoc.xml

BSIMCMGDEPFILES =bsimcmg_binning_parameters.include \
	bsimcmg_body.include \
	bsimcmg_cfringe.include \
	bsimcmg_quasi_static_cv.include \
	bsimcmg_rdsmod.include \
	common_defs.include 

CLEANFILES = N_DEV_ADMSbsimcmg_110.h N_DEV_ADMSbsimcmg_110.C disciplines.vams constants.vams *.orig *.rej *.tex

all-source: N_DEV_ADMSbsimcmg_110.C

clean:
	rm -f $(CLEANFILES)

N_DEV_ADMSbsimcmg_110.h N_DEV_ADMSbsimcmg_110.C: bsimcmg.va $(XyceADMSFILES) $(BSIMCMGDEPFILES)
	admsXml -D__XYCE__=1 -x -e $(ADMSDIR)/adms.implicit.xml -e $(ADMSDIR)/xyceVersion_nosac.xml -e $(ADMSDIR)/xyceBasicTemplates_nosac.xml -e $(ADMSDIR)/xyceAnalogFunction_nosac.xml -e $(ADMSDIR)/xyceImplementationFile_nosac.xml -e $(ADMSDIR)/xyceHeaderFile_nosac.xml -e $(ADMSDIR)/xyceOutVarsDoc.xml bsimcmg.va
	mv M_110_OutputVars.tex M_110_OutputVars.tex_borken
	sed -e 's|\([TLA]\)_|\1\\_|g' < M_110_OutputVars.tex_borken > M_110_OutputVars.tex
	rm M_110_OutputVars.tex_borken
#	emacs N_DEV_ADMSbsimcmg_110.C --batch --eval="(require 'cc-mode)" --eval="(c-set-offset 'substatement-open 0)" --eval="(c-set-offset 'arglist-intro 3)" --eval="(c-set-offset 'innamespace -2)" --eval="(setq-default indent-tabs-mode nil)" --eval='(indent-region (point-min) (point-max) nil)' -f save-buffer
#	emacs N_DEV_ADMSbsimcmg_110.h --batch --eval="(require 'cc-mode)" --eval="(c-set-offset 'substatement-open 0)" --eval="(c-set-offset 'arglist-intro 3)" --eval="(c-set-offset 'innamespace -2)" --eval="(setq-default indent-tabs-mode nil)" --eval='(indent-region (point-min) (point-max) nil)' -f save-buffer
