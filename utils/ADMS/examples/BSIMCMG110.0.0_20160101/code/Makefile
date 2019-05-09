ADMSDIR = ../../..

XyceADMSFILES = $(ADMSDIR)/xyceVersion.xml \
	$(ADMSDIR)/xyceBasicTemplates.xml \
	$(ADMSDIR)/xyceHeaderFile.xml \
	$(ADMSDIR)/xyceImplementationFile.xml

BSIMCMGDEPFILES =bsimcmg_binning_parameters.include \
	bsimcmg_body.include \
	bsimcmg_cfringe.include \
	bsimcmg_quasi_static_cv.include \
	bsimcmg_rdsmod.include \
	common_defs.include 

CLEANFILES = N_DEV_ADMSbsimcmg_110.h N_DEV_ADMSbsimcmg_110.C disciplines.vams constants.vams

all-source: N_DEV_ADMSbsimcmg_110.C

clean:
	rm -f $(CLEANFILES)

N_DEV_ADMSbsimcmg_110.h N_DEV_ADMSbsimcmg_110.C: bsimcmg.va $(XyceADMSFILES) $(BSIMCMGDEPFILES)
	admsXml -D__XYCE__=1 -e $(ADMSDIR)/xyceVersion.xml -e $(ADMSDIR)/xyceBasicTemplates.xml -e $(ADMSDIR)/xyceImplementationFile.xml -e $(ADMSDIR)/xyceHeaderFile.xml bsimcmg.va
	emacs N_DEV_ADMSbsimcmg_110.C --batch --eval="(require 'cc-mode)" --eval="(c-set-offset 'substatement-open 0)" --eval="(c-set-offset 'arglist-intro 3)" --eval="(c-set-offset 'innamespace -2)" --eval="(setq-default indent-tabs-mode nil)" --eval='(indent-region (point-min) (point-max) nil)' -f save-buffer
	emacs N_DEV_ADMSbsimcmg_110.h --batch --eval="(require 'cc-mode)" --eval="(c-set-offset 'substatement-open 0)" --eval="(c-set-offset 'arglist-intro 3)" --eval="(c-set-offset 'innamespace -2)" --eval="(setq-default indent-tabs-mode nil)" --eval='(indent-region (point-min) (point-max) nil)' -f save-buffer
	patch < make_bsimcmg_usable.diff
