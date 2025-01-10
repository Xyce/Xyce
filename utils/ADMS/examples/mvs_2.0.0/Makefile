ADMSDIR=../..

XyceADMSFILES = $(ADMSDIR)/adms.implicit.xml \
	$(ADMSDIR)/xyceVersion_nosac.xml \
	$(ADMSDIR)/xyceBasicTemplates_nosac.xml \
	$(ADMSDIR)/xyceAnalogFunction_nosac.xml \
        $(ADMSDIR)/xyceHeaderFile_nosac.xml \
        $(ADMSDIR)/xyceImplementationFile_nosac.xml


all-source: N_DEV_ADMSmvs_2_0_0_hemt.C N_DEV_ADMSmvs_2_0_0_etsoi.C

clean:
	rm -f N_DEV_ADMSmvs_2_0_0_hemt.C N_DEV_ADMSmvs_2_0_0_etsoi.C N_DEV_ADMSmvs_2_0_0_hemt.h N_DEV_ADMSmvs_2_0_0_etsoi.h *.o *.so*

N_DEV_ADMSmvs_2_0_0_hemt.C: mvs_2_0_0_hemt.va  $(XyceADMSFILES)
	admsXml -x -e $(ADMSDIR)/adms.implicit.xml  -e $(ADMSDIR)/xyceVersion_nosac.xml -e $(ADMSDIR)/xyceBasicTemplates_nosac.xml -e $(ADMSDIR)/xyceAnalogFunction_nosac.xml -e $(ADMSDIR)/xyceHeaderFile_nosac.xml -e $(ADMSDIR)/xyceImplementationFile_nosac.xml mvs_2_0_0_hemt.va
#	emacs N_DEV_ADMSmvs_2_0_0_hemt.C --batch --eval="(require 'cc-mode)" --eval="(c-set-offset 'substatement-open 0)" --eval="(c-set-offset 'arglist-intro 3)" --eval="(c-set-offset 'innamespace -2)" --eval="(setq-default indent-tabs-mode nil)" --eval='(indent-region (point-min) (point-max) nil)' -f save-buffer

N_DEV_ADMSmvs_2_0_0_etsoi.C: mvs_2_0_0_etsoi.va  $(XyceADMSFILES)
	admsXml -x -e $(ADMSDIR)/adms.implicit.xml  -e $(ADMSDIR)/xyceVersion_nosac.xml -e $(ADMSDIR)/xyceBasicTemplates_nosac.xml -e $(ADMSDIR)/xyceAnalogFunction_nosac.xml -e $(ADMSDIR)/xyceHeaderFile_nosac.xml -e $(ADMSDIR)/xyceImplementationFile_nosac.xml mvs_2_0_0_etsoi.va
#	emacs N_DEV_ADMSmvs_2_0_0_etsoi.C --batch --eval="(require 'cc-mode)" --eval="(c-set-offset 'substatement-open 0)" --eval="(c-set-offset 'arglist-intro 3)" --eval="(c-set-offset 'innamespace -2)" --eval="(setq-default indent-tabs-mode nil)" --eval='(indent-region (point-min) (point-max) nil)' -f save-buffer

