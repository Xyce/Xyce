ADMSDIR=../..

XyceADMSFILES = $(ADMSDIR)/adms.implicit.xml \
	$(ADMSDIR)/xyceVersion_nosac.xml \
        $(ADMSDIR)/xyceBasicTemplates_nosac.xml \
        $(ADMSDIR)/xyceAnalogFunction_nosac.xml \
        $(ADMSDIR)/xyceHeaderFile_nosac.xml \
        $(ADMSDIR)/xyceImplementationFile_nosac.xml


all-source: N_DEV_ADMShic0_full.C

clean:
	rm -f N_DEV_ADMShic0_full.C N_DEV_ADMShic0_full.h *.o *.so*

N_DEV_ADMShic0_full.C: hicumL0V1p32.va  $(XyceADMSFILES)
	admsXml -x -e $(ADMSDIR)/adms.implicit.xml  -e $(ADMSDIR)/xyceVersion_nosac.xml -e $(ADMSDIR)/xyceBasicTemplates_nosac.xml -e $(ADMSDIR)/xyceAnalogFunction_nosac.xml -e $(ADMSDIR)/xyceHeaderFile_nosac.xml -e $(ADMSDIR)/xyceImplementationFile_nosac.xml hicumL0V1p32.va
#	emacs N_DEV_ADMShic0_full.C --batch --eval="(require 'cc-mode)" --eval="(c-set-offset 'substatement-open 0)" --eval="(c-set-offset 'arglist-intro 3)" --eval="(c-set-offset 'innamespace -2)" --eval="(setq-default indent-tabs-mode nil)" --eval='(indent-region (point-min) (point-max) nil)' -f save-buffer
#	emacs N_DEV_ADMShic0_full.h --batch --eval="(require 'cc-mode)" --eval="(c-set-offset 'substatement-open 0)" --eval="(c-set-offset 'arglist-intro 3)" --eval="(c-set-offset 'innamespace -2)" --eval="(setq-default indent-tabs-mode nil)" --eval='(indent-region (point-min) (point-max) nil)' -f save-buffer

