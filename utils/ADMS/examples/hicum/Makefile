ADMSDIR=../..

XyceADMSFILES = $(ADMSDIR)/xyceVersion_nosac.xml \
        $(ADMSDIR)/xyceBasicTemplates_nosac.xml \
        $(ADMSDIR)/xyceHeaderFile.xml \
        $(ADMSDIR)/xyceImplementationFile_nosac.xml


all-source: N_DEV_ADMShicumL2va.C

clean:
	rm -f N_DEV_ADMShicumL2va.C *.o *.so*

N_DEV_ADMShicumL2va.C: hicumL2V2p4p0.va  $(XyceADMSFILES)
	admsXml  -e $(ADMSDIR)/xyceVersion_nosac.xml -e $(ADMSDIR)/xyceBasicTemplates_nosac.xml -e $(ADMSDIR)/xyceHeaderFile.xml -e $(ADMSDIR)/xyceImplementationFile_nosac.xml hicumL2V2p4p0.va
	emacs N_DEV_ADMShicumL2va.C --batch --eval="(require 'cc-mode)" --eval="(c-set-offset 'substatement-open 0)" --eval="(c-set-offset 'arglist-intro 3)" --eval="(c-set-offset 'innamespace -2)" --eval="(setq-default indent-tabs-mode nil)" --eval='(indent-region (point-min) (point-max) nil)' -f save-buffer
	emacs N_DEV_ADMShicumL2va.h --batch --eval="(require 'cc-mode)" --eval="(c-set-offset 'substatement-open 0)" --eval="(c-set-offset 'arglist-intro 3)" --eval="(c-set-offset 'innamespace -2)" --eval="(setq-default indent-tabs-mode nil)" --eval='(indent-region (point-min) (point-max) nil)' -f save-buffer

