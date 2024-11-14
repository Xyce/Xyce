ADMSDIR=../..

XyceADMSFILES = $(ADMSDIR)/adms.implicit.xml \
	$(ADMSDIR)/xyceVersion_nosac.xml \
        $(ADMSDIR)/xyceBasicTemplates_nosac.xml \
        $(ADMSDIR)/xyceAnalogFunction_nosac.xml \
        $(ADMSDIR)/xyceHeaderFile_nosac.xml \
        $(ADMSDIR)/xyceImplementationFile_nosac.xml \
        $(ADMSDIR)/xyceOutVarsDoc.xml


all-source: N_DEV_ADMShicumL2va.C

clean:
	rm -f N_DEV_ADMShicumL2va.C N_DEV_ADMShicumL2va.h *.o *.so* *.tex

N_DEV_ADMShicumL2va.C: hicumL2V2p4p0.va  $(XyceADMSFILES)
	admsXml -x -e $(ADMSDIR)/adms.implicit.xml -e $(ADMSDIR)/xyceVersion_nosac.xml -e $(ADMSDIR)/xyceBasicTemplates_nosac.xml -e $(ADMSDIR)/xyceAnalogFunction_nosac.xml -e $(ADMSDIR)/xyceImplementationFile_nosac.xml  -e $(ADMSDIR)/xyceHeaderFile_nosac.xml -e $(ADMSDIR)/xyceOutVarsDoc.xml hicumL2V2p4p0.va
	mv Q_234_OutputVars.tex Q_234_OutputVars.tex_borken
	sed -e 's|\(r[a-z]*\)_|\1\\_|g' < Q_234_OutputVars.tex_borken > Q_234_OutputVars.tex
	rm Q_234_OutputVars.tex_borken
#	emacs N_DEV_ADMShicumL2va.C --batch --eval="(require 'cc-mode)" --eval="(c-set-offset 'substatement-open 0)" --eval="(c-set-offset 'arglist-intro 3)" --eval="(c-set-offset 'innamespace -2)" --eval="(setq-default indent-tabs-mode nil)" --eval='(indent-region (point-min) (point-max) nil)' -f save-buffer
#	emacs N_DEV_ADMShicumL2va.h --batch --eval="(require 'cc-mode)" --eval="(c-set-offset 'substatement-open 0)" --eval="(c-set-offset 'arglist-intro 3)" --eval="(c-set-offset 'innamespace -2)" --eval="(setq-default indent-tabs-mode nil)" --eval='(indent-region (point-min) (point-max) nil)' -f save-buffer

