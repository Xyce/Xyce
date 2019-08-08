ADMSDIR = ../..

XyceADMSDeps = $(ADMSDIR)/adms.implicit.xml \
	$(ADMSDIR)/xyceVersion_nosac.xml \
	$(ADMSDIR)/xyceBasicTemplates_nosac.xml \
	$(ADMSDIR)/xyceHeaderFile.xml \
	$(ADMSDIR)/xyceImplementationFile_nosac.xml
all-source: N_DEV_ADMSvbic13_4t.C N_DEV_ADMSvbic13.C

clean:
	rm -f N_DEV_ADMSvbic13.C N_DEV_ADMSvbic13_4t.C N_DEV_ADMSvbic13.h N_DEV_ADMSvbic13_4t.h *.o *.so*

N_DEV_ADMSvbic13_4t.C: vbic_1p3.va cmcGeneralMacrosAndDefines.va cmcStandardModelMacros.va make_vbic13_usable.diff make_vbic13_headers_usable.diff $(XyceADMSDeps)
	admsXml -D__XYCE_COMPACT_MODELING__=1 -x -e $(ADMSDIR)/adms.implicit.xml -e $(ADMSDIR)/xyceVersion_nosac.xml -e $(ADMSDIR)/xyceBasicTemplates_nosac.xml -e $(ADMSDIR)/xyceImplementationFile_nosac.xml  -e $(ADMSDIR)/xyceHeaderFile.xml vbic_1p3.va
	emacs N_DEV_ADMSvbic13_4t.C --batch --eval="(require 'cc-mode)" --eval="(c-set-offset 'substatement-open 0)" --eval="(c-set-offset 'arglist-intro 3)" --eval="(c-set-offset 'innamespace -2)" --eval="(setq-default indent-tabs-mode nil)" --eval='(indent-region (point-min) (point-max) nil)' -f save-buffer
	emacs N_DEV_ADMSvbic13_4t.h --batch --eval="(require 'cc-mode)" --eval="(c-set-offset 'substatement-open 0)" --eval="(c-set-offset 'arglist-intro 3)" --eval="(c-set-offset 'innamespace -2)" --eval="(setq-default indent-tabs-mode nil)" --eval='(indent-region (point-min) (point-max) nil)' -f save-buffer
	patch N_DEV_ADMSvbic13_4t.C make_vbic13_usable.diff
	patch N_DEV_ADMSvbic13_4t.h make_vbic13_headers_usable.diff

N_DEV_ADMSvbic13.C: vbic_1p3.va cmcGeneralMacrosAndDefines.va cmcStandardModelMacros.va make_vbic13_usable.diff make_vbic13_3t_header_usable.diff $(XyceADMSDeps)
	admsXml -D__XYCE_COMPACT_MODELING__=1 -DthreeTerminal=1  -x -e $(ADMSDIR)/adms.implicit.xml -e $(ADMSDIR)/xyceVersion_nosac.xml -e $(ADMSDIR)/xyceBasicTemplates_nosac.xml -e $(ADMSDIR)/xyceImplementationFile_nosac.xml -e $(ADMSDIR)/xyceHeaderFile.xml vbic_1p3.va
	emacs N_DEV_ADMSvbic13.C --batch --eval="(require 'cc-mode)" --eval="(c-set-offset 'substatement-open 0)" --eval="(c-set-offset 'arglist-intro 3)" --eval="(c-set-offset 'innamespace -2)" --eval="(setq-default indent-tabs-mode nil)" --eval='(indent-region (point-min) (point-max) nil)' -f save-buffer
	emacs N_DEV_ADMSvbic13.h --batch --eval="(require 'cc-mode)" --eval="(c-set-offset 'substatement-open 0)" --eval="(c-set-offset 'arglist-intro 3)" --eval="(c-set-offset 'innamespace -2)" --eval="(setq-default indent-tabs-mode nil)" --eval='(indent-region (point-min) (point-max) nil)' -f save-buffer
	patch N_DEV_ADMSvbic13.C make_vbic13_usable.diff
	patch N_DEV_ADMSvbic13.h make_vbic13_3t_header_usable.diff

