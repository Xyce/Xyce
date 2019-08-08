ADMSDIR = ../..

XyceADMSFILES = $(ADMSDIR)/adms.implicit.xml \
	$(ADMSDIR)/xyceVersion.xml \
	$(ADMSDIR)/xyceBasicTemplates.xml \
	$(ADMSDIR)/xyceHeaderFile.xml \
	$(ADMSDIR)/xyceImplementationFile.xml

CLEANFILES = N_DEV_ADMSHBT_X.C N_DEV_ADMSHBT_X.h

all-source:  N_DEV_ADMSHBT_X.C 

clean:
	rm -f $(CLEANFILES)

N_DEV_ADMSHBT_X.h N_DEV_ADMSHBT_X.C: fbhhbt-2.1_nonoise_limited_inductors_typed.va $(XyceADMSFILES) limited_typed_source.patch limited_typed_header.patch
	admsXml -x -e $(ADMSDIR)/adms.implicit.xml -e $(ADMSDIR)/xyceVersion.xml -e $(ADMSDIR)/xyceBasicTemplates.xml -e $(ADMSDIR)/xyceImplementationFile.xml -e $(ADMSDIR)/xyceHeaderFile.xml fbhhbt-2.1_nonoise_limited_inductors_typed.va
	emacs N_DEV_ADMSHBT_X.C --batch --eval="(require 'cc-mode)" --eval="(c-set-offset 'substatement-open 0)" --eval="(c-set-offset 'arglist-intro 3)" --eval="(c-set-offset 'innamespace -2)" --eval="(setq-default indent-tabs-mode nil)" --eval='(indent-region (point-min) (point-max) nil)' -f save-buffer
	emacs N_DEV_ADMSHBT_X.h --batch --eval="(require 'cc-mode)" --eval="(c-set-offset 'substatement-open 0)" --eval="(c-set-offset 'arglist-intro 3)" --eval="(c-set-offset 'innamespace -2)" --eval="(setq-default indent-tabs-mode nil)" --eval='(indent-region (point-min) (point-max) nil)' -f save-buffer
	patch < limited_typed_source.patch
	patch < limited_typed_header.patch
