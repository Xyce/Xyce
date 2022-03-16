# NOTE:  This is stupidly assuming your git repos are all cloned
# in the same way that mine are, and is not flexible.
ADMSDIR = ../../../Xyce/utils/ADMS

XyceADMSDeps = $(ADMSDIR)/adms.implicit.xml \
	$(ADMSDIR)/xyceVersion_nosac.xml \
	$(ADMSDIR)/xyceBasicTemplates_nosac.xml \
	$(ADMSDIR)/xyceAnalogFunction_nosac.xml \
	$(ADMSDIR)/xyceHeaderFile_nosac.xml \
	$(ADMSDIR)/xyceImplementationFile_nosac.xml \
	$(ADMSDIR)/xyceOutVarsDoc.xml
all-source: N_DEV_ADMSl_utsoi.C N_DEV_ADMSl_utsoi.h 

clean:
	rm -f N_DEV_ADMSl_utsoi.C N_DEV_ADMSl_utsoi.h *.o *.so* *.orig *.rej *.tex

N_DEV_ADMSl_utsoi.C: L_UTSOI_102.va $(XyceADMSDeps)
	admsXml -x -e $(ADMSDIR)/adms.implicit.xml -e $(ADMSDIR)/xyceVersion_nosac.xml -e $(ADMSDIR)/xyceBasicTemplates_nosac.xml -e $(ADMSDIR)/xyceAnalogFunction_nosac.xml -e $(ADMSDIR)/xyceImplementationFile_nosac.xml  -e $(ADMSDIR)/xyceHeaderFile_nosac.xml -e $(ADMSDIR)/xyceOutVarsDoc.xml L_UTSOI_102.va
	mv M_10240_OutputVars.tex M_10240_OutputVars.tex_borken
	sed -e 's|\([thf]\)_|\1\\_|g' < M_10240_OutputVars.tex_borken > M_10240_OutputVars.tex
	rm M_10240_OutputVars.tex_borken
