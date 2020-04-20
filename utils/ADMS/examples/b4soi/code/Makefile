ADMSDIR = ../../..

XyceADMSFILES = $(ADMSDIR)/adms.implicit.xml \
	$(ADMSDIR)/xyceVersion_nosac.xml \
	$(ADMSDIR)/xyceBasicTemplates_nosac.xml \
	$(ADMSDIR)/xyceAnalogFunction_nosac.xml \
	$(ADMSDIR)/xyceHeaderFile_nosac.xml \
	$(ADMSDIR)/xyceImplementationFile_nosac.xml \
	$(ADMSDIR)/xyceOutVarsDoc.xml

CLEANFILES = N_DEV_ADMSbsimsoi.h N_DEV_ADMSbsimsoi.C disciplines.vams constants.vams *.tex

all-source: N_DEV_ADMSbsimsoi.C

clean:
	rm -f $(CLEANFILES)

N_DEV_ADMSbsim6.h N_DEV_ADMSbsimsoi.C: bsimsoi.va $(XyceADMSFILES)
	admsXml -x -e $(ADMSDIR)/adms.implicit.xml -e $(ADMSDIR)/xyceVersion_nosac.xml -e $(ADMSDIR)/xyceBasicTemplates_nosac.xml -e $(ADMSDIR)/xyceAnalogFunction_nosac.xml -e $(ADMSDIR)/xyceImplementationFile_nosac.xml -e $(ADMSDIR)/xyceHeaderFile_nosac.xml -e $(ADMSDIR)/xyceOutVarsDoc.xml bsimsoi.va
	mv M_70_OutputVars.tex M_70_OutputVars.tex_borken
	sed -e 's|T\_|T\\_|g' < M_70_OutputVars.tex_borken > M_70_OutputVars.tex
	rm M_70_OutputVars.tex_borken
