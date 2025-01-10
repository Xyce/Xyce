ADMSDIR = ../../..

XyceADMSDeps = $(ADMSDIR)/adms.implicit.xml \
	$(ADMSDIR)/xyceVersion_nosac.xml \
	$(ADMSDIR)/xyceBasicTemplates_nosac.xml \
	$(ADMSDIR)/xyceAnalogFunction_nosac.xml \
	$(ADMSDIR)/xyceHeaderFile_nosac.xml \
	$(ADMSDIR)/xyceImplementationFile_nosac.xml
all-source: N_DEV_ADMSekv_va.C

clean:
	rm -f N_DEV_ADMSekv_va.C N_DEV_ADMSekv_va.h *.o *.so* *.orig *.rej

N_DEV_ADMSekv_va.C: ekv26_SDext_Verilog-A.va $(XyceADMSDeps)
	admsXml -x -e $(ADMSDIR)/adms.implicit.xml -e $(ADMSDIR)/xyceVersion_nosac.xml -e $(ADMSDIR)/xyceBasicTemplates_nosac.xml -e $(ADMSDIR)/xyceAnalogFunction_nosac.xml -e $(ADMSDIR)/xyceImplementationFile_nosac.xml  -e $(ADMSDIR)/xyceHeaderFile_nosac.xml ekv26_SDext_Verilog-A.va
