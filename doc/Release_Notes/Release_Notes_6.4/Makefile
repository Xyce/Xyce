# Makefile for the Xyce release notes.
#

# Set the TEXINPUTS environment variable (not MAKE variable!)
# The final colon on this variable is mandatory.  It indicates that the given paths are
# PREPENDED to the defaults, and do not replace the defaults.
export TEXINPUTS:=../../Common_Guide_Files/:../../Common_Guide_Files/SANDreport:


BASE=Release_Notes_6.4
all: $(BASE).pdf

$(BASE).pdf: $(BASE).tex knownDefectsTbl.tex fixedDefectsTbl.tex changesInputTbl.tex 

	pdflatex $(BASE)
	pdflatex $(BASE)

clean:
	@rm -f $(BASE).aux $(BASE).toc $(BASE).idx $(BASE).ind $(BASE).ilg $(BASE).log $(BASE).out \
           $(BASE).bbl $(BASE).blg $(BASE).lof $(BASE).lot

realclean:	clean
	@rm -f $(BASE).ps $(BASE).pdf $(BASE).dvi
