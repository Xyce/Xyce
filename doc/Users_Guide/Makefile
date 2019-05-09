# From the latexmk Users' Guide (or 'man latexmk')

export TEXINPUTS=:../Common_Guide_Files//:
export BIBINPUTS=:../Common_Guide_Files/

.PHONY : FORCE_MAKE
all : Xyce_UG.pdf
%.pdf : %.tex FORCE_MAKE
	latexmk -pdf $<

clean:
	latexmk -c -bibtex

realclean:
	latexmk -C -bibtex

