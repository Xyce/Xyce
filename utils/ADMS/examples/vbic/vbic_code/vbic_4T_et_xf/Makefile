
.SUFFIXES:	.f

.f.o:
	g77 -c $(FFLAGS) $<

all:	vbic_4T_et_xf_solverv vbic_4T_et_xf_solveri

vbic_4T_et_xf_solverv: vbic_4T_et_xf_t.o vbic_4T_et_xf_elim_fj.o vbic_4T_et_xf_solverv.o
	g77 -o vbic_4T_et_xf_solverv vbic_4T_et_xf_t.o vbic_4T_et_xf_elim_fj.o vbic_4T_et_xf_solverv.o

vbic_4T_et_xf_solveri: vbic_4T_et_xf_t.o vbic_4T_et_xf_elim_fj.o vbic_4T_et_xf_solveri.o
	g77 -o vbic_4T_et_xf_solveri vbic_4T_et_xf_t.o vbic_4T_et_xf_elim_fj.o vbic_4T_et_xf_solveri.o

