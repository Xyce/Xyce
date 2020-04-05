
.SUFFIXES:	.f

.f.o:
	g77 -c $(FFLAGS) $<

all:	vbic_3T_it_xf_solverv vbic_3T_it_xf_solveri

vbic_3T_it_xf_solverv: vbic_3T_it_xf_t.o vbic_3T_it_xf_elim_fj.o vbic_3T_it_xf_solverv.o
	g77 -o vbic_3T_it_xf_solverv vbic_3T_it_xf_t.o vbic_3T_it_xf_elim_fj.o vbic_3T_it_xf_solverv.o

vbic_3T_it_xf_solveri: vbic_3T_it_xf_t.o vbic_3T_it_xf_elim_fj.o vbic_3T_it_xf_solveri.o
	g77 -o vbic_3T_it_xf_solveri vbic_3T_it_xf_t.o vbic_3T_it_xf_elim_fj.o vbic_3T_it_xf_solveri.o

