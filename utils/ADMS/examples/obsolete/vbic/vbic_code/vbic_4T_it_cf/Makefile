
.SUFFIXES:	.f

.f.o:
	g77 -c $(FFLAGS) $<

all:	vbic_4T_it_cf_solverv vbic_4T_it_cf_solveri

vbic_4T_it_cf_solverv: vbic_4T_it_cf_t.o vbic_4T_it_cf_elim_fj.o vbic_4T_it_cf_solverv.o
	g77 -o vbic_4T_it_cf_solverv vbic_4T_it_cf_t.o vbic_4T_it_cf_elim_fj.o vbic_4T_it_cf_solverv.o

vbic_4T_it_cf_solveri: vbic_4T_it_cf_t.o vbic_4T_it_cf_elim_fj.o vbic_4T_it_cf_solveri.o
	g77 -o vbic_4T_it_cf_solveri vbic_4T_it_cf_t.o vbic_4T_it_cf_elim_fj.o vbic_4T_it_cf_solveri.o

