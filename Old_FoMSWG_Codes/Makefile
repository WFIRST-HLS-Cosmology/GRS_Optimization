clean: *.c *.h Makefile
	-rm wl-fisher.x
	-rm bao-fisher.x
test: test.c fishercalc.c ps.c utils.c nr_utils.c fom.c
	gcc test.c ps.c utils.c nr_utils.c fom.c -lm -Wall -O3\
          -DSWG36 -DGROWTH_GAMMA -DGAUSS_CONVERGENCE -DPLANCK_PAR\
          -DWL_VERSION -o wl-fisher.x
	gcc test.c ps.c utils.c nr_utils.c fom.c -lm -Wall -O3\
          -DSWG36 -DGROWTH_GAMMA -DGAUSS_CONVERGENCE -DPLANCK_PAR\
          -DBAO_VERSION -o bao-fisher.x
