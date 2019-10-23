CC=icc

stencil: stencil.c
	icc stencil.c -Ofast -fast -march=native -qopt-report=1 -qopt-report-phase=vec -o stencil
