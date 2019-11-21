CC=mpicc

stencil: stencil.c
	mpiicc stencil.c -Ofast -fast -xHOST -qopt-report=5 -qopt-report-phase=vec -o stencil
