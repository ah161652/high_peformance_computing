CC=mpicc

stencil: stencil.c
	mpiicc -std=c99 stencil.c -Ofast -xHOST -qopt-report=5 -qopt-report-phase=vec -o stencil
