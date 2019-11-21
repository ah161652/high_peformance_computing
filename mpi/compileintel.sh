rm stencil
rm stencil.out
rm stencil.pgm
mpiicc -std=c99 stencil.c -Ofast -xHOST -qopt-report=5 -qopt-report-phase=vec -o stencil
sbatch stencil.job
sleep 20
python check.py --ref-stencil-file stencil_1024_1024_100.pgm --stencil-file stencil.pgm
cat stencil.out
