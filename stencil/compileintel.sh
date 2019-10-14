rm stencil
rm stencil.out
rm stencil.pgm
icc stencil.c -Ofast -march=native -o stencil
sbatch stencil.job
sleep 10
python check.py --ref-stencil-file stencil_1024_1024_100.pgm --stencil-file stencil.pgm
cat stencil.out
