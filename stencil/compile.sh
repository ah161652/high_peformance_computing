rm stencil
rm stencil.out
rm stencil.pgm
gcc -std=c99 -Wall stencil.c -Ofast -march=native -o stencil
sbatch stencil.job
sleep 10
python check.py --ref-stencil-file stencil_1024_1024_100.pgm --stencil-file stencil.pgm
cat stencil.out
