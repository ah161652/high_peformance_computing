make
sbatch stencil.job
wait 10
python check.py --ref-stencil-file stencil_1024_1024_100.pgm --stencil-file stencil.pgm
cat stencil.out
