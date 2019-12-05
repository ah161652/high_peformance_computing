#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <mpi.h>
// Define output file name
#define OUTPUT_FILE "stencil.pgm"

void stencil(const int nx, const int ny, const int width, const int height,
             float* image, float* tmp_image);
void init_image(const int nx, const int ny, const int width, const int height,
                float* image, float* tmp_image);
void output_image(const char* file_name, const int nx, const int ny,
                  const int width, const int height, float* image);
double wtime(void);

int main(int argc, char* argv[])
{
  // Check usage
  if (argc != 4) {
    fprintf(stderr, "Usage: %s nx ny niters\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  // Initiliase problem dimensions from command line arguments
  int nx = atoi(argv[1]);
  int ny = atoi(argv[2]);
  int niters = atoi(argv[3]);

  // we pad the outer edge of the image to avoid out of range address issues in
  // stencil
  int width = nx + 2;
  int height = ny + 2;

  // Allocate the image
  float* image = malloc(sizeof(float) * width * height);
  float* tmp_image = malloc(sizeof(float) * width * height);

  //MPI setup
  MPI_Init(&argc, &argv);
  int nprocs, rank, flag;

  //Check if init worked
  MPI_Initialized(&flag);
  if ( flag != 1 ) {
    MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
  }

  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(rank ==0){
  // Set the input image
  init_image(nx, ny, width, height, image, tmp_image);
}

  int num_cols;
  int num_cols_with_halo;

  // if (nx % nprocs == 0 ){
  //   num_cols = width/(nprocs);
  // }
  // else{
  //   if (rank != nprocs -1){
  //     num_cols = width/(nprocs);
  //   }
  //   else{
  //     num_cols = (width % nprocs) + (width/nprocs);
  //   }
  // }
  //
  // if (rank == 0){
  //   int num_cols_with_halo = num_cols + 1;
  // }
  // else if (rank = nprocs - 1){
  //   num_cols_with_halo - num_cols + 1;
  // }
  // else{
  //   num_cols_with_halo = num_cols +2;
  // }
  //
  // int working_size = num_cols * height;
  // int working_size_with_halo = num_cols_with_halo * height;

  int working_size = height * width/nprocs;
  int remainder_size = (width%nprocs)*height;


  float* buffer;
  float* tmp_buffer;
  if (rank == nprocs -1 ){
    buffer = (float*)malloc(sizeof(float) * (working_size + remainder_size));
    tmp_buffer = (float*)malloc(sizeof(float) * (working_size + remainder_size));
  }
  else {
  buffer = (float*)malloc(sizeof(float) * working_size);
  tmp_buffer = (float*)malloc(sizeof(float) * working_size);
}


  MPI_Scatter(image, working_size_with_halo, MPI_FLOAT, buffer, working_size_with_halo, MPI_FLOAT, 0, MPI_COMM_WORLD);

  int num_cols = (rank == nprocs-1) ? nx/nprocs + nx%nprocs : nx/nprocs;






  // Call the stencil kernel
  double tic = wtime();
  for (int t = 0; t < niters; ++t) {
    stencil(num_cols, ny, width, height, buffer, tmp_buffer);
    stencil(num_cols, ny, width, height, buffer, tmp_buffer);
  }
  double toc = wtime();




float* final_image = malloc(sizeof(float)*width*height);

MPI_Gather(tmp_buffer, working_size_with_halo, MPI_FLOAT,final_image ,working_size_with_halo, MPI_FLOAT,0, MPI_COMM_WORLD);





  if(rank == 0){

  output_image(OUTPUT_FILE, nx, ny, width, height, image);
}
MPI_Finalize();

// Output
printf("------------------------------------\n");
printf(" runtime: %lf s\n", toc - tic);
printf("------------------------------------\n");


  free(image);
  free(tmp_image);

}

void stencil(const int nx, const int ny, const int width, const int height,
             float* image, float* tmp_image)
{
  float calc = 3.0/5.0;
  float calc2 = 0.5/5.0;
  for (int i = 1; i < ny + 1; ++i) {
    for (int j = 1; j < nx + 1; ++j) {

      tmp_image[j + i * height] =  (image[j     + i       * height] * calc) + (image[j     + (i - 1) * height] * calc2) + (image[j     + (i + 1) * height] * calc2) + (image[j - 1 + i       * height] * calc2) + (image[j + 1 + i       * height] * calc2) ;
      // tmp_image[j + i * height] += (image[j     + (i - 1) * height] * calc2) + (image[j     + (i + 1) * height] * calc2) + (image[j - 1 + i       * height] * calc2) + (image[j + 1 + i       * height] * calc2);
      // tmp_image[j + i * height] += image[j     + (i + 1) * height] * calc2;
      // tmp_image[j + i * height] += image[j - 1 + i       * height] * calc2;
      // tmp_image[j + i * height] += image[j + 1 + i       * height] * calc2;
    }
  }
}



// Create the input image
void init_image(const int nx, const int ny, const int width, const int height,
                float* image, float* tmp_image)
{
  // Zero everything
  for (int j = 0; j < ny + 2; ++j) {
    for (int i = 0; i < nx + 2; ++i) {
      image[j + i * height] = 0.0;
      tmp_image[j + i * height] = 0.0;
    }
  }

  const int tile_size = 64;
  // checkerboard pattern
  for (int jb = 0; jb < ny; jb += tile_size) {
    for (int ib = 0; ib < nx; ib += tile_size) {
      if ((ib + jb) % (tile_size * 2)) {
        const int jlim = (jb + tile_size > ny) ? ny : jb + tile_size;
        const int ilim = (ib + tile_size > nx) ? nx : ib + tile_size;
        for (int j = jb + 1; j < jlim + 1; ++j) {
          for (int i = ib + 1; i < ilim + 1; ++i) {
            image[j + i * height] = 100.0;
          }
        }
      }
    }
  }
}

// Routine to output the image in Netpbm grayscale binary image format
void output_image(const char* file_name, const int nx, const int ny,
                  const int width, const int height, float* image)
{
  // Open output file
  FILE* fp = fopen(file_name, "w");
  if (!fp) {
    fprintf(stderr, "Error: Could not open %s\n", OUTPUT_FILE);
    exit(EXIT_FAILURE);
  }

  // Ouptut image header
  fprintf(fp, "P5 %d %d 255\n", nx, ny);

  // Calculate maximum value of image
  // This is used to rescale the values
  // to a range of 0-255 for output
  float maximum = 0.0;
  for (int j = 1; j < ny + 1; ++j) {
    for (int i = 1; i < nx + 1; ++i) {
      if (image[j + i * height] > maximum) maximum = image[j + i * height];
    }
  }

  // Output image, converting to numbers 0-255
  for (int j = 1; j < ny + 1; ++j) {
    for (int i = 1; i < nx + 1; ++i) {
      fputc((char)(255.0 * image[j + i * height] / maximum), fp);
    }
  }

  // Close the file
  fclose(fp);
}

// Get the current time in seconds since the Epoch.
double wtime(void)
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + tv.tv_usec * 1e-6;
}
