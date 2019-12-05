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
void halo(int rank, int nprocs, float* image, int startcol, int endcol, int working_cols);

int main(int argc, char* argv[])
{
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

  // divide up work into columns
  int working_cols = nx/(nprocs);

  //find starting col and ending col dependant on rank
  int endcol;
  int startcol;

  startcol = (rank * working_cols) + 1;

  if (rank == nprocs - 1) {
        endcol = nx + 1;
    }
  else {
       endcol = startcol + working_cols ;
    }

  // Allocate the image
  float* image = malloc(sizeof(float) * width * height);
  float* tmp_image = malloc(sizeof(float) * width * height);

  // Set the input image
  init_image(nx, ny, width, height, image, tmp_image);

  // Call the stencil kernel
  double tic = wtime();
  for (int t = 0; t < niters; ++t) {
    //halo(rank, nprocs, image, startcol, endcol, working_cols);
    stencil(nx, ny, width, height, image, tmp_image);
    // halo(rank, nprocs, tmp_image, startcol, endcol, working_cols);
    stencil(nx, ny, width, height, image, tmp_image);
  }
  double toc = wtime();



  // stitch it back together
  float* final_image = malloc(sizeof(float) * width * height);

  final_image  = image;


  // Output
  printf("------------------------------------\n");
  printf(" runtime: %lf s\n", toc - tic);
  printf("------------------------------------\n");

  output_image(OUTPUT_FILE, nx, ny, width, height, image);
  free(image);
  free(final_image);
  free(tmp_image);

  MPI_Finalize();
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


void halo(int rank, int nprocs, float* image, int startcol, int endcol, int working_cols){
    if (rank == 0){
      // send right, receive right
      //store in array
      //input into image



    }
    else if (rank == nprocs - 1){
      //send left receive left
      //store in array
      //input into image
    }
    else{
      //send and receive both ways
      //store in arrays
      //input into image
    }

    return;




}
