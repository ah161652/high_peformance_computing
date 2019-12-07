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
void stencil_mpi(const int nx, const int ny, const int width, const int height,
             float* image, float* tmp_image, int rank, int nprocs);
void halo(int rank, float* image, int height, float* buff, int start_pxl, int nx_mpi, int section_size, int nprocs);

int main(int argc, char* argv[])
{
  // Check usage
  if (argc != 4) {
    fprintf(stderr, "Usage: %s nx ny niters\n", argv[0]);
    exit(EXIT_FAILURE);
  }

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


  init_image(nx, ny, width, height, image, tmp_image);


  // if 1 processor then run sequentially
  if (nprocs == 1){

double tic = wtime();
 for (int t = 0; t < niters; ++t) {
   stencil(nx, ny, width, height, image, tmp_image);
   stencil(nx, ny, width, height, tmp_image, image);
 }
 double toc = wtime();

 // Output
 printf("------------------------------------\n");
 printf(" runtime: %lf s\n", toc - tic);
 printf("------------------------------------\n");

 output_image(OUTPUT_FILE, nx, ny, width, height, image);
 free(image);
 free(tmp_image);

 return 0;
  }


  // Split up columns
  int nx_mpi;
 if(rank == nprocs - 1){
   nx_mpi = (nx/nprocs) + (nx%nprocs);
 }
 else{
   nx_mpi = nx/nprocs;
 }

 printf("%s\n",nx_mpi);

 int section_size = height*nx_mpi;

 //define starting pixel
 int start_pxl = ((rank*nx_mpi)+1)*height;

 float* buff = malloc(height*sizeof(float));

  MPI_Barrier(MPI_COMM_WORLD);
   // Call the stencil kernel under mpi
  double tic = wtime();
  for (int t = 0; t < niters; ++t) {
    stencil_mpi(nx_mpi, ny, width, height, image, tmp_image, rank, nprocs);
    halo(rank, tmp_image, height, buff, start_pxl, nx_mpi, section_size, nprocs);
    stencil_mpi(nx_mpi, ny, width, height, tmp_image, image, rank , nprocs);
    halo(rank, image, height, buff, start_pxl, nx_mpi, section_size, nprocs);
  }
  double toc = wtime();


float* final_buff = malloc(sizeof(float)*section_size);


if (rank == 0){
  for(int i = 1; i < nprocs; ++i){
    MPI_Recv(final_buff, section_size, MPI_FLOAT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (int j = 0; j< section_size; ++j){
      image[(((nx_mpi*i)+1)*height) + j] = final_buff[j];
    }
  }
}

else{
  for (int i = 0; i < section_size; i++) {
    final_buff[i] = image[start_pxl + i];
  }
  MPI_Send(final_buff,section_size, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
}

// Output
printf("------------------------------------\n");
printf(" runtime: %lf s\n", toc - tic);
printf("------------------------------------\n");


if (rank == 0){
  output_image(OUTPUT_FILE, nx, ny, width, height, image);
}





  free(image);
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

void stencil_mpi(const int nx, const int ny, const int width, const int height,
             float* image, float* tmp_image, int rank, int nprocs)
{
  float calc = 3.0/5.0;
  float calc2 = 0.5/5.0;

  int start = 1 + (rank * nx);
  int end = start + nx;


  if (rank ==0){
  for (int i = 1; i < end; ++i) {
    for (int j = 1; j < ny + 1; ++j) {

      tmp_image[j + i * height] =  (image[j     + i       * height] * calc) + (image[j     + (i - 1) * height] * calc2) + (image[j     + (i + 1) * height] * calc2) + (image[j - 1 + i       * height] * calc2) + (image[j + 1 + i       * height] * calc2) ;

    }
  }
}

else if (rank == nprocs -1){

  for (int i = start; i < end; ++i) {
    for (int j = 1; j <  ny + 1; ++j) {

      tmp_image[j + i * height] =  (image[j     + i       * height] * calc) + (image[j     + (i - 1) * height] * calc2) + (image[j     + (i + 1) * height] * calc2) + (image[j - 1 + i       * height] * calc2) + (image[j + 1 + i       * height] * calc2) ;

    }
  }

}

else{

  for (int i = start; i < end; ++i) {
    for (int j = 1; j < ny + 1; ++j) {

      tmp_image[j + i * height] =  (image[j     + i       * height] * calc) + (image[j     + (i - 1) * height] * calc2) + (image[j     + (i + 1) * height] * calc2) + (image[j - 1 + i       * height] * calc2) + (image[j + 1 + i       * height] * calc2) ;

    }
  }

}


}

void halo(int rank, float* image, int height, float* buff, int start_pxl, int nx_mpi, int section_size, int nprocs)
{
if(rank == 0){
  MPI_Sendrecv(&image[start_pxl + (nx_mpi-1)*height], height,  MPI_FLOAT, rank + 1, 0,
               buff, height, MPI_FLOAT, rank+1, 0,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                for (int i = 0; i < height; ++i) {
                  image[start_pxl + section_size + i] = buff[i];
                }

}

else if ( rank == nprocs - 1){
  MPI_Sendrecv(&image[start_pxl], height,  MPI_FLOAT, rank - 1, 0,
               buff, height, MPI_FLOAT, rank-1, 0,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);

               for (int i = 0; i < height; ++i) {
                 image[start_pxl - height + i] = buff[i];
}
}

else {
  MPI_Sendrecv(&image[start_pxl+ ((nx_mpi - 1)*height)], height,  MPI_FLOAT, rank + 1, 0,
               buff, height, MPI_FLOAT, rank-1, 0,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);

               for (int i = 0; i < height; ++i) {
                 image[start_pxl - height + i] = buff[i];
               }

 MPI_Sendrecv(&image[start_pxl], height,  MPI_FLOAT, rank - 1, 0,
              buff, height, MPI_FLOAT, rank+1, 0,
              MPI_COMM_WORLD, MPI_STATUS_IGNORE);

              for (int i = 0; i < height; ++i) {
                image[start_pxl + section_size + i] = buff[i];

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
