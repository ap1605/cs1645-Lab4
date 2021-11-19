#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <sys/unistd.h>
#include <mpi.h>

#define XMAX 2.0
#define XMIN -2.0
#define YMAX 2.0
#define YMIN -2.0
#define N 10000                //number of divisions for the grid
#define ITER 50                //number of iterations for each point

#define TEST_RESULTS

#define 	NUM_THREADS 8


int pixels[N][N];

struct timeval startTime;
struct timeval finishTime;

//The printing is only for fun :)
void printMandelBrot();

void verifyMatrixSum();

double timeIntervalLength;

int main() {
    int i, j, count;

    double x, y;            //(x,y) point on the complex plane
    double x0, y0, tempx;
    double dx, dy;

    //Added
    int my_rank;    // processes id
    int comm_sz;    // total number of processes
    int local_n;    // each processes end position
    int local_i;    // each processes starting position
    int chunk;      // amount of work for each process
    double avg;     
    int totalSum;
    int localSum;

    //increments in the real and imaginary directions
    dx = (XMAX - XMIN) / N;
    dy = (YMAX - YMIN) / N;

    MPI_Init(NULL, NULL);

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

    chunk = N/comm_sz;
    local_i = chunk*my_rank;
    local_n = local_i+chunk;


    //Get the start time
    gettimeofday(&startTime, NULL); /* START TIME */
    //calculations for mandelbrot
    for (i = local_i+1; i <= local_n; i++) {
        for (j = 1; j <= N; j++) {
            // c_real
            x0 = XMAX -
                 (dx * (i));  //scaled i coordinate of pixel (scaled to lie in the Mandelbrot X scale (-2.00, 2.00))
            // c_imaginary
            y0 = YMAX -
                 (dy * (j));  //scaled j coordinate of pixel (scaled to lie in the Mandelbrot X scale (-2.00, 2.00))

            // z_real
            x = 0;

            // z_imaginary
            y = 0;
            count = 0;

            while ((x * x + y * y < 4) && (count < ITER)) {
                tempx = (x * x) - (y * y) + x0;
                y = (2 * x * y) + y0;
                x = tempx;
                count++;
            }
            pixels[i - 1][j - 1] = count;
        }
    }



    // Normalize the result based on the avg value
    localSum = 0.0;
    for (i = local_i; i < local_n; i++) {
        for (j = 0; j < N; j++) {
            localSum += (double) pixels[i][j];
        }
    }
    totalSum = 0.0;

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allreduce(&localSum, &totalSum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    avg = totalSum/(N*N);

    for (i = local_i; i < local_n; i++) {
        for (j = 0; j < N; j++) {
            pixels[i][j]=pixels[i][j]/avg;
        }
    }


    //Get the end time
    gettimeofday(&finishTime, NULL);  /* END TIME */

    //Calculate the interval length
    timeIntervalLength = (double) (finishTime.tv_sec - startTime.tv_sec) * 1000000
                         + (double) (finishTime.tv_usec - startTime.tv_usec);
    timeIntervalLength = timeIntervalLength / 1000;
    //Print the interval length
    printf("Interval length: %g msec.\n", timeIntervalLength);

    MPI_Finalize();


#ifdef TEST_RESULTS
    verifyMatrixSum(pixels);
//    printMandelBrot(pixels);
#endif
    return 0;
}

// Helper function to verify if the sum from parallel and serial versions match
void verifyMatrixSum() {
    int i, j;

    double totalSum;
    totalSum = 0;
    //
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            totalSum += (double) pixels[i][j];
        }
    }

    printf("\nTotal Sum = %g\n", totalSum);
}

void printMandelBrot() {
    int i, j;
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            printf("%d ", pixels[i][j]);
        }
        printf("\n");
    }
}
