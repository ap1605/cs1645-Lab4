#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <mpi.h>

#define		NROW	1024
#define		NCOL	NROW

#define 	NUM_THREADS 8

#define TEST_RESULTS


//Matrix A
int matrixA  [NROW][NCOL];
//Matrix B
int matrixB  [NROW][NCOL];
//Matrix C
int matrixC [NROW][NCOL];

//Temp array
int tempMatrix [NROW][NCOL];

//Output Array C
int outputMatrix [NROW][NCOL];

struct timeval startTime;
struct timeval finishTime;
double timeIntervalLength;

void verifyMatrixSum();

int main(int argc, char* argv[])
{
    int i,j,k;
    // Matrix initialization. Just filling with arbitrary numbers.
    for(i=0;i<NROW;i++)
    {
        for(j=0;j<NCOL;j++)
        {
            matrixA[i][j]= (i + j)/128;
            matrixB[i][j]= (j + j)/128;
            matrixC[i][j]= (i + j)/128;
            tempMatrix[i][j] = 0;
            outputMatrix[i][j]= 0;
        }
    }


    //Added
    int my_rank;    // processes id
    int comm_sz;    // total number of processes
    int chunk;      // amount of work for each process



    //MPI_Status sts;

    MPI_Init(NULL, NULL);

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

    chunk = NROW/comm_sz;


    int local_i = chunk*my_rank;
    //int local_n = chunk*(my_rank+1);

    int local_temp[chunk][NROW];
    int local_out[chunk][NROW];
    int local_a[chunk][NROW];

    //Get the start time
    gettimeofday(&startTime, NULL); /* START TIME */

    MPI_Scatter(&matrixA, chunk*NROW, MPI_INT, &local_a, chunk*NROW, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&matrixB, NROW*NCOL, MPI_INT, 0, MPI_COMM_WORLD);

    //A*B*C
    for(i=0;i<chunk;i++){
        for(j=0; j<NCOL; j++){
            for(k=0;k<NROW;k++){
                // K iterates rows on matrixA and columns in matrixB
                local_temp[i][j] += local_a[i][k] * matrixB[k][j];
            }
        }
        for(j=0; j<NCOL; j++) {
            for (k = 0; k < NROW; k++) {
                // K iterates rows on matrixA and columns in matrixB
                local_out[i][j] += local_temp[i][k] * matrixC[k][j];
            }
        }
    }
    if(my_rank == 0) printf("\n MAIN THREAD \n");


    MPI_Gather(&local_out, chunk*NROW, MPI_INT, &outputMatrix[local_i][0], chunk*NROW, MPI_INT, 0, MPI_COMM_WORLD);


    MPI_Finalize();

    



    //Get the end time
    gettimeofday(&finishTime, NULL);  /* END TIME */



    //Calculate the interval length
    timeIntervalLength = (double)(finishTime.tv_sec-startTime.tv_sec) * 1000000
                         + (double)(finishTime.tv_usec-startTime.tv_usec);
    timeIntervalLength=timeIntervalLength/1000;

    #ifdef TEST_RESULTS
    //[Verifying the matrix summation]
    if(my_rank == 0) verifyMatrixSum();
    #endif

    //Print the interval length
    printf("Interval length: %g msec.\n", timeIntervalLength);

    return 0;
}


// Helper function to verify if the sum from parallel and serial versions match
void verifyMatrixSum() {
    int i, j;

    double totalSum;
    totalSum=0;
    //
    for(i=0;i<NROW;i++){
        for(j=0;j<NCOL;j++)
        {
            totalSum+=(double)outputMatrix[i][j];
        }
    }
    printf("\nTotal Sum = %g\n",totalSum);
}
