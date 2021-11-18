#include <stdio.h>
#include <time.h>
#include <sys/time.h>

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

    //Get the start time
    gettimeofday(&startTime, NULL); /* START TIME */

    //A*B*C
    for(i=0;i<NROW;i++){
        for(j=0; j<NCOL; j++){
            for(k=0;k<NROW;k++){
                // K iterates rows on matrixA and columns in matrixB
                tempMatrix[i][j] += matrixA[i][k] * matrixB[k][j];
            }
        }

        for(j=0; j<NCOL; j++) {
            for (k = 0; k < NROW; k++) {
                // K iterates rows on matrixA and columns in matrixB
                outputMatrix[i][j] += tempMatrix[i][k] * matrixC[k][j];
            }
        }
    }

    //Get the end time
    gettimeofday(&finishTime, NULL);  /* END TIME */

    //Calculate the interval length
    timeIntervalLength = (double)(finishTime.tv_sec-startTime.tv_sec) * 1000000
                         + (double)(finishTime.tv_usec-startTime.tv_usec);
    timeIntervalLength=timeIntervalLength/1000;

    #ifdef TEST_RESULTS
    //[Verifying the matrix summation]
    verifyMatrixSum();
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
