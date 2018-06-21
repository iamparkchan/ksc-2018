#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include "gram_schmidt.h"

static void initialize(double **vector);
static double dot_product(double *A, double *B);
static double rms_error(double **vector);

int main(int argc, char *argv[])
{
    int i, rank;
    double *vector[M];
    double time_interval = 0.;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    if (0 == rank) {
        for (i = 0; i < M; i++) {
            vector[i] = malloc(sizeof(double) * N);
        }
        
        initialize(vector);
        
        time_interval = -MPI_Wtime();
    }
    
    gram_schmidt(vector);
    
    if (0 == rank) {
        time_interval += MPI_Wtime();
        
        printf("error: %e\n", rms_error(vector));
        printf("elapsed time: %f sec\n", time_interval);
        
        for (i = 0; i < M; i++) {
            free(vector[i]);
        }
    }
    
    MPI_Finalize();
    
    return 0;
}

void initialize(double **vector)
{
    int i, k;
    
    srand((unsigned)time(NULL));
    for (i = 0; i < M; i++) {
        for (k = 0; k < N; k++) {
            vector[i][k] = 2. * rand() / RAND_MAX - 1.;
        }
    }
}

double dot_product(double *A, double *B)
{
    int k;
    double sum = 0.;
    
    for (k = 0; k < N; k++) {
        sum += A[k] * B[k];
    }
    
    return sum;
}

double rms_error(double **vector)
{
    int i, j;
    double sum = 0., err;
    
    for (i = 0; i < M; i++) {
        for (j = 0; j < M; j++) {
            err = dot_product(vector[i], vector[j]);
            if (i == j) err -= 1.;
            err /= N;
            sum += err * err;
        }
    }
    
    return sqrt(sum/M/M);
}
