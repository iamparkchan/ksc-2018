#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include "gram_schmidt.h"

static void initialize(double **vector);
static void export(double **vector);

int main(int argc, char *argv[])
{
    int i, rank;
    double *vector[N];
    double time_interval = 0.;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    if (0 == rank) {
        for (i = 0; i < N; i++) {
            vector[i] = malloc(sizeof(double) * M);
        }
        
        initialize(vector);
        
        time_interval = -MPI_Wtime();
    }

    gram_schmidt(vector);
    
    if (0 == rank) {
        time_interval += MPI_Wtime();
        
        printf("elapsed time: %f sec\n", time_interval);
        
        export(vector);
        
        for (i = 0; i < N; i++) {
            free(vector[i]);
        }
    }
    
    MPI_Finalize();
    
    return 0;
}

static void initialize(double **vector)
{
    int i, k;
    
    srand((unsigned)time(NULL));
    for (i = 0; i < N; i++) {
        for (k = 0; k < M; k++) {
            vector[i][k] = 2. * rand() / RAND_MAX - 1.;
        }
    }
}

static void export(double **vector)
{
    int i, k;
    FILE *fp;
    
    fp = fopen("result", "w");
    for (i = 0; i < N; i++) {
        fwrite(vector[i], sizeof(double), M, fp);
    }
    fclose(fp);
}