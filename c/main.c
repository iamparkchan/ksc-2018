#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include <mpi.h>
#include "gram_schmidt.h"

static void initialize(double *vector[]);
static double orthogonality_test(double *basis[]);
static double span_test(double *vector[], double *basis[]);

int main(int argc, char *argv[])
{
    int i, j, rank;
    double *vector[N], *basis[N];
    double time_interval = 0.;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    if (0 == rank) {
        for (i = 0; i < N; i++) {
            vector[i] = malloc(sizeof(double) * M);
            basis [i] = malloc(sizeof(double) * M);
        }
        initialize(vector);

        for (i = 0; i < N; i++) {
            for (j = 0; j < M; j++) {
                basis[i][j] = vector[i][j];
            }
        }
        
        time_interval = -MPI_Wtime();
    }

    gram_schmidt(basis);
    
    if (0 == rank) {
        time_interval += MPI_Wtime();
        
        printf("Elapsed Time: %f sec\n", time_interval);
        
        printf("Orthogonality Test: %e\n", orthogonality_test(basis));
        printf("Span Test: %e\n", span_test(vector, basis));
        
        for (i = 0; i < N; i++) {
            free(vector[i]);
            free(basis[i]);
        }
    }
    
    MPI_Finalize();
    
    return 0;
}

static void initialize(double *vector[])
{
    int i, k;
    
    srand((unsigned)time(NULL));
    for (i = 0; i < N; i++) {
        for (k = 0; k < M; k++) {
            vector[i][k] = 2. * rand() / RAND_MAX - 1.;
        }
    }
}

static double orthogonality_test(double *basis[])
{
    double max = 0., sum = 0.;
    int i, j, k;

    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            sum = 0.;
            for (k = 0; k < M; k++) {
                sum += basis[i][k]*basis[j][k];
            }
            sum += (i == j) ? -1. : 0.;
            sum = fabs(sum);
            if (sum > max) max = sum;
        }
    }
    
    return max;
}

static double span_test(double *vector[], double *basis[])
{
    double max = 0., sum = 0.;
    double *v;
    int i, j, k;
    
    v = malloc(sizeof(double) * M);
    for (i = 0; i < N; i++) {
        for (k = 0; k < M; k++) {
            v[k] = vector[i][k];
        }
        for (j = 0; j < N; j++) {
            sum = 0.;
            for (k = 0; k < M; k++) {
                sum += v[k] * basis[j][k];
            }
            for (k = 0; k < M; k++) {
                v[k] -= sum * basis[j][k];
            }
        }
        for (k = 0; k < M; k++) {
            if (max < fabs(v[k])) max = fabs(v[k]);
        }
    }
    free(v);
    
    return max;
}