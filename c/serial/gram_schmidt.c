#include "gram_schmidt.h"
#include <mpi.h>
#include <math.h>

static double dot_product(double *A, double *B)
{
    int k;
    double sum = 0.;
    
    for (k = 0; k < M; k++) {
        sum += A[k] * B[k];
    }
    
    return sum;
}


static void multiply_add(double *A, double c, double *B)
{
    int k;

    for (k = 0; k < M; k++) {
        A[k] += c * B[k];
    }
}

static void multiply(double *A, double c)
{
    int k;
    
    for (k = 0; k < M; k++) {
        A[k] *= c; 
    }
}

void gram_schmidt(double *vector[])
{
    int rank;
    int i, j;
    double coef;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    if (0 == rank) {
        for (j = 0; j < N; j++) {
            for (i = 0; i < j; i++) {
                coef = -dot_product(vector[i], vector[j]);
                multiply_add(vector[j], coef, vector[i]);
            }
            coef = 1. / sqrt(dot_product(vector[j], vector[j]));
            multiply(vector[j], coef);
        }
    }
}
