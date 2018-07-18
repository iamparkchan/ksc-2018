#include "gram_schmidt.h"
#include <mpi.h>
#include <math.h>
#include <stdlib.h>

static double dot_product(double *A, double *B, int n)
{
    int k;
    double sum = 0.;

    #pragma omp parallel for reduction(+:sum)
    for (k = 0; k < n; k++) {
        sum += A[k] * B[k];
    }

    MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return sum;
}


static void multiply_add(double *A, double c, double *B, int n)
{
    int k;

    #pragma omp parallel for
    for (k = 0; k < n; k++) {
        A[k] += c * B[k];
    }
}

static void multiply(double *A, double c, int n)
{
    int k;

    #pragma omp parallel for
    for (k = 0; k < n; k++) {
        A[k] *= c; 
    }
}

void gram_schmidt(double *vector[])
{
    int i, j, r;
    double coef;
    double *v[N];
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int counts[size], displs[size];

    int unit = (M + size - 1) / size;
    for (r = 0; r < size; r++) {
        if (unit*(r + 1) <= M) {
            counts[r] = unit;
        } else if (unit*r < M) {
            counts[r] = M - unit*r;
        } else {
            counts[r] = 0;
        }
        displs[r] = unit*r;
    }
    int chunk = counts[rank];

    if (chunk > 0) {
        for (i = 0; i < N; i++) {
            v[i] = malloc(sizeof(double) * chunk);
        }
    }

    for (i = 0; i < N; i++) {
        MPI_Scatterv(vector[i], counts, displs, MPI_DOUBLE, v[i], chunk, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    if (chunk > 0) {
        for (j = 0; j < N; j++) {
            for (i = 0; i < j; i++) {
                coef = -dot_product(v[i], v[j], chunk);
                multiply_add(v[j], coef, v[i], chunk);
            }
            coef = 1. / sqrt(dot_product(v[j], v[j], chunk));
            multiply(v[j], coef, chunk);
        }
    }

    for (i = 0; i < N; i++) {
        MPI_Gatherv(v[i], chunk, MPI_DOUBLE, vector[i], counts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    if (chunk > 0) {
        for (i = 0; i < N; i++) {
            free(v[i]);
        }
    }
}