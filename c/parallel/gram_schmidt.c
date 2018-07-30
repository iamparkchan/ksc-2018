#include "gram_schmidt.h"
#include <mpi.h>
#include <math.h>
#include <stdlib.h>

static double dot_product(double *A, double *B, int n, MPI_Comm comm)
{
    int k;
    double sum = 0.;

    #pragma omp parallel for reduction(+:sum)
    for (k = 0; k < n; k++) {
        sum += A[k] * B[k];
    }

    MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, comm);

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

#define CEIL_DIV(x, y) (((x) + (y) - 1) / (y))

void gram_schmidt(double *vector[])
{
    int i, j, r;
    double coef;
    double *v[N];
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    int unit = CEIL_DIV(M, size);
    int nprocs = CEIL_DIV(M, unit);

    MPI_Group world_group;
    MPI_Comm_group(MPI_COMM_WORLD, &world_group);
    MPI_Group work_group;
    int range[3] = {0, nprocs - 1, 1};
    MPI_Group_range_incl(world_group, 1, &range, &work_group);
    MPI_Comm comm;
    MPI_Comm_create(MPI_COMM_WORLD, work_group, &comm);
    MPI_Group_free(&world_group);
    MPI_Group_free(&work_group);
    
    if (rank >= nprocs) return;
    MPI_Comm_rank(comm, &rank);

    int counts[nprocs], displs[nprocs];

    for (r = 0; r < nprocs; r++) {
        if (r != nprocs - 1) {
            counts[r] = unit;
        } else {
            counts[r] = M - unit*r;
        }
        displs[r] = unit*r;
    }
    
    int chunk = counts[rank];

    for (i = 0; i < N; i++) {
        v[i] = malloc(sizeof(double) * chunk);
    }

    for (i = 0; i < N; i++) {
        MPI_Scatterv(vector[i], counts, displs, MPI_DOUBLE, v[i], chunk, MPI_DOUBLE, 0, comm);
    }

    for (j = 0; j < N; j++) {
        for (i = 0; i < j; i++) {
            coef = -dot_product(v[i], v[j], chunk, comm);
            multiply_add(v[j], coef, v[i], chunk);
        }
        coef = 1. / sqrt(dot_product(v[j], v[j], chunk, comm));
        multiply(v[j], coef, chunk);
    }

    for (i = 0; i < N; i++) {
        MPI_Gatherv(v[i], chunk, MPI_DOUBLE, vector[i], counts, displs, MPI_DOUBLE, 0, comm);
    }

    for (i = 0; i < N; i++) {
        free(v[i]);
    }
    
    MPI_Comm_free(&comm);
}
