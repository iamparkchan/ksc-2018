/*  serial code */
#include "gram_schmidt.h"
#include <mpi.h>
#include <math.h>

static double dot_product(double *A, double *B)
{
    int k;
    double sum = 0.;
    
#pragma omp parallel for reduction(+:sum)
    for (k = 0; k < N; k++) {
        sum += A[k] * B[k];
    }
    
    return sum;
}

static void multiply_add(double *A, double c, double *B)
{
    int k;

#pragma omp parallel for
    for (k = 0; k < N; k++) {
        A[k] += c * B[k];
    }
}

static void multiply(double *A, double c)
{
    int k;
    
#pragma omp parallel for
    for (k = 0; k < N; k++) {
        A[k] *= c; 
    }
}

void gram_schmidt(double **vector)
{
    int rank;
    int i, j;
    double coef;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    if (0 == rank) {
        for (j = 0; j < M; j++) {
            for (i = 0; i < j; i++) {
                coef = -dot_product(vector[i], vector[j]);
                multiply_add(vector[j], coef, vector[i]);
            }
            coef = 1. / sqrt(dot_product(vector[j], vector[j]));
            multiply(vector[j], coef);
        }
    }
}




/* OpenMP code */
/*
#include "gram_schmidt.h"
#include <mpi.h>
#include <omp.h>

#define CEIL_OF_DIVISION(a,b) (((a) + (b) - 1)/(b))

static double dot_product(double *A, double *B, int size)
{
  int k;
  double sum = 0.;

  for (k = 0; k < size; k++) {
    sum += A[k] * B[k];
  }

  return sum;
}

static void multiply_add(double *A, double c, double *B, int size)
{
  int k;

  for (k = 0; k < size; k++) {
    A[k] += c * B[k];
  }
}

void gram_schmidt(double **vector)
{
  int rank;
  int i, j, k;
  double dot[N] = {0}, partial_dot;
  double ***slice;
  int num_threads, chunk, tid;
  double coef;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (0 == rank) {
    #pragma omp parallel
    {
      num_threads = omp_get_num_threads();
    }

    slice = malloc(sizeof(double) * num_threads);
    chunk = CEIL_OF_DIVISION(N, num_threads);
    for (k = 0; k < num_threads; k++) {
      slice[k] = malloc(sizeof(double) * N);
      if (k == num_threads - 1) {
        chunk = N - chunk*(num_threads - 1);
      }
      for (i = 0; i < N; i++) {
        slice[k][i] = malloc(sizeof(double) * chunk);
        for (j = 0; j < chunk; j++) {
          slice[k][i][j] = vector[i][chunk*(k - 1) + j];
        }
      }
    }

    #pragma omp parallel private(tid, chunk, partial_dot)
    {
      tid = omp_get_thread_num();
      if (tid == num_threads - 1) {
        chunk = N - CEIL_OF_DIVISION(N, num_threads)*(num_threads - 1);
      } else {
        chunk = CEIL_OF_DIVISION(N, num_threads);
      }

      for (j = 1; j < N; j++) {
        partial_dot = dot_product(slice[tid][j - 1], slice[tid][j - 1], chunk);
        #pragma omp atomic
          dot[j - 1] += partial_dot;

        #pragma omp barrier
        for (i = 0; i < j; i++) {
          partial_dot = - dot_product(slice[tid][i], slice[tid][j], chunk) / dot[i];
          coef = 0.;
          #pragma omp barrier
          #pragma omp atomic
            coef += partial_dot;
          #pragma omp barrier
          // multiply_add(slice[tid][j], coef, slice[tid][i], chunk);
        }
      }
    }

    chunk = CEIL_OF_DIVISION(N, num_threads);
    for (k = 0; k < num_threads; k++) {
      if (k == num_threads - 1) {
        chunk = N - chunk*(num_threads - 1);
      }
      for (i = 0; i < N; i++) {
        for (j = 0; j < chunk; j++) {
          vector[i][chunk*(k - 1) + j] = slice[k][i][j];
        }
      }
    }


    for (k = 0; k < num_threads; k++) {
      for (i = 0; i < N; i++) {
        free(slice[k][i]);
      }
      free(slice[k]);
    }
    free(slice);
  }
}
*/
