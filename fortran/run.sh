#!/bin/bash

export OMP_NUM_THREADS=1
mpirun -np 2 ./gs
