#!/bin/bash
#$ -V
#$ -cwd
#$ -N ksc2018_job
#$ -pe mpi_8cpu 32
#$ -q ksc2018@tachyon,ksc2018@tachyon,ksc2018@tachyon,ksc2018@tachyon
#$ -R yes
#$ -l h_rt=01:00:00
#$ -l OMP_NUM_THREADS=1
export OMP_NUM_THREADS=1
mpirun -machinefile $TMPDIR/machines -np $NSLOTS ./gs