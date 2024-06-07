#!/bin/bash
#PBS -N il_time_2
#PBS -A UMIC0093
#PBS -l walltime=1:00:00
#PBS -q main
#PBS -j oe
#PBS -k eod
#PBS -m abe
#PBS -M cygnari@umich.edu
#PBS -l select=2:ncpus=128:mpiprocs=128
#PBS -l place=group=rack

export TMPDIR=/glade/derecho/scratch/$USER/temp
mkdir -p $TMPDIR

mpiexec ../build/executables/inverse_laplacian > $TMPDIR/run_out2.txt
