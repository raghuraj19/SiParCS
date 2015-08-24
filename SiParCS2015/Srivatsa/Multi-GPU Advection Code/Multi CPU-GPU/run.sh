#!/bin/tcsh


#BSUB -W 00:05

#BSUB -J TauGB

#BSUB -o TauGB.out

#BSUB -e TauGB.err

#BSUB -n 8

#BSUB -R "span[ptile=2]"

#BSUB -q gpgpu

#BSUB -P P86850055
module load cuda/6.5
module load intel/15.0.3

setenv LD_LIBRARY_PATH /ncar/opt/cuda/6.5/lib64:$LD_LIBRARY_PATH

mpirun.lsf ./BlosseyDurranMPI_CUDA partitions_401_8.txt neighbors_401.txt

