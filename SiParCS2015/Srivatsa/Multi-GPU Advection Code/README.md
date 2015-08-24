# blosseydurran

Serial Code

1. Change directory to 'Serial C'
2. Copy the input text files from one of the folders (51x51, 101x101,....) to the current directory
3. Run 'make'
4. Execute the output using ./BlosseyDurranSerial

Single GPU Code

1. Change directory to 'Single GPU'
2. Copy the input text files from one of the folders (51x51, 101x101,....) to the current directory
3. Run 'make'
4. Execute the output using ./BlosseyDurranCUDA

Multi CPU Code

1. Change directory to 'Multi CPU'
2. Copy the input text files from one of the folders (51x51, 101x101,....) to the current directory
3. Copy the neighbor and partition input files from respective folders to the current directory
4. Change NUM_NODES and NUM_PARTITIONS in blosseydurran.c accordingly
5. Run 'make'
6. Execute the output using mpirun -n # ./BlosseyDurranMPI partitions.txt neighbors.txt (OR) run.sh script

Multi CPU-GPU Code

1. Change directory to 'Multi CPU-GPU'
2. Copy the input text files from one of the folders (51x51, 101x101,....) to the current directory
3. Copy the neighbor and partition input files from respective folders to the current directory
4. Change NUM_NODES and NUM_PARTITIONS in blosseydurran.c accordingly
5. Run 'make'
6. Execute the output using mpirun -n # ./BlosseyDurranMPI_CUDA partitions.txt neighbors.txt (OR) run.sh script
