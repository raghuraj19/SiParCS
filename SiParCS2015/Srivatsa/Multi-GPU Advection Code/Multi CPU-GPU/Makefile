# compiler option
#CC = gcc
CC = mpicc
NVCC = nvcc 

CUDA_PATH = /ncar/opt/cuda/6.5

INCLUDES := -I$(CUDA_PATH)/include
LIBS := -L$(CUDA_PATH)/lib64/

INC = ./include

SRCS = timer.c blosseydurran.c loadinputdata.c geninputdata.c vector.c

LIB = 

EXE = BlosseyDurranMPI_CUDA

all:
	$(NVCC)-I$(INC) $(INCLUDES) bdcuda.cu -c
	$(CC) -I$(INC) $(SRCS) $(LIB) -c  -lm -w
	$(CC) *.o -o $(EXE) $(LIBS) -lcudart -lm 
clean: 
	rm ./BlosseyDurranMPI_CUDA
