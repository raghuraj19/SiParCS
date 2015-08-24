# compiler option
#CC = gcc
CC = icc
NVCC = nvcc 

CUDA_PATH = /usr/local/cuda-6.5

INCLUDES := -I$(CUDA_PATH)/include
LIBS := -L$(CUDA_PATH)/lib64/

INC = ./include

SRCS = timer.c blosseydurran.c loadinputdata.c geninputdata.c 

LIB = 

EXE = BlosseyDurranCUDA

all:
	$(NVCC) -O3 -arch=sm_35 -Xptxas -dlcm=cg -I$(INC) $(INCLUDES) bdcuda.cu -c
	$(CC) -I$(INC) $(SRCS) $(LIB) -c  -lm
	$(CC) *.o -o $(EXE) $(LIBS) -lcudart -lm 
clean: 
	rm ./BlosseyDurranCUDA
