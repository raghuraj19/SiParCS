# compiler option
#CC = gcc
CC = mpicc

INC = ./include

SRCS = timer.c blosseydurran.c loadinputdata.c geninputdata.c vector.c

LIB = 

EXE = BlosseyDurranMPI

all:
	$(CC) -I$(INC) $(SRCS) $(LIB) -o $(EXE) -lm -w

clean: 
	rm ./BlosseyDurranMPI
