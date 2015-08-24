# compiler option
#CC = gcc
CC = icc

INC = ./include

SRCS = timer.c blosseydurran.c loadinputdata.c geninputdata.c 

LIB = 

EXE = BlosseyDurranSerial

all:
	$(CC) -O3 -no-prec-div -ipo -xHost -I$(INC) $(SRCS) $(LIB) -o $(EXE)

clean: 
	rm ./BlosseyDurranSerial
