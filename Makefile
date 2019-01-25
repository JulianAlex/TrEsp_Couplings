#================================================================
#
#   Makefile Couplings, Julian Adolphs 2015
#
#================================================================
#
CC_FLAGS = -Wall

LD_FLAGS = -lm   

OBJ = couplings_main.o couplings_calc.o couplings_input.o couplings_output.o\

SRC = couplings_main.c couplings_calc.c couplings_input.c couplings_output.c\

HDR = couplings.h couplings_func.h

CC = gcc  #-ggdb -Wall    # Compiler

EXECUTABLE = couplings



$(EXECUTABLE) : $(SRC) $(HDR)                                        # $(OBJ)
	$(CC) $(CC_FLAGS) -o couplings $(SRC) $(LD_FLAGS)              

