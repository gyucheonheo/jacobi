OPSYS = $(shell uname -s)
CC=gcc
CFLAGS=-Wall -fopenmp
LFLAGS=-lm

TARGET=jacobi jacobi_mp

all: ${TARGET}

help:; @echo "  "
	@echo "Operating System Detected: $(OPSYS)"
	@echo "  "
	@echo "USAGE: "
	@echo "make help	To get this listing"
	@echo "make		To compile all jacobi program in current environment"
	@echo "make jacobi	To compile seriel-version jacobi program"
	@echo "make jacobi_mp	To compile openMP-version jacobi program"
	@echo "make clean	Remove *.o and executable files"
	@echo "  "

jacobi: jacobi.c
	$(CC) jacobi.c -o jacobi -lm
clean:
	rm -f ${TARGET} ${TARGET:=.o}

.PHONY: all clean
