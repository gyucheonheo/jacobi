OPSYS = $(shell uname -s)
CC=gcc
CFLAGS=-Wall -fopenmp
LIBS=-lm

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
	$(CC) $(CFLAGS) jacobi.c -o jacobi $(LIBS)
jacobi_mp: jacobi_mp.c
	$(CC) $(CFLAGS) jacobi_mp.c -o jacobi_mp $(LIBS)
test:
	./jacobi -n 500 -i 25 -c 0.0001
	./jacobi_mp -t 8 -n 500 -i 25 -c 0.0001
test1:
	./jacobi -n 5000 -i 25 -c 0.0001
	./jacobi_mp -t 8 -n 5000 -i 25 -c 0.0001
test2:
	./jacobi -n 10000 -i 25 -c 0.0001
	./jacobi_mp -t t -n 10000 -i 25 -c 0.0001

clean:
	rm -f ${TARGET} ${TARGET:=.o}

.PHONY: all clean
