CC=gcc

all: jacobi

jacobi: jacobi.c
	$(CC) jacobi.c -o jacobi -lm

.PHONY: all	
