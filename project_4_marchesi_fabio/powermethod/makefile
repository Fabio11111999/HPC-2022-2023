CC=mpicc
CFLAGS=-O3 -lm

all: powermethod

powermethod: hpc-power.c powermethod.c
	$(CC) $(CFLAGS) $^ -o $@

clean:
	rm -f powermethod
