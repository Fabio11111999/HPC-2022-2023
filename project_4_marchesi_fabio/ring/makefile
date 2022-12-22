CC=mpicc
CFLAGS=-O3

all: max_ring

sum_ring: max_ring.c
	$(CC) $(CFLAGS) $^ -o $@

clean:
	rm -f max_ring

