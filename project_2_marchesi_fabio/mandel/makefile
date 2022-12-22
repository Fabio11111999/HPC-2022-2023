all: mandel_seq mandel_parallel

mandel_seq: mandel_seq.c pngwriter.c
	gcc -o $@ -I. -O3 $^ -lpng

mandel_parallel: mandel_parallel.c pngwriter.c
	gcc -o $@ -I. -O3 $^ -lpng -fopenmp
clean:
	rm -rf mandel_seq

