/* 
    Please include compiler name below (you may also include any other modules you would like to be loaded)

COMPILER= gnu

    Please include All compiler flags and libraries as you want them run. You can simply copy this over from the Makefile's first few lines
 
CC = cc
OPT = -O3
CFLAGS = -Wall -std=gnu99 $(OPT)
MKLROOT = /opt/intel/composer_xe_2013.1.117/mkl
LDLIBS = -lrt -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm

*/

#include<stdlib.h>
#include <stdio.h>

const char* dgemm_desc = "Naive, three-loop dgemm.";

/* This routine performs a dgemm operation
 *  C := C + A * B
 * where A, B, and C are lda-by-lda matrices stored in column-major format.
 * On exit, A and B maintain their input values. */    
void square_dgemm(int n, double* A, double *B, double *C) {
	int sz = 32;
	int padded_n = sz * ((n + sz - 1) / sz);
	double *padded_A = (double*)aligned_alloc(64, padded_n * padded_n * sizeof(double));
	double *padded_B = (double*)aligned_alloc(64, padded_n * padded_n * sizeof(double));
	double *padded_C = (double*)aligned_alloc(64, padded_n * padded_n * sizeof(double));

	for (int i = 0; i < padded_n; i++) {
		for (int j = 0; j < padded_n; j++) {
			if (i >= n || j >= n) {
				padded_A[i * padded_n + j] = 0.0;
				padded_B[i * padded_n + j] = 0.0;
			} else {
				padded_A[i * padded_n + j] = A[j * n + i];
				padded_C[i * padded_n + j] = 0.0;
			}
		}
	}
	for (int j = 0; j < n; j++) {
		for (int i = 0; i < n; i++) {
				padded_B[j * padded_n + i] = B[j * n + i];
		}
	}
	for (int i_block = 0; i_block < padded_n; i_block += sz) {
		for (int j_block = 0; j_block < padded_n; j_block += sz) {
			for (int k_block = 0; k_block < padded_n; k_block += sz) {
				for (int i = i_block; i < i_block + sz; i++) {
					for (int j = j_block; j < j_block + sz; j++) {
						for (int k = k_block; k < k_block + sz; k++) {
							padded_C[padded_n * j + i] += padded_A[i * padded_n + k] * padded_B[j * padded_n + k];
						}
					}
				}
			}
		}
	}

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			C[i * n + j] = padded_C[i * padded_n + j];
		}
	}
	free(padded_A);
	free(padded_B);
	free(padded_C);
}
