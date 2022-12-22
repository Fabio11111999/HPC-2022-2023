/****************************************************************
 *                                                              *
 * This file has been written as a sample solution to an        *
 * exercise in a course given at the CSCS-USI Summer School.    *
 * It is made freely available with the understanding that      *
 * every copy of this file must include this header and that    *
 * CSCS/USI take no responsibility for the use of the enclosed  *
 * teaching material.                                           *
 *                                                              *
 * Purpose: : Parallel matrix-vector multiplication and the     *
 *            and power method                                  *
 *                                                              *
 * Contents: C-Source                                           *
 *                                                              *
 ****************************************************************/


#include <stdio.h>
#include <mpi.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include "hpc-power.h"


double* generateMatrix(int n, int numrows) {
	srand((int)time(NULL));
	double *A = (double*)calloc(n * numrows, sizeof(double));
	for (int i = 0; i < n * numrows; i++) {
		A[i] = (double) rand() / ((double) RAND_MAX + 1);
	}
	return A;
}

double norm(int n, double* x) {
	double ans = 0.0;
	for (int i = 0; i < n; i++) {
		ans += x[i] * x[i];
	}
	return sqrt(ans);
}

double* matVec(double* A, double *x, int n, int numrows) {
	double *ans = (double*)calloc(numrows, sizeof(double));
	for (int i = 0; i < numrows; i++) {
		for (int j = 0; j < n; j++) {
			ans[i] += A[i * n + j] * x[j];
		}
	}
	return ans;
}

int main (int argc, char *argv[]) {
	/*
	This code uses the provided matrix generator and does 1000 iterations.
	The size of the side of the matrix is hard coded to 1024.
	The size of the side of the matrix needs to be a multiple of the number of processes.
	*/

    int my_rank, size;
    int snd_buf, rcv_buf;
    int right, left;
    int sum, i;
	int iterations = 1000;
    MPI_Status  status;
    MPI_Request request;


    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    MPI_Comm_size(MPI_COMM_WORLD, &size);

    /* This subproject is about to write a parallel program to
       multiply a matrix A by a vector x, and to use this routine in
       an implementation of the power method to find the absolute
       value of the largest eigenvalue of the matrix. Your code will
       call routines that we supply to generate matrices, record
       timings, and validate the answer.
    */
	int N = 1024;
	int startrow = N / size * my_rank;
	int endrow = N / size * (my_rank + 1) - 1;
	int sz = endrow - startrow + 1;
	// Uncomment the following line and comment the next one to generate the Matrix using my own method (the check on the result doesn't work on the random matrix)
	// double* A = generateMatrix(N, sz);
	double* A = hpc_generateMatrix(N, startrow, sz);
	double* x = generateMatrix(N, 1);
	for (int i = 0; i < N; i++) {
		x[i] = 1.0;
	}
	double start_time = hpc_timer();
	for (int i = 0; i < iterations; i++) {
		// Process 0 computes the norm
		// And divides x by it
		if (my_rank == 0) {
			double current_norm = norm(N, x);
			for (int j = 0; j < N; j++) {
				x[j] /= current_norm;
			}
		}
		// Need to share new vector x
		int ierr = MPI_Bcast(x, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		// Each process compute the part of the product that uses its rows
		double* prod = matVec(A, x, N, sz);

		ierr = MPI_Gather(prod, sz, MPI_DOUBLE, x, sz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
	if (my_rank == 0) {
		double total_time = hpc_timer() - start_time;
		printf("Size of the matrix: %d, Number of processes: %d, Time required: %fs\n", N, size, total_time);
		// This check can be used only if we used hpc_generateMatrix()
		if (abs(norm(N, x) - (double)N) > 0.0001) {
			printf("Wrong computation\n");
		}
	}
    MPI_Finalize();
    return 0;
}
