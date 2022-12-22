/****************************************************************
 *                                                              *
 * This file has been written as a sample solution to an        *
 * exercise in a course given at the CSCS-USI Summer School.    *
 * It is made freely available with the understanding that      *
 * every copy of this file must include this header and that    *
 * CSCS/USI take no responsibility for the use of the enclosed  *
 * teaching material.                                           *
 *                                                              *
 * Purpose: Parallel sum using a ping-pong                      *
 *                                                              *
 * Contents: C-Source                                           *
 *                                                              *
 ****************************************************************/


#include <stdio.h>
#include <mpi.h>


int main (int argc, char *argv[])
{
    int my_rank, size;
    int snd_buf, rcv_buf;
    int right, left;
    int max = 0, i;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    MPI_Comm_size(MPI_COMM_WORLD, &size);

	MPI_Status  status;
	MPI_Request request;

	right = (my_rank + 1) % size;
	left = (my_rank - 1 + size) % size;
	max = (3 * my_rank) % (2 * size);
	for (int i = 0, snd_buf = max; i < size - 1; i++, snd_buf = max) {
		MPI_Issend(&snd_buf, 1, MPI_INT, right, 1, MPI_COMM_WORLD, &request);
		MPI_Recv(&rcv_buf, 1, MPI_INT, left, 1, MPI_COMM_WORLD, &status);
		MPI_Wait(&request, &status);
		if (rcv_buf > max) {
			max = rcv_buf;
		}
	}

    printf ("Process %i:\tMax = %i\n", my_rank, max);

    MPI_Finalize();
    return 0;
}
