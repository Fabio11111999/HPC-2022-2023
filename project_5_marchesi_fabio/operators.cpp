//******************************************
// operators.cpp
// based on min-app code written by Oliver Fuhrer, MeteoSwiss
// modified by Ben Cumming, CSCS
// *****************************************

// Description: Contains simple operators which can be used on 3d-meshes

#include <iostream>

#include <mpi.h>

#include "data.h"
#include "operators.h"
#include "stats.h"

namespace operators {

void diffusion(const data::Field &U, data::Field &S)
{
    using data::options;
    using data::domain;

    using data::bndE;
    using data::bndW;
    using data::bndN;
    using data::bndS;

    using data::buffE;
    using data::buffW;
    using data::buffN;
    using data::buffS;

    using data::x_old;

    MPI_Comm comm_cart = domain.comm_cart;

    double dxs = 1000. * (options.dx * options.dx);
    double alpha = options.alpha;
    int nx = domain.nx;
    int ny = domain.ny;
    int iend  = nx - 1;
    int jend  = ny - 1;

	// Requests are associated to immediate send/receive operations in order to being able to wait for their completion
    MPI_Request requests[8];


	// Exchanging information with the neighbours
    if(domain.neighbour_north>=0) {
        MPI_Irecv(&bndN[0], nx, MPI_DOUBLE, domain.neighbour_north, domain.neighbour_north, comm_cart, requests);
        for(int i=0; i<nx; i++) {
			buffN[i] = U(i,ny-1);
		}
        MPI_Isend(&buffN[0], nx, MPI_DOUBLE, domain.neighbour_north, domain.rank, comm_cart, requests + 1);
    }

    if(domain.neighbour_south>=0) {
		MPI_Irecv(&bndS[0], nx, MPI_DOUBLE, domain.neighbour_south, domain.neighbour_south, comm_cart, requests + 2);
		for (int i = 0; i < nx; i++) {
			buffS[i] = U(i, 0);
		}
		MPI_Isend(&buffS[0], nx, MPI_DOUBLE, domain.neighbour_south, domain.rank, comm_cart, requests + 3);
    }

    if(domain.neighbour_east>=0) {
		MPI_Irecv(&bndE[0], ny, MPI_DOUBLE, domain.neighbour_east, domain.neighbour_east, comm_cart, requests + 4);
		for (int i = 0; i < ny; i++) {
			buffE[i] = U(nx - 1, i);
		}
		MPI_Isend(&buffE[0], ny, MPI_DOUBLE, domain.neighbour_east, domain.rank, comm_cart, requests + 5);
    }

    if(domain.neighbour_west>=0) {
		MPI_Irecv(&bndW[0], ny, MPI_DOUBLE, domain.neighbour_west, domain.neighbour_west, comm_cart, requests + 6);
		for (int i = 0; i < ny; i++) {
			buffW[i] = U(0, i);
		}
		MPI_Isend(&buffW[0], ny, MPI_DOUBLE, domain.neighbour_west, domain.rank, comm_cart, requests + 7);
    }


    // Compute the inner points
    #pragma omp parallel for
    for (int j=1; j < jend; j++) {
        for (int i=1; i < iend; i++) {
            S(i,j) = -(4. + alpha) * U(i,j)               // central point
                                    + U(i-1,j) + U(i+1,j) // east and west
                                    + U(i,j-1) + U(i,j+1) // north and south
                                    + alpha * x_old(i,j)
                                    + dxs * U(i,j) * (1.0 - U(i,j));
        }
    }

	// Waiting for the eastbound
	if (domain.neighbour_east >= 0) {
		MPI_Wait(requests + 4, MPI_STATUS_IGNORE);
	}
	// Computing the east outer points
    {
        int i = nx - 1;
        for (int j = 1; j < jend; j++)
        {
            S(i,j) = -(4. + alpha) * U(i,j)
                        + U(i-1,j) + U(i,j-1) + U(i,j+1)
                        + alpha*x_old(i,j) + bndE[j]
                        + dxs * U(i,j) * (1.0 - U(i,j));
        }
    }

    // Waiting for the westbound
	if (domain.neighbour_west >= 0) {
		MPI_Wait(requests + 6, MPI_STATUS_IGNORE);
	}
	// Compute the west outer points
    {
        int i = 0;
        for (int j = 1; j < jend; j++)
        {
            S(i,j) = -(4. + alpha) * U(i,j)
                        + U(i+1,j) + U(i,j-1) + U(i,j+1)
                        + alpha * x_old(i,j) + bndW[j]
                        + dxs * U(i,j) * (1.0 - U(i,j));
        }
    }

	// Waiting for the northbound
	if (domain.neighbour_north >= 0) {
		MPI_Wait(requests, MPI_STATUS_IGNORE);
	}
	// Computing the north outer points and the top corners
    {
        int j = ny - 1;

        {
            int i = 0; // NW corner
            S(i,j) = -(4. + alpha) * U(i,j)
                        + U(i+1,j) + U(i,j-1)
                        + alpha * x_old(i,j) + bndW[j] + bndN[i]
                        + dxs * U(i,j) * (1.0 - U(i,j));
        }

        for (int i = 1; i < iend; i++)
        {
            S(i,j) = -(4. + alpha) * U(i,j)
                        + U(i-1,j) + U(i+1,j) + U(i,j-1)
                        + alpha*x_old(i,j) + bndN[i]
                        + dxs * U(i,j) * (1.0 - U(i,j));
        }

        {
            int i = nx-1; // NE corner
            S(i,j) = -(4. + alpha) * U(i,j)
                        + U(i-1,j) + U(i,j-1)
                        + alpha * x_old(i,j) + bndE[j] + bndN[i]
                        + dxs * U(i,j) * (1.0 - U(i,j));
        }
    }

	// Waiting for the southbound
	if (domain.neighbour_south >= 0) {
		MPI_Wait(requests + 2, MPI_STATUS_IGNORE);
	}
	// Computing the south outer points and the bottom corners
    {
        int j = 0;

        {
            int i = 0; // SW corner
            S(i,j) = -(4. + alpha) * U(i,j)
                        + U(i+1,j) + U(i,j+1)
                        + alpha * x_old(i,j) + bndW[j] + bndS[i]
                        + dxs * U(i,j) * (1.0 - U(i,j));
        }

        for (int i = 1; i < iend; i++)
        {
            S(i,j) = -(4. + alpha) * U(i,j)
                        + U(i-1,j) + U(i+1,j) + U(i,j+1)
                        + alpha * x_old(i,j) + bndS[i]
                        + dxs * U(i,j) * (1.0 - U(i,j));
        }

        {
            int i = nx - 1; // SE corner
            S(i,j) = -(4. + alpha) * U(i,j)
                        + U(i-1,j) + U(i,j+1)
                        + alpha * x_old(i,j) + bndE[j] + bndS[i]
                        + dxs * U(i,j) * (1.0 - U(i,j));
        }
    }

	// Wait for the send operations to finish
        if (domain.neighbour_north >= 0) {
            MPI_Wait(requests + 1, MPI_STATUS_IGNORE);
        }
        if (domain.neighbour_south >= 0) {
            MPI_Wait(requests + 3, MPI_STATUS_IGNORE);
        }
        if (domain.neighbour_east >= 0) {
            MPI_Wait(requests + 5, MPI_STATUS_IGNORE);
        }
        if (domain.neighbour_west >= 0) {
            MPI_Wait(requests + 7, MPI_STATUS_IGNORE);
        }



    // Accumulate the flop counts
    // 8 ops total per point
    stats::flops_diff +=
        + 12 * (nx - 2) * (ny - 2) // interior points
        + 11 * (nx - 2  +  ny - 2) // NESW boundary points
        + 11 * 4;                                  // corner points
}

} // namespace operators
