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

// input: s, gives updated solution for f
// only handles interior grid points, as boundary points are fixed
// those inner grid points neighbouring a boundary point, will in the following
// be referred to as boundary points and only those grid points
// only neighbouring non-boundary points are called inner grid points
void diffusion(const data::Field &s, data::Field &f)
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

    using data::y_old;

	MPI_Comm comm_cart = domain.comm_cart;

    double alpha = options.alpha;
    double beta = options.beta;

    int nx = domain.nx;
    int ny = domain.ny;
    int iend  = nx - 1;
    int jend  = ny - 1;

    MPI_Request requests[8];


    // Exchanging information with the neighbours
    if(domain.neighbour_north>=0) {
        MPI_Irecv(&bndN[0], nx, MPI_DOUBLE, domain.neighbour_north, domain.neighbour_north, comm_cart, requests);
        for(int i=0; i<nx; i++) {
            buffN[i] = s(i,ny-1);
        }
        MPI_Isend(&buffN[0], nx, MPI_DOUBLE, domain.neighbour_north, domain.rank, comm_cart, requests + 1);
    }

    if(domain.neighbour_south>=0) {
        MPI_Irecv(&bndS[0], nx, MPI_DOUBLE, domain.neighbour_south, domain.neighbour_south, comm_cart, requests + 2);
        for (int i = 0; i < nx; i++) {
            buffS[i] = s(i, 0);
        }
        MPI_Isend(&buffS[0], nx, MPI_DOUBLE, domain.neighbour_south, domain.rank, comm_cart, requests + 3);
    }

    if(domain.neighbour_east>=0) {
        MPI_Irecv(&bndE[0], ny, MPI_DOUBLE, domain.neighbour_east, domain.neighbour_east, comm_cart, requests + 4);
        for (int i = 0; i < ny; i++) {
            buffE[i] = s(nx - 1, i);
        }
        MPI_Isend(&buffE[0], ny, MPI_DOUBLE, domain.neighbour_east, domain.rank, comm_cart, requests + 5);
    }

    if(domain.neighbour_west>=0) {
        MPI_Irecv(&bndW[0], ny, MPI_DOUBLE, domain.neighbour_west, domain.neighbour_west, comm_cart, requests + 6);
        for (int i = 0; i < ny; i++) {
            buffW[i] = s(0, i);
        }
        MPI_Isend(&buffW[0], ny, MPI_DOUBLE, domain.neighbour_west, domain.rank, comm_cart, requests + 7);
    }

    // the interior grid points
    for (int j=1; j < jend; j++) {
        for (int i=1; i < iend; i++) {
            f(i,j) = -(4. + alpha) * s(i,j)               // central point
                                    + s(i-1,j) + s(i+1,j) // east and west
                                    + s(i,j-1) + s(i,j+1) // north and south
                                    + alpha * y_old(i,j)
                                    + beta * s(i,j) * (1.0 - s(i,j));
        }
    }

    if (domain.neighbour_east >= 0) {
        MPI_Wait(requests + 4, MPI_STATUS_IGNORE);
    }
    // the east boundary
    {
        int i = nx - 1;
        for (int j = 1; j < jend; j++)
        {
            f(i,j) = -(4. + alpha) * s(i,j)
                        + s(i-1,j) + s(i,j-1) + s(i,j+1)
                        + alpha*y_old(i,j) + bndE[j]
                        + beta * s(i,j) * (1.0 - s(i,j));
        }
    }

    if (domain.neighbour_west >= 0) {
        MPI_Wait(requests + 6, MPI_STATUS_IGNORE);
    }
    {
        int i = 0;
        for (int j = 1; j < jend; j++)
        {
            f(i,j) = -(4. + alpha) * s(i,j)
                        + s(i+1,j) + s(i,j-1) + s(i,j+1)
                        + alpha * y_old(i,j) + bndW[j]
                        + beta * s(i,j) * (1.0 - s(i,j));
        }
    }

    if (domain.neighbour_north >= 0) {
        MPI_Wait(requests, MPI_STATUS_IGNORE);
    }
    {
        int j = ny - 1;

        {
            int i = 0; // NW corner
            f(i,j) = -(4. + alpha) * s(i,j)
                        + s(i+1,j) + s(i,j-1)
                        + alpha * y_old(i,j) + bndW[j] + bndN[i]
                        + beta * s(i,j) * (1.0 - s(i,j));
        }

        // north boundary
        for (int i = 1; i < iend; i++)
        {
            f(i,j) = -(4. + alpha) * s(i,j)
                        + s(i-1,j) + s(i+1,j) + s(i,j-1)
                        + alpha*y_old(i,j) + bndN[i]
                        + beta * s(i,j) * (1.0 - s(i,j));
        }

        {
            int i = nx-1; // NE corner
            f(i,j) = -(4. + alpha) * s(i,j)
                        + s(i-1,j) + s(i,j-1)
                        + alpha * y_old(i,j) + bndE[j] + bndN[i]
                        + beta * s(i,j) * (1.0 - s(i,j));
        }
    }

    if (domain.neighbour_south >= 0) {
        MPI_Wait(requests + 2, MPI_STATUS_IGNORE);
    }
    {
        int j = 0;

        {
            int i = 0; // SW corner
            f(i,j) = -(4. + alpha) * s(i,j)
                        + s(i+1,j) + s(i,j+1)
                        + alpha * y_old(i,j) + bndW[j] + bndS[i]
                        + beta * s(i,j) * (1.0 - s(i,j));
        }

        // south boundary
        for (int i = 1; i < iend; i++)
        {
            f(i,j) = -(4. + alpha) * s(i,j)
                        + s(i-1,j) + s(i+1,j) + s(i,j+1)
                        + alpha * y_old(i,j) + bndS[i]
                        + beta * s(i,j) * (1.0 - s(i,j));
        }

        {
            int i = nx - 1; // SE corner
            f(i,j) = -(4. + alpha) * s(i,j)
                        + s(i-1,j) + s(i,j+1)
                        + alpha * y_old(i,j) + bndE[j] + bndS[i]
                        + beta * s(i,j) * (1.0 - s(i,j));
        }
    }
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
