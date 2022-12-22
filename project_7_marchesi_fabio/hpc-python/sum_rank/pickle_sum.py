from mpi4py import MPI
import numpy as np

# get comm, size & rank
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
total_sum = comm.allreduce(rank, op=MPI.SUM)
print(f"Rank: {rank}. Sum: {total_sum}")


