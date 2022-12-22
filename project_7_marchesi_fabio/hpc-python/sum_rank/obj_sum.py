from mpi4py import MPI
import numpy as np

# get comm, size & rank
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = np.array([comm.Get_rank()])
total_sum = np.array([0])
comm.Allreduce(rank, total_sum, op=MPI.SUM)
print(f"Rank: {rank[0]}. Sum: {total_sum[0]}")


