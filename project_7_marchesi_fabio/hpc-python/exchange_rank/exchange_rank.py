from mpi4py import MPI

# get comm, size & rank
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
dims = [0, 0]
dims = MPI.Compute_dims(size, dims)
periods = [True, True]
comm_cart = comm.Create_cart(dims, periods=periods)
coords = comm_cart.Get_coords(rank)
neigh_west, neigh_east = comm_cart.Shift(1, -1)
neigh_south, neigh_north = comm_cart.Shift(0, -1)

west_rank = comm_cart.sendrecv(rank, neigh_east, neigh_west)
east_rank = comm_cart.sendrecv(rank, neigh_west, neigh_east)
south_rank = comm_cart.sendrecv(rank, neigh_north, neigh_south)
north_rank = comm_cart.sendrecv(rank, neigh_south, neigh_north)

print(f'Process: {rank}. Coordinates: {coords}. West: {neigh_west}. East: {neigh_east}. South: {neigh_south}. North: {neigh_north}')

print(f'Process: {rank}. West neighbour rank: {west_rank}')
print(f'Process: {rank}. East neighbour rank: {east_rank}')
print(f'Process: {rank}. South neighbour rank: {south_rank}')
print(f'Process: {rank}. North neighbour rank: {north_rank}')

