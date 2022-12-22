addpaths_GP;
A = blkdiag(ones(5), ones(5));
A(1, 10) = 1;
A(10, 1) = 1;
A(5, 6) = 1;
A(6, 5) = 1;
[p1, p2] = bisection_metis(sparse(A), 0, 0)
