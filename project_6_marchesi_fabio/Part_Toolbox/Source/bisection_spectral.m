function [part1,part2] = bisection_spectral(A,xy,picture)
% bisection_spectral : Spectral partition of a graph.
%
% D.P & O.S for the "HPC Course" at USI and
%                   "HPC Lab for CSE" at ETH Zuric
%
% [part1,part2] = bisection_spectral(A) returns a partition of the n vertices
%                 of A into two lists part1 and part2 according to the
%                 spectral bisection algorithm of Simon et al:
%                 Label the vertices with the components of the Fiedler vector
%                 (the second eigenvector of the Laplacian matrix) and partition
%                 them around the median value or 0.



% Steps
% 1. Construct the Laplacian.
rows = size(A, 1);
cols = size(A, 2);
L = sparse(rows, cols);
[ii, jj, ss] = find(A);
for k = 1:size(ii,1)
    i = ii(k);
    j = jj(k);
    L(i, j) = L(i, j) - ss(k);
    L(i, i) = L(i, i) + ss(k);
end
% 2. Calculate its eigensdecomposition.
[eig_vec, eig_val] = eigs(L, 2, 1e-14);
eig_vec = eig_vec(:, 2);

eig_vec_median = median(eig_vec);
% 3. Label the vertices with the components of the Fiedler vector.
part1 = [];
part2 = [];
% 4. Partition them around their median value, or 0.
for i = 1 : rows
    if eig_vec(i) < eig_vec_median
        part1 = [part1, i];
    else
        part2 = [part2, i];
    end
end


if picture == 1
    gplotpart(A,xy,part1);
    title('Spectral bisection using the Fiedler Eigenvector');
end


end
