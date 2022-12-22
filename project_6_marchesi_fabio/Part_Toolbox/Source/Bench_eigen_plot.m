% Visualize information from the eigenspectrum of the graph Laplacian
%
% D.P & O.S for the "HPC Course" at USI and
%                   "HPC Lab for CSE" at ETH Zurich

% add necessary paths
addpaths_GP;

% Graphical output at bisection level
picture = 0;
cases = {
    'airfoil1.mat';
    '3elt.mat';
    'barth4.mat';
    'mesh3e1.mat';
    'crack.mat';
    };

nc = length(cases);

for c=1:nc
    load(cases{c});
    W      = Problem.A;
    coords = Problem.aux.coord;
    rows = size(W, 1);
    cols = size(W, 2);
    L = sparse(rows, cols);
    [ii, jj, ss] = find(W);

    for k = 1:size(ii,1)
        i = ii(k);
        j = jj(k);
        L(i, j) = L(i, j) - ss(k);
        L(i, i) = L(i, i) + ss(k);
    end
    [part1,part2] = bisection_spectral(W,coords,0);
    % 2. Compute eigenvectors associated with the smallest eigenvalues.
    [eig_vec, eig_val] = eigs(L, 3, 1e-15);
    % Printing the eigenvectors for airfoil1
    if c == 1
        disp('Eigenvector associated to the smallest Eigenvalue');
        plot(eig_vec(:, 1));
        pause();
        disp('Eigenvector associated to the second smallest Eigenvalue');
        plot(eig_vec(:, 2));
        pause();  
    end
    % Plotting the partition using the second smallest eigenvector as the
    % third dimension
    v2 = eig_vec(:, 2);
    %verify that %axis equal; is commented on line 45 of gplotg.m
    gplotpart3(W,[coords v2], part1);
    pause;
    coords1 = coords;
    v3 = eig_vec(:, 3);
    for i=1:rows
        coords1(i, 1) = v2(i);
        coords1(i, 2) = v3(i);
    end
    gplotpart(W, coords1, part1);
    pause;

end
