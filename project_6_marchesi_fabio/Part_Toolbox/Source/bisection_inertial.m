%
% D.P & O.S for the "HPC Course" at USI and
%                   "HPC Lab for CSE" at ETH Zurich

function [part1,part2] = bisection_inertial(A,xy,picture)
% bisection_inertial : Inertial partition of a graph.
%
% [p1,p2] = bisection_inertial(A,xy) returns a list of the vertices on one side of a partition
%     obtained by bisection with a line or plane normal to a moment of inertia
%     of the vertices, considered as points in Euclidean space.
%     Input A is the adjacency matrix of the mesh (used only for the picture!);
%     each row of xy is the coordinates of a point in d-space.
%
% bisection_inertial(A,xy,1) also draws a picture.

% Steps
% 1. Calculate the center of mass.
x_center = mean(xy(:,1));
y_center = mean(xy(:,2));
% 2. Construct the matrix M.
s_xx = 0;
s_yy = 0;
s_xy = 0;

for i=(1:size(xy(:,1)))
    x_diff = xy(i,1) - x_center;
    y_diff = xy(i,2) - y_center;
    s_xx = s_xx + x_diff^2;
    s_yy = s_yy + y_diff^2;
    s_xy = s_xy + x_diff * y_diff;
end
M = [s_yy s_xy; s_xy s_xx];
% 3. Calculate the smallest eigenvector of M.  
[V, D] = eig(M);
eig_vec = V(:, 1);
% 4. Find the line L on which the center of mass lies.
L = eig_vec;
% 5. Partition the points around the line L.
%   (you may use the function partition.m)
[part1, part2] = partition(xy, L);

if picture == 1
    gplotpart(A,xy,part1);
    title('Inertial bisection using the Fiedler Eigenvector');
end

% Dummy implementation to generate a partitioning >>>>

end
