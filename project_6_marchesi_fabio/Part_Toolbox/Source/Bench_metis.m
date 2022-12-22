function [cut_recursive,cut_kway] = Bench_metis(picture)
% Compare recursive bisection and direct k-way partitioning,
% as implemented in the Metis 5.0.2 library.

%  Add necessary paths
addpaths_GP;
visualize = false;

cases = {
    '../Datasets/Roads/luxembourg_osm.mat';
    '../Datasets/Roads/usroads.mat';
    };

nc = length(cases);


for c = 1:nc
    load (cases{c});
    W = Problem.A;
    coords = Problem.aux.coord;
    if visualize == true
        gplotg(W,coords);
        pause;
    end
    [map_rec_16, edgecut_rec_16] = metismex('PartGraphRecursive', W, 16);
    fprintf('Map: %s, Partitions: %d, Method: bisection, Cutsize: %i\n',cases{c},16,edgecut_rec_16);
    if visualize == true
        gplotmap(W,coords, map_rec_16);
        pause;
    end
    [map_k_16, edgecut_k_16] = metismex('PartGraphKway', W, 16);
    fprintf('Map: %s, Partitions: %d, Method: K-way, Cutsize: %i\n',cases{c},16,edgecut_k_16)
    if visualize == true
        gplotmap(W,coords, map_k_16);
        pause;
    end
    [map_rec_32, edgecut_rec_32] = metismex('PartGraphRecursive', W, 32);
    fprintf('Map: %s, Partitions: %d, Method: bisection, Cutsize: %i\n',cases{c},32,edgecut_rec_32);
    if visualize == true
        gplotmap(W,coords, map_rec_32);
        pause;
    end
    [map_k_32, edgecut_k_32] = metismex('PartGraphKway', W, 32);
    fprintf('Map: %s, Partitions: %d, Method: K-way, Cutsize: %i\n',cases{c},32,edgecut_k_32)
    if visualize == false
        gplotmap(W,coords, map_k_32);
        pause;
    end
end

cases = {
    '../Datasets/Countries_Meshes/mat/gr.mat';
    '../Datasets/Countries_Meshes/mat/ch.mat';
    '../Datasets/Countries_Meshes/mat/vn.mat';
    '../Datasets/Countries_Meshes/mat/no.mat';
    '../Datasets/Countries_Meshes/mat/ru.mat';
};

nc = length(cases);
for c = nc:nc
    load (cases{c});
    if visualize == true
        gplotg(adj,positions);
        pause;
    end
    [map_rec_16, edgecut_rec_16] = metismex('PartGraphRecursive', adj + transpose(adj), 16);
    fprintf('Map: %s, Partitions: %d, Method: bisection, Cutsize: %i\n',cases{c},16,edgecut_rec_16);
    if visualize == true
        gplotmap(adj,positions, map_rec_16);
        pause;
    end
    [map_k_16, edgecut_k_16] = metismex('PartGraphKway', adj, 16);
    fprintf('Map: %s, Partitions: %d, Method: K-way, Cutsize: %i\n',cases{c},16,edgecut_k_16)
    if visualize == true
        gplotmap(adj,positions, map_k_16);
        pause;
    end
    [map_rec_32, edgecut_rec_32] = metismex('PartGraphRecursive', adj, 32);
    fprintf('Map: %s, Partitions: %d, Method: bisection, Cutsize: %i\n',cases{c},32,edgecut_rec_32);
    if visualize == true
        gplotmap(adj,positions, map_rec_32);
        pause;
    end
    [map_k_32, edgecut_k_32] = metismex('PartGraphKway', adj, 32);
    fprintf('Map: %s, Partitions: %d, Method: K-way, Cutsize: %i\n',cases{c},32,edgecut_k_32)
    if visualize == true
        gplotmap(adj,positions, map_k_32);
        pause;
    end
end



% Graphs in question
% load usroads;
% load luxembourg_osm;

% Steps
% 1. Initialize the cases
% 2. Call metismex to
%     a) Recursively partition the graphs in 16 and 32 subsets.
%     b) Perform direct k-way partitioning of the graphs in 16 and 32 subsets.
% 3. Visualize the results for 32 partitions.


end