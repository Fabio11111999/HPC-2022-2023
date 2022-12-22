% Script to load .csv lists of adjacency matrices and the corresponding 
% coordinates. 
% The resulting graphs should be visualized and saved in a .csv file.
%
% D.P & O.S for the "HPC Course" at USI and
%                   "HPC Lab for CSE" at ETH Zurich

addpaths_GP;

% Steps
% 1. Load the .csv files
adjs = {
    '..\Datasets\Countries_Meshes\csv/CH-4468-adj.csv';
    '..\Datasets\Countries_Meshes\csv/CL-13042-adj.csv';
    '..\Datasets\Countries_Meshes\csv/GB-5946-adj.csv';
    '..\Datasets\Countries_Meshes\csv/GR-3117-adj.csv';
    '..\Datasets\Countries_Meshes\csv/NO-9935-adj.csv';
    '..\Datasets\Countries_Meshes\csv/RU-40527-adj.csv';
    '..\Datasets\Countries_Meshes\csv/VN-4031-adj.csv';
    };
pts = {
    '..\Datasets\Countries_Meshes\csv/CH-4468-pts.csv';
    '..\Datasets\Countries_Meshes\csv/CL-13042-pts.csv';
    '..\Datasets\Countries_Meshes\csv/GB-5946-pts.csv';
    '..\Datasets\Countries_Meshes\csv/GR-3117-pts.csv';
    '..\Datasets\Countries_Meshes\csv/NO-9935-pts.csv';
    '..\Datasets\Countries_Meshes\csv/RU-40527-pts.csv';
    '..\Datasets\Countries_Meshes\csv/VN-4031-pts.csv';
    };
flnames = {
    '..\Datasets\Countries_Meshes\mat\ch.mat'; 
    '..\Datasets\Countries_Meshes\mat\cl.mat';
    '..\Datasets\Countries_Meshes\mat\gb.mat';
    '..\Datasets\Countries_Meshes\mat\gr.mat';
    '..\Datasets\Countries_Meshes\mat\no.mat';
    '..\Datasets\Countries_Meshes\mat\ru.mat';
    '..\Datasets\Countries_Meshes\mat\vn.mat';
    };
nc = length(adjs);

for c = 1:nc
    edges = csvread(adjs{c},1);
    positions = csvread(pts{c}, 1);
    % 2. Construct the adjaceny matrix (NxN). There are multiple ways to do so
    adj = accumarray(edges, 1, [], [], [], true);
    adj = adj + transpose(adj);
    % 3. Visualize the resulting graphs

    %verify that %axis equal; is not commented on line 45 of gplotg.m
    gplotg(adj, positions)
    pause;
    % 4. Save the resulting graphs
    save(flnames{c}, "adj", "positions")
end
