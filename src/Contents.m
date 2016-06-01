% Commands for manipulating graph density of states.
% Version 0.1 01-Jun-2016
%
% Loading data
%   load_koblenz     - Form weighted adjacency from a Koblenz TSV file
%   load_koblenz_raw - Raw data from a Koblenz TSV file
%
% Forming associated matrices and matrix functions
%   matrix_adjacency - Strip weights from a matrix and return just adjacency
%   matrix_laplacian - Form (weighted) graph Laplacian from an adjacency matrix
%   matrix_normalize - Normalize a weighted adjacency
%   mfunc_laplacian  - Form a function to apply a weighted Laplacian
%   mfunc_modularity - Form a function to apply the modularity matrix
%   mfunc_normalize  - Form a normalized version of a matrix apply function
