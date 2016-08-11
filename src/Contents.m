% Commands for manipulating graph density of states.
% Version 0.1 01-Jun-2016
%
% Loading data
%   load_koblenz     - Form weighted adjacency from a Koblenz TSV file
%   load_koblenz_raw - Raw data from a Koblenz TSV file
%
% Saving data
%   save_gexf        - Save graph+st in Gephi's GEXF format
%
% Forming associated matrices and matrix functions
%   matrix_adjacency - Strip weights from a matrix and return just adjacency
%   matrix_laplacian - Form (weighted) graph Laplacian from an adjacency matrix
%   matrix_normalize - Normalize a weighted adjacency
%   mfunc_laplacian  - Form a function to apply a weighted Laplacian
%   mfunc_modularity - Form a function to apply the modularity matrix
%   mfunc_normalize  - Form a normalized version of a matrix apply function
%
% Estimate moments of DoS
%   rescale_matrix     - Apply linear scaling to put spectrum in [-1,1] (matrix)
%   rescale_mfunc      - Like rescale_matrix, but returns a function
%   moments_cheb       - Compute Chebyshev moments for given start vector(s)
%   moments_cheb_dos   - Estimate Chebyshev moments of DoS
%   moments_cheb_ldos  - Estimate Chebyshev moments of local DoS
%   moments_cheb_ldosx - Exact Chebyshev moments of local DoS
%   moments_lanczos    - Quadrature rule via Lanczos
%   moments_lan_dos    - Moments of DoS via randomized Lanczos
%   moments_lan_ldos   - Moments of LDoS for select nodes via Lanczos
%   moments_quad2cheb  - Convert weighted deltas to Chebyshev moments
%
% Filtering Chebyshev moments
%   filter_jackson - Apply Jackson filter
%   filter_lorentz - Apply Lorentz filter
