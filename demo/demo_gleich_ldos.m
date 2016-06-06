% demo_gleich(dname, Nprobe, Ncheb, Nbin)
%
% Loads the graph in the 'gleich' collection (with reference eigenvalues)
% and plots histograms computed with true eigenvalue distribution adn with
% the filtered Chebyshev estimate.
%
% Input:
%   dname: Data set name
%   Nprobe: Number of probe vectors for moment estimation
%   Ncheb: Number of Chebyshev moments
%   Nbin:  Number of histogram bins
%
function demo_gleich_ldos(dname, Nprobe, Ncheb)


  if nargin < 2, Nprobe = 20;  end
  if nargin < 3, Ncheb = 1000; end
  if nargin < 4, Nbin  = 50;   end

  tic;
  A = load_graph('gleich', dname);
  A = matrix_normalize(A, 's');
  fprintf('Time to load and convert: %g\n', toc);

  tic;
  c = moments_cheb_ldos(A, Nprobe, Ncheb);
  cf = filter_jackson(c);
  fprintf('Time to compute filtered moments: %g\n', toc);

  tic;
  plot_cheb_ldos(cf);
  fprintf('Form/sort/plot local DoS: %g\n', toc);

end
