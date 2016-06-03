% test_gleich(dname, Ncheb, Nbin)
%
% Loads the graph in the 'gleich' collection (with reference eigenvalues)
% and plots histograms computed with true eigenvalue distribution adn with
% the filtered Chebyshev estimate.
%
% Input:
%   dname: Data set name
%   Ncheb: Number of Chebyshev moments
%   Nbin:  Number of histogram bins
%
function test_gleich(dname, Ncheb, Nbin)

  if nargin < 2, Ncheb = 1000; end
  if nargin < 3, Nbin  = 50;   end

  tic;
  [A,lambda] = load_graph('gleich', dname);
  lambda = 1-lambda;
  N = matrix_normalize(A, 's');
  fprintf('Time to load and convert: %g\n', toc);

  tic;
  c = moments_cheb_dos(N, 10, Ncheb);
  cf = filter_jackson(c);
  fprintf('Time to compute filtered moments: %g\n', toc);

  clf
  hold on
  lmin = max(min(lambda),-1);
  lmax = min(max(lambda), 1);
  x = linspace(lmin,lmax,Nbin+1);
  xm = (x(1:end-1)+x(2:end))/2;
  hist(lambda, Nbin)
  plot(xm, plot_chebhist(cf,x), 'r.', 'markersize', 20);
  hold off;
