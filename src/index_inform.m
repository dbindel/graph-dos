% IC = index_inform(A, nZ, N)
% IC = index_inform(Afun, n, nZ, N)
%
% Compute the information centrality by the Chebyshev expansion of 
% 1/(ab(1)*x_ab(2)) and moments of (L+ee^T/n).
%
% Inputs:
%    A: Matrix or function to apply matrix (to multiple RHS)
%    n: Dimension of the space (if A is a function)
%    nZ: Number of probe vectors with which we want to compute moments
%    N: Number of moments to compute
%
% Output:
%    IC: Information centrality
%
function IC = index_inform(varargin)

    defaults = {'Afun', NaN, 'n', NaN, ...
              'nZ', 100, 'N', 10, 'kind', 1};
    [Afun, n, nZ, N] = mfuncify(defaults, varargin{:});
    if N < 2, N = 2; end
    
    % Compute the eigenrange of Laplacian
    Lfun = mfunc_laplacian(Afun,n);
    opts = [];
    opts.isreal = 1;
    opts.issym = 1;
    range = sort(eigs(Lfun, n, 4, 'be', opts));
    range(1) = 1;
    range = sort(range);
    range = range([1,end]);
    
    % Linear rescaling of (L+ee^T/n)
    Bfun = @(x) (bsxfun(@plus,mean(x,1),Lfun(x)));
    [Bfuns,ab] = rescale_mfunc(Bfun,n,range);
    
    % Compute the ldos moments
    [c,~] = moments_cheb_ldos(Bfuns,n,nZ,N);
    
    % Compute the Chebyshev expansion of 1/(ab(1)x+ab(2))
    z = -ab(2)/ab(1);
    w = 1./(sqrt(z^2-1)*(z-sqrt(z^2-1)).^(0:N-1))/ab(1);
    w(2:N)=2*w(2:N);
    
    % Compute information centrality
    d = c'*w'-(n-1)/n^2;
    IC = 1./(d+mean(d)-2/n^2);
    
end