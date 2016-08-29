% IC = index_inform(A, nZ, N)
%
% Compute the information centrality by the Chebyshev expansion of 
% 1/(ab(1)*x_ab(2)) and moments of (L+ee^T/n).
%
%    A: Adjacency matrix
%    nZ: Number of probe vectors with which we want to compute moments
%    N: Number of moments to compute
%
% Output:
%    IC: Information centrality
%
function IC = index_inform(A, nZ, N)
    if nargin<3
        N = 150;
    end
    
    if nargin<2
        nZ=100;
    end
    
    n = size(A,1);
   
    % Compute the eigenrange of Laplacian
    L = matrix_laplacian(A);
    opts = [];
    opts.isreal = 1;
    opts.issym = 1;
    range = sort(eigs(L, 2, 'sm', opts));
    range(3) = eigs(L, 1, 'lm', opts);
    range(1) = 1;
    range = sort(range);
    range = range([1,end]);
    
    % Linear rescaling of (L+ee^T/n)
    Lfun = @(X) L*X;
    Bfun = @(X) (bsxfun(@plus,mean(X,1),Lfun(X)));
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