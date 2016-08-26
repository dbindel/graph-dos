% [c] = moments_resolvent(N)
% [c] = moments_resolvent(N, alpha)
%
% Compute Chebyshev moments 0 through N-1 of the function
% f(x)=1/(1-alpha*x) on [-1,1]
%
% Input:
%    N: Number of moments
%    alpha: Resolvent coefficient in (0,1)
%
% Output:
%    c: A column vector of N moments
%
function c = moments_resolvent(N,alpha)
    if nargin<2
        alpha = 0.9;
    end
    
    z = 1/alpha;
    c = z./(sqrt(z^2-1)*(z+sqrt(z^2-1)).^(0:N-1));
    c(2:N) = c(2:N)*2;
    c=c';
    
end