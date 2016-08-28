% [c] = moments_resolvent(N)
% [c] = moments_resolvent(N, alpha)
%
% Compute Chebyshev moments 0 through N-1 of the function
% f(x)=1/(1-alpha*x) on [-1,1]
%
% Input:
%    N: Number of moments
%    alpha: Resolvent coefficient in (0,1)
%    ab: Rescaling parameters. Original matrix should be ab(1)*A+ab(2)
%
% Output:
%    c: A column vector of N moments
%
function c = moments_resolvent(N,alpha,ab)
    if nargin<3
        ab = [1,0];
    end
    if nargin<2
        alpha = 0.9;
    end
    
    ab = alpha*ab;
    
    z = (1-ab(2))/ab(1);
    c = 1./(sqrt(z^2-1)*(z+sqrt(z^2-1)).^(0:N-1))/ab(1);
    c(2:N) = c(2:N)*2;
    c=c';
    
end