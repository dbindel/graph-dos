% [c] = moments_exponential(N)
% [c] = moments_exponential(N, beta)
%
% Compute Chebyshev moments 0 through N-1 of exponential function
% f(x)=exp(beta*x) on [-1,1]
%
% Input:
%    N: Number of moments
%    beta: Exponential coefficient
%
% Output:
%    c: A column vector of N moments
%
function c = moments_exponential(N,beta)
    if nargin<2
        beta = 1;
    end
    
    c = besseli(0:N-1,beta);
    c(2:N) = c(2:N)*2;
    c=c';
end