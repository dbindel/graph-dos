% [c] = moments_exponential(N)
%
% Compute Chebyshev moments 0 through N-1 of exponential function
% f(x)=exp(x) on [-1,1]
%
function c = moments_exponential(N)
    c = besseli(0:N-1,1);
    c(2:N) = c(2:N)*2;
    c=c';
end