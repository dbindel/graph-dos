function [c] = filter_test(N,x)
%FILTER_TEST
%

%Calculate Coeefcients 
P=zeros(N,1);
P(1) = 1;
P(2) = x;

for n = 3:N
    P(n) = 2*x*P(n-1)-P(n-2);  
end


c=P;


