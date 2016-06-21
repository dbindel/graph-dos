function [Y] = polyomial_filter(c,A)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


N = length(c);
P0 = eye(length(A));
P1 = A;
Y = (c(1)/2)*P0 + c(2)*A;
for np = 3:N
    Pn = 2*A*P1 - P0;
    Y = Y + c(np)*Pn;
    P0 = P1;
    P1 = Pn;
end

end

