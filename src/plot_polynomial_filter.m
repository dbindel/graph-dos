function plot_polynomial_filter(Y,N)
%PLOT_POLYNOMIAL_FILTER
%   Y is p(A)
%   V is a random/specified vector of size nx2
%   N is a cell of conneted components of A

%if V is random
n=length(Y);
V= rand(n,2);

%compute p(A)*V
PV=Y*V;

%Plot connected components using coordinates given in P
n=length(N);

%in each component

fprintf('case specific test \n');
for k=1:1
   %plot each point
   
   mems=N{k};
   m=length(mems);
   xs=zeros(1,m);
   ys=xs;
   for j=1:m
       xs(j)=PV(mems(j),1);
       ys(j)=PV(mems(j),2);
   end 
   figure();
   plot(xs,ys,'*');
end




end

