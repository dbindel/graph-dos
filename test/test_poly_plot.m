%A = load_graph('zachary');
A = load_graph('ufsparse', 'Gleich/minnesota');
n = length(A);
L = matrix_laplacian(A);
%lambdas = eig(L);

%eigenvalue soecific
[Ls,ab] = rescale_matrix(L);
c = moments_delta(0.5,100);
cf = filter_jackson(c);
Y = polyomial_filter(cf,Ls);
[nComponents,sizes,members] = networkComponents(A);
%plot_polynomial_filter(Y,members)

%{
%Normal
c1 = moments_cheb_dos(Ls,n,100,100);
cf1= filter_jackson(c1);
Y1= polyomial_filter(cf1,Ls);
plot_polynomial_filter(Y1,members)
%}
