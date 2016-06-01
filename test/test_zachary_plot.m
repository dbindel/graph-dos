% Test plotting on Zachary example
A = load_zachary();
n = length(A);
L = matrix_laplacian(A);
lambdas = eig(L);

[Ls,ab] = rescale_matrix(L);
c = moments_cheb_dos(Ls,n,100,100);

xx = ab(1)*linspace(-1+1e-8,1-1e-8,1001)+ab(2);

figure(1); clf; hold on;
plot_chebint(filter_jackson(c),xx,ab);
plot(lambdas,0:length(lambdas)-1, '*');
title('Eigenvalue CDF vs KPM estimate');

figure(2); clf; hold on;
plot_cheb(filter_jackson(c),xx,ab,1);
plot(lambdas,0*lambdas,'*');
title('Eigenvalue PDF vs KPM estimate');
