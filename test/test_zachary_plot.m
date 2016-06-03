% Test plotting on Zachary example
A = load_graph('zachary');
n = length(A);
L = matrix_laplacian(A);
lambdas = eig(L);

[Ls,ab] = rescale_matrix(L);
c = moments_cheb_dos(Ls,n,100,100);
cf = filter_jackson(c);

figure(1); clf; hold on;
plot_chebint(cf,ab);
plot(lambdas,0:length(lambdas)-1, '*');
title('Eigenvalue CDF vs KPM estimate');

figure(2); clf; hold on;
plot_cheb(cf,ab);
plot(lambdas,0*lambdas,'*');
title('Eigenvalue PDF vs KPM estimate');

figure(3); clf; hold on
xx = ab(1)*linspace(-1+1e-8,1-1e-8,11)+ab(2);
xm = (xx(1:end-1)+xx(2:end))/2;
nlam  = hist(lambdas,xm);
nkpm  = plot_chebhist(c,xx,ab);
nkpmj = plot_chebhist(cf,xx,ab);
bar(xm', [nlam; nkpm; nkpmj]');
legend('True', 'KPM', 'KPMJ');
title('Eigenvalue histogram vs KPM variants');
