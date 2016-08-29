clear
clc

addpath([fileparts(pwd),'/src']);
addpath([fileparts(pwd),'/data/gleich']);

pass = 'Test passed.';
fail = 'Test failed.';
lb = '------------------------------------------------------------\n';

% test moments_exponential, which computes Chebyshev coefficient of
% exp(beta*x) bewteen [-1,1]

beta = 0.9;
N = 20;
fprintf(['Test momemts_expoential with beta = ',num2str(beta),...
        ' and N = ',num2str(N),'.\n']);
w = moments_exponential(N,beta);
xx = linspace(-1+1e-8,1-1e-8,1001);
yy1 = plot_chebp(w,xx);
yy2 = exp(beta*xx);
relerr = norm(yy1-yy2)/norm(yy2);
if relerr < 1e-12
    fprintf(['The relative error is ',num2str(relerr),'. ',pass,'\n']);
else
    fprintf(['The relative error is ',num2str(relerr),'. ',fail,'\n']);
end
fprintf(lb);

% test moments_resolvent, which computes Chebyshev coefficient of 
% 1/(1-alpha*x) on [-1,1].

alpha = 0.9;
N = 70;
fprintf(['Test momemts_resolvent with alpha = ',num2str(alpha),...
        ' and N = ',num2str(N),'.\n']);
w = moments_resolvent(N,alpha);
xx = linspace(-1+1e-8,1-1e-8,1001);
yy1 = plot_chebp(w,xx);
yy2 = 1./(1-alpha*xx);
relerr = norm(yy1-yy2)/norm(yy2);
if relerr < 1e-12
    fprintf(['The relative error is ',num2str(relerr),'. ',pass,'\n']);
else
    fprintf(['The relative error is ',num2str(relerr),'. ',fail,'\n']);
end
fprintf(lb);

% Computing eigenvalue of Erdos02-cc for test purposes
fprintf(['Loading Erdos02-cc and computing its eigen-decomposition \n'...
        'for testing purposes. Takes 50s.\n']);

A = load_graph('Erdos02-cc');
L = matrix_normalize(A);
[VA,eA] = eig(full(A),'vector');
[VL,eL] = eig(full(L),'vector');
fprintf(lb)

% Test on index_estrada, which computes the Estrada index of matrix
fprintf('Test index_estrada on normalized adjacency. No rescaling.\n');
EE2 = sum(exp(eL));
N = 10;
for nZ =10:10:50
    fprintf(['With nZ = ',num2str(nZ),' and N = ',num2str(N),':\n']);
    EE1 = index_estrada(L,nZ,N);
    relerr = abs(EE1-EE2)/EE2;
    if relerr < 0.01
        fprintf(['The relative error is ',num2str(relerr),'. '...
                pass,'\n']);
    else
        fprintf(['The relative error is ',num2str(relerr),'. '...
                fail,'\n']);
    end
end
fprintf(lb)

fprintf(['Test index_estrada on normalized adjacency with '... 
    'precomputed dos. No rescaling.\n']);
N = 10;
for nZ =10:10:50
    fprintf(['With nZ = ',num2str(nZ),' and N = ',num2str(N),':\n']);
    c = moments_cheb_dos(L,nZ,N);
    EE1 = index_estrada(c);
    relerr = abs(EE1-EE2)/EE2;
    if relerr < 0.01
        fprintf(['The relative error is ',num2str(relerr),'. '...
                pass,'\n']);
    else
        fprintf(['The relative error is ',num2str(relerr),'. '...
                fail,'\n']);
    end
end
fprintf(lb)

fprintf('Test index_estrada on adjacency. With rescaling.\n');
EE2 = sum(exp(eA));
N = 30;
[As,ab] = rescale_matrix(A);
for nZ =200:200:1000
    fprintf(['With nZ = ',num2str(nZ),' and N = ',num2str(N),':\n']);
    EE1 = index_estrada(As,nZ,N,ab);
    relerr = abs(EE1-EE2)/EE2;
    if relerr < 0.05
        fprintf(['The relative error is ',num2str(relerr),'. '...
                pass,'\n']);
    else
        fprintf(['The relative error is ',num2str(relerr),'. '...
                fail,'\n']);
    end
end
fprintf(lb)

% Test on index_sub_exp, which computes the exponential subgraph
% centrality

fprintf(['Test index_sub_exp on normalized adjacency. No rescaling.'...
        'beta = 0.9\n']);
beta = 0.9;
ESC2 = VL.^2*exp(beta*eL);
N = 10;
for nZ =200:200:1000
    fprintf(['With nZ = ',num2str(nZ),' and N = ',num2str(N),':\n']);
    ESC1 = index_sub_exp(beta,L,nZ,N);
    relerr = norm(ESC1-ESC2)/norm(ESC2);
    if relerr < 0.05
        fprintf(['The relative error is ',num2str(relerr),'. '...
                pass,'\n']);
    else
        fprintf(['The relative error is ',num2str(relerr),'. '...
                fail,'\n']);
    end
end
fprintf(lb)

fprintf(['Test index_sub_exp on normalized adjacency with '... 
    'precomputed ldos. No rescaling. beta = 0.9\n']);
beta = 0.9;
N = 10;
for nZ =200:200:1000
    fprintf(['With nZ = ',num2str(nZ),' and N = ',num2str(N),':\n']);
    c = moments_cheb_ldos(L,nZ,N);
    ESC1 = index_sub_exp(beta,c);
    relerr = norm(ESC1-ESC2)/norm(ESC2);
    if relerr < 0.05
        fprintf(['The relative error is ',num2str(relerr),'. '...
                pass,'\n']);
    else
        fprintf(['The relative error is ',num2str(relerr),'. '...
                fail,'\n']);
    end
end
fprintf(lb)

fprintf('Test index_sub_exp on adjacency. With rescaling. beta = 0.9\n');
beta = 0.9;
ESC2 = VA.^2*exp(beta*eA);
N = 20;
[As,ab] = rescale_matrix(A);
for nZ =400:400:2000
    fprintf(['With nZ = ',num2str(nZ),' and N = ',num2str(N),':\n']);
    ESC1 = index_sub_exp(beta,As,nZ,N,ab);
    relerr = norm(ESC1-ESC2)/norm(ESC2);
    if relerr < 0.05
        fprintf(['The relative error is ',num2str(relerr),'. '...
                pass,'\n']);
    else
        fprintf(['The relative error is ',num2str(relerr),'. '...
                fail,'\n']);
    end
end
fprintf(lb)

% Test on index_sub_res, which computes the resolvent subgraph
% centrality

fprintf(['Test index_sub_res on normalized adjacency. No rescaling.'...
        ' alpha = 0.9\n']);
alpha = 0.9;
RSC2 = diag((eye(5534)-alpha*L)^(-1));
N = 10;
for nZ =300:300:1500
    fprintf(['With nZ = ',num2str(nZ),' and N = ',num2str(N),':\n']);
    RSC1 = index_sub_res(alpha,L,nZ,N);
    relerr = norm(RSC1-RSC2)/norm(RSC2);
    if relerr < 0.05
        fprintf(['The relative error is ',num2str(relerr),'. '...
                pass,'\n']);
    else
        fprintf(['The relative error is ',num2str(relerr),'. '...
                fail,'\n']);
    end
end
fprintf(lb)

fprintf(['Test index_sub_res on normalized adjacency with '... 
    'precomputed ldos. No rescaling. alpha = 0.9\n']);
alpha = 0.9;
N = 10;
for nZ =300:300:1500
    fprintf(['With nZ = ',num2str(nZ),' and N = ',num2str(N),':\n']);
    c = moments_cheb_ldos(L,nZ,N);
    RSC1 = index_sub_res(alpha,c);
    relerr = norm(RSC1-RSC2)/norm(RSC2);
    if relerr < 0.05
        fprintf(['The relative error is ',num2str(relerr),'. '...
                pass,'\n']);
    else
        fprintf(['The relative error is ',num2str(relerr),'. '...
                fail,'\n']);
    end
end
fprintf(lb)

[As,ab] = rescale_matrix(A);
alpha = 0.9/sum(ab);
fprintf(['Test index_sub_res on adjacency. With rescaling. alpha = ',...
    num2str(alpha),'\n']);
RSC2 = diag((eye(5534)-alpha*A)^(-1));
N = 20;
for nZ =300:300:1500
    fprintf(['With nZ = ',num2str(nZ),' and N = ',num2str(N),':\n']);
    RSC1 = index_sub_res(alpha,As,nZ,N,ab);
    relerr = norm(RSC1-RSC2)/norm(RSC2);
    if relerr < 0.05
        fprintf(['The relative error is ',num2str(relerr),'. '...
                pass,'\n']);
    else
        fprintf(['The relative error is ',num2str(relerr),'. '...
                fail,'\n']);
    end
end
fprintf(lb)

% Test on index_inform, which computes the information centrality

fprintf('Test index_inform on adjacency.\n');
B = matrix_laplacian(A)+1;
K = B^(-1);
I = repmat(diag(K),[1,5534])+repmat(diag(K)',[5534,1])-2*K; 
IC2 = 1./(mean(I,2));
N = 150;
for nZ =100:100:500
    fprintf(['With nZ = ',num2str(nZ),' and N = ',num2str(N),':\n']);
    IC1 = index_inform(A,nZ,N);
    relerr = norm(IC1-IC2)/norm(IC2);
    if relerr < 0.05
        fprintf(['The relative error is ',num2str(relerr),'. '...
                pass,'\n']);
    else
        fprintf(['The relative error is ',num2str(relerr),'. '...
                fail,'\n']);
    end
end
fprintf(lb)
fprintf('Test Over\n')