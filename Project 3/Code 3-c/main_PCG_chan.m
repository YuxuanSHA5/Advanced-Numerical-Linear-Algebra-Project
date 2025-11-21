% main_PCG_chan.m
clear; clc;

n = 128;

%% construct Toeplitz matrix Tn
k = (0:n-1)';
c = zeros(n,1);
c(1)     = 4;
c(2:end) = -4 ./ (4*k(2:end).^2 - 1);      % c_k = -4/(4k^2-1), k>=1

Tn = toeplitz(c);


b = ones(n,1);

%% Chan'scirculant  c(Tn)
cc = chan_circulant_column(c);      
Cn = gallery('circul', cc);         

%% PCG  Tnx = b Chan circulant
tol   = 1e-8;
maxit = 1000;

Mfun = @(z) Cn \ z;                 

[x, flag, relres, iter, resvec] = pcg(@(x) Tn*x, b, tol, maxit, Mfun);

fprintf('Chan-PCG 结束：flag=%d, 迭代次数=%d, 残差范数=%.3e\n', ...
        flag, iter, norm(b - Tn*x));

%% 4. convergence curve
figure;
semilogy(0:iter, resvec, 'o-');
xlabel('Number of Iterations k');
ylabel('||r_k||_2');
title('PCG with Chan circulant preconditioner c(T_n)');
grid on;
