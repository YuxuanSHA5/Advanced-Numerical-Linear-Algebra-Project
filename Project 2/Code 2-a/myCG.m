function [x, flag, iter, resvec] = myCG(A, b, x, tol, maxit)

if nargin < 3 || isempty(x)
    x = zeros(size(b));
end
if nargin < 4 || isempty(tol)
    tol = 1e-8;
end
if nargin < 5 || isempty(maxit)
    maxit = length(b);
end

r = b - A*x;            
p = r;                  
rho = r' * r;

res0 = sqrt(rho);
resvec = zeros(maxit+1,1);
resvec(1) = res0;

if res0 < tol
    flag = 0;
    iter = 0;
    resvec = resvec(1);
    return;
end

for k = 1:maxit
    Ap = A * p;
    alpha = rho / (p' * Ap);

    x = x + alpha * p;
    r = r - alpha * Ap;

    rho_new = r' * r;
    resvec(k+1) = sqrt(rho_new);

    if resvec(k+1) < tol
        flag = 0;
        iter = k;
        resvec = resvec(1:k+1);
        return;
    end

    beta = rho_new / rho;
    p = r + beta * p;

    rho = rho_new;
end


flag = 1;          
iter = maxit;
resvec = resvec(1:maxit+1);
end
