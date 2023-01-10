function [Pis, XFW, cost, constraintResidual] = QFW(Q, N, alg, maxit, beta0)

if (~exist('alg','var'))
    alg = 'CGAL';
end

if (~exist('maxit','var'))
    maxit = 10000;
end

if (~exist('beta0','var'))
    beta0 = 1;
end

b = ones(size(Q,1), 1);

% % frank-wolfe (Alp's edits)
Aoperator1 = @(x) (A*x).^2;
ATranspose1 = @(y)  A' * bsxfun(@times, A, y);
b1 = b.^2; % this shouldn't matter if b is always 0 or 1

m = size(q,1);
Aoperator2 = @(x) x.^2 - x(1).*x;
ATranspose2 = @(y) diag(y) - [y';zeros(m-1,m)] - [y,zeros(m,m-1)];
b2 = zeros(size(q));

Aoperator = @(x) [Aoperator1(x);Aoperator2(x)];
ATranspose = @(y) ATranspose1(y(1:m)) + ATranspose2(y(m+1:end));
bb = [b1;b2];

if (strcmp(alg, 'CGAL')==0 || strcmp(alg, 'cgal')==0)
    [XFW, constraintResidual] = CGAL(-Q,Aoperator,ATranspose,bb,beta0,maxit);
else
    [XFW, constraintResidual] = HCGM(Q,A,ATranspose,b,beta0,maxit);
end

% post-processs
[u,~] = svd(XFW);
u = u(:,1);
u = - u;

% re-arrange into solutions
% costMat = u*u';
Pis = perms_q_to_cell(u, N);
costs = zeros(n,1);
for i=1:n
    costMat = Pis{i};
    [assignment, cost] = munkres(costMat);
    Pis{i} = assignment;
    costs(i) = cost;
end

end