function [X] = CGAL(Q,A,ATranspose,b,beta0,maxit)
%HCGM2 Summary of this function goes here
%   Detailed explanation goes here

n = size(Q,1);
if size(Q,2) ~= n, error('Q must be square.'); end
X = zeros(n,n);
y = zeros(size(b));
constraintResidual = -b;
beta = beta0;

for t = 1:maxit
    
    % parameters
    beta = beta0 * sqrt(t+1); 
%     beta = beta * 1.01;
    eta = 2 / (t+1);
    
    % gradient evaluation
    grad = Q + ATranspose(y + beta .* constraintResidual);
    
    % linear minimization oracle
    x = solve_exhaustive_binary(grad);

    % convex combination step
    X =  X + eta .* (x*x' - X);
    constraintResidual = constraintResidual + eta * (A(x) - b - constraintResidual);
    
    % dual update
    y = y + beta0 * constraintResidual;

    fprintf("Constraint: %2e, Cost: %2e \n",[norm(constraintResidual), Q(:)'*X(:)]);
% 
%     if 2^floor(log2(t)) == t
%         [uu,zz] = eig(X);
%         zz = diag(zz);
%         [~,ind] = max(abs(zz));
%         
%         AAA = perms_q_to_cell(uu(:,ind),2);
%         AAA{:}
%     end

end


end

function out = LMO(in)
    [~, out] = solve_exhaustive(grad, nCams, N);
end
