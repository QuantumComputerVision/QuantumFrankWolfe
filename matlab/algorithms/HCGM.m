function [X, constraintResidual] = HCGM(Q,A,ATranspose,b,beta0,maxit)

n = size(Q,1);
if size(Q,2) ~= n, error('Q must be square.'); end
if (~exist('verbose', 'var')) verbose = 1; end
X = zeros(n,n);
constraintResidual = -b;
beta = beta0;

for t = 1:maxit

    % parameters
    beta = beta0 * sqrt(t+1);
    % beta = beta * 1.04;
    eta = 2 / (t+1);

    % gradient evaluation
    grad = Q + beta .* ATranspose(constraintResidual);

    % symmetrize for numerical reasons
    % grad = 0.5.*(grad + grad');

    % linear minimization oracle
    x = solve_exhaustive_binary(grad);

    % convex combination step
    X =  X + eta .* (x*x' - X);
    constraintResidual = constraintResidual + eta * (A(x) - b - constraintResidual);

    if (verbose)
        fprintf("Constraint: %2e, Cost: %2e \n",[norm(constraintResidual), Q(:)'*X(:)]);
    end
end

end

function out = LMO(in)
[~, out] = solve_exhaustive(grad, nCams, N);
end
