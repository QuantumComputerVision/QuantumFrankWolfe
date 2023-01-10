function [data] = HCGM_with_data(folder,data,A,ATranspose,b,beta0,maxit)

Q = data.Q;

n = size(Q,1);
if size(Q,2) ~= n, error('Q must be square.'); end
X = zeros(n,n);
constraintResidual = -b;
beta = beta0;

data.maxit = maxit;
data.beta0 = beta;
data.Qorg = Q;
for t = 1:maxit

    % parameters
    beta = beta0 * sqrt(t+1);
    % beta = beta * 1.04;
    eta = 2 / (t+1);

    % gradient evaluation
    grad = Q + beta .* ATranspose(constraintResidual);

    % symmetrize for numerical reasons
    grad = 0.5.*(grad + grad');

    % linear minimization oracle
    x = solve_exhaustive_binary(grad);
    data.Q = grad;
    data.xExh = x;
    data.t = t;
    data.beta = beta;
    fileName = [folder '/FW_Q_t=' num2str(t) '.mat'];
    save(fileName, 'data');

    % convex combination step
    X =  X + eta .* (x*x' - X);
    constraintResidual = constraintResidual + eta * (A(x) - b - constraintResidual);

    fprintf("Constraint: %2e, Cost: %2e \n",[norm(constraintResidual), Q(:)'*X(:)]);
end

end
