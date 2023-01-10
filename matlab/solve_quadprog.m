function [Pis] = solve_quadprog(Q, nCams, nPoints)
% for now it's a chance that this works. it doesn't.
% the reason is: Q needs to be invertible. ours is not -> nonconvex.

n = size(Q,1);
s = zeros(n, 1);

[Ai, ~] = compute_constraint_system(nPoints);
A = kron(eye(nCams), Ai);
b = ones(size(A,1), 1);

x = quadprog(-2*Q, s, [], [], A, b);

% project onto permutations (uses Hungarian algorithm)
x = abs(x);
x = x./max(x(:));
Pis = perms_q_to_cell(x, nPoints);

for i=1:length(Pis)
    I = matchpairs(1-Pis{i}, 9999);
    ind = sub2ind([nPoints, nPoints], I(:,1), I(:,2));
    P = Pis{i};
    P(:) = 0;
    P(ind) = 1;
    Pis{i} = P;
end


end