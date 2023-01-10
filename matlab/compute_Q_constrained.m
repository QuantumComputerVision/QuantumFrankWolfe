% computes the big Q matrix given a cell array of relative permutations and
% the corresponding edges
function [Q, s] = compute_Q_constrained(Pijs, I, numCams, lambda)

Q = compute_Q(Pijs, I, numCams);

n = size(Pijs{1},1);

[Ai, ~] = compute_constraint_system(n);

A = kron(eye(numCams), Ai);
b = ones(size(A,1), 1);

% these are for the argmax not argmin, so put a minus if you like to
% minimize.
Q = Q - lambda*(A'*A);
s = 2*lambda*A'*b;

end