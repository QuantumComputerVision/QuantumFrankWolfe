% generates an NxN random permutation matrix
function [P] = perm_rand(N)

P = eye(N);
P = P(randperm(N),:);

end