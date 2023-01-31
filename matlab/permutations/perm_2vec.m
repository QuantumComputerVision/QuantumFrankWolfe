% perm_2vec converts a permutation matrix into a vector
function [v] = perm_2vec(P)

v = P*(1:length(P))';

end