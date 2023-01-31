% use the hungariam algorithm to project onto permutations
% X has to be a non-negative cost matrix 
% this is equivalent to solving the linear program:
% \sum \sum c_{ij} x_{ij} with constraints:
% \sum_i x_{ij} = 1
% \sum_j x_{ij} = 1
% x_{ij}=0 or 1
function [P] = perm_project(X)

[P,cost] = munkres(1-X);
P = full(perm_2matrix(P));

end