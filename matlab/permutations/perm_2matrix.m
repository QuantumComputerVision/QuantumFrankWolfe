% V2P converts a vector representing a permutation into a permutation matrix
% n is the size of the permutation matrix
function [P] = perm_2matrix(v, n)

m = length(v);
if nargin < 2
    n = m;
end
i = (1:m);
i(v==0) = [];
v(v==0) = [];
P = full(sparse(i,v,1,m,n));

% double stochastic
%assert(all(sum(P,1)<=1) && all(sum(P,2)<=1))

end