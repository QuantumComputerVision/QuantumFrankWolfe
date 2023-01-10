% computes the big Q matrix given a cell array of relative permutations and
% the corresponding edges
function [Q] = compute_Q(Pijs, I, numCams)

K = numCams;
numEdges = length(I);

% compute the size of each permutation
dimP = zeros(K,1);
dimPKron = zeros(K,1);
for i=1:K
    dimP(i) = size(Pijs{i},1);
    dimPKron(i) = dimP(i)*dimP(i);
end

% retrieve the indices for each block (pairwise permutation)
m = [0;cumsum(dimPKron(1:end-1))];
blk = @(k) 1+m(k):m(k)+dimPKron(k);

N = dimPKron(end) + m(end);
Q = sparse(N,N);
Id = eye(dimP(1));
for i=1:numEdges
    edge = I(:,i);
    ei = edge(1); ej = edge(2);
    Q(blk(ei),blk(ej)) = kron(Id, Pijs{i});
end

% also set the diagonal to identity
for i=1:K
    Q(blk(i),blk(i)) = kron(Id, Id);
end

end