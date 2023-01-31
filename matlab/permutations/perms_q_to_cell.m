% convert a big q-vector to a cell array of permutation matrices
function [Ps] = perms_q_to_cell(q, N)

numAll = length(q);
N2 = (N*N);
K = numAll/N2;

Ps = cell(K, 1);

for i=1:K
    qi = q((i-1)*N2+1:i*N2);
    Pi = reshape(qi,N,N);
    Ps{i} = Pi;
end

end