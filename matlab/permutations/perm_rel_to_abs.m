% convert a matrix of relative permutations to absolute permutations
% as devised in:
function [Ps] = perm_rel_to_abs(Pijs, I, n)

d = size(Pijs{1},1);
dimPerm = d*ones(n,1);
allN = sum(dimPerm);
Z = zeros(allN, allN);
e = size(I,2);

% compose Z
for k=1:e
    i = I(1, k); j = I(2, k);
    Pij = Pijs{k};
    zstarti = (i-1)*d+1;
    zstartj = (j-1)*d+1;
    zendi = zstarti + d - 1;
    zendj = zstartj + d - 1;
    Z(zstarti:zendi,zstartj:zendj) = Pij;
end

Ps = perm_rel_to_abs_Z(Z, I, n);

end