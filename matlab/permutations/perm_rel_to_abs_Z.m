% convert a matrix of relative permutations to absolute permutations
% as devised in:
%
function [Ps] = perm_rel_to_abs_Z(Z, d)

[U,Ds] = eigs(Z,d,'lm');
U = real(U)*sqrt(abs(Ds));

n = size(U,1)./size(U,2);
d = size(U,2);
Ps={};
for i=1:n
    zstarti = (i-1)*d+1;
    zendi = zstarti + d - 1;
    Pi1 = U(zstarti:zendi, :)*U(1:d, 1:d)';
    %X = full(matrix2perm(Pi1));
    X = full(perm_2matrix( munkres(-Pi1)));
    Ps{i} = X;
end

end