% convert a a cell array of permutation matrices to a big q-vector
function [q] = perms_cell_to_q(Ps)

K = length(Ps);

q = [];
for i=1:K
    P = Ps{i};
    q = [q; P(:)];
end

end