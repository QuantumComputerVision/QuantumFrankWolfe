% if second argument is avoided, it's expected to be all identity
function [acc] = get_synth_solution_acc(Psol, Pgnd)

K = length(Psol);
[M,N]=size(Psol{1});
if (~exist('Pgnd', 'var'))
    Pgnd = cell(K,1);
    for i=1:K
        Pgnd{i} = eye(N);
    end
end

numTotalElements = K*N;

acc = 0;
for i=1:K
    compare = bitand(Psol{i}, Pgnd{i});
    acc = acc + sum(sum(compare));
end

acc = acc/(numTotalElements);

end