% if second argument is avoided, it's expected to be all identity
function [precision, recall, F, acc] = get_synth_solution_PR(Psol, Pgnd)

K = length(Psol);
[M,N]=size(Psol{1});
if (~exist('Pgnd', 'var'))
    Pgnd = cell(K,1);
    for i=1:K
        Pgnd{i} = eye(N);
    end
end

TP = 0;
FP = 0;
FN = 0;
acc = 0;

for i=1:K
    
    % total number of bits guessed correctly
    compare = Psol{i}==Pgnd{i};
    acc = acc + sum(compare(:))./(size(Psol{i},1)*size(Psol{i},2));
    
    % count how many 1s we identify correctly
    compare = bitand(Psol{i}, Pgnd{i}) & Pgnd{i};
    TP = TP + sum(sum(compare));
    
    % how many 1s we identify that are not in ground truth
    compare = Psol{i}==1 & Pgnd{i}==0;
    FP = FP + sum(sum(compare));
    
    % how many 1s are in ground truth and we donot identify
    compare = Psol{i}==0 & Pgnd{i}==1;
    FN = FN + sum(sum(compare));
end

precision = (TP)./(TP+FP);
recall = TP./(TP+FN);

F = 2* (precision.*recall)./(precision+recall);

acc = acc./K;

end
