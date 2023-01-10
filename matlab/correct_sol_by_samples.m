function [Pis] = correct_sol_by_samples(PisCells, numSolCheck)

Pis = cell(size(PisCells,1),1);
M = size(PisCells{1,1},1);
I = eye(M);
for i=1:size(PisCells,1)
    PisCur = PisCells(i, :);
    PisSol = PisCells{i, 1};
    mask = PisSol~=I;
    
    totalP = min(size(PisCells,2), numSolCheck);
    M = size(PisCur{1},1);
    N = size(PisCur{1},2);
    P = zeros(M,N, totalP);
    
    for j=1:totalP
        P(:,:,j) = PisCur{j};
    end
    Pmode = mode(P,3);
    PisI = PisSol;
    PisI(mask) = Pmode(mask);
    Pis{i} = PisI;
end

end