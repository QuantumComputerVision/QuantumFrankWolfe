% exhaustively solves the qubo problem
function [Pis, qResult] = solve_exhaustive(Q, nCams, nPoints)

permIndicesP = uint16(perms(1:nPoints));

numPerms = size(permIndicesP,1);
%allPossibleSolutions = nchoosek(1:numPerms, nCams);
allPossibleSolutions =nmultichoosek(1:numPerms, nCams);

numSol = length(allPossibleSolutions);

minErr = 999999999;
solInd = 1;
n2 = nPoints*nPoints;
sizeP = [nPoints,nPoints];
q = zeros(nPoints*nPoints*nCams,1);
indBasis = 1:nPoints;
for i=1:numSol
    curSolInd = allPossibleSolutions(i, :);
    selectedPerms = permIndicesP(curSolInd,:);
    q(:) = 0;
    for j=1:size(selectedPerms, 1)
        q((j-1)*n2 + sub2ind(sizeP, indBasis, selectedPerms(j,:))) = 1;
    end
    
    err = -q'*Q*q;
    if (err<minErr)
        qResult = q;
        minErr = err;
        solInd = i;
    end
    
end

Pis = perms_q_to_cell(qResult, nPoints);
Pis = perms_transform_first_to_eye(Pis);

end