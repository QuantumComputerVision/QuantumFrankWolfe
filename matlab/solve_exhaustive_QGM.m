% exhaustively solves the qubo problem
function [Pis, qResult] = solve_exhaustive_QGM(Q, nPoints)

permIndicesP = uint16(perms(1:nPoints));
numPerms = size(permIndicesP,1);

minErr = 999999999;
for i=1:numPerms
    selectedPerms = permIndicesP(i,:);
    qM = zeros(nPoints, nPoints);
    for j=1:nPoints
        qM(j, selectedPerms(j)) = 1;
    end

    q = qM(:);
    
    err = q'*Q*q;
    if (err<minErr)
        qResult = q;
        minErr = err;
    end
    
end

Pis = perms_q_to_cell(qResult, nPoints);

end