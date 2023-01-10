% exhaustively solves the qubo problem
function [qResult, minErr] = solve_exhaustive_binary(Q)

n = size(Q,2);
N = 2^n;

de2bif = @(x)  2.^[(floor(log2(x)):-1:1),0];

q = zeros(n, 1);

minErr = inf;
qResult = nan;
for i = 0:(N-1)
    % q = str2num(dec2bin(i,n)'); %#ok
    de2bix = de2bif(i);
    qi = rem(i,2*de2bix)>(de2bix-1);

    q(:) = 0;
    % [~,l]=log2(i); % length of qi
    q(1:length(qi)) = qi;
    err = q'*Q*q;
    if (err<minErr)
        qResult = q;
        minErr = err;
    end
    
end

end
