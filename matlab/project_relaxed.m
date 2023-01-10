function [qx, energyMunkres, u, PsX] = project_relaxed(Q, X)

% post-processs
XFW = X(2:end, 2:end);
[u,~] = svd(XFW);
u = u(:,1);
if (sum(u>0)<length(u)/2)
    u = - u;
end

PsX = perms_q_to_cell(u, 3);
[PsX,~] = round_munkres(PsX);
qx = perms_cell_to_q(PsX);
energyMunkres = qx'*Q(2:end,2:end)*qx;
end