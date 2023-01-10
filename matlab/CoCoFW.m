function [Pis, XFW, cost, constraintResidual] = CoCoFW(Q, N, A, b, alg, maxit, beta0, verbose, dWaveMode, isSave,isVisualize)

m = size(Q,1);
mb = size(b,1);

if (~exist('alg','var'))            alg = 'FWAL';       end
if (~exist('maxit','var'))          maxit = 10000;      end
if (~exist('beta0','var'))          beta0 = 1;          end
if (~exist('b','var'))              b = ones(m, 1);     end
if (~exist('verbose', 'var'))       verbose = 1;        end
if (~exist('dWaveMode', 'var'))     dWaveMode = 1;      end
if (~exist('isSave', 'var'))        isSave = 0;         end
if (~exist('isVisualize', 'var'))   isVisualize = 0;    end

% % frank-wolfe (Alp's edits)
% xVec = @(W) 0.5 .* ( W(1,2:end)' + W(2:end,1) );
% xDiag = @(W) diag(W(2:end,2:end));


Aoperator0 = @(x) x(1)^2;
ATranspose0 = @(y) [y, zeros(1,m); zeros(m,1), zeros(m,m)];
b0 = 1;

Aoperator1 = @(x) (A*x(2:end)).^2;
ATranspose1 = @(y)  [0, zeros(1,m); zeros(m,1), A' * bsxfun(@times, A, y)];
b1 = b.^2; % this shouldn't matter if b is always 0 or 1

% m = size(q,1);
Aoperator2 = @(x) x(2:end).^2 - x(1).*x(2:end);
% ATranspose2 = @(y) diag(y) - [y';zeros(m-1,m)] - [y,zeros(m,m-1)];
ATranspose2 = @(y) diag([0;y]) -  [0,y';y,zeros(m,m)]; 
b2 = zeros(m, 1);


% m = size(q,1);
Aoperator3 = @(x) A*(x(1)*x(2:end));
% ATranspose2 = @(y) diag(y) - [y';zeros(m-1,m)] - [y,zeros(m,m-1)];
ATranspose3 = @(y) [0,y'*A;A'*y,zeros(m,m)]; 
b3 = b;


% Aoperator = @(x) [Aoperator1(x);Aoperator2(x)];
% ATranspose = @(y) ATranspose1(y(1:mb)) + ATranspose2(y(mb+1:end));
% bb = [b1;b2];
% Aoperator = @(x) [Aoperator0(x);Aoperator1(x);Aoperator2(x)];
% ATranspose = @(y) ATranspose0(y(1)) + ATranspose1(y(2:mb+1)) + ATranspose2(y(mb+2:end));
% bb = [b0;b1;b2];

Aoperator = @(x) [Aoperator0(x);Aoperator1(x);Aoperator2(x);Aoperator3(x)];
ATranspose = @(y) ATranspose0(y(1)) + ATranspose1(y(2:mb+1)) + ATranspose2(y(mb+2:mb+1+m)) + ATranspose3(y(mb+2+m:end));
bb = [b0;b1;b2;b3];


% symmetrize for numerical reasons
Q = 0.5.*(Q + Q');
Qc = [0,zeros(1,m);zeros(m,1),Q];


if (strcmp(alg, 'FWAL')==1 || strcmp(alg, 'fwal')==1)
    [XFW, constraintResidual] = FWAL(-Qc,Aoperator,ATranspose,bb,beta0,maxit,verbose,0,isSave,isVisualize);
elseif (strcmp(alg, 'dwaveFWAL')==1 || strcmp(alg, 'dwavefwal')==1)
    [XFW, constraintResidual] = FWAL(-Qc,Aoperator,ATranspose,bb,beta0,maxit,verbose,dWaveMode,isSave,isVisualize);
elseif (strcmp(alg, 'HCGM')==1 || strcmp(alg, 'hcgm')==1)
    [XFW, constraintResidual] = HCGM(-Qc,Aoperator,ATranspose,bb,beta0,maxit,verbose,dWaveMode,isSave,isVisualize);
elseif (strcmp(alg, 'dwaveHCGM')==1 || strcmp(alg, 'dwavehcgm')==1)
    [XFW, constraintResidual] = HCGM_with_data(-Qc,Aoperator,ATranspose,bb,beta0,maxit,verbose,dWaveMode,isSave,isVisualize);
end

% post-processs
XFW = XFW(2:end, 2:end);
[u,~] = svd(XFW);
u = u(:,1);
if (sum(u>0)<length(u)/2)
    u = - u;
end

% re-arrange into solutions
% costMat = u*u';
Pis = perms_q_to_cell(u, N);
costs = zeros(length(Pis),1);
for i=1:length(Pis)
    costMat = 1 - Pis{i};
%     assignment = Pis{i};
%     thr = 0.5*(mean(assignment(:))+min(assignment(:)));
%     assignment(assignment<thr) = 0;
%     assignment(assignment>=thr) = 1;
    [assignment, cost] = munkres(costMat);
    Pis{i} = assignment;
    costs(i) = cost;
end

end