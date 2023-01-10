clear all;
rng(0,'twister');

addpath('./permutations/');
addpath('./dataset/');
addpath('./helpers/');
addpath('./algorithms/');

close all;

N = 3; % the number of the grid points. permutations are NxN
nCams = 3; % number of cameras (# nodes)
completeness = 1.0; % 1 is fully connected, 0 has no connectivity
swapRatio = 0.0; % like to add some noise? 0.01-0.15
lambda = 1.25; % regularizer (permutation constratins)

% set up the synth. problem
% Pijs are noisy relative permutations
% PijsGnd are non-noisy (true) relative permutations
% Pis are the solution whereas Pijs are the data.
[Pis, Pijs, PijsGnd, G, I, Iupper] = gen_synth_multi_graph_2d_no_grid(N, nCams, completeness, swapRatio);
% [Pis, Pijs, PijsGnd, G, I, Iupper, points] = gen_synth_multi_graph_2d(N, nCams, completeness, swapRatio);

% compute the unconstrained Q matrix from ground truth relative permutations
% without permutation constrints (unconstrained)
% QGnd induces the quadratic energy: -q'*QGnd*q where qs are vectorized.
QGnd = compute_Q(PijsGnd, I, nCams);

% verify that QGnd is maximized at q=vec(pis)
q = perms_cell_to_q(Pis);
-q'*QGnd*q

% noisy Q should lead to a higher cost
Q = compute_Q(Pijs, I, nCams);
-q'*Q*q

% compute the exhaustive solution - this should align with the gnd-truth
[PisExhaustive, qResultExhaustive] = solve_exhaustive(Q, nCams, N);
-qResultExhaustive'*Q*qResultExhaustive
XResultExhasutive = qResultExhaustive*qResultExhaustive';

% we can also try a continuous relaxation and use matlab's quadprog
% NOTE: Below call doesn't work at all. (good for us)
[PisQP] = solve_quadprog(Q, nCams, N);

% now let's play with the constrained version
[Ai, ~] = compute_constraint_system(N);
A = kron(eye(nCams), Ai);
b = ones(size(A,1), 1);
% let's verify the constraints:
disp(['Cosntraint debugging: sum(A*q-b) = ', num2str(sum(A*q-b))]);

%bsquared = b.^2; % this shouldn't matter if b is always 0 or 1
maxit = 200;
beta0 = 1;
verbose = false;
tic;
[Pis, XFW, cost, constraintResidual] = CoCoFW(Q, N, A, b, 'FWAL', maxit, beta0, verbose, 0, 1);
timeCoco = toc;

qsol = perms_cell_to_q(Pis);    
costCoco = -qsol'*Q*qsol;
costExh = -qResultExhaustive'*Q*qResultExhaustive;
cocoObj = XFW(:)'*Q(:);
disp(['cost CoCo: ' num2str(costCoco) ', CoCo-obj: ' num2str(cocoObj) ', cost exh: ' num2str(costExh)]);

if swapRatio == 0
save(['test_synt_data_sample-results-N=',num2str(N),'-nCams=',num2str(nCams),'-noiseless']);
else
save(['test_synt_data_sample-results-N=',num2str(N),'-nCams=',num2str(nCams),'-noisy']);
end

return;

[Qcons, s] = compute_Q_constrained(Pijs, I, nCams, lambda);
% [Qcons, s] = compute_Q_constrained(PijsGnd, I, nCams, lambda);
q = perms_cell_to_q(Pis);
-q'*Qcons*q - s'*q

[PisConsExh, qConsExh] = solve_exhaustive_constrained(Qcons, s, nCams, N);
-qConsExh'*Qcons*qConsExh - s'*qConsExh

figure, spy(QGnd);
figure, spy(Q);


%%%% IF YOU WANT TO SAVE %%%%
data.QGnd = QGnd;
data.Q = Q;
data.Qcons = sparse(Qcons);
data.s = s;
data.PijsGnd = PijsGnd;
data.Pijs = Pijs;
data.Pis = Pis;
data.N = N; % # of points
data.G = G; % # of points
data.I = I; % # of points
data.Iupper = Iupper; % # of points
data.nCams = nCams;
data.lambda = lambda;
data.swapRatio = swapRatio;
data.completeness = completeness;

[data.QconsRows, data.QconsCols, data.QconsValues] = find(Qcons);
[data.QRows, data.QCols, data.QValues] = find(Q);
[data.QGndRows, data.QGndCols, data.QGndValues] = find(QGnd);


fileName = './dataset/synth2d_dataN3NC3';
fileNameQ = './dataset/QMatrixN3NC3';
fileNameAllData = './dataset/synth_data_outputN3NC3.mat';
saveData = 0;

if (saveData)
    save(fileNameAllData, 'data');
end

