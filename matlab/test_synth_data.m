
addpath('./permutations/');
addpath('./dataset/');
addpath('./helpers/');

close all;

fileName = './dataset/synth2d_dataN3NC3';
fileNameQ = './dataset/QMatrixN3NC3';
fileNameAllData = './dataset/synth_data_outputN3NC3.mat';
saveData = 0;

N = 4; % the number of the grid points. permutations are NxN
nCams = 4;
completeness = 1.0; % 1 is fully connected, 0 has no connectivity
swapRatio = 0; % like to add some noise?
lambda = 1.25;

t = tiledlayout(1,3);
t.Padding = 'none';
t.TileSpacing = 'compact';
h1=nexttile; beautify_plot; axis off;
h2=nexttile; beautify_plot; axis off;
h3=nexttile; beautify_plot; axis off;
set(gcf, 'currentaxes', h1);
% [Pis, Pijs, PijsGnd, G, I, Iupper, points] = gen_synth_multi_graph_2d(N, nCams, completeness, swapRatio);
[Pis, Pijs, PijsGnd, G, I, Iupper] = gen_synth_multi_graph_2d_no_grid(N, nCams, completeness, swapRatio);

if (saveData)
    save(fileName);
end

colors = distinguishable_colors(length(I));

visualize_synth_match_data_all(Pis, Pijs, I);
th=title('(a) Ground truth', 'FontSize', 12);
titlePos = get( th , 'position');
titlePos(2) = 8.5;
set( th , 'position' , titlePos);

set(gcf, 'currentaxes', h2);
swapRatio = 0.115;
[Pis, Pijs, PijsGnd, G, I, Iupper] = gen_synth_multi_graph_2d_no_grid(N, nCams, completeness, swapRatio);
visualize_synth_match_data_all(Pis, Pijs, I);
th=title('(b) \sigma=0.1', 'FontSize', 12);
titlePos = get( th , 'position');
titlePos(2) = 8.5;
set( th , 'position' , titlePos);

set(gcf, 'currentaxes', h3);
swapRatio = 0.2;
[Pis, Pijs, PijsGnd, G, I, Iupper] = gen_synth_multi_graph_2d_no_grid(N, nCams, completeness, swapRatio);
visualize_synth_match_data_all(Pis, Pijs, I);
th=title('(c) \sigma=0.2', 'FontSize', 12);
titlePos = get( th , 'position');
titlePos(2) = 8.5;
set( th , 'position' , titlePos);
return;
% get the position of the title


%visualize_synth_match_data(Pis, Pijs, points, I, [N N]);
%beautify_plot();
%title('visualizing just some of the pairwise matrices - the graph is fully connected');
return;
% check how much the noise corruped the data:
% note that if 10% swaps are introduced 20% is expected to be wrong
[acc] = get_synth_solution_acc(Pijs)

% compute the unconstrained Q matrix from ground truth relative permutations
QGnd = compute_Q(PijsGnd, I, nCams);

% verify that QGnd is maximized at q=vec(pis)
q = perms_cell_to_q(Pis);
-q'*QGnd*q

% noisy Q should lead to a higher cost
Q = compute_Q(Pijs, I, nCams);
-q'*Q*q

if (saveData)
    save(fileNameQ, 'Q');
end

figure, spy(Q); title('Q matrix');

% compute the exhaustive solution - this should align with the gnd-truth
[PisExhaustive, qResultExhaustive] = solve_exhaustive(Q, nCams, N);
-qResultExhaustive'*QGnd*qResultExhaustive

% we can also try a continuous relaxation and use matlab's quadprog
% NOTE: Below call doesn't work at all. (good for us)
[PisQP] = solve_quadprog(Q, nCams, N);

% now let's play with the constrained version
[Ai, ~] = compute_constraint_system(N);
A = kron(eye(nCams), Ai);
b = ones(size(A,1), 1);
% let's verify the constraints:
disp(['Cosntraint debugging: sum(A*q-b) = ', num2str(sum(A*q-b))]);

% What happpens if we just satisfy constraints? (doesn't work)
% [Pis] = solve_quadprog(A'*A, nCams, N); 

[Qcons, s] = compute_Q_constrained(Pijs, I, nCams, lambda);
% [Qcons, s] = compute_Q_constrained(PijsGnd, I, nCams, lambda);
q = perms_cell_to_q(Pis);
-q'*Qcons*q - s'*q

[PisConsExh, qConsExh] = solve_exhaustive_constrained(Qcons, s, nCams, N);
[QconsGnd, sGnd] = compute_Q_constrained(PijsGnd, I, nCams, lambda);
-qConsExh'*QconsGnd*qConsExh - sGnd'*qConsExh

figure, spy(Qcons); title('Constrained noisy Q matrix');
figure, spy(QconsGnd); title('Constrained ground truth Q matrix');

% now this does something but still really bad.
[PisQPCons] = solve_quadprog(QconsGnd, nCams, N); 

%% a tidier way to save 
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

if (saveData)
    save(fileNameAllData, 'data');
end
