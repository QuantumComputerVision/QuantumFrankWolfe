
addpath('./permutations/');
addpath('./dataset/');
addpath('./helpers/');

close all;

fileName = './dataset/synth2d_dataN3NC3';
fileNameQ = './dataset/QMatrixN3NC3';
fileNameAllData = './dataset/synth_data_outputN3NC3.mat';
folderToSave = './dataset/synth_data_multi_random_exp/';
saveData = 1;

% Ns = [3,4,5];
% nCamss = [3,4,5];
% comps = [0.5 0.6 0.75 0.9 1];
% swapRatios = [0, 0.025 0.05 0.1 0.15 0.2 0.25];
% lambdas = [0, 0.5 1 1.25 1.5 2 2.5];

Ns = 4;
nCamss = [4];
comps = [0.5 0.6 0.75 0.9 1];
swapRatios = [0, 0.025 0.05 0.1 0.15 0.2 0.25];
%lambdas = [0, 0.5 1 1.25 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8];
lambdas = [0, 0.5 1 1.25 1.5 2 2.5];
numExperiments = 10;

for nsi=1:length(Ns)
    N = Ns(nsi); % the number of the grid points. permutations are NxN
    for nCamssi=1:length(nCamss)
        nCams = nCamss(nCamssi);
        for compsi=1:length(comps)
            completeness = comps(compsi); % 1 is fully connected, 0 has no connectivity
            for swapRatiosi=1:length(swapRatios)
                swapRatio = swapRatios(swapRatiosi); % like to add some noise?
                for iexp = 1:numExperiments
                    % [Pis, Pijs, PijsGnd, G, I, Iupper, points] = gen_synth_multi_graph_2d(N, nCams, completeness, swapRatio);
                    [Pis, Pijs, PijsGnd, G, I, Iupper] = gen_synth_multi_graph_2d_no_grid(N, nCams, completeness, swapRatio);
                    % compute the unconstrained Q matrix from ground truth relative permutations
                    QGnd = compute_Q(PijsGnd, I, nCams);
                    % noisy Q should lead to a higher cost
                    Q = compute_Q(Pijs, I, nCams);
                    for lambdasi=1:length(lambdas)
                        lambda = lambdas(lambdasi);
                        
                        % What happpens if we just satisfy constraints? (doesn't work)
                        % [Pis] = solve_quadprog(A'*A, nCams, N);
                        
                        [Qcons, s] = compute_Q_constrained(Pijs, I, nCams, lambda);
                        
                        
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
                            fileNameAllData = [folderToSave 'Q_N' num2str(N) '_NC' num2str(nCams) '_C' num2str(completeness) '_SR' num2str(swapRatio) '_L' num2str(lambda) '_' num2str(iexp) '.mat'];
                            save(fileNameAllData, 'data');
                            disp(fileNameAllData)
                        end
                    end
                end
            end
        end
    end
end
