clear all;
clc;
close all;

addpath('./permutations/');
addpath('./dataset/');
addpath('./helpers/');


N = 3; % the number of the grid points. permutations are NxN
nCams = 3;
completeness = 1.0; % 1 is fully connected, 0 has no connectivity
swapRatio = 0.2; % like to add some noise?

genereatePRPlots = 0;

% 18 bits preceeded by 3x3 identity matrices
allPossibleSolutions = makebits_augmented(N*N*(nCams-1), N);

swapRatios = [0, 0.01 0.02 0.03 0.04 0.05, 0.075 0.1, 0.0125 0.15 0.175 0.2, 0.225 0.25 0.275 0.3 0.325 0.35 0.375 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75];

numTotalExperiments = 7;
% clrs = distinguishable_colors(numTotalExperiments);
%clrs = colormap;
clrs = lines(numTotalExperiments);
markers = {'o', '+', '*', 'x', 'square', 'diamond', 'pentagram', 'hexagram', 'v'};

% debug lambda now
%lambdas = [0, 0.01, 0.015 0.025, 0.04 0.05, 0.075 0.1 0.125 0.15 0.2 0.25 0.3 0.35 0.4 0.425 0.45 0.5 0.525 0.55 0.6 0.75 1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 10 11];
% lambdas = [0, 0.01, 0.015 0.025, 0.04 0.05, 0.075 0.1 0.125 0.15 0.2 0.25 0.3 0.31 0.32 0.33 0.34 0.35 0.36 0.37 0.38 0.39 0.4 0.425 0.45 0.5 0.525 0.55 0.6 0.75 1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 6.5 7 7.5];
lambdas = [0, 0.01, 0.015 0.025, 0.04 0.05, 0.075 0.1 0.125 0.15 0.2 0.25 0.3 0.31 0.32 0.33 0.34 0.35 0.36 0.37 0.38 0.39 0.4 0.425 0.45 0.5 0.525 0.55 0.6 0.75 1 1.5 2 2.5 3 3.5 4 4.5 5];

lambdaLegendNames = {};
for i=1:length(lambdas)
    lambdaLegendNames = [lambdaLegendNames ['\lambda=' num2str(lambdas(i))]];
end
plotHandles = zeros(2, 1);

if (genereatePRPlots)
    % figure, subplot(2,2,1); beautify_plot;
    % subplot(2,2,2); beautify_plot;
    % subplot(2,2,3); beautify_plot;
    % subplot(2,2,4); beautify_plot;
    figure;
    t = tiledlayout(2,2);
    t.Padding = 'none';
    t.TileSpacing = 'compact';
    h1=nexttile; beautify_plot;
    h2=nexttile; beautify_plot;
    h3=nexttile; beautify_plot;
    h4=nexttile; beautify_plot;
    for experiment=1:numTotalExperiments
        
        disp (['Experiment ' num2str(experiment)]);
        
        [Pis, Pijs, PijsGnd, G, I, Iupper] = gen_synth_multi_graph_2d_no_grid(N, nCams, completeness, swapRatio);
        
        [Qcons0, s0] = compute_Q_constrained(Pijs, I, nCams, 0);
        
        [PisCons, qCons] = solve_exhaustive_constrained(Qcons0, s0, nCams, N); % solves only permutations
        currEnergy = -qCons'*Qcons0*qCons - s0'*qCons;
        
        energies = []; energies0 = [];
        differences = [];
        accuracies = []; accuraciesCons = [];
        PR = []; PRCons = [];
        minDiffEnergy = 9999;
        bestInd = 1;
        lambdaSol = lambdas(bestInd);
        bestq = zeros(length(qCons),1);
        bestPis = PisCons;
        for i=1:length(lambdas)
            lambda = lambdas(i);
            [Qcons, s] = compute_Q_constrained(Pijs, I, nCams, lambda);
            % [Qcons, s] = compute_Q_constrained(PijsGnd, I, nCams, lambda);
            [PisCons, qCons] = solve_exhaustive_constrained(Qcons, s, nCams, N);
            currEnergy = -qCons'*Qcons*qCons - s'*qCons;
            
            [PisExh, qConsExh] = solve_exhaustive_unconstrained(Qcons, s, nCams, N, allPossibleSolutions);
            
            currEnergyExh = -qConsExh'*Qcons*qConsExh - s'*qConsExh;
            
            diffEnergy = abs(currEnergyExh-currEnergy);
            disp(['lambda: ' num2str(lambda) ', diff: ', num2str(diffEnergy)]);
            
            if (diffEnergy<minDiffEnergy)
                lambdaSol = lambda;
                bestInd = i;
                bestq = qConsExh;
                bestPis = PisExh;
            end
            
            % accCons = get_synth_solution_acc(PisCons);
            % acc = get_synth_solution_acc(PisExh);
            [precision, recall, F, acc] = get_synth_solution_PR(PisExh);
            [precisionCons, recallCons, FCons, accCons] = get_synth_solution_PR(PisCons);
            
            energies0 = [energies0, currEnergy];
            energies = [energies, currEnergyExh];
            differences = [differences, diffEnergy];
            accuracies = [accuracies, acc];
            accuraciesCons = [accuraciesCons, accCons];
            PR = [PR [precision; recall]];
            PRCons = [PRCons [precisionCons; recallCons]];
        end
        
        set(gcf, 'currentaxes', h1); %subplot(2,2,1),
        hold on; plotHandles(1) = plot(lambdas, energies, 'LineWidth', 1, 'Marker','o', 'Color', clrs(experiment,:));
        hold on; plotHandles(2) = plot(lambdas, energies0, 'LineWidth', 1, 'Color', clrs(experiment,:));
        
        set(gcf, 'currentaxes', h2);  % subplot(2,2,2);
        hold on, plot(lambdas, differences, 'LineWidth', 1, 'Color', clrs(experiment,:), 'DisplayName', ['Exp - ', num2str(experiment)], 'Marker', markers{experiment});
        
        set(gcf, 'currentaxes', h3); %subplot(2,2,3);
        hold on, plot(lambdas, accuracies, 'LineWidth', 1, 'Color', clrs(experiment,:), 'DisplayName', ['Exp - ', num2str(experiment)], 'Marker', markers{experiment});
        % hold on, plot(PR(1,:), PR(2,:), 'LineWidth', 1, 'Color', clrs(experiment,:), 'DisplayName', ['Exp - ', num2str(experiment)]);
        
        set(gcf, 'currentaxes', h4); % subplot(2,2,4);
        hold on, plot(lambdas, accuraciesCons, 'LineWidth', 1, 'Color', clrs(experiment,:), 'DisplayName', ['Exp - ', num2str(experiment)], 'Marker', markers{experiment});
        % hold on, plot(PRCons(1,:), PRCons(2,:), 'LineWidth', 1, 'Color', clrs(experiment,:), 'DisplayName', ['Exp - ', num2str(experiment)]);
        
        drawnow;
        
    end
    
    set(gcf, 'currentaxes', h1); %subplot(2,2,1)
    lh = legend(plotHandles, {'binary', 'permutation'});
    set(lh,'FontSize',10);
    title ('(a) Energies from binary vs. permutation variables', 'FontSize', 12);
    xlabel('\lambda','FontSize',12); ylabel('$\mathbf{q}^\top\mathbf{Q}\mathbf{q} + \mathbf{s}^{\top}\mathbf{q}$','Interpreter','latex','FontSize',12);
    set(gcf, 'currentaxes', h2); % subplot(2,2,2);
    lh = legend('show'); set(lh,'FontSize',10);
    title ('(b) Gap between lowest energies', 'FontSize', 12);
    xlabel('\lambda','FontSize',12); ylabel('$| E_{\mathrm{binary}} - E_{\mathrm{perm.}} |$','Interpreter','latex','FontSize',12);
    set(gcf, 'currentaxes', h3); % subplot(2,2,3);
    title ('(c) Accuracy with binary constraints w.r.t. \lambda', 'FontSize', 12);
    xlabel('\lambda','FontSize',12); ylabel('% correct bits','FontSize',12,'Interpreter','none');
    %title ('(c) PR curve of the unconstrained solution w.r.t. \lambda', 'FontSize', 12);
    % xlabel('precision','FontSize',12); ylabel('recall','FontSize',12,'Interpreter','none');
    % lh = legend('show'); set(lh,'FontSize',10);
    set(gcf, 'currentaxes', h4); % subplot(2,2,4);
    title ('(d) Accuracy with perm. constraints w.r.t. \lambda', 'FontSize', 12);
    xlabel('\lambda','FontSize',12); ylabel('% correct bits','FontSize',12,'Interpreter','none');
    %title ('(d) PR curve of the constrained solution w.r.t. \lambda', 'FontSize', 12);
    % xlabel('precision', 'FontSize', 12); ylabel('recall','FontSize',12,'Interpreter','none');
    % lh = legend('show'); set(lh,'FontSize',10);
    
    supertitle = ['Over ' num2str(numTotalExperiments) ' Different Experiments using N = ' num2str(N) ', #Nodes = ', num2str(nCams), ', \sigma = ', num2str(swapRatio)];
    % suptitle(supertitle);
    title(t, supertitle);
    
end

legendNames = {};
lambda = 2; % some optimum lambda
energiesNoise = []; energies0Noise = [];
differencesNoise = []; 
accuraciesNoise = []; accuraciesConsNoise = [];
PRNoise = []; PRConsNoise = [];
figure;
t2 = tiledlayout(1,2);
t2.Padding = 'none';
t2.TileSpacing = 'compact';
h21=nexttile; beautify_plot;
h22=nexttile; beautify_plot;
for i=1:length(swapRatios)
    swapRatio = swapRatios(i);
    currEnergy=0;
    currEnergyExh=0;
    acc=0;
    accCons=0;
    diffEnergy=0;
    for experiment=1:numTotalExperiments
        [Pis, Pijs, PijsGnd, G, I, Iupper] = gen_synth_multi_graph_2d_no_grid(N, nCams, completeness, swapRatio);
        
        [Qcons, s] = compute_Q_constrained(Pijs, I, nCams, lambda);
        % [Qcons, s] = compute_Q_constrained(PijsGnd, I, nCams, lambda);
        [PisCons, qCons] = solve_exhaustive_constrained(Qcons, s, nCams, N);
        currEnergy = currEnergy + (-qCons'*Qcons*qCons - s'*qCons);
        
        [PisExh, qConsExh] = solve_exhaustive_unconstrained(Qcons, s, nCams, N, allPossibleSolutions);
        
        currEnergyExh = currEnergyExh + (-qConsExh'*Qcons*qConsExh - s'*qConsExh);
        
        diffEnergy = diffEnergy + abs(currEnergyExh-currEnergy);
        disp(['lambda: ' num2str(lambda) ', swapRatio: ' num2str(swapRatio) ', diff: ', num2str(diffEnergy)]);
        
        [precision, recall, F, accCur] = get_synth_solution_PR(PisExh);
        [precisionCons, recallCons, FCons, accConsCur] = get_synth_solution_PR(PisCons);
        acc = acc + accCur;
        accCons = accCons + accConsCur;
    end
    accCons = accCons./numTotalExperiments;
    acc = acc./numTotalExperiments;
    currEnergyExh = currEnergyExh./numTotalExperiments;
    currEnergy = currEnergy./numTotalExperiments;
    diffEnergy = diffEnergy./numTotalExperiments;
        
    energies0Noise = [energies0Noise, currEnergy];
    energiesNoise = [energiesNoise, currEnergyExh];
    differencesNoise = [differencesNoise, diffEnergy];
    accuraciesNoise = [accuraciesNoise, acc];
    accuraciesConsNoise = [accuraciesConsNoise, accCons];
    
    legendNames = [legendNames, ['$\sigma=$' num2str(swapRatio)] ];
end

set(gcf, 'currentaxes', h21);
hold on; plotHandles(2) = plot(swapRatios, energies0Noise, 'LineWidth', 2, 'Color', clrs(2,:), 'Marker', markers{1});
hold on; plotHandles(1) = plot(swapRatios, energiesNoise, 'LineWidth', 2, 'Color', clrs(1,:), 'Marker', markers{2});
lh = legend(plotHandles, {'binary constraints', 'permutations'}); set(lh,'FontSize',10);
title ('(a) Effect of Constraints on Energy', 'FontSize', 12);
xlabel('swap ratios (\sigma)','FontSize',12); ylabel('$\mathbf{q}^\top\mathbf{Q}\mathbf{q} + \mathbf{s}^{\top}\mathbf{q}$','Interpreter','latex','FontSize',12);
%hold on, yline(currEnergy);

set(gcf, 'currentaxes', h22);
%hold on, plot(swapRatios, differencesNoise, 'LineWidth', 2);
%title ('(b) Gap between constrained and unconstrained energies', 'FontSize', 12);
%xlabel('swap ratios (\sigma)','FontSize',12); ylabel('$| E_{\mathrm{cons}} - E_{\mathrm{uncons}} |$','Interpreter','latex','FontSize',12);
hold on, plot(swapRatios, accuraciesNoise, 'LineWidth', 2, 'Color', clrs(1,:), 'Marker', markers{1});
hold on, plot(swapRatios, accuraciesConsNoise, 'LineWidth', 2, 'Color', clrs(2,:), 'Marker', markers{2});
lh = legend(plotHandles, {'binary constraints', 'permutations'}); set(lh,'FontSize',10);
title ('(b) Accuracy vs. Noise', 'FontSize', 12);
xlabel('swap ratios (\sigma)','FontSize',12); ylabel('% correct bits','FontSize',12);
supertitle2 = ['Behavior under noise. \lambda = ' num2str(lambda) ', N = ' num2str(N) ', #Nodes = ', num2str(nCams), ', \sigma = ', num2str(swapRatio)];
% suptitle(supertitle);
title(t2, supertitle2);

drawnow;
