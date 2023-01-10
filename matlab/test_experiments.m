clear all;
clc;
close all;

addpath('./permutations/');
addpath('./dataset/');
addpath('./helpers/');

% N, NC, #occurances out of 200 or #occurances/total exp., #qubits, chain length. 
% NC {N =3 N = 4 N = 5, number of physical qubits} {N =3 N = 4 N = 5, maximum chain length} 
dataFile = 'C:\Users\tolga\Box Sync\gitforks\QuantumSync\output\dwave\n_physical_qubits_table.txt';

D = load(dataFile);

figure;
t = tiledlayout(1,2);
t.Padding = 'none';
t.TileSpacing = 'compact';
h1=nexttile; beautify_plot;
h2=nexttile; beautify_plot;

set(gcf, 'currentaxes', h1); %subplot(2,2,1),
hold on; plotHandles(1) = plot(D(:,1), D(:,2), 'LineWidth', 2, 'Marker','o', 'DisplayName', 'N=3');
hold on; plotHandles(2) = plot(D(1:end-1,1), D(1:end-1,3), 'LineWidth', 2, 'Marker','+', 'DisplayName', 'N=4');
hold on; plotHandles(3) = plot(D(1:3,1), D(1:3,4), 'LineWidth', 2, 'Marker','diamond', 'DisplayName', 'N=5');
title ('(a) #qubits vs. #nodes', 'FontSize', 12);


set(gcf, 'currentaxes', h2);  % subplot(2,2,2);
hold on; plotHandles(1) = plot(D(:,1), D(:,5), 'LineWidth', 2, 'Marker','o', 'DisplayName', 'N=3');
hold on; plotHandles(2) = plot(D(1:end-1,1), D(1:end-1,6), 'LineWidth', 2, 'Marker','+', 'DisplayName', 'N=4');
hold on; plotHandles(3) = plot(D(1:3,1), D(1:3,7), 'LineWidth', 2, 'Marker','diamond', 'DisplayName', 'N=5');
title ('(b) chain strength vs. #nodes', 'FontSize', 12);
lh=legend('show');
set(lh,'FontSize',10);

export_fig('C:\Users\tolga\Box Sync\gitforks\QuantumSync\visualizations\dwave_experiment_NCvsQubits.pdf');

return;

set(gcf, 'currentaxes', h1); %subplot(2,2,1)
lh = legend(plotHandles, {'binary', 'permutation'});
set(lh,'FontSize',10);
title ('(a) Energies from binary vs. permutation variables', 'FontSize', 12);
xlabel('\lambda','FontSize',12); ylabel('$\mathbf{q}^\top\mathbf{Q}\mathbf{q} + \mathbf{s}^{\top}\mathbf{q}$','Interpreter','latex','FontSize',12);