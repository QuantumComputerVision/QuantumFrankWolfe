% clear all;
% clc;
% close all;

addpath('./permutations/');
addpath('./dataset/');
addpath('./helpers/');

% N, NC, #occurances out of 200 or #occurances/total exp., #qubits, chain length. 
% NC {N =3 N = 4 N = 5, number of physical qubits} {N =3 N = 4 N = 5, maximum chain length} 
dataFile = 'C:\Users\tolga\Box Sync\gitforks\QuantumSync\output\dwave\n_physical_qubits_table_NEW.txt';

D = load(dataFile);

figure;
t = tiledlayout(1,3);
t.Padding = 'none';
t.TileSpacing = 'none';
h1=nexttile; beautify_plot;
h2=nexttile; beautify_plot;
h3=nexttile; beautify_plot;

set(gcf, 'currentaxes', h1); %subplot(2,2,1),
hold on; plotHandles(1) = plot(D(:,1), D(:,2), 'LineWidth', 2, 'Marker','o', 'DisplayName', 'N=3');
hold on; plotHandles(2) = plot(D(1:end-1,1), D(1:end-1,3), 'LineWidth', 2, 'Marker','+', 'DisplayName', 'N=4');
hold on; plotHandles(3) = plot(D(1:3,1), D(1:3,4), 'LineWidth', 2, 'Marker','diamond', 'DisplayName', 'N=5');
% title ('(a) #qubits vs. #nodes', 'FontSize', 12);
title ('(a)', 'FontSize', 11);

set(gcf, 'currentaxes', h2);  % subplot(2,2,2);
hold on; plotHandles(1) = plot(D(:,1), D(:,5), 'LineWidth', 2, 'Marker','o', 'DisplayName', 'N=3');
hold on; plotHandles(2) = plot(D(1:end-1,1), D(1:end-1,6), 'LineWidth', 2, 'Marker','+', 'DisplayName', 'N=4');
hold on; plotHandles(3) = plot(D(1:3,1), D(1:3,7), 'LineWidth', 2, 'Marker','diamond', 'DisplayName', 'N=5');
%h = title ('(b) max. chain length vs. #nodes', 'FontSize', 12);
h = title ('(b)', 'FontSize', 11);
%pos = get ( h, 'position' );
%set(h,'position', pos-[0 29 0]);
lh=legend('show');
set(lh,'FontSize',9);

dataFile2 = 'C:\Users\tolga\Box Sync\gitforks\QuantumSync\output\dwave\PROBLEM_SIZE_EXPERIMENT.txt';
D2 = dlmread(dataFile2);

 Ns = [3,4,5];
 Ms = [3 4 5 6 7 8];
 
 [Xm,Ym] = ndgrid(Ns,Ms);
 Zm = reshape(D2(:,3),6,3)';
 Zm(isnan(Zm))=0;
 set(gcf, 'currentaxes', h3);  % subplot(2,2,2);
 surf(Xm,Ym,Zm);
 shading interp;
 xlabel('n');
 ylabel('m');
 zlabel('# occurances');
 h = title ('(c)', 'FontSize', 11);
 
export_fig('C:\Users\tolga\Box Sync\gitforks\QuantumSync\visualizations\dwave_experiment_NCvsQubits.pdf');

return;
