clear all;
clc;
close all;

addpath('./permutations/');
addpath('./dataset/');
addpath('./helpers/');
addpath('./sota/');

N = 4;
nCams = 4;
Qmatrices = {};
t = tiledlayout(2,9);
t.Padding = 'none';
t.TileSpacing = 'compact';
for h=1:18
    hs(h)=nexttile; beautify_plot;
end

h = 1;
for lambda=[0, 1.5]
    for swapRatio=0:0.1:0.2
        for completeness=1:-0.25:0.5
            [Pis, Pijs, PijsGnd, G, I, Iupper] = gen_synth_multi_graph_2d_no_grid(N, nCams, completeness, swapRatio);
            [Qcons, s] = compute_Q_constrained(Pijs, I, nCams, lambda);
            
            Qmatrices = [Qmatrices Qcons];
            
            set(gcf, 'currentaxes', hs(h)); %subplot(2,2,1),
            spy(Qcons);
            title(['C=',num2str(completeness), ', \sigma=' num2str(swapRatio)], 'FontSize', 12);
            axis off;
            drawnow;
            h=h+1;
        end
    end
end

