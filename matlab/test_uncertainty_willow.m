

addpath('./permutations/');
addpath('./dataset/willow');
addpath('./dataset/willow/utils');
addpath('./dataset/willow/utils/rrwm');

%% Load data
%clear
%startup
imgsets = {'Car','Duck','Motorbike','Winebottle'};
plotTitles = {'(a) Car','(b) Duck','(c) Motorbike','(d) Winebottle'};
% if show match
showmatch = false;

baseFolder = './dataset/real_data/willow/';
folderNameDWaveOut = ['./output/dwave/willow/'];
saveData = 0;

N = 4; NC = 4;
lambda = 2.5;
completeness = 1;
swapRatio = 0 ;
nCams = NC;


figure('Renderer', 'painters', 'Position', [272 102 1048 834]);
beautify_plot;
% t = tiledlayout(2,2);
% t.Padding = 'none';
% t.TileSpacing = 'none';
% hs(1)=nexttile; beautify_plot;
% hs(2)=nexttile; beautify_plot;
% hs(3)=nexttile; beautify_plot;
% hs(4)=nexttile; beautify_plot;
clrs = lines(7);
markers = {'o', '+', '*', 'x', 'square', 'diamond', 'pentagram', 'hexagram', 'v'};
%dispNames = {'exhaustive', 'matchEIG', 'matchLift', 'matchQuantum'};
% dispNames = {'single-best', '2-best', '4-best', '6-best', '8-best'};
dispNames = {'single-best', '2-best', '3-best', '4-best', '5-best', '6-best'};

meanAccuracies = [];
numSolChecks = 2:6;
numData = 35;

dispNames = {'single-best'};
for i=1:length(numSolChecks)
    dispNames = [dispNames [num2str(numSolChecks(i)) '-best'] ];
end

meanAccuracies = [];
stdAccuracies = [];

accuracies = zeros(1+length(numSolChecks), numData);
for imgInd = 1:length(imgsets)
    imgset = imgsets{imgInd};
    folderToSave = ['./dataset/real_data/willow/' imgset '/'];
    
    datapath = sprintf('./dataset/willow/WILLOW-ObjectClass/%s/',imgset);
    savepath = sprintf('./dataset/willow/WILLOW-ObjectClass/%s',imgset);
    savefile = sprintf('%s/match_kpts.mat',savepath);
    imgList = dir([datapath,'*.png']);
    imgList = {imgList.name};
    viewList = dir([datapath,'*hypercols_kpts.mat']);
    viewList = {viewList.name};
    
    folderNameDWave = folderNameDWaveOut;
    folderNameDWave = [folderNameDWave imgset];
    
    xValues = [];
    
    for startIndex=1:numData
        startIndex
        
        fileNameDWave = [folderNameDWave '/Q_N4_NC4_C1_SR' num2str(swapRatio) '_L2.5_' num2str(startIndex)];
        fileNameDWave = ['.' erase(fileNameDWave,'.')  '.txt'];
        [PisCells, energies, qubits] = load_dwave_text_output(fileNameDWave, N, NC);
        PisDWave = PisCells(:,1);
        qDWave = perms_cell_to_q(PisDWave);
        
        
        pMatch = runGraphMatchBatch_N_NC(datapath,viewList,N,NC,startIndex,'all',[],'wEdge',0.0,'thscore',0);
        if (isempty(pMatch))
            continue;
        end
        % save(savefile,'pMatch');
        %         if (mod(startIndex,7)==0)
        %            nexttile; beautify_plot;
        %            visualize_willow_data(datapath,viewList,N,NC,startIndex, pMatch,0);
        %            %visPMatch(datapath,pMatch);
        %         end
        %         continue;
        %% pairwise matching
        % load(savefile,'pMatch');
        
        %     C = [];
        %     %for i = 1:length(viewList)
        %     for i = 1:length(viewList)
        %         views(i) = load(sprintf('%s/%s',datapath,viewList{i}));
        %         cnt(:,i) = sum(views(i).frame,2)/double(views(i).nfeature);
        %         C = [C,views(i).frame - repmat(cnt(:,i),1,views(i).nfeature)];
        %     end
        
        % X1 = pMatch2perm(pMatch); % pairwise matching result
        
        jMatch = runJointMatch(pMatch, [],'Method','prepare','univsize',N,...
            'rank',3,'lambda',lambda, 'initonly', 0, 'excludegeom', 0);
        
        Qcons = jMatch.Q ;
        s = jMatch.s ;
        Pijs = jMatch.Pijs ;
        I = jMatch.I ;
        
        [precision, recall, ~, accDWave] = get_synth_solution_PR(PisDWave);
        accuracies(1, startIndex)=accDWave;
        k=2;
        for numSolCheck = numSolChecks
            PisDWaveImprove = correct_sol_by_samples(PisCells, numSolCheck);
            [precision, recall, ~, accDWaveImprove] = get_synth_solution_PR(PisDWaveImprove);
            
            %accuracies = [accuracies, [accCons, accMatchEIG accMatchLift accDWave]'];
            accuracies(k, startIndex)=accDWaveImprove;
            k=k+1;
        end
        xValues = [xValues startIndex];
    end
%     
%     set(gcf, 'currentaxes', hs(imgInd)); %subplot(2,2,1),
%     plot_multiple(xValues, accuracies, clrs, markers, dispNames, 1);
%     beautify_plot;
%     title(plotTitles{imgInd}, 'FontSize', 12);
%     drawnow;
       
    meanAccuracy = mean(accuracies, 2);
    stdAccuracy = std(accuracies, [], 2);
    meanAccuracies = [meanAccuracies meanAccuracy];
    stdAccuracies = [stdAccuracies stdAccuracy];
end

supertitle('Effect of Incorporating Multiple Solutions');

figure; beautify_plot;
%hold on, bar(xValues, yMatrix(experiment,:), 'LineWidth', 2, 'Color', clrs(experiment,:), 'DisplayName', dispNames{experiment}, 'Marker', markers{experiment});
hBar = bar(categorical(imgsets(1:4)), meanAccuracies, 'group');

legend(dispNames);
ylim([0.75 1.0]);
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'XColor'      , [.0 .0 .0], ...
  'YColor'      , [.0 .0 .0], ...
  'GridLineStyle','--', ...
  'LineWidth'   , 1.0         );

if (plotLegend)
    lh=legend('show');
    set(lh,'FontSize',10);
end


