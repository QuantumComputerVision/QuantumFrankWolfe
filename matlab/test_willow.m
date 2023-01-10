

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
%folderNameDWaveOut = ['./output/dwave/willow/'];
folderNameDWaveOut = ['./output/dwave/willow/2000Q/'];
saveData = 0;

N = 4; NC = 4;
lambda = 2.5;
completeness = 1;
swapRatio = 0 ;
nCams = NC;

figure('Renderer', 'painters', 'Position', [272 102 1048 834]);
t = tiledlayout(2,2);
t.Padding = 'none';
t.TileSpacing = 'none';
hs(1)=nexttile; beautify_plot;
hs(2)=nexttile; beautify_plot;
hs(3)=nexttile; beautify_plot;
hs(4)=nexttile; beautify_plot;
clrs = lines(7);
markers = {'o', '+', '*', 'x', 'square', 'diamond', 'pentagram', 'hexagram', 'v'};
%dispNames = {'exhaustive', 'matchEIG', 'matchLift', 'matchQuantum'};
dispNames = {'exhaustive', 'matchEIG', 'matchALS', 'matchLift', 'matchBirkhoff', 'matchQuantum (ours)'};

meanAccuracies = [];
stdAccuracies = [];
%%

for imgInd = 4:4 %length(imgsets)
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
    accuracies = [];
    
    %figure('Renderer', 'painters', 'Position', [272 102 1048 834]);
%     t = tiledlayout(5,1);
%     t.Padding = 'none';
%     t.TileSpacing = 'none';
    
    for startIndex=1:35
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
        % run state of the art        
        [PisMatchEIG, ~] = perm_sync_sota(Pijs, I, N, NC, 'matchEIG'); % MatchEIG
        [PisMatchALS, ~] = perm_sync_sota(Pijs, I, N, NC, 'matchALS'); % MatchALS
        [PisMatchLift, ~] = perm_sync_sota(Pijs, I, N, NC, 'matchLift'); % MatchLift
        [PisBirdal, ~] = perm_sync_sota(Pijs, I, N, NC, 'matchBirdal'); % Birdal'19
        
        [PisCons, qCons, currEnergyCons] = solve_exhaustive_constrained(Qcons, s, NC, N);
        
        [precisionCons, recallCons, ~, accCons] = get_synth_solution_PR(PisCons);
        [precision, recall, ~, accDWave] = get_synth_solution_PR(PisDWave);
        [precision, recall, ~, accMatchEIG] = get_synth_solution_PR(PisMatchEIG);
        [precision, recall, ~, accMatchALS] = get_synth_solution_PR(PisMatchALS);
        [precision, recall, ~, accMatchLift] = get_synth_solution_PR(PisMatchLift);
        [precision, recall, ~, accBirdal] = get_synth_solution_PR(PisBirdal);
        
        %accuracies = [accuracies, [accCons, accMatchEIG accMatchLift accDWave]'];
        accuracies = [accuracies, [accCons, accMatchEIG, accMatchALS, accMatchLift, accBirdal, accDWave]'];
        
        % save also ?
        if (saveData)
            data.Qcons = sparse(Qcons);
            data.s = s;
            data.Pijs = Pijs;
            data.N = N; % # of points
            data.I = I; % # of points
            data.nCams = nCams;
            data.lambda = lambda;
            data.swapRatio = swapRatio;
            data.completeness = completeness;
            [data.QconsRows, data.QconsCols, data.QconsValues] = find(Qcons);
            
            fileNameAllData = [folderToSave 'Q_N' num2str(N) '_NC' num2str(nCams) '_C' num2str(completeness) '_SR' num2str(swapRatio) '_L' num2str(lambda) '_' num2str(startIndex) '.mat'];
            save(fileNameAllData, 'data');
            disp(fileNameAllData)
        end
        
        xValues = [xValues startIndex];
    end
    
    set(gcf, 'currentaxes', hs(imgInd)); %subplot(2,2,1),
    plot_multiple(xValues, accuracies, clrs, markers, dispNames, 1); 
    beautify_plot;
    title(plotTitles{imgInd}, 'FontSize', 12);
    drawnow;
    
    meanAccuracy = mean(accuracies, 2);
    stdAccuracy = std(accuracies, [], 2);
    meanAccuracies = [meanAccuracies meanAccuracy];
    stdAccuracies = [stdAccuracies stdAccuracy];
end
meanAccuracies
stdAccuracies
figure; beautify_plot;
%hold on, bar(xValues, yMatrix(experiment,:), 'LineWidth', 2, 'Color', clrs(experiment,:), 'DisplayName', dispNames{experiment}, 'Marker', markers{experiment});
hBar = bar(categorical(imgsets(1:4)), meanAccuracies, 'group');

%hBar = bar(meanAccuracies, 0.8); 
for k1 = 1:size(meanAccuracies,2)
    ctr(k1,:) = bsxfun(@plus, k1, hBar(k1).XOffset');       % Note: ‘XOffset’ Is An Undocumented Feature, This Selects The ‘bar’ Centres
    ydt(k1,:) = hBar(k1).YData;                                         % Individual Bar Heights
end
hold on, errorbar(ctr', ydt',stdAccuracies ,'.');
set(gca,'xticklabel',categorical(imgsets(1:4)))
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



