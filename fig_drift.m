%% fig_drift.m
% Simon Frew | NNL | BCCHRI
% overview and export of head motion data 
% generate high/med/low groups within cohort

load hm_analysis.mat

%% R1, R2, Movie DM Drift in time across groups

rest_volumes = 354; movie_volumes = 729; tr = 0.8;
gap = 120; gap_volumes = 120 / 0.8; % number of volumes in gap 

time_vol = 1:(rest_volumes*2 + gap_volumes*2 + movie_volumes);
time = time_vol*tr / 60; % time in minutes

condition = ["Rest1", "Rest2", "MovieDM"];
groups = ["low","med","high"]; 
dataFields = ["r1_fd", "r2_fd", "mDM_full_fd"]; 

meanPlotData = zeros(length(time), length(condition)); 
sdPlotData = zeros(length(time), length(condition)); 

figure
hold on 
for conditionIdx = 1:3 
    for grpIdx = 1:3
        fd_data(conditionIdx).motionGrouping = condition(1);
        fd_data(conditionIdx).condition = condition(conditionIdx);
        fd_data(conditionIdx).(groups(grpIdx)) = vertcat( hm_data( motionIdx(1).(groups(grpIdx)) ).( dataFields(conditionIdx) ) );
    end
end
for grpIdx = 1:3 
    meanPlotData(:, grpIdx) = [mean( fd_data(1).(groups(grpIdx)) ), ...
                           [1:gap_volumes]*0, mean( fd_data(2).(groups(grpIdx)) ), ...
                           [1:gap_volumes]*0, mean( fd_data(3).(groups(grpIdx)) )];
%    sdPlotData(:, grpIdx) = [std( fd_data(1).(groups(grpIdx)) ), ...
%                           [1:gap_volumes]*0, std( fd_data(2).(groups(grpIdx)) ), ...
%                           [1:gap_volumes]*0, std( fd_data(3).(groups(grpIdx)) )];
                       
    sdPlotData(:, grpIdx) = [mad( fd_data(1).(groups(grpIdx)) ), ...
                           [1:gap_volumes]*0, mad( fd_data(2).(groups(grpIdx)) ), ...
                           [1:gap_volumes]*0, mad( fd_data(3).(groups(grpIdx)) )];
                       
    patch([time, fliplr(time)],...
        [meanPlotData(:, grpIdx)' + sdPlotData(:, grpIdx)', fliplr(meanPlotData(:, grpIdx)' - sdPlotData(:, grpIdx)')],...
        'k', 'FaceColor', [grpIdx/3,0.5,0.5], 'FaceAlpha', 0.25, 'EdgeColor', 'none')
    
end



plot(time, meanPlotData)
legend(["MAD " + groups, groups])

xlabel("Scan time (min)")
ylabel("Within-volume mean FD (mm)")
title(["HBN-1388: Head motion drift over scan session", "Rest1 Motion Groups"])

%% %%%% Drift Slope 
clearvars -except hm_data motionIdx


%% calculate fd slopes for each volume 
index = 1:354'; 
indexFull = 1:729';

for i = 1:length(hm_data)
    slopes(i).Rest1 = polyfit(index, [hm_data(i).r1_fd], 1);
    slopes(i).Rest2 = polyfit(index, [hm_data(i).r2_fd], 1);
    slopes(i).MovieDM = polyfit(index, [hm_data(i).mDM_fd], 1);
    slopes(i).MovieDMFull = polyfit(indexFull, [hm_data(i).mDM_full_fd], 1);
end

condition = ["Rest1", "Rest2", "MovieDM", "MovieDMFull"];
conditionNames = ["Rest1", "Rest2", "Movie-H", "Movie-F"];
fd = ["r1_mean_fd", "r2_mean_fd", "mDM_mean_fd", "mDM_full_mean_fd"];
CBCLIndex = find(arrayfun(@(hm_data) ~isempty(hm_data.CBCL_Total_T), hm_data));

%% Plot 
figure
for conditionIdx = 1:length(condition)
    Y = cat(1, slopes.(condition(conditionIdx))); Y = Y(:, 1); 

    % mean FD x slope
        subplot(4, 3, (conditionIdx - 1)*3 + 1)

        X = [hm_data.(fd(conditionIdx))]';
        plot(X, Y, 'x')
        [r, p] = corr(X, Y)
        title(condition(conditionIdx) + ": abs(Drift Slope) vs. Mean FD")
        xlabel("Mean FD (mm)")
        ylabel("Drift slope (mm/volume)")
        ylim([0, 0.04])
    
    % age x slope 
        subplot(4, 3, (conditionIdx - 1)*3 + 2)

        X = [hm_data.age]';
        plot(X, Y, 'x')
        [r, p] = corr(X, Y)
        title(condition(conditionIdx) + ": abs(Drift Slope) vs. Age")
        xlabel("Age (y)")
        ylabel("Drift slope (mm/volume)")
        ylim([0, 0.04])

        
    % cbcl x slope 
        subplot(4, 3, (conditionIdx - 1)*3 + 3)

        X = [hm_data(CBCLIndex).CBCL_Total_T]';
        Y = Y(CBCLIndex); 
        plot(X, Y, 'x')
        [r, p] = corr(X, Y)
        title(condition(conditionIdx) + ": abs(Drift Slope) vs. CBCL")
        xlabel("CBCL Total T")
        ylabel("Drift slope (mm/volume)")
        ylim([0, 0.04])

    
end

%% Group vs. Condition Mean Drift Slope 

motionGrp = ["low", "med", "high"]; 
figure
%sgtitle(["Drift slope by motion group by condition", ""])

for conditionIdx = 1:length(condition)
    
    subplot(1,4,conditionIdx)
    
    Y = cat(1, slopes.(condition(conditionIdx))); Y = Y(:, 1); 
    % within-condition goruping
%     idx = [motionIdx(conditionIdx).low] + [motionIdx(conditionIdx).med]*2 + [motionIdx(conditionIdx).high]*3;
    % rest1 grouping 
    idx = [motionIdx(1).low] + [motionIdx(1).med]*2 + [motionIdx(1).high]*3;
    
    
    boxplot(Y, idx)
    
    xticklabels(motionGrp)
    ylabel("Drift slope (mm/volume)")
    xlabel("Motion Group")
    title(condition(conditionIdx))
    ylim([0, 0.04])
    
end 

%% export Group vs. Condition Mean Drift Slope 

for conditionIdx = 1:length(condition)
    Y = cat(1, slopes.(condition(conditionIdx))); Y = Y(:, 1); 
    outTable = array2table(nan(length(Y), 3), "VariableNames", motionGrp);
    
    for grpIdx = 1:length(motionGrp)
        tmpIdx = motionIdx(conditionIdx).( motionGrp(grpIdx) );
        outTable{1:sum(tmpIdx), grpIdx} = Y(tmpIdx);
    end
    % export
    writetable(outTable, fullfile("out", sprintf("fig_drift_slopes-%s.csv", condition(conditionIdx))));
end

%% export motion groups for all conditions 
for grpIdx = 1:length(motionGrp)
    outTable = array2table(nan(length(slopes), 4), "VariableNames", conditionNames);
    
    for conditionIdx = 1:length(condition)
        Y = cat(1, slopes.(condition(conditionIdx))); Y = Y(:, 1); 
        tmpIdx = motionIdx(conditionIdx).( motionGrp(grpIdx) );
        outTable{1:sum(tmpIdx), conditionIdx} = Y(tmpIdx);
    end
    % export
    writetable(outTable, fullfile("out", sprintf("fig_drift_slopes-%s.csv", motionGrp(grpIdx))));
end















