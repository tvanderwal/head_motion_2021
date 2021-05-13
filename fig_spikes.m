%% fig_spikes.m
% Simon Frew | NNL | BCCHRI
% determination and visualization of spikes

load hm_analysis.mat

condition = ["Rest1", "Rest2", "MovieDM"];
groups = ["low","med","high"]; 
thresh = 0.3;

%% HBN-1388 FD Spike Composition Breakdown by Condition

dataFields = ["r1_fd_components", "r2_fd_components", "mDM_fd_components"]; 
axesName = ["x", "y", "z", "pitch", "roll", "yaw"]; 

spikeIdx.group = "all";
spikeIdx.Rest1 = arrayfun(@(sub) sub.r1_fd > thresh, hm_data, 'UniformOutput', false);
spikeIdx.Rest2 = arrayfun(@(sub) sub.r2_fd > thresh, hm_data, 'UniformOutput', false);
spikeIdx.MovieDM = arrayfun(@(sub) sub.mDM_fd > thresh, hm_data, 'UniformOutput', false);

for conditionIdx = 1:length(condition)
    for grpIdx = 1:length(groups)
        
        fd = cat(3, hm_data(motionIdx(conditionIdx).(groups(grpIdx))).(dataFields(conditionIdx))); 
        
        idx = spikeIdx.(condition(conditionIdx)) (motionIdx(conditionIdx).(groups(grpIdx)) );
        idx = cat(1, idx{:,:});
        
        % remove subjects with no spikes
        noSpikeIdx = find(sum(idx, 2) == 0); 
        idx(noSpikeIdx, :) = []; 
        fd(:, :, noSpikeIdx) = [];
        
        subOut = zeros(size(idx, 1), 6);
        
        for subj = 1:size(idx, 1)
            % extract 
            tmp = fd(idx(subj, :), :, subj);
%             % normalize to percentage
%             tmpNorm = tmp ./ sum(tmp, 2); 
            % avg 
            subOut(subj, :) = mean(tmp);
        end
        
        plotData(conditionIdx, grpIdx, :) = mean(subOut); 
        
        % plot
        subplot(3, 3, (conditionIdx - 1)*3 + grpIdx)
        bar(squeeze(plotData(conditionIdx, grpIdx, :))')
        title(condition(conditionIdx) + " : " + groups(grpIdx) + ", FD Spike Composition (mm)")
        ylabel("Mean FD (mm)")
        xticklabels(["pitch", "roll", "yaw", "x", "y", "z"])

    end
end



%% FD Spike Composition Breakdown by conditionExport for Rest1 and MovieDM, with CONDITION SPECIFIC GROUPINGS
condition = ["Rest1", "MovieDM"];
conditionMap = [1,3];
groups = ["low","med","high"]; 
dataFields = ["r1_fd_components", "mDM_fd_components"]; 
axesName = ["x", "y", "z", "pitch", "roll", "yaw"]; 

for conditionIdx = 1:length(condition)
    for grpIdx = 1:length(groups)
        
        fd = cat(3, hm_data(motionIdx(conditionIdx).(groups(grpIdx))).(dataFields(conditionIdx))); 
        
        idx = spikeIdx.(condition(conditionIdx)) (motionIdx(conditionIdx).(groups(grpIdx)) );
        idx = cat(1, idx{:,:});
        
        % remove subjects with no spikes
        noSpikeIdx = find(sum(idx, 2) == 0); 
        idx(noSpikeIdx, :) = []; 
        fd(:, :, noSpikeIdx) = [];
        
        subOut = zeros(size(idx, 1), 6);
        
        for subj = 1:size(idx, 1)
            % extract 
            tmp = fd(idx(subj, :), :, subj);
%             % normalize to percentage
%             tmpNorm = tmp ./ sum(tmp, 2); 
            % avg 
            subOut(subj, :) = mean(tmp);
        end
        
        % permute columns for visualization
        subOut = subOut(:, [4,5,6,1,2,3]);
        T = array2table(subOut, 'VariableNames', axes);
        writetable(T, fullfile('out', 'fig_spikes-HBN1388-' + condition(conditionIdx) + '-' + groups(grpIdx) + '-RD.csv')) 
    end
end

%% # of spikes per subject vs. age, cbcl
thresh = 0.3; % 0.3mm thresh for spikes

spikeN.Rest1 = arrayfun(@(sub) sum(sub.r1_fd > thresh), hm_data);
spikeN.Rest2 = arrayfun(@(sub) sum(sub.r2_fd > thresh), hm_data);
spikeN.MovieDM = arrayfun(@(sub) sum(sub.mDM_fd > thresh), hm_data);

condition = ["Rest1", "Rest2", "MovieDM"];
groups = ["low","med","high"]; 
CBCLIndex = find(arrayfun(@(hm_data) ~isempty(hm_data.CBCL_Total_T), hm_data));

figure
for conditionIdx = 1:length(condition)
    Y = spikeN.(condition(conditionIdx)); 

    % Age x SpikeN
        subplot(3, 2, (conditionIdx - 1)*2 + 1)

        X = [hm_data.age]';
        plot(X, Y, 'x')
        title(condition(conditionIdx) + ": # of spikes vs. Age")
        xlabel("Age (y)")
        ylabel("# of spikes")
        
    % CBCL x SpikeN
        subplot(3, 2, (conditionIdx - 1)*2 + 2)

        X = [hm_data(CBCLIndex).CBCL_Total_T]';
        Y = Y(CBCLIndex); 
        
        plot(X, Y, 'x')
        title(condition(conditionIdx) + ": # of spikes vs. CBCL")
        xlabel("CBCL Total T")
        ylabel("# of spikes")
end

legend(["pitch / x-rot", "roll / y-rot", "yaw / z-rot", "x", "y", "z"])

%% Spikes over time 
rest_volumes = 354; movie_volumes = 729; tr = 0.8;
gap = 120; gap_volumes = 120 / 0.8; % number of volumes in gap 

time_vol = 1:(rest_volumes*2 + gap_volumes*2 + movie_volumes);
time = time_vol*tr / 60; % time in minutes

condition = ["Rest1", "Rest2", "MovieDMFull"];
groups = ["low","med","high"]; 
dataFields = ["r1_fd", "r2_fd", "mDM_full_fd"]; 

meanPlotData = zeros(length(time), length(condition)); 
sdPlotData = zeros(length(time), length(condition)); 

figure
hold on 


spikeIdx.group = "all";
tmp = arrayfun(@(sub) sub.r1_fd > thresh, hm_data, 'UniformOutput', false);
spikeIdx.Rest1 = cat(1, tmp{:,:});
tmp = arrayfun(@(sub) sub.r2_fd > thresh, hm_data, 'UniformOutput', false);
spikeIdx.Rest2 = cat(1, tmp{:,:});
tmp = arrayfun(@(sub) sub.mDM_full_fd > thresh, hm_data, 'UniformOutput', false);
spikeIdx.MovieDMFull = cat(1, tmp{:,:});

for conditionIdx = 1:length(condition)
    for grpIdx = 1:length(groups)
        grpSpike(conditionIdx).condition = condition(conditionIdx);
        grpSpike(conditionIdx).(groups(grpIdx)) = spikeIdx.(condition(conditionIdx)) (motionIdx(conditionIdx).(groups(grpIdx)), :); 
    end
end

for grpIdx = 1:3 
    meanPlotData(:, grpIdx) = [mean( grpSpike(1).(groups(grpIdx)) ), ...
                           [1:gap_volumes]*0, mean( grpSpike(2).(groups(grpIdx)) ), ...
                           [1:gap_volumes]*0, mean( grpSpike(3).(groups(grpIdx)) )];
%     sdPlotData(:, grpIdx) = [std( grpSpike(1).(groups(grpIdx)) ), ...
%                            [1:gap_volumes]*0, std( grpSpike(2).(groups(grpIdx)) ), ...
%                            [1:gap_volumes]*0, std( grpSpike(3).(groups(grpIdx)) )];
%                        
%     sdPlotData(:, grpIdx) = [mad( grpSpike(1).(groups(grpIdx)) ), ...
%                            [1:gap_volumes]*0, mad( grpSpike(2).(groups(grpIdx)) ), ...
%                            [1:gap_volumes]*0, mad( grpSpike(3).(groups(grpIdx)) )];
%                        
%     patch([time, fliplr(time)],...
%         [meanPlotData(:, grpIdx)' + sdPlotData(:, grpIdx)', fliplr(meanPlotData(:, grpIdx)' - sdPlotData(:, grpIdx)')],...
%         'k', 'FaceColor', [0.5,0.5,0.5], 'FaceAlpha', 0.25, 'EdgeColor', 'none')
%     
end

plot(time, meanPlotData)
% legend(["MAD " + groups, groups])
legend(groups)
ylim([0, 1])

xlabel("Scan time (min)")
ylabel("% of subjects with spike")
title(["HBN-1388: Head motion spikes over scan session", "Condition-specific Motion Groups"])

 
%% %% R1, R2, Movie DM Drift in time across groups

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
        fd_data(conditionIdx).condition = condition(conditionIdx);
        fd_data(conditionIdx).(groups(grpIdx)) = vertcat( hm_data( motionIdx(conditionIdx).(groups(grpIdx)) ).( dataFields(conditionIdx) ) );
    end
end




for grpIdx = 1:3 
    meanPlotData(:, grpIdx) = [mean( fd_data(1).(groups(grpIdx)) ), ...
                           [1:gap_volumes]*0, mean( fd_data(2).(groups(grpIdx)) ), ...
                           [1:gap_volumes]*0, mean( fd_data(3).(groups(grpIdx)) )];
%     sdPlotData(:, grpIdx) = [std( fd_data(1).(groups(grpIdx)) ), ...
%                            [1:gap_volumes]*0, std( fd_data(2).(groups(grpIdx)) ), ...
%                            [1:gap_volumes]*0, std( fd_data(3).(groups(grpIdx)) )];
                       
    sdPlotData(:, grpIdx) = [mad( fd_data(1).(groups(grpIdx)) ), ...
                           [1:gap_volumes]*0, mad( fd_data(2).(groups(grpIdx)) ), ...
                           [1:gap_volumes]*0, mad( fd_data(3).(groups(grpIdx)) )];
                       
    patch([time, fliplr(time)],...
        [meanPlotData(:, grpIdx)' + sdPlotData(:, grpIdx)', fliplr(meanPlotData(:, grpIdx)' - sdPlotData(:, grpIdx)')],...
        'k', 'FaceColor', [0.5,0.5,0.5], 'FaceAlpha', 0.25, 'EdgeColor', 'none')
    
end



plot(time, meanPlotData)
legend(["MAD " + groups, groups])

xlabel("Scan time (min)")
ylabel("Within-volume mean FD (mm)")
title(["HBN-1388: Head motion drift over scan session", "Condition-specific Motion Groups"])

%% Spike % Composition by condition and group
clearvars -except hm_data motionIdx thresh groups condition 

condition = ["Rest1", "Rest2", "MovieDM"];
groups = ["low","med","high"]; 
dataFields = ["r1_fd_components", "r2_fd_components", "mDM_fd_components"]; 

spikeIdx.group = "all";
spikeIdx.Rest1 = arrayfun(@(sub) sub.r1_fd > thresh, hm_data, 'UniformOutput', false);
spikeIdx.Rest2 = arrayfun(@(sub) sub.r2_fd > thresh, hm_data, 'UniformOutput', false);
spikeIdx.MovieDM = arrayfun(@(sub) sub.mDM_fd > thresh, hm_data, 'UniformOutput', false);

plotData = zeros(length(condition), length(groups), 6);

figure

for conditionIdx = 1:length(condition)
    for grpIdx = 1:length(groups)
        
        fd = cat(3, hm_data(motionIdx(conditionIdx).(groups(grpIdx))).(dataFields(conditionIdx))); 
        
        idx = spikeIdx.(condition(conditionIdx)) (motionIdx(conditionIdx).(groups(grpIdx)) );
        idx = cat(1, idx{:,:});
        
        % remove subjects with no spikes
        noSpikeIdx = find(sum(idx, 2) == 0); 
        idx(noSpikeIdx, :) = []; 
        fd(:, :, noSpikeIdx) = [];
        
        subOut = zeros(size(idx, 1), 6);
        
        for subj = 1:size(idx, 1)
            % extract 
            tmp = fd(idx(subj, :), :, subj);
            % normalize to percentage
            tmpNorm = tmp ./ sum(tmp, 2); 
            % avg 
            subOut(subj, :) = mean(tmpNorm);
        end
        
        plotData(conditionIdx, grpIdx, :) = mean(subOut); 
        
        % plot
        subplot(3, 3, (conditionIdx - 1)*3 + grpIdx)
        pie(plotData(conditionIdx, grpIdx, :))
        title(condition(conditionIdx) + " : " + groups(grpIdx) + ", FD Composition of spikes (%)")
    end
end

legend(["pitch", "roll", "yaw", "x", "y", "z"])





