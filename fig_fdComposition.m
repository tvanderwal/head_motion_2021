%% fig_fdComposition.m
% Simon Frew | NNL | BCCHRI
% FD composition by axis
load hm_analysis.mat

%% HBN-1388 Normalized RD Components vs. Mean FD by condition
% moving average of percentages plotted against mean fd
figure

condition = ["Rest1", "Rest2", "MovieDM"];
dataFields = ["r1_mean_fd_components", "r2_mean_fd_components", "mDM_mean_fd_components"]; 
for i = 1:3
    subplot(3, 1, i)
    rd_data = table2array( vertcat(hm_data.( dataFields(i) )) );

    % sort by mean fd so that moving average can happen
    [~, order] = sort(sum(rd_data, 2));
    
    rd_data = rd_data(order, :); 
    %normalize
    norm_rd_data = rd_data ./ sum(rd_data, 2);

    semilogx(sum(rd_data, 2), smoothdata(norm_rd_data, 'movmean', 200), 'x');

    title(["Normalized RD Components vs. Mean FD", "HBN-1388, "+condition(i)])
    ylabel("Moving Average of RD Components (%)")
    xlabel("Mean FD (mm), log scale")
    legend(["pitch / x-rot", "roll / y-rot", "yaw / z-rot", "x", "y", "z"])

    xline(0.3, "--k", "0.3 mm", "LineWidth", 2);
    xline(0.15, "--k", "0.15 mm", "LineWidth", 2);
    xlim([0.04, 12]);
    
end

% Rest1 Mean FD export for prism visualization
    rd_data = table2array( vertcat(hm_data.( dataFields(i) )) );
    % sort by mean fd so that moving average can happen
    [~, order] = sort(sum(rd_data, 2));
    rd_data = rd_data(order, :); 
    %normalize and smooth
    norm_rd_data = rd_data ./ sum(rd_data, 2);
    smooth_norm_rd_data = smoothdata(norm_rd_data, 'movmean', 200);

T = array2table([[hm_data(order).id]', sum(rd_data, 2), smooth_norm_rd_data], 'VariableNames', ["id", "Mean FD", "pitch", "roll", "yaw", "x", "y", "z"]);
writetable(T, fullfile('out', 'fig_fdComposition-HBN1388-Rest1-SmoothNormalizedRD.csv'))



%% HBN-1388 FD Composition Breakdown by Condition

condition = ["Rest1", "Rest2", "MovieDM"];
groups = ["low","med","high"]; 
dataFields = ["r1_mean_fd_components", "r2_mean_fd_components", "mDM_mean_fd_components"]; 
axesName = ["x", "y", "z", "pitch", "roll", "yaw"]; 

figure
for conditionIdx = 1:3 

    for grpIdx = 1:3
        
        rd_data = table2array( vertcat( hm_data( motionIdx(conditionIdx).(groups(grpIdx)) ).( dataFields(conditionIdx) ) ) );
        
        % permute columns for visualization
        rd_data = rd_data(:, [4,5,6,1,2,3]);
        
        subplot(3, 3, (conditionIdx - 1)*3 +  grpIdx)
        boxplot(rd_data, axesName)
        xlabel(condition(conditionIdx) + " - " + groups(grpIdx))
        ylabel("Mean FD (mm)")
        

    end
end
sgtitle("FD Composition Breakdown by Condition")

% Export for Rest1 and MovieDM, with CONDITION SPECIFIC GROUPINGS
condition = ["Rest1", "MovieDM"];
conditionMap = [1,3];
groups = ["low","med","high"]; 
dataFields = ["r1_mean_fd_components", "mDM_mean_fd_components"]; 
axesName = ["x", "y", "z", "pitch", "roll", "yaw"]; 

for conditionIdx = 1:length(condition) 

    for grpIdx = 1:length(groups)
        
        rd_data = table2array( vertcat( hm_data( motionIdx(conditionMap(conditionIdx)).(groups(grpIdx)) ).( dataFields(conditionIdx) ) ) );
        % permute columns for visualization
        rd_data = rd_data(:, [4,5,6,1,2,3]);
        
        T = array2table([[hm_data( motionIdx(conditionMap(conditionIdx)).(groups(grpIdx)) ).id]', rd_data], 'VariableNames', ["id", "x", "y", "z", "pitch", "roll", "yaw"]);
        writetable(T, fullfile('out', 'fig_fdComposition-HBN1388-' + condition(conditionIdx) + '-' + groups(grpIdx) + '-RD.csv'))
        
        
    end
end

%% HBN-1388 Normalized FD Composition Breakdown by Condition

condition = ["Rest1", "Rest2", "MovieDM"];
groups = ["low","med","high"]; 
dataFields = ["r1_mean_fd_components", "r2_mean_fd_components", "mDM_mean_fd_components"]; 
axesName = ["x", "y", "z", "pitch", "roll", "yaw"]; 

figure
for conditionIdx = 1:3 

    for grpIdx = 1:3
        
        rd_data = table2array( vertcat( hm_data( motionIdx(conditionIdx).(groups(grpIdx)) ).( dataFields(conditionIdx) ) ) );
        
        % permute columns for visualization
        rd_data = rd_data(:, [4,5,6,1,2,3]) ./ sum(rd_data, 2);
        
        subplot(3, 3, (conditionIdx - 1)*3 +  grpIdx)
        boxplot(rd_data, axesName)
        xlabel(condition(conditionIdx) + " - " + groups(grpIdx))
        ylabel("Mean FD (mm)")
        

    end
end
sgtitle("FD Composition Breakdown by Condition")

% Export for Rest1 and MovieDM, with CONDITION SPECIFIC GROUPINGS
condition = ["Rest1", "MovieDM"];
conditionMap = [1,3];
groups = ["low","med","high"]; 
dataFields = ["r1_mean_fd_components", "mDM_mean_fd_components"]; 
axesName = ["x", "y", "z", "pitch", "roll", "yaw"]; 

for conditionIdx = 1:length(condition) 

    for grpIdx = 1:length(groups)
        
        rd_data = table2array( vertcat( hm_data( motionIdx(conditionMap(conditionIdx)).(groups(grpIdx)) ).( dataFields(conditionIdx) ) ) );
        % permute columns for visualization
        rd_data = rd_data(:, [4,5,6,1,2,3]) ./ sum(rd_data, 2);
        
        T = array2table([[hm_data( motionIdx(conditionMap(conditionIdx)).(groups(grpIdx)) ).id]', rd_data], 'VariableNames', ["id", "x", "y", "z", "pitch", "roll", "yaw"]);
        writetable(T, fullfile('out', 'fig_fdComposition-HBN1388-' + condition(conditionIdx) + '-' + groups(grpIdx) + '-Normalized-RD.csv'))
        
        
    end
end
%% Spike FD Composition Breakdown by Condition
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
        T = array2table(subOut, 'VariableNames', axesName);
        writetable(T, fullfile('out', 'fig_fdComposition-HBN1388-' + condition(conditionIdx) + '-' + groups(grpIdx) + '-RD-spikes.csv')) 
    end
end

%% Spike Percent Composition

%% Spike % Composition by condition and group
clearvars -except hm_data motionIdx thresh groups condition 

condition = ["Rest1", "Rest2", "MovieDM"];
groups = ["low","med","high"]; 
dataFields = ["r1_fd_components", "r2_fd_components", "mDM_fd_components"]; 
thresh = 0.3;

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

%%
% for export to prism: 
% permute to ["x", "y", "z", "pitch", "roll", "yaw"];
% only export rest1 and movieDM
exportData = plotData([1,3], :, [4,5,6,1,2,3]);
axesName = ["x", "y", "z", "pitch", "roll", "yaw"]; 
condition = ["Rest1", "MovieDM"];

for conditionIdx = 1:2

    T = array2table(permute(exportData(conditionIdx, :, :), [2,3,1]), 'VariableNames', ["x", "y", "z", "pitch", "roll", "yaw"], 'RowNames', ["Low", "Med", "High"]);
    writetable(T, fullfile('out', 'fig_fdComposition-HBN1388-SpikeComposition-' + condition(conditionIdx) + '.csv'), 'WriteRowNames', true)



end
