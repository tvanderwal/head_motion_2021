% fig_supp.m
% Simon Frew | NNL | BCCHRI
% supplementary analyses on motion spikes and new cohort for psychometrics
% export of datasets for use in fig_supp.R

load(fullfile("..", "hm_analysis.mat"))
hm_table = struct2table(hm_data);
hm_table = hm_table(:, [1:5, 10, 14, 18, 22]);
hm_table{:, 1} = erase(hm_table{:,1}, "sub-");
%writetable(hm_table, "hm_table.csv")
%% motion spikes 
% mean and SD of motion spikes for each run in HBN 1388


fd_spikes = table2array(cat(1, hm_data.fd_spikes));
varnames = hm_data(1).fd_spikes.Properties.VariableNames;
condition = ["Rest1", "Rest2", "MovieDM", "MovieDM_lasthalf"];
groups = ["low","med","high"]; 

for conditionIdx = 1:4
    for grpIdx = 1:3
        idx = motionIdx(conditionIdx).(groups(grpIdx));
        
        outSpikes.(condition(conditionIdx))(grpIdx, :) = [ mean(fd_spikes(idx, conditionIdx)) ; std(fd_spikes(idx, conditionIdx)) ];
    end
end
% key: first column: mean, second column: std
outSpikesTable = struct2table(outSpikes, 'RowNames', groups);

outSpikesTable


% calculate mean and SD spikes per run 
mean(fd_spikes)
std(fd_spikes)

%% psychometric cohort --> see R script
%% binned motion --> see R script



