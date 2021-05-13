%% fig_hmOverview.m
% Simon Frew | NNL | BCCHRI
% overview and export of head motion data 
% generate high/med/low groups within cohort

load hm_analysis.mat


%% HBN-1388 Mean FD by condition
clearvars -except hm_data motionIdx

fd = [[hm_data.r1_mean_fd]', [hm_data.r2_mean_fd]', [hm_data.mDM_mean_fd]', [hm_data.mDM_full_mean_fd]'];

boxplot(fd)

% export for PRISM
T = array2table([[hm_data.id]', fd], 'VariableNames', ["id", "Rest1", "Rest2", "MovieDM", "MovieDM Full"]);
writetable(T, fullfile("out", "fig_hmOverview-HBN1388-MeanFDbyCondition.csv"))


%% HBN-1388 Mean FD vs. Age by condition
clearvars -except hm_data motionIdx

age = [hm_data.age];

subplot(1,3,1)
scatter(age, [hm_data.r1_mean_fd], 'x') 
ylabel("Mean FD (mm)")
xlabel("Age (y)")
title("Rest1")
ylim([0, 12])

subplot(1,3,2)
scatter(age, [hm_data.r2_mean_fd], 'x') 
ylabel("Mean FD (mm)")
xlabel("Age (y)")
title("Rest2")
ylim([0, 12])

subplot(1,3,3)
scatter(age, [hm_data.mDM_mean_fd], 'x') 
ylabel("Mean FD (mm)")
xlabel("Age (y)")
title("MovieDM")
ylim([0, 12])

% export for PRISM
out = [ age',  [hm_data.r1_mean_fd]',  [hm_data.r2_mean_fd]',  [hm_data.mDM_mean_fd]', [hm_data.mDM_full_mean_fd]'];
T = array2table([[hm_data.id]', out], 'VariableNames', ["id", "Age", "Rest1", "Rest2", "Movie-H", "Movie-F"]);
writetable(T, fullfile("out", "fig_hmOverview-HBN1388-MeanFDvsAgeByCondition.csv"))


%% HBN-865 Mean FD vs. CBCL by condition
clearvars -except hm_data motionIdx

cbclIndex = find(arrayfun(@(hm_data) ~isempty(hm_data.CBCL_Total_T), hm_data));
cbcl = [hm_data(cbclIndex).CBCL_Total_T];

subplot(1,3,1)
scatter(cbcl, [hm_data(cbclIndex).r1_mean_fd], 'x') 
ylabel("Mean FD (mm)")
xlabel("CBCL Total T")
title("Rest1")
ylim([0, 12])

subplot(1,3,2)
scatter(cbcl, [hm_data(cbclIndex).r2_mean_fd], 'x') 
ylabel("Mean FD (mm)")
xlabel("CBCL Total T")
title("Rest2")
ylim([0, 12])

subplot(1,3,3)
scatter(cbcl, [hm_data(cbclIndex).mDM_mean_fd], 'x') 
ylabel("Mean FD (mm)")
xlabel("CBCL Total T")
title("MovieDM")
ylim([0, 12])


% export for PRISM
out = [ cbcl',  [hm_data(cbclIndex).r1_mean_fd]',  [hm_data(cbclIndex).r2_mean_fd]',  [hm_data(cbclIndex).mDM_mean_fd]', [hm_data(cbclIndex).mDM_full_mean_fd]' ];
T = array2table([[hm_data(cbclIndex).id]', out], 'VariableNames', ["id", "CBCL Total T", "Rest1", "Rest2", "Movie-H", "Movie-F"]);
writetable(T, fullfile("out", "fig_hmOverview-HBN865-MeanFDvsCBCLByCondition.csv"))


%% HBN 1388/865 Groupings
clearvars -except hm_data motionIdx

% determine % overlap and export to console 
motionIdx(5).condition = "group % overlap"; 
motionIdx(5).nhigh = all(cat(1, motionIdx(1:3).high), 1) / any(cat(1, motionIdx(1:3).high), 1);
motionIdx(5).nmed = all(cat(1, motionIdx(1:3).med), 1) / any(cat(1, motionIdx(1:3).med), 1);
motionIdx(5).nlow = all(cat(1, motionIdx(1:3).low), 1) / any(cat(1, motionIdx(1:3).low), 1);
T = struct2table(motionIdx); disp(T(:, [1,5,6,7])) % visualize

for i = 1:3
    grpMembership(:, i) = [[motionIdx(i).high]'*3 + [motionIdx(i).med]'*2 + [motionIdx(i).low]'*1];
end 
T = array2table([[hm_data.id]', grpMembership], 'VariableNames', ["id", "Rest1", "Rest2", "MovieDM"]);
writetable(T, fullfile("out", "fig_hmOverview-HBN1388-MotionGroupByCondition.csv"))
% high = 3, med = 2, low = 1





%% Motion grouping numbers 
% for each group: % female, mean age, mean CBCL
conditionList = ["Rest1", "Rest2", "MovieDM", "MovieDM_full"];
groupList = ["low", "med", "high"];


for conditionIdx = 1:length(conditionList) % loop over conditions in order of ConditionList
    age = array2table(nan(1388, 3), 'VariableNames', groupList);
    cbcl = array2table(nan(1388, 3), 'VariableNames', groupList);
    
    for groupIdx = 1:length(groupList)
        idx = motionIdx(conditionIdx).( groupList(groupIdx) ); 
        
        % age 
        age{1:length([hm_data(idx).age]), groupIdx} = [hm_data(idx).age]';
        % cbcl
        cbcl{1:length([hm_data(idx).CBCL_Total_T]), groupIdx} = [hm_data(idx).CBCL_Total_T]';

    end
    writetable(cbcl, fullfile("out", "fig_hmOverview-" + conditionList(conditionIdx) + "-cbcl.csv"))
    writetable(age, fullfile("out", "fig_hmOverview-" + conditionList(conditionIdx) + "-age.csv"))
    
end

% writetable(sex, fullfile("out", "fig_demographics_motionGroups-sex.csv"), 'WriteRowNames', true)
% writetable(age, fullfile("out", "fig_demographics_motionGroups-age.csv"), 'WriteRowNames', true)

