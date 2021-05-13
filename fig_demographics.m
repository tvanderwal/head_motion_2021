%% fig_demographics.m
% Simon Frew | NNL | BCCHRI
% setup and export demographic data

load hm_analysis.mat

%% overlap 

% low 
low = sum(all(cat(1, motionIdx(1:3).low), 1)) / sum(any(cat(1, motionIdx(1:3).low), 1))
% med 
med = sum(all(cat(1, motionIdx(1:3).med), 1)) / sum(any(cat(1, motionIdx(1:3).med), 1))
% high 
high = sum(all(cat(1, motionIdx(1:3).high), 1)) / sum(any(cat(1, motionIdx(1:3).high), 1))



%% HBN-1388 Age/Sex Distribution
clearvars -except hm_data

% age histogram with sex 
ageData = [hm_data.age];
sexData = [hm_data.sex]; % 0 = male, 1 = female

% bin from ages 5-22
edges = 5:22; 

histFemale = histcounts( ageData(logical(sexData)), edges);
histMale = histcounts( ageData(~logical(sexData)), edges); 

% previsualize

bar(edges(1:end-1), [histFemale; histMale]', 'stacked')
legend(["Female", "Male"])
ylabel("# of subjects")
xlabel("age (y)")
title("HBN-1388: Age/Sex Distribution")

% export for prism visualization
T = array2table([edges(1:end-1); histFemale; histMale], 'RowNames', ["age", "Female", "Male"]);
writetable(T, fullfile('out', 'fig_demographics_hbn1388_ageSexHistogram.csv'))

%% HBN-865 Age/Sex Distribution
clearvars -except hm_data

% index subject with CBCL values
cbclIndex = find(arrayfun(@(hm_data) ~isempty(hm_data.CBCL_Total_T), hm_data));


% age histogram with sex 
ageData = [hm_data(cbclIndex).age];
sexData = [hm_data(cbclIndex).sex]; % 0 = male, 1 = female

% bin from ages 5-22
edges = 5:22; 

histFemale = histcounts( ageData(logical(sexData)), edges);
histMale = histcounts( ageData(~logical(sexData)), edges); 

% previsualize

bar(edges(1:end-1), [histFemale; histMale]', 'stacked')
legend(["Female", "Male"])
ylabel("# of subjects")
xlabel("age (y)")
title("HBN-865: Age/Sex Distribution")

% export for prism visualization
T = array2table([edges(1:end-1); histFemale; histMale], 'RowNames', ["age", "Female", "Male"]);
writetable(T, fullfile('out', 'fig_demographics_hbn865_ageSexHistogram.csv'))

%% HBN-865 CBCL/Sex Distribution
clearvars -except hm_data

% index subject with CBCL values
cbclIndex = find(arrayfun(@(hm_data) ~isempty(hm_data.CBCL_Total_T), hm_data));

cbclData = [hm_data(cbclIndex).CBCL_Total_T]; 
sexData = [hm_data(cbclIndex).sex]; % 0 = male, 1 = female
cbclID = [hm_data(cbclIndex).id];

edges = [20:5:90] - 2.5; 
centers = edges(1:end-1) + 2.5;

histFemale = histcounts( cbclData(logical(sexData)), edges);
histMale = histcounts( cbclData(~logical(sexData)), edges); 

bar(centers, [histFemale; histMale]', 'stacked')
xline(65, "--k");
legend(["Female", "Male"])

ylabel("# of subjects")
xlabel("CBCL Total T")
title("HBN-865: CBCL Total T/Sex Distribution")
% export for prism visualization
T = table(centers', histFemale', histMale', 'VariableNames', ["CBCL Total T", "Female", "Male"]);
writetable(T, fullfile('out', 'fig_demographics_hbn865_cbclTotalTscoreSexHistogram.csv'))

%% Motion grouping numbers 
% for each group: % female, mean age, mean CBCL
conditionList = ["Rest1", "Rest2", "MovieDM", "MovieDM_full"];
groupList = ["low", "med", "high"];

sex = array2table(zeros(4, 3), 'VariableNames', groupList, 'RowNames', conditionList);
age = array2table(zeros(4, 3), 'VariableNames', groupList, 'RowNames', conditionList);
ageSD = array2table(zeros(4, 3), 'VariableNames', groupList, 'RowNames', conditionList);
cbcl = array2table(zeros(4, 3), 'VariableNames', groupList, 'RowNames', conditionList);
cbclSD = array2table(zeros(4, 3), 'VariableNames', groupList, 'RowNames', conditionList);
n = array2table(zeros(4, 3), 'VariableNames', groupList, 'RowNames', conditionList);
formatTbl = array2table(string(zeros(4, 3)), 'VariableNames', groupList, 'RowNames', conditionList);

for conditionIdx = 1:length(conditionList) % loop over conditions in order of ConditionList
    for groupIdx = 1:length(groupList)
        idx = motionIdx(conditionIdx).( groupList(groupIdx) ); 
        
        % sex
        sex{conditionIdx, groupIdx} = sum([hm_data(idx).sex]); % 1 = female, sum returns n female
        % age 
        age{conditionIdx, groupIdx} = mean([hm_data(idx).age]);
        ageSD{conditionIdx, groupIdx} = std([hm_data(idx).age]);
        % cbcl
        cbcl{conditionIdx, groupIdx} = mean([hm_data(idx).CBCL_Total_T], 'omitnan');
        cbclSD{conditionIdx, groupIdx} = std([hm_data(idx).CBCL_Total_T], 'omitnan');
        % n 
        n{conditionIdx, groupIdx} = sum(idx);
        % formatTbl
        formatTbl{conditionIdx, groupIdx} = sprintf("N %i (%i F)\nAge %.1f ± %.1f\nCBCL %.1f ± %.1f,",...
                                                        n{conditionIdx, groupIdx},...
                                                        sex{conditionIdx, groupIdx},...
                                                        age{conditionIdx, groupIdx},...
                                                        ageSD{conditionIdx, groupIdx},...
                                                        cbcl{conditionIdx, groupIdx},...
                                                        cbclSD{conditionIdx, groupIdx});
    end
end


sex
ageSD
cbcl
cbclSD
n
formatTbl


writetable(sex, fullfile("out", "fig_demographics_motionGroups-sex.csv"), 'WriteRowNames', true)
writetable(age, fullfile("out", "fig_demographics_motionGroups-age.csv"), 'WriteRowNames', true)
writetable(cbcl, fullfile("out", "fig_demographics_motionGroups-cbcl.csv"), 'WriteRowNames', true)
writetable(n, fullfile("out", "fig_demographics_motionGroups-n.csv"), 'WriteRowNames', true)
writetable(formatTbl, fullfile("out", "fig_demographics_motionGroups-formatTbl.csv"), 'WriteRowNames', true)

    



