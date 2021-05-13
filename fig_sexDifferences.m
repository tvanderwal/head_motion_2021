%% fig_sexDifferences.m
% Simon Frew | NNL | BCCHRI
% sex differences
load hm_analysis.mat


%% % 0 = male, 1 = female

% mean fd grouped by sex 
figure 
sgtitle("Mean FD grouped by sex (HBN-1388)")
subplot(1,3,1)
boxplot([hm_data.r1_mean_fd], {[hm_data.sex]})
xticklabels(["male", "female"])
ylabel("mean fd")
xlabel("Rest1")

subplot(1,3,2)
boxplot([hm_data.r2_mean_fd], {[hm_data.sex]})
xticklabels(["male", "female"])
ylabel("mean fd")
xlabel("Rest2")

subplot(1,3,3)
boxplot([hm_data.mDM_mean_fd], {[hm_data.sex]})
xticklabels(["male", "female"])
ylabel("mean fd")
xlabel("MovieDM")

%% Export age, cbcl grouped by sex
age = [hm_data.age]; 
sex = logical([hm_data.sex]); 
cbcl = cat(1, hm_data.CBCL_Total_T); 
cbclIndex = find(arrayfun(@(hm_data) ~isempty(hm_data.CBCL_Total_T), hm_data));

out = nan(length(age), 4);

%age
out(1:sum(sex), 1) = age(sex); 
out(1:sum(~sex), 2) = age(~sex); 
%cbcl
out(1:sum(sex(cbclIndex)), 3) = cbcl(sex(cbclIndex)); 
out(1:sum(~sex(cbclIndex)), 4) = cbcl(~sex(cbclIndex)); 

out = array2table(out, 'VariableNames', ["Age-F", "Age-M", "CBCL-F", "CBCL-M"]);

writetable(out, fullfile("out", "fig_sexDifferences-demos.csv"))


%% Export mean FD grouped by sex 
condList = ["r1_mean_fd", "r2_mean_fd", "mDM_mean_fd", "mDM_full_mean_fd"];




for condIdx = 1:length(condList)
    out = array2table(nan(max([ sum([hm_data.sex]), sum(~[hm_data.sex])]), 2), 'VariableNames', ["female", "male"]);
    
    out{ 1:sum([hm_data.sex]) , 1} = [ hm_data(logical([hm_data.sex])).(condList(condIdx)) ]'; % female
    out{ 1:sum(~[hm_data.sex]) , 2} = [ hm_data(logical(~[hm_data.sex])).(condList(condIdx)) ]'; % male
    
    
    writetable(out, fullfile("out", "fig_sexDifferences-" + condList(condIdx) + ".csv")); 
end


%% Overall FD composition by sex
exportAxes = ["x", "y", "z", "pitch", "roll", "yaw"]; 

condTitle = ["Rest1", "Rest2", "MovieDM", "MovieDM Full"];
condList = ["r1_mean_fd_components", "r2_mean_fd_components", "mDM_mean_fd_components", "mDM_full_mean_fd_components"];
sexList = ["male", "female"];

outTable = array2table(zeros(length(condList), 6), 'RowNames', condTitle, 'VariableNames', exportAxes);


for idx = 0:1 % 0 = male, 1 = female
    sexIdx = logical([hm_data.sex]) == idx; 
    
    for condIdx = 1:length(condList)

        out = table2array(cat(1, hm_data(logical([hm_data.sex])).(condList(condIdx)))); 
        out = mean(out ./ sum(out, 2));
        out = out(:, [4,5,6,1,2,3]);

        outTable{condIdx, :} = out;     

    end
    
    writetable(outTable, fullfile("out", sprintf("fig_sexDifferences-composition-full-%s.csv", sexList(idx+1))), 'WriteRowNames', true)
end 


%% Overall FD spike composition by sex
exportAxes = ["x", "y", "z", "pitch", "roll", "yaw"]; 

condTitle = ["Rest1", "Rest2", "MovieDM", "MovieDM Full"];
condList = ["r1_mean_fd_components", "r2_mean_fd_components", "mDM_mean_fd_components", "mDM_full_mean_fd_components"];
sexList = ["male", "female"];

outTable = array2table(zeros(length(condList), 6), 'RowNames', condTitle, 'VariableNames', exportAxes);


for idx = 0:1 % 0 = male, 1 = female
    sexIdx = logical([hm_data.sex]) == idx; 
    
    for condIdx = 1:length(condList)

        out = table2array(cat(1, hm_data(logical([hm_data.sex])).(condList(condIdx)))); 
        out = mean(out ./ sum(out, 2));
        out = out(:, [4,5,6,1,2,3]);

        outTable{condIdx, :} = out;     

    end
    
    writetable(outTable, fullfile("out", sprintf("fig_sexDifferences-composition-full-%s.csv", sexList(idx+1))), 'WriteRowNames', true)
end 
%%
thresh = 0.3;
spikeIdx.group = "all";
spikeIdx.Rest1 = arrayfun(@(sub) sub.r1_fd > thresh, hm_data, 'UniformOutput', false);
spikeIdx.Rest2 = arrayfun(@(sub) sub.r2_fd > thresh, hm_data, 'UniformOutput', false);
spikeIdx.MovieDM = arrayfun(@(sub) sub.mDM_fd > thresh, hm_data, 'UniformOutput', false);

condition = ["Rest1", "MovieDM"];
conditionMap = [1,3];
dataFields = ["r1_fd_components", "mDM_fd_components"]; 
axesName = ["x", "y", "z", "pitch", "roll", "yaw"]; 
sexList = ["male", "female"];

outTable = array2table(zeros(length(condition), 6), 'RowNames', condition, 'VariableNames', axesName);

for sexIdx = 0:1 % 0 = male, 1 = female
    for conditionIdx = 1:length(condition)
        
        fd = cat(3, hm_data( logical([hm_data.sex]) == sexIdx ).(dataFields(conditionIdx))); 
        
        idx = spikeIdx.(condition(conditionIdx)) ( logical([hm_data.sex]) == sexIdx );
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
        
        % permute columns for visualization
        subOut = subOut(:, [4,5,6,1,2,3]);
        
        % calculate mean composition
        subOut = mean(subOut); 
        
        outTable{conditionIdx, :} = subOut; 
        
    end
    
    writetable(outTable, fullfile("out", sprintf("fig_sexDifferences-composition-spike-%s.csv", sexList(sexIdx+1))), 'WriteRowNames', true)

end

%% export subject level mean spike and overall FD composition for t-testing (Fig 7C)
thresh = 0.3;
spikeIdx.group = "all";
spikeIdx.Rest1 = arrayfun(@(sub) sub.r1_fd > thresh, hm_data, 'UniformOutput', false);
spikeIdx.Rest2 = arrayfun(@(sub) sub.r2_fd > thresh, hm_data, 'UniformOutput', false);
spikeIdx.MovieDM = arrayfun(@(sub) sub.mDM_fd > thresh, hm_data, 'UniformOutput', false);

condition = ["Rest1", "MovieDM"];
conditionMap = [1,3];
dataFields = ["r1_fd_components", "mDM_fd_components"]; 
condList = ["r1_mean_fd_components", "mDM_mean_fd_components"];
axesName = ["x", "y", "z", "pitch", "roll", "yaw"]; 
sexList = ["male", "female"];

% outTable = array2table(zeros(length(condition), 6), 'RowNames', condition, 'VariableNames', axesName);

for sexIdx = 0:1 % 0 = male, 1 = female
    for conditionIdx = 1:length(condition)
        
        fd = cat(3, hm_data( logical([hm_data.sex]) == sexIdx ).(dataFields(conditionIdx))); 
        
        idx = spikeIdx.(condition(conditionIdx)) ( logical([hm_data.sex]) == sexIdx );
        idx = cat(1, idx{:,:});
        
        % non spike composition
        
        out = table2array(cat(1, hm_data(logical([hm_data.sex] == sexIdx )).(condList(conditionIdx)))); 
        out = out ./ sum(out, 2);
        out = out(:, [4,5,6,1,2,3]);
        writetable(array2table(out, 'VariableNames', axesName), fullfile("out", sprintf("fig_sexDifferences-composition-subj-overall-%s-%s.csv", sexList(sexIdx+1), condition(conditionIdx))))

        
        
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
        
        % permute columns for visualization
        subOut = subOut(:, [4,5,6,1,2,3]);
        
        writetable(array2table(subOut, 'VariableNames', axesName), fullfile("out", sprintf("fig_sexDifferences-composition-subj-spike-%s-%s.csv", sexList(sexIdx+1), condition(conditionIdx))))

    end
    

end


%% 
% mean FD grouped by sex and age <12
age = categorical(floor([hm_data.age])); 
sex = renamecats(categorical([hm_data.sex]), ["0", "1"], ["M", "F"]);

figure
sgtitle("Mean FD sex differences by year (HBN-1388)")

subplot(2,1,1)
boxplot([hm_data.r1_mean_fd], {age, sex})
xlabel("Rest1")
ylabel("mean fd")

subplot(2,1,2)
boxplot([hm_data.mDM_mean_fd], {age, sex})
xlabel("MovieDM")
ylabel("mean fd")

%% anova with age (cts), sex 

[~, tbl.r1, stats] = anovan([hm_data.r1_mean_fd], {[hm_data.age], [hm_data.sex]}, 'continuous', [1], 'model','interaction','varnames',{'age','sex'});
table(stats.coeffnames, stats.coeffs);

[~, tbl.r2, stats] = anovan([hm_data.r2_mean_fd], {[hm_data.age], [hm_data.sex]}, 'continuous', [1], 'model','interaction','varnames',{'age','sex'});
table(stats.coeffnames, stats.coeffs);

[~, tbl.mDM, stats] = anovan([hm_data.mDM_mean_fd], {[hm_data.age], [hm_data.sex]}, 'continuous', [1], 'model','interaction','varnames',{'age','sex'});
table(stats.coeffnames, stats.coeffs);

[~, tbl.mDM_full, stats] = anovan([hm_data.mDM_full_mean_fd], {[hm_data.age], [hm_data.sex]}, 'continuous', [1], 'model','interaction','varnames',{'age','sex'});
table(stats.coeffnames, stats.coeffs);

outTbl = [cell2table(tbl.r1);cell2table(tbl.r2);cell2table(tbl.mDM);cell2table(tbl.mDM_full);];

writetable(outTbl, fullfile("out", "fig_sexDifferences_ANOVA-r1-r2-mDMH-mDMF.csv"))

%% export for anova 

writetable(hm_table(:, [2:4, 9, 13, 17, 21]), fullfile("out", "demos-for-R.csv"))


%% anova with age, sex, psychopathology 
cbclIndex = find(arrayfun(@(hm_data) ~isempty(hm_data.CBCL_Total_T), hm_data));
fd = [hm_data(cbclIndex).r1_mean_fd];
age = [hm_data(cbclIndex).age];
ageCat = age < 12; 
sex = [hm_data(cbclIndex).sex];
cbcl = [hm_data(cbclIndex).CBCL_Total_T] > 65; % 1: clinical, 0: non-clinical

anovan(fd, {age, sex, cbcl}, 'continuous', [1], 'model', 'interaction', 'varnames',{'age','sex', 'cbcl'})

boxplot(fd, {floor(age), sex, cbcl})
vartestn(fd, {ageCat, sex, cbcl})

%% anova with age, sex, psychopathology 
cbclIndex = find(arrayfun(@(hm_data) ~isempty(hm_data.CBCL_Total_T), hm_data));
fd = [hm_data(cbclIndex).mDM_mean_fd];
age = [hm_data(cbclIndex).age];
sex = [hm_data(cbclIndex).sex];
cbcl = [hm_data(cbclIndex).CBCL_Total_T]; % 1: clinical, 0: non-clinical

anovan(fd, {age, sex, cbcl}, 'continuous', [1, 3], 'model', 1, 'varnames',{'age','sex', 'cbcl'})

