%% hm_preprocess.m
% Simon Frew | NNL | BCCHRI
% calculate FD, mean FD, and component FD according to Power et al. 2012
% generate motion groupings
fprintf("hm_preprocess.m\n")

%% Load hm_data (output from hm_import.m)
fprintf("\t 1/4\tloading hm_data.mat\n")
load hm_data.mat

%% calculate Powers: fd, mean fd, mean component fd for 355 volumes
fprintf("\t 2/4\tcalculating fd, mean fd, components\n")

% ONLY FIRST 375 VOLUMES FOR MOVIEDM
% omit first and last 10 volumes (e.g., only volumes 11-365)
for i = 1:length(hm_data)
    if class(hm_data(i).rest1) == 'table'
        [hm_data(i).r1_fd, hm_data(i).r1_mean_fd, hm_data(i).r1_fd_components, hm_data(i).r1_mean_fd_components] = fd_power(hm_data(i).rest1(11:365, :));
    end   
    if class(hm_data(i).rest2) == 'table'
        [hm_data(i).r2_fd, hm_data(i).r2_mean_fd, hm_data(i).r2_fd_components, hm_data(i).r2_mean_fd_components] = fd_power(hm_data(i).rest2(11:365, :));
    end
    if class(hm_data(i).movieDM) == 'table'
        [hm_data(i).mDM_fd, hm_data(i).mDM_mean_fd, hm_data(i).mDM_fd_components, hm_data(i).mDM_mean_fd_components] = fd_power(hm_data(i).movieDM(11:365, :));
    end
    if class(hm_data(i).movieDM) == 'table'
        [hm_data(i).mDM_full_fd, hm_data(i).mDM_full_mean_fd, hm_data(i).mDM_full_fd_components, hm_data(i).mDM_full_mean_fd_components] = fd_power(hm_data(i).movieDM(11:740, :));
    end
end

% calculate number of FD spikes
table_vars = ["r1", "r2", "mDM", "mDM_lasthalf"];
fd_thresh = 0.3;
for i = 1:length(hm_data)
    hm_data(i).fd_spikes = array2table([sum([hm_data(i).r1_fd] > fd_thresh), sum([hm_data(i).r2_fd] > fd_thresh), sum([hm_data(i).mDM_fd] > fd_thresh), sum([hm_data(i).mDM_full_fd(356:729)] > fd_thresh)], 'VariableNames', table_vars);
end

%% Generate motion grouping structure 
fprintf("\t 3/4\tgenerating motion groups...\n")

motionIdx(1).condition = "Rest1"; 
motionIdx(1).high = [hm_data.r1_mean_fd] > 0.3;
motionIdx(1).med = [hm_data.r1_mean_fd] < 0.3 & 0.15 <= [hm_data.r1_mean_fd];
motionIdx(1).low = [hm_data.r1_mean_fd] < 0.15;
motionIdx(1).nhigh = sum(motionIdx(1).high); 
motionIdx(1).nmed = sum(motionIdx(1).med); 
motionIdx(1).nlow = sum(motionIdx(1).low); 

motionIdx(2).condition = "Rest2"; 
motionIdx(2).high = [hm_data.r2_mean_fd] > 0.3;
motionIdx(2).med = [hm_data.r2_mean_fd] < 0.3 & 0.15 <= [hm_data.r2_mean_fd];
motionIdx(2).low = [hm_data.r2_mean_fd] < 0.15;
motionIdx(2).nhigh = sum(motionIdx(2).high); 
motionIdx(2).nmed = sum(motionIdx(2).med); 
motionIdx(2).nlow = sum(motionIdx(2).low); 

motionIdx(3).condition = "MovieDM"; 
motionIdx(3).high = [hm_data.mDM_mean_fd] > 0.3;
motionIdx(3).med = [hm_data.mDM_mean_fd] < 0.3 & 0.15 <= [hm_data.mDM_mean_fd];
motionIdx(3).low = [hm_data.mDM_mean_fd] < 0.15;
motionIdx(3).nhigh = sum(motionIdx(3).high); 
motionIdx(3).nmed = sum(motionIdx(3).med); 
motionIdx(3).nlow = sum(motionIdx(3).low); 

motionIdx(4).condition = "MovieDM_full"; 
motionIdx(4).high = [hm_data.mDM_full_mean_fd] > 0.3;
motionIdx(4).med = [hm_data.mDM_full_mean_fd] < 0.3 & 0.15 <= [hm_data.mDM_full_mean_fd];
motionIdx(4).low = [hm_data.mDM_full_mean_fd] < 0.15;
motionIdx(4).nhigh = sum(motionIdx(4).high); 
motionIdx(4).nmed = sum(motionIdx(4).med); 
motionIdx(4).nlow = sum(motionIdx(4).low); 





%% save as hm_analysis.mat
clearvars -except hm_data motionIdx
fprintf("\t 4/4\tsaving hm_analysis.mat\n")

save("hm_analysis.mat")


%% helper functions

function [fd, mean_fd, fd_components, mean_fd_components] = fd_power(data)
% [fd, mean_fd, fd_components, mean_fd_components] = fd_power(data)
% calculate fd and mean fd according to Power et al. 2012
% 
% Inputs: 
%    data: table with 6 axis of motion data (x_rad, y_rad, z_rad, x_mm, y_mm, z_mm)
% Outputs:
%    fd: array of fd values of length(data) - 1
%    mean_fd: the mean of the absolute values of fd
%    fd_components: components along 6 axes of motion data
%    mean_fd_components: table of mean of components along 6 axes of motion data
    data = data{:,:}; % convert to array
    table_vars = ["x_rot_mm", "y_rot_mm", "z_rot_mm", "x_mm", "y_mm", "z_mm"];
    radius = 50; % sphere of radius 50mm

    motion = [rad2km(data(:,1),radius), rad2km(data(:,2),radius), rad2km(data(:,3),radius), data(:, 4:6)];

    n = size(motion, 1); 
    fd = zeros(1, n-1);
    fd_components = zeros(n-1, 6);

    % loop through each data point and compute fd (based on prev point)
    % starting at the last volume and working backwards
    while n > 1
        fd_components(n-1, :)= abs(motion(n-1, :) - motion(n, :));
        fd(n-1) = sum(fd_components(n-1, :));
        n = n-1;
    end

    mean_fd = mean(fd);
    mean_fd_components = array2table(mean(fd_components, 1), 'VariableNames', table_vars);
end