%% fig_frequency.m
% Simon Frew | NNL | BCCHRI
% frequency analysis -> fair et al. 

load hm_analysis.mat

%% Rest 1 FFT
% build index to sort by r1 mean fd
condition_array = [[1:length(hm_data)]', [hm_data.r1_mean_fd]'];
condition_index = sortrows(condition_array, 2, 'ascend'); % low to high index
names = {'pitch', 'roll', 'yaw', 'x', 'y', 'z'};
order = [4, 5, 6, 1, 2, 3]; % order to plot

hp = 0.0; % no high pass filter
fftdata = zeros(257, 1388, 6); % initialize matrix

% calculate fft using pwelch
for i = 1:length(condition_index)
    [t, fftdata(:, i, :)] = hm_fft(hm_data(condition_index(i, 1)).rest1{11:355, :}, 800, hp, 0);
end

sgtitle(["Power spectra of absolute displacement in Rest1" "HBN-1388, sorted by mean FD"])

for n = 1:6
    subplot(1, 8, order(n)+1)

    imagesc(t, [], normcdf(fftdata(:, :, n))')
    
    title(string(names(n)))
    xlabel("Hz")
    set(gca,'ytick',[])
    yticklabels({})
end

colormap parula

% plot mean FD
subplot(1, 8, 1)
plot(condition_index(:, 2), 1:length(condition_index), '-k')
set(gca, 'Ydir', 'reverse')
set(gca, 'Xdir', 'reverse')

yticklabels({})

yline(427, '-r', '', 'LineWidth', 2);
text(4, 428/2,' Low');
yline(427+432, '-r', '', 'LineWidth', 2);
text(4, 428+(432)/2,' Med');
text(4, 530/2+427+432,' High');


title("Rest1 Mean FD")
xlabel("Mean FD (mm)")
ylabel(["HBN-1388: Sorted low to high movers", ""])

% plot colorbar
sp7 = get(subplot(1,8,7), 'position');
sp6 = get(subplot(1,8,6), 'position');
sp_space = sp7(1) - sp6(1) - sp7(3);
colorbar_location = [sp7(1) + sp7(3) + sp_space, sp7(2), 0.0100, sp7(4)];
%text_location = colorbar_location + [colorbar_location(1) + colorbar_location(3), [sp_7(1) + 2*sp_7(3) + sp_space, sp_7(2), 0.0100, sp_7(4)];

sp8 = subplot(1,8,8);

imagesc([], [1, 0], linspace(1, 0, 1024)');
set(gca, 'Ydir', 'normal')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
ylabel("Relative Spectral Power (%)")

sp8.Position = colorbar_location;

%% MovieDM FFT
figure
% build index to sort by r1 mean fd
condition_array = [[1:length(hm_data)]', [hm_data.mDM_mean_fd]'];
condition_index = sortrows(condition_array, 2, 'ascend'); % low to high index
names = {'pitch', 'roll', 'yaw', 'x', 'y', 'z'};
order = [4, 5, 6, 1, 2, 3]; % order to plot

hp = 0.0; % no high pass filter
fftdata = zeros(257, 1388, 6); % initialize matrix

% calculate fft using pwelch
for i = 1:length(condition_index)
    [t, fftdata(:, i, :)] = hm_fft(hm_data(condition_index(i, 1)).movieDM{11:355, :}, 800, hp, 0);
end

sgtitle(["Power spectra of absolute displacement in MovieDM" "HBN-1388, sorted by mean FD"])

for n = 1:6
    subplot(1, 8, order(n)+1)

    imagesc(t, [], normcdf(fftdata(:, :, n))')
    
    title(string(names(n)))
    xlabel("Hz")
    set(gca,'ytick',[])
    yticklabels({})
end

colormap parula

% plot mean FD
subplot(1, 8, 1)
plot(condition_index(:, 2), 1:length(condition_index), '-k')
set(gca, 'Ydir', 'reverse')
set(gca, 'Xdir', 'reverse')

yticklabels({})

yline(427, '-r', '', 'LineWidth', 2);
text(4, 428/2,' Low');
yline(427+432, '-r', '', 'LineWidth', 2);
text(4, 428+(432)/2,' Med');
text(4, 530/2+427+432,' High');


title("MovieDM Mean FD")
xlabel("Mean FD (mm)")
ylabel(["HBN-1388: Sorted low to high movers", ""])

% plot colorbar
sp7 = get(subplot(1,8,7), 'position');
sp6 = get(subplot(1,8,6), 'position');
sp_space = sp7(1) - sp6(1) - sp7(3);
colorbar_location = [sp7(1) + sp7(3) + sp_space, sp7(2), 0.0100, sp7(4)];
%text_location = colorbar_location + [colorbar_location(1) + colorbar_location(3), [sp_7(1) + 2*sp_7(3) + sp_space, sp_7(2), 0.0100, sp_7(4)];

sp8 = subplot(1,8,8);

imagesc([], [1, 0], linspace(1, 0, 1024)');
set(gca, 'Ydir', 'normal')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
ylabel("Relative Spectral Power (%)")

sp8.Position = colorbar_location;
