%% fig_stimulusCorrelations.m
% Simon Frew | NNL | BCCHRI
% FD composition by axis

%% load / generate dataset 
% load hm_analysis.mat
% 
% % hm_fdTimeSeries.mat
% %   Contents: 
% %       data: 1x1 struct variable containing 4 fields
% %           data.r1_fd      : 1388 x 354 (subjects x volumes) matrix with raw FD values for Rest1 
% %           data.r2_fd      : 1388 x 354 (subjects x volumes) matrix with raw FD values for Rest2 
% %           data.mDM_fd     : 1388 x 354 (subjects x volumes) matrix with raw FD values for MovieDM 
% %           data.mDM_full_fd: 1388 x 729 (subjects x volumes) matrix with raw FD values for the untruncated MovieDM, referred to as MovieDM Full 
% 
% data.r1_fd = cat(1, hm_data.r1_fd);
% data.r2_fd = cat(1, hm_data.r2_fd);
% data.mDM_fd = cat(1, hm_data.mDM_fd);
% data.mDM_full_fd = cat(1, hm_data.mDM_full_fd);
% 
% save("hm_fdTimeSeries.mat", data)
% load("hm_fdTimeSeries.mat")

%% generate 
isc.vol15 = hm_isc(15); 
save('hm_fdTimeSeries_multipleVols.mat')

%% load
load('hm_fdTimeSeries_multipleVols.mat')

%% visualize


hm_stimFig(isc.vol15, data, 1, 0, 0, 1)



%% video vis 
hm_stimFig(vol15, data, 1, 1)


%% figure generation
% select 15 volumes for final window 
% keep all data in fisher Z values
out = rmfield(isc.vol15, "mDM_full_fd"); 
outVols = "vol15"; 
windowSize = 15;
numberOfWindows = length(out.mDM_fd);
% export ISC data by condition for PRISM 

writetable(struct2table(out), fullfile("out", sprintf("fig_stimulusCorrelations-rawISCbyCondition-%s.csv", outVols)))


% save figure 

% hm_stimFig(out, data, 0, 1, 1, 1, 1)

% export ISC data vs. windowed mean FD for scatterplot 
tmpFDMeasure = zeros(numberOfWindows, 1); 
for windowIdx = 1:numberOfWindows
    volIdx = windowIdx:(windowIdx + windowSize - 1); 
    tmpWindowFD = data.mDM_fd(:, volIdx); 
    tmpFDMeasure(windowIdx) = mean(tmpWindowFD, "all");
end

writetable(array2table([tmpFDMeasure, out.mDM_fd], 'VariableNames', ["Windowed Mean FD", "ISC"]), fullfile("out", sprintf("fig_stimulusCorrelations-MovieHiscVsFD-%s.csv", outVols)))

% max/min values 
fieldList = string(fields(out)); 
for i = 1:length(fields(out))
    fprintf("%s: \n\tMean, \t\tMax, \t\tMin\n", fieldList(i))
     disp([mean(out.(fieldList(i))), max(out.(fieldList(i))), min(out.(fieldList(i)))])
end


%% correlation function 

hm_corrISCFD(isc, data, 0, 0)

hm_stimFig(out, data, 0, 0, 0, 0)


%% functions 


function [] = hm_stimFig(out, data, r, saveFig, timepoints, playvideo, saveFrame)
    % INPUTS
    %   out         :   structure containing ISC values as generated by hm_isc
    %   data        :   structure containing mean FD values 
    %   r           :   logical flag, if 1, convert ISCs to R values with tanh
    %   saveFig     :   logical flag, if 1, save figure as .png in subdirectory 'out'
    %   timepoints  :   logical flag, if 1, print movie time points into console
    %   playvideo   :   logical flag, if 1, load movie from dropbox and play top 10 window ISCs
    %   saveFrame   :   logical flag, if 1, save first frame of each window
    
    % input checking
    if ~exist('saveFig', 'var') % 
        saveFig = 0; 
    end
    if ~exist('timepoints', 'var')
        timepoints = 0; 
    end
    if ~exist('playvideo', 'var')
        playvideo = 0; 
    end
    if ~exist('saveFrame', 'var')
        saveFrame = 0; 
    end
    
    % set convert from fisher Z to R by default
    if r == 1
        hm_convertZtoR(out);
        yLabels = "R Value";
    else
        yLabels = "Z Value";
    end
    
    % plot figure 
    windowSize = size(data.r1_fd, 2) - size(out.r1_fd, 1) + 1; 
    figure
    subplot(2,1,1)
        hold on 
        plot([(out.r1_fd), (out.r2_fd)]); 
        findpeaks(out.mDM_fd, 'SortStr', 'descend', 'npeaks', 10);
        title(sprintf("R Value Timeseries\n %i volumes per window", windowSize))
        legend(["Rest1", "Rest2", "Movie-H", "Peaks"])
        ylim([0, 0.04])
        xlim([1, 354])
        ylabel(yLabels)
        xlabel("Windows")
        
    [~, locs] = findpeaks(out.mDM_fd, 'SortStr', 'descend', 'npeaks', 10);
    volLocs = [locs, locs + windowSize];

    subplot(2,1,2)
        plot(seconds(([1:354]+10).*0.8), mean(data.mDM_fd), 'DurationTickFormat', 'mm:ss')
        ylabel("Mean FD (mm)")
        xlabel("Movie Time (s)")
        title("Movie-H Volumewise Mean FD")
        xlim(seconds(([1,354]+10).*0.8))
        yline(mean(data.mDM_fd, 'all'), '-', 'mean fd');

        for i = 1:length(locs)
            
            
            tmplocFD = data.mDM_fd(:, volLocs(i, 1):volLocs(i, 2)-1);
            tmplocMeanFD = mean(data.mDM_fd(:, volLocs(i, 1):volLocs(i, 2)-1), 'all');
            
%             % sort by slope
%                 tmpSlopes = zeros(size(tmplocFD, 1), 2);
%                 for j = 1:size(tmplocFD, 1)
%                     tmpSlopes(j, :) = polyfit(1:windowSize, tmplocFD(j, :), 1);
%                 end
%                 tmpMeanSlopes = mean(tmpSlopes(:, 1));
%                 if tmpMeanSlopes < 0
%                     patchColour = 'b';
%                 else 
%                     patchColour = 'r';
%                 end   
                
            % sort by mean FD             
                if tmplocMeanFD < mean(data.mDM_fd, 'all')
                    patchColour = 'b';
                else 
                    patchColour = 'r';
                end
            % create patch
            patch( ([volLocs(i, 1), volLocs(i, 1), volLocs(i, 2), volLocs(i, 2)] + 10).*0.8, [flip(ylim), ylim] , patchColour, 'FaceAlpha', 0.2, 'EdgeAlpha', 0)
        end
    
    % optional: save out figure 
    if saveFig
        saveas(gcf, fullfile("out", sprintf("fig_stimulusCorrelations-previs-vol%i.pdf", windowSize)))
    end
            
    % optional: calculate timepoints 
    if timepoints
        timepoints = table('Size', [length(locs), 3], 'variableTypes', ["double", "string", "string"], 'variableNames', ["WindowIdx", "Start", "End"]);
        timepoints{:, 1} = locs;

        for i = 1:length(locs)
            timepoints{i, 2:3} = string(datestr((seconds(volLocs(i, :) + 10 )* 0.8), 'MM:ss.fff'))'; % +10 for truncation
        end

        sortrows(timepoints)
    end
    
    % optional: play video at specified timepoints
    if playvideo
        vidObj = VideoReader("C:\Users\simon\Dropbox\Vanderwal Lab\Studies\HBN\stimuli\DespicableMe_10min_Eng.avi");
        fps = round(get(vidObj, 'FrameRate'));
        sortLocs = sort(volLocs);

        for i = 1:length(locs)
            implay(read(vidObj, fps * (sortLocs(i, :) + 10 )* 0.8), fps) % +10 for truncation
            
            if saveFrame
                imwrite(read(vidObj, fps * (sortLocs(i, 1) + 10 )* 0.8), fullfile("out", sprintf("fig_stimulusCorrelations_videoFrame-%i.png", fps * (sortLocs(i, 1) + 10 )* 0.8)))
            end
            
        end        
    end
    
    
end


function [] = hm_corrISCFD(isc, data, r, slope)
%     correlate ISC with windowed mean FD 
%     INPUT: 
%         isc     :   structure of outputs from the function hm_isc
%         data    :   raw FD data 
%         r       :   logical flag, if 1, iscs will be converted to R value
%         slope   :   logical flag, if 1, correlations will be computed against window slope, otherwise window mean FD
    
    figure
    iscFields = string(fieldnames(isc));
    for iscIdx = 1:length(iscFields)
        tmpField = iscFields(iscIdx);
        
        windowSize = size(data.r1_fd, 2) - size(isc.(tmpField).r1_fd, 1) + 1; 
        numberOfWindows = size(isc.(tmpField).r1_fd, 1);
        
        if r == 1
            tmpIsc = tanh(isc.(tmpField).mDM_fd); 
            yLabel = "ISC (R values)";
        else
            tmpIsc = isc.(tmpField).mDM_fd;
            yLabel = "ISC (Z values)";
        end
        
        tmpFDMeasure = zeros(numberOfWindows, 1); 

        % for each window... 
        for windowIdx = 1:numberOfWindows
            volIdx = windowIdx:(windowIdx + windowSize - 1); 
            tmpWindowFD = data.mDM_fd(:, volIdx); 
            
            if slope == 1
                tmpSlopes = zeros(size(tmplocFD, 1), 2);
                for j = 1:size(tmplocFD, 1)
                    tmpSlopes(j, :) = polyfit(1:windowSize, tmpWindowFD, 1);
                end
                tmpFDMeasure(windowIdx) = mean(tmpSlopes(:, 1));
            else
                tmpFDMeasure(windowIdx) = mean(data.mDM_fd(:, volIdx), "all");
            end
        end
        
        
        subplot(4, 3, iscIdx)
        scatter(tmpFDMeasure, tmpIsc)
        title(tmpField)
        ylabel(yLabel)
        xlabel("windowed mean FD")
        xlim([0.3, 0.7])
        ylim([0, 0.03])
    end
    
end



function out = hm_convertZtoR(out)
    
    fields = string(fieldnames(out));
    for i = 1:length(fields)
        out.(fields(i)) = tanh(out.(fields(i))); 
    end

end

function out = hm_isc(windowLength)
    if ~exist('windowLength','var')
       windowLength = 5; 
    end
    fprintf('Using sliding window of %d timepoints\n',windowLength);

    load('hm_fdTimeSeries.mat'); % data
    fields = fieldnames(data);
    % loop over each field (task)
    for i=1:numel(fields) -1 % to skip movie-F
       fprintf('\t%d/%d\t%s\n',i,numel(fields),fields{i});
       mat = data.(fields{i});
       nOut = size(mat,2)-windowLength+1;
       out.(fields{i}) = nan(nOut, 1);
       % loop over each volume (timepoint)
       for n=1:nOut
           r = corr(mat(:,n:n+windowLength-1)');
           % select lower triangular 
           r = r( logical( tril( ones(size(r, 1)), -1) ) ); 
           % z-transform average r-values
           rz = atanh(r);
           rzavg = mean(rz(rz~=Inf),'all');
           rzstd = std(rz(rz~=Inf), 0, 'all');
           out.(fields{i})(n) = rzavg;
       end
    end
end