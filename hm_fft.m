function [tt, y] = hm_fft(data, TR, hp, fig)
% [tt, y] = hm_fft(data, TR, fig)
% calculates fft of basic head motion to identify underlying frequencies (or noise)
% returns Hz for plotting 
% 
% Replicates: 
%	Fair, D. A., Miranda-Dominguez, O., Snyder, A. Z., Perrone, A., Earl, E. A., Van, A. N., Koller, J. M., Feczko, E., Tisdall, M. D., van der Kouwe, A., Klein, R. L., Mirro, A. E., Hampton, J. M., Adeyemo, B., Laumann, T. O., Gratton, C., Greene, D. J., Schlaggar, B. L., Hagler, D. J., … Dosenbach, N. U. F. (2020). Correction of respiratory artifacts in MRI head motion estimates. NeuroImage, 208, 116400. https://doi.org/10.1016/j.neuroimage.2019.116400
% 
%
% inputs:
%   data = n by 1 column vector of head motion data with n volumes and m axes.
%   TR = MRI TR sampling rate (duration between volumes) in milliseconds
%   hp = Value for high pass filter in Hz % depreciated
%   fig = 1: plot values, 0: do not plot
% returns: 
%   tt = freq domain scale to plot points (Hz)
%   y = Fourier transform of input data (either vector or columns)

TR = TR / 1000;         % MRI TR: convert from ms to s
Fs = (TR)^-1;           % sampling frequency (Hz)

[y, tt] = pmtm(data, 8, 512, Fs);

% recreating paper: 10log10, then zscore
y = 10*log10(y);
y = zscore(y);

%y(tt < hp, :) = 0;
if fig
    subplot(3,1,1)
    plot(t, data)%, 10, 'MarkerFaceColor', 'b', 'MarkerFaceAlpha', 0.075, 'MarkerEdgeAlpha', 0)
    subplot(3,1,2)
    plot(f, P1)%, 10,'MarkerFaceColor', 'b', 'MarkerFaceAlpha', 0.075, 'MarkerEdgeAlpha', 0)
    subplot(3,1,3)
    plot(tt, Y)%, 10,'MarkerFaceColor', 'b', 'MarkerFaceAlpha', 0.075, 'MarkerEdgeAlpha', 0)
end

end


