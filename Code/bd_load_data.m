% bd_load_data.m
% Script for Multi-Site Frequency Drift 3T Phantom Project
% Richard Edden's Lab
% Department of Radiology and Radiological Science, The Johns Hopkins University School of Medicine 2020
%
% USAGE:
% Rename data file 
%   - GE Example     : G01a_day1_post.7; G01a_day1_pre.7; G01a_day2_post.7; G01a_day2_post.7
%   - Philips Example: P01a_days1_post_raw_act.sdat, P01a_days1_post_raw_ref.spar, 
%                      P01a_days1_pre_raw_act.sdat, P01a_days1_pre_raw_ref.spar
%   - Siemens Exmple : S01a_Day1_PRESS_Post.dat and S01a_Day1_PRESS_Pre.dat
% Run this script at the data directory
%
% DESCRIPTION:
% This code simply extract water peak from each individual transient

%% Load and plot data from Big Drift project

clear; clc;
outDir        = '~/Desktop';
showModelFits = 0;
printPlot     = 0;

[~,siteID] = fileparts(pwd); % assumes the name of the current directory is the site ID, e.g., 'S01a'

% Automatically elect all the MRS files in the current directory
switch siteID(1)
    case 'G'
        fnames = SelectMRSFiles('7');
    case 'P'
%        fnames      = SelectMRSFiles('data');
         fnames      = SelectMRSFiles('sdat');
        spar_fnames = dir('*.spar');
        spar_fnames = {spar_fnames.name};
    case 'S'
        fnames = SelectMRSFiles('dat');
end

%%%%%%%%%%%%%% CHANGE THIS FOR EACH SITE %%%%%%%%%%%%%%
if length(fnames) == 2
    fnames = fnames([2 1]); % reorder the filenames so they are pre/post (for visualization purposes)
else
    fnames = fnames([2 1 4 3]); % reorder the filenames so they are day 1/2 and pre/post (for visualization purposes)
end

% For Philips data only
% Also reorder the spar filenames so they match the order of the .data files
if strcmp(siteID(1),'P')
    if length(fnames) == 2
        spar_fnames = spar_fnames([2 1]);
    else
        spar_fnames = spar_fnames([2 1 4 3]);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load data

data = cell(1, length(fnames));
for ii = 1:length(fnames)
    switch siteID(1)
        case 'G'
            data{ii} = bd_GERead(fnames{ii});
        case 'P'
            [~,~,ext] = fileparts(fnames{ii});
            switch lower(ext)
                case '.data'
                    data{ii} = bd_PhilipsReadRaw(fnames{ii}, spar_fnames{ii});
                case '.sdat'
                    data{ii} = bd_PhilipsRead(fnames{ii}, spar_fnames{ii});
            end
        case 'S'
            data{ii} = bd_SiemensTwixRead(fnames{ii});
    end
end

%% ECC, zero-fill and FFT

for ii = 1:length(fnames)
    data{ii}.fids = bd_EddyCurrentCorrection(data{ii}.fids, data{ii}.fids(:,1));
    zeroFillTo    = round(32768 / 2000 * data{ii}.p.sw);
    freqRange     = data{ii}.p.sw / data{ii}.p.LarmorFreq;
    data{ii}.freq = ((zeroFillTo + 1 - (1:zeroFillTo)) / zeroFillTo * freqRange + 4.8 - freqRange/2)';
    data{ii}.spec = fftshift(fft(data{ii}.fids, zeroFillTo, 1), 1);
end

%% Modeling

warning('off','stats:nlinfit:ModelConstantWRTParam');
warning('off','stats:nlinfit:IllConditionedJacobian');
warning('off','MATLAB:rankDeficientMatrix');

lsqopts  = optimset('lsqcurvefit');
lsqopts  = optimset(lsqopts,'MaxIter',800,'TolX',1e-4,'TolFun',1e-4,'Display','off');
nlinopts = statset('nlinfit');
nlinopts = statset(nlinopts,'MaxIter',400,'TolX',1e-7,'TolFun',1e-7,'FunValCheck','off');
a0 = max(real(data{1}.spec(:,1)));
x0 = [0   a0   4.8   -0.01   0   0];
lb = [0   0    4.55 -10    -10 -10];
ub = [Inf 2*a0 5.05   0     10  10];

pars = cell(1,length(fnames));
FWHM = cell(1,length(fnames));

for ii = 1:length(fnames)    
    freqLim = data{ii}.freq >= 4.8-0.35 & data{ii}.freq <= 4.8+0.35;    
    for jj = 1:size(data{ii}.spec,2)
        if jj == 1
            x0 = [0 a0 4.8 -0.01 0 0];
        end
        x0 = lsqcurvefit(@Voigt, x0, data{ii}.freq(freqLim), real(data{ii}.spec(freqLim,jj)), lb, ub, lsqopts);
        pars{ii}(jj,:) = nlinfit(data{ii}.freq(freqLim), real(data{ii}.spec(freqLim,jj)), @Voigt, x0, nlinopts);
        x0 = pars{ii}(jj,:);
        
        V   = Voigt([pars{ii}(jj,1:4) 0 0], data{ii}.freq(freqLim));
        V   = V/max(V);
        ind = find(V >= 0.5);
        f   = data{ii}.freq(freqLim);
        W   = abs(f(ind(1)) - f(ind(end)));
        FWHM{ii}(jj) = W * data{ii}.p.LarmorFreq;
        
        if showModelFits == 1
            cla;
            plot(data{ii}.freq(freqLim), real(data{ii}.spec(freqLim,jj)), ...
                data{ii}.freq(freqLim), Voigt(pars{ii}(jj,:), data{ii}.freq(freqLim)));
            set(gca, 'XDir', 'reverse', 'XLim', [4.8-0.25 4.8+0.25]);
            drawnow;
        end
    end
end

close all;

warning('on','stats:nlinfit:ModelConstantWRTParam');
warning('on','stats:nlinfit:IllConditionedJacobian');
warning('on','MATLAB:rankDeficientMatrix');

%% Plot spectra and frequency drift

H = figure(100);
clf;
d.w = 1;
d.h = 1;
d.l = (1-d.w)/2;
d.b = (1-d.h)/2;
set(H, 'Color', 'w', 'Units', 'Normalized', 'OuterPosition', [d.l d.b d.w d.h]);

xlabs          = string(0:0.025:10);
xlabs(2:2:end) = "";
titles         = [", pre-fMRI" ", post-fMRI" ...
                  ", pre-fMRI (day 2)" ", post-fMRI (day 2)"];
titles = siteID + titles;
c = 1;

% Spectra
for ii = 1:length(fnames)
    h1 = subplot(2, length(fnames), c);
    p = plot(data{ii}.freq, real(data{ii}.spec), 'LineWidth', 0.7);
    title(titles(ii), 'FontSize', 18);
    set(gca, 'XDir', 'reverse', 'XLim', [4.8-0.25 4.8+0.25], 'XTick', 0:0.025:10, 'XTickLabel', xlabs, ...
        'Box', 'off', 'TickDir', 'out', 'LineWidth', 1, 'FontSize', 14);
    h1.YAxis.Visible = 'off';
    ylim([min(real(data{ii}.spec),[],'all') - 0.1 * abs(min(real(data{ii}.spec),[],'all')) ...
        max(real(data{ii}.spec),[],'all') + 0.1 * abs(max(real(data{ii}.spec),[],'all'))]);
    xlabel('ppm', 'FontWeight', 'bold', 'FontSize', 16);
    set(p, {'color'}, num2cell(ones(size(data{ii}.spec,2),3) .* linspace(0, 0.7, size(data{ii}.spec,2))', 2));
    c = c + 1;
end

%xticks         = [1 100:100:2048];
xticks         = [1 20:20:360]; % scnh
xlabs          = string(xticks);
xlabs(2:2:end) = "";

% Frequency shifts
if length(fnames) == 2
    h2 = subplot(2, length(fnames), c);
else
    h2 = subplot(2, length(fnames), [c c + 1]);
end
hold on;
for ii = 1:length(fnames)
    plot(pars{ii}(:,3), 'LineWidth', 1);
    xlabel('Average', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('ppm (\Delta{\itF}_{0} [Hz])', 'FontSize', 16, 'FontWeight', 'bold');
    set(gca, 'XTick', xticks, 'XLim', [1 360], 'XTickLabel', xlabs, ... % scnh, changed 2048 to 64
        'Box', 'off', 'TickDir', 'out', 'LineWidth', 1, 'TickLength', [0.005 0.005], 'FontSize', 14);
end
hold off;
legend(titles(1:length(fnames)), 'Box', 'off', 'Location', 'best', 'FontSize', 16);
title('Water frequency shifts over time', 'FontSize', 18);

drawnow;
for ii = 1:length(h2.YTickLabel)
    a = (str2double(h2.YTickLabel(ii)) - 4.8) * data{1}.p.LarmorFreq;
    h2.YTickLabel{ii} = sprintf([h2.YTickLabel{ii} ' (%.2f)'], a);
end

% Linewidth
if length(fnames) == 2
    subplot(2, length(fnames), c + 1);
else
    subplot(2, length(fnames), [c + 2 c + 3]);
end
hold on;
for ii = 1:length(fnames)
%     plot(pi * abs(pars{ii}(:,4)) * data{ii}.p.LarmorFreq, 'LineWidth', 1);
    plot(FWHM{ii}, 'LineWidth', 1);
    set(gca, 'XTick', [1 20:20:360], 'XLim', [1 360], 'XTickLabel', xlabs, 'YLim', [0 7], ...
        'Box', 'off', 'TickDir', 'out', 'LineWidth', 1, 'TickLength', [0.005 0.005], 'FontSize', 14);
    xlabel('Average', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('FWHM (Hz)', 'FontSize', 16, 'FontWeight', 'bold');
end
hold off;
legend(titles(1:length(fnames)), 'Box', 'off', 'Location', 'best', 'FontSize', 16);
title('Water linewidth changes over time', 'FontSize', 18);

set(findall(H, '-property', 'FontName'), 'FontName', 'Arial');

%% Print plots to desktop

if printPlot == 1
    PrintPlot(fullfile(outDir, [siteID '_big_drift_data']));
    close all;
end







