% run script for broadening simluated spectrum with acquired spectrum trace
clc
clear

% reload bd_load_data_SH_all_vendors.m for a new big_drift_data.mat if analyzing a new dataset
inputSpec = load('/Users/steve/Documents/My_Studies/Big_Drift/Code/sim_spectrum_for_big drift/in_vivo_like_spectrum.mat');
load('/Users/steve/Documents/My_Studies/Big_Drift/data/big_drift_data_final_all.mat');
%load('/Users/steve/Documents/My_Studies/Big_Drift/data/big_drift_data_day1_repro.mat');

vendor = {'GE','Philips','Siemens'};
%vendor = {'GE_Day2','Philips_Day2','Siemens_Day2'};

site{1} = {'G01a','G02a','G03a','G04a','G05a','G06a','G07a','G08a','G09a','G10a',...
          'G11a','G12a','G15a','G16a','G17a','G18a','G19a','G20a','G23a','G24a','G24b'}; % GE (21)
site{2} = {'P01a','P01b','P02a','P02b','P03a','P03b','P04a','P05a','P06a','P06b',...
          'P06c','P07a','P08a','P09a','P10a','P11a','P12a','P13a','P15a','P16a',...
          'P17a','P18a','P18b','P18c','P19a','P20a','P22a','P23a','P24a','P25a'};        % Philips (30)
site{3} = {'S01a','S03a','S03b','S03c','S04a','S05a','S06a','S07a','S07b','S08a',...
          'S09a','S10a','S11a','S11c','S12a','S13a','S13b','S14a','S15a','S15b',...
          'S16a','S17a','S18a','S18b','S19a','S20a','S21a','S22a','S23a','S24a',...
          'S25a','S26a','S27a','S27b','S28a','S29a','S30a','S31a','S33a','S37a',...
          'S38a','S40a','S40b','S99a'};                                                   % Siemens (44)

%Day2 data
% site{1} = {'G01a','G02a','G03a','G04a','G09a','G10a','G11a','G15a','G16a','G17a',...
%            'G18a','G19a','G23a','G24a','G24b'};                                         % GE (15)
% site{2} = {'P01a','P01b','P02a','P02b','P03a','P03b','P05a','P07a','P08a','P11a',...
%            'P12a','P13a','P15a','P17a','P18a','P18b','P18c','P19a','P20a','P23a'};             % Philips (20)
% site{3} = {'S01a','S03a','S03b','S03c','S04a','S05a','S07a','S07b','S08a','S09a',...
%            'S11a','S11c','S15a','S15b','S20a','S21a','S24a','S25a','S26a','S27a',...
%            'S27b','S29a','S30a','S37a','S38a'};                                         % Siemens (25)

% Removed: P25a, S02a
% SP21a has been changed to S99a

x=0;
freqs=zeros(64,(length(site{1})+length(site{2})+length(site{3})));
freqs_post=zeros(64,(length(site{1})+length(site{2})+length(site{3})));
%freqs_post=zeros(360,(length(site{1})+length(site{2})+length(site{3})));
% Loop through vendor
for ii = 1:length(vendor)
    
    % Loop through site
    for jj = 1:length(site{ii})
        x=x+1;
%        if size(deltaF0.(vendor{ii}).(site{ii}{jj}){2},1) < 64
%            %deltaF0.(vendor{ii}).(site{ii}{jj}){1}((end+1):64)=0;
%            deltaF0.(vendor{ii}).(site{ii}{jj}){2}((end+1):64)=0;
%        end        
        freqs(1:64,x)=deltaF0.(vendor{ii}).(site{ii}{jj}){1}(1:64);
        freqs_post(1:64,x)=deltaF0.(vendor{ii}).(site{ii}{jj}){2}(1:64);
       %freqs_post(1:360,x)=deltaF0.(vendor{ii}).(site{ii}{jj}){2};

    end
end
%freqs=ones(1,128).*(0:1:127)/10;

input=hilbert(inputSpec.out.specs);
for ii=(1:x)
    output(:,ii) = broaden_spectrum(input,freqs(:,ii));
    output_post(:,ii) = broaden_spectrum(input,freqs_post(:,ii));
end

% Automatically select max./median/min. spectrum and plot on a 3x3 figure
GE=zeros(95,1);
GE(1:length(site{1}))=1;
Philips=zeros(95,1);
Philips((length(site{1})+1):(length(site{1})+length(site{2})))=1;
Siemens=zeros(95,1);
Siemens((length(site{1})+1+length(site{2})):end)=1;

% Pre-fMRI PRESS
drift_pre_std=std(freqs);
drift_post_std=std(freqs_post(1:64,:));
drift_pre_mean=mean(abs(freqs));
drift_post_mean=mean(abs(freqs_post(1:64,:)));
%drift_post_std=std(freqs_post(1:360,:));
[value_GE_pre_min, G_index_pre_min]=min(drift_pre_mean(GE==1));
[value_GE_pre_med, G_index_pre_median]=median_SH(drift_pre_mean(GE==1));
[value_GE_pre_max, G_index_pre_max]=max(drift_pre_mean(GE==1));
[value_Ph_pre_min, P_index_pre_min]=min(drift_pre_mean(Philips==1));
[value_Ph_pre_med, P_index_pre_median]=median_SH(drift_pre_mean(Philips==1));
[value_Ph_pre_max, P_index_pre_max]=max(drift_pre_mean(Philips==1));
[value_Si_pre_min, S_index_pre_min]=min(drift_pre_mean(Siemens==1));
[value_Si_pre_med, S_index_pre_median]=median_SH(drift_pre_mean(Siemens==1));
[value_Si_pre_max, S_index_pre_max]=max(drift_pre_mean(Siemens==1));

P_index_pre_min = P_index_pre_min + length(site{1});
P_index_pre_median = P_index_pre_median + length(site{1});
P_index_pre_max = P_index_pre_max + length(site{1});
S_index_pre_min = S_index_pre_min + length(site{1})+length(site{2});
S_index_pre_median = S_index_pre_median + length(site{1})+length(site{2});
S_index_pre_max = S_index_pre_max + length(site{1})+length(site{2});

% Post-fMRI PRESS
[value_GE_post_min G_index_post_min]=min(drift_post_mean(GE==1));
[value_GE_post_med G_index_post_median]=median_SH(drift_post_mean(GE==1));
[value_GE_post_max G_index_post_max]=max(drift_post_mean(GE==1));
[value_Ph_post_min P_index_post_min]=min(drift_post_mean(Philips==1));
[value_Ph_post_med P_index_post_median]=median_SH(drift_post_mean(Philips==1));
[value_Ph_post_max P_index_post_max]=max(drift_post_mean(Philips==1));
[value_Si_post_min S_index_post_min]=min(drift_post_mean(Siemens==1));
[value_Si_post_med S_index_post_median]=median_SH(drift_post_mean(Siemens==1));
[value_Si_post_max S_index_post_max]=max(drift_post_mean(Siemens==1));
        
P_index_post_min = P_index_post_min + length(site{1});
P_index_post_median = P_index_post_median + length(site{1});
P_index_post_max = P_index_post_max + length(site{1});
S_index_post_min = S_index_post_min + length(site{1})+length(site{2});
S_index_post_median = S_index_post_median + length(site{1})+length(site{2});
S_index_post_max = S_index_post_max + length(site{1})+length(site{2});

%%
ppm_axis=inputSpec.out.ppm;

% Plot Pre-PRESS
figure (1)
subplot(3,3,1)
plot(ppm_axis,real(output(:,G_index_pre_min)),'r')
set(gca, 'XDir', 'reverse'),title('GE pre-PRESS Drift'),legend('Min. Drift','Location','northwest'),xlabel('PPM');
xlim([1 5]),ylim([0 2000])
subplot(3,3,4)
plot(ppm_axis,real(output(:,G_index_pre_median)),'r')
set(gca, 'XDir', 'reverse'),title('GE pre-PRESS Drift'),legend('Median Drift','Location','northwest'),xlabel('PPM');
xlim([1 5]),ylim([0 2000])
subplot(3,3,7)
plot(ppm_axis,real(output(:,G_index_pre_max)),'r')
set(gca, 'XDir', 'reverse'),title('GE pre-PRESS Drift'),legend('Max. Drift','Location','northwest'),xlabel('PPM');
xlim([1 5]),ylim([0 2000])

subplot(3,3,2)
plot(ppm_axis,real(output(:,P_index_pre_min)),'b')
set(gca, 'XDir', 'reverse'),title('Philips pre-PRESS Drift'),legend('Min.Drift','Location','northwest'),xlabel('PPM');
xlim([1 5]),ylim([0 2000])
subplot(3,3,5)
plot(ppm_axis,real(output(:,P_index_pre_median)),'b')
set(gca, 'XDir', 'reverse'),title('Philips pre-PRESS Drift'),legend('Median Drift','Location','northwest'),xlabel('PPM');
xlim([1 5]),ylim([0 2000])
subplot(3,3,8)
plot(ppm_axis,real(output(:,P_index_pre_max)),'b')
set(gca, 'XDir', 'reverse'),title('Philips pre-PRESS Drift'),legend('Max. Drift','Location','northwest'),xlabel('PPM');
xlim([1 5]),ylim([0 2000])

subplot(3,3,3)
plot(ppm_axis,real(output(:,S_index_pre_min)),'k')
set(gca, 'XDir', 'reverse'),title('Siemens pre-PRESS Drift'),legend('Min. Drift','Location','northwest'),xlabel('PPM');
xlim([1 5]),ylim([0 2000])
subplot(3,3,6)
plot(ppm_axis,real(output(:,S_index_pre_median)),'k')
set(gca, 'XDir', 'reverse'),title('Siemens pre-PRESS Drift'),legend('Median Drift','Location','northwest'),xlabel('PPM');
xlim([1 5]),ylim([0 2000])
subplot(3,3,9)
plot(ppm_axis,real(output(:,S_index_pre_max)),'k')
set(gca, 'XDir', 'reverse'),title('Siemens pre-PRESS Drift'),legend('Max. Drift','Location','northwest'),xlabel('PPM');
xlim([1 5]),ylim([0 2000])

% Plot Post-PRESS
figure (2)
subplot(3,3,1)
plot(ppm_axis,real(output_post(:,G_index_post_min)),'r')
set(gca, 'XDir', 'reverse'),title('GE post-PRESS Drift'),legend('Min. Drift','Location','northwest'),xlabel('PPM');
xlim([1 5]),ylim([0 2000])
subplot(3,3,4)
plot(ppm_axis,real(output_post(:,G_index_post_median)),'r')
set(gca, 'XDir', 'reverse'),title('GE post-PRESS Drift'),legend('Median Drift','Location','northwest'),xlabel('PPM');
xlim([1 5]),ylim([0 2000])
subplot(3,3,7)
plot(ppm_axis,real(output_post(:,G_index_post_max)),'r')
set(gca, 'XDir', 'reverse'),title('GE post-PRESS Drift'),legend('Max. Drift','Location','northwest'),xlabel('PPM');
xlim([1 5]),ylim([0 2000])

subplot(3,3,2)
plot(ppm_axis,real(output_post(:,P_index_post_min)),'b')
set(gca, 'XDir', 'reverse'),title('Philips post-PRESS Drift'),legend('Min.Drift','Location','northwest'),xlabel('PPM');
xlim([1 5]),ylim([0 2000])
subplot(3,3,5)
plot(ppm_axis,real(output_post(:,P_index_post_median)),'b')
set(gca, 'XDir', 'reverse'),title('Philips post-PRESS Drift'),legend('Median Drift','Location','northwest'),xlabel('PPM');
xlim([1 5]),ylim([0 2000])
subplot(3,3,8)
plot(ppm_axis,real(output_post(:,P_index_post_max)),'b')
set(gca, 'XDir', 'reverse'),title('Philips post-PRESS Drift'),legend('Max. Drift','Location','northwest'),xlabel('PPM');
xlim([1 5]),ylim([0 2000])

subplot(3,3,3)
plot(ppm_axis,real(output_post(:,S_index_post_min)),'k')
set(gca, 'XDir', 'reverse'),title('Siemens post-PRESS Drift'),legend('Min. Drift','Location','northwest'),xlabel('PPM');
xlim([1 5]),ylim([0 2000])
subplot(3,3,6)
plot(ppm_axis,real(output_post(:,S_index_post_median)),'k')
set(gca, 'XDir', 'reverse'),title('Siemens post-PRESS Drift'),legend('Median Drift','Location','northwest'),xlabel('PPM');
xlim([1 5]),ylim([0 2000])
subplot(3,3,9)
plot(ppm_axis,real(output_post(:,S_index_post_max)),'k')
set(gca, 'XDir', 'reverse'),title('Siemens post-PRESS Drift'),legend('Max. Drift','Location','northwest'),xlabel('PPM');
xlim([1 5]),ylim([0 2000])

