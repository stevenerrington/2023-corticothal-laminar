clear all; clc

% Directories
if ispc()
    dirs.root = 'C:\Users\Steven\Documents\GitHub\2023-corticothal-laminar';
    dirs.data = 'Y:\Steven\2023-acc-vlpfc-thal\data';
else
    dirs.root = '/Users/stevenerrington/Desktop/Projects/2023-corticothal-laminar';
    dirs.data = '/Volumes/Share2/Steven/2023-acc-vlpfc-thal/data';
end

% Add paths
addpath(dirs.root);

% Session
session_name = 'LemmyKim-08162023-001_MYInfoPavChoice';

%% Load data
area_label = 'vlpfc';

clear data_in
% Load neural data (local field potentials)
fprintf(['Loading ' area_label ' data... \n']);
clear data_in; load_data = load(fullfile(dirs.data,[session_name '_' area_label '_lfp.mat'])); % Load LFP data
data_in = load_data.([area_label '_data']); clear load_data

% Load behavioral data
clear beh
load(fullfile(dirs.data,[session_name '_pds.mat'])); % Load behavioral data

%% Set parameters
% Get experiment parameters
n_electrodes = size(data_in,1);
n_trials = size(data_in,2);

% Set filter
filterFreq = [1 120];

% Set time windows for extraction
time_win = [-500:1000];
align_event = 'timetargeton'; %timefpon, timetargeton, timereward

%% Data structure

% Restructure data for CSD/spectrolaminar analysis
% - for each electrode
for electrode_i = 1:n_electrodes
    fprintf(['Analysing ' area_label ' electrode %i of %i \n'],...    
        electrode_i, n_electrodes);
    
    % - across each trial
    for trial_i = 1:n_trials

        try % if there is an alignment time
            % - get the time to align on
            onset_time = round(beh.(align_event)(trial_i)*1000) ;

            % - get the LFP data for the given 
            signal_in = data_in{electrode_i, trial_i} (onset_time+time_win,1)';
            signal_in_flitered = filter_signal(signal_in, filterFreq(1) , filterFreq(2), 1000);

            % - save the relevant data in an output array for future use
            csd_signal_out(electrode_i,:,trial_i) = signal_in_flitered;
            lfp_contact_out{electrode_i}(trial_i,:) = signal_in_flitered;

        catch % NaN out, if no event time
            csd_signal_out(electrode_i,:,trial_i) = nan(1,length(time_win));
            lfp_contact_out{electrode_i}(trial_i,:) = nan(1,length(time_win));
        end

    end
end

%% Test: check alignment

line_colors = magma(n_electrodes); % Generate a colormap

figure; hold on % Open a figure window

for electrode_i = 1:n_electrodes % Loop through the electrodes
    plot(time_win,nanmean(lfp_contact_out{electrode_i}),'color',line_colors(electrode_i,:)) % Plot the average LFP
end

xlabel(['Time from ' align_event ' (ms)']) % X-axis label
ylabel(['Voltage (\muV)']) % Y-axis label

%% Curate: clean faulty channels

fault_ch = [];

if ~isempty(fault_ch)
    % Average channels before and after together    
    for trial_i = 1:n_trials
        csd_signal_out(fault_ch,:,trial_i) = nanmean(csd_signal_out(fault_ch+[-1 1],:,trial_i));
        lfp_contact_out{fault_ch}(trial_i,:) = nanmean(csd_signal_out(fault_ch+[-1 1],:,trial_i));
    end
end

%% Analysis: Run laminar toolbox (Westerberg)
laminar_analysis = SUITE_LAM(csd_signal_out, 'spc', 0.05, 'times', time_win); % Run laminar toolbox

%% Analysis: Run spectrolaminar suite (Sajad)
% Define parameters
freq_in = [3:120]; freq_step = 3; fs = 1000;
fSteps = freq_in(1):freq_step:freq_in(end);

% Determine the array index for gamma and alpha values
clear gamma_index alpha_index
gamma_freq = [40:120]; gamma_index = find(ismember(fSteps,gamma_freq));
alpha_freq = [8:25]; alpha_index = find(ismember(fSteps,alpha_freq));

spectrolaminar_analysis = getSpectroLaminar(nanmean(csd_signal_out,3), fs, freq_in, freq_step);
plot_psd = H_SMOOTHD1_50um(spectrolaminar_analysis.norm_peak_1);

%% Analysis: Estimate power using matlab bandpower function
clear power_depth*

for electrode_i = 1:n_electrodes % For each electrode
    power_depth_alpha(electrode_i,1) = bandpower(nanmean(squeeze(csd_signal_out(electrode_i,:,:))'),1000,[8 30]); % Get alpha/beta power
    power_depth_gamma(electrode_i,1) = bandpower(nanmean(squeeze(csd_signal_out(electrode_i,:,:))'),1000,[40 120]); % Get gamma power
end

% Normalize power of each contact to maximum power across electrode
power_depth_alpha = power_depth_alpha./max(power_depth_alpha); % for alpha
power_depth_gamma = power_depth_gamma./max(power_depth_gamma); % for gamma

%% Figure: generate a summary figure

f_h = figure('Renderer', 'painters', 'Position', [100 300 1400 450]); hold on; % open figure window

% Power spectral density
ax1 = subplot(2, 3, 1);
P_PSD_BASIC(laminar_analysis.PSD_NORM, laminar_analysis.PSD_F, f_h, ax1)
xlabel('Frequency (Hz)'); ylabel('Contact'); 
title('Power spectral density')

ax1 = subplot(2, 3, 4);
imagesc('XData', fSteps, 'YData', [1:size(plot_psd,1)],'CData', plot_psd)
xlim([fSteps(1) fSteps(end)]); ylim([1 size(plot_psd,1)])
set(gca,'YDir','Reverse'); ylim([1 size(plot_psd,1)])

% Alpha/beta & gamma power x depth
ax2 = subplot(2, 3, 2); hold on
plot(power_depth_alpha, 1:n_electrodes,'LineWidth',2)
plot(power_depth_gamma, 1:n_electrodes,'LineWidth',2)
legend({'alpha','gamma'})
set(gca,'XLim',[0 1],'YLim',[1 n_electrodes],'YDir','reverse')
xlabel('Normalized power'); ylabel('Electrode contact')
title('Power x depth')

% PSD derived 
subplot(2,3,5); hold on
plot(nanmean(spectrolaminar_analysis.norm_peak_1(:,gamma_index),2),1:64)
plot(nanmean(spectrolaminar_analysis.norm_peak_1(:,alpha_index),2),1:64)
set(gca,'YDir','Reverse'); ylim([1 64])
legend({'Gamma','Alpha'})

% Cross-contact correlation
ax3 = subplot(2, 3, 3);
P_CORRE_BASIC(laminar_analysis.CORRE, 1:64, f_h, ax3)
title('Cross-contact correlation')

sgtitle(session_name,'Interpreter','None')  % Whole figure title

