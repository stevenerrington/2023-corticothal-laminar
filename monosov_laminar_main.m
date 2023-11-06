clear all; clc

addpath('/Users/stevenerrington/Desktop/Projects/2023-corticothal-laminar');

area_label = 'acc';
session_name = 'LemmyKim-08112023-003_MYInfoPavChoice_NovelReward';

%% Load data

clear data_in
% Load neural data (local field potentials)
switch area_label
    case 'vlpfc'
        load('/Users/stevenerrington/Desktop/monosov_laminar/data/vlpfc_lfp.mat') % Load VLPFC LFP data
        data_in = vlpfc_data; clear vlpfc_data;
    case 'acc'
        load('/Users/stevenerrington/Desktop/monosov_laminar/data/acc_lfp.mat') % Load ACC LFP data
        data_in = acc_data; clear acc_data;
     case 'thalamus'
        load('/Users/stevenerrington/Desktop/monosov_laminar/data/thalamus_lfp.mat') % Load ACC LFP data       
        data_in = thalamus_data; clear thalamus_data;
end

% Load behavioral data
load(['/Users/stevenerrington/Desktop/monosov_laminar/data/' session_name '_cl1_PDS.mat']) % Beh

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

clear signal_out test

% Restructure data for CSD/spectrolaminar analysis
% - for each electrode
for electrode_i = 1:n_electrodes
    % - across each trial
    for trial_i = 1:n_trials

        try % if there is an alignment time
            % - get the time to align on
            onset_time = round(PDS.(align_event)(trial_i)*1000) ;

            % - get the LFP data for the given 
            signal_in = data_in{electrode_i, trial_i} (onset_time+time_win,1)';
            signal_in_flitered = filter_signal(signal_in, filterFreq(1) , filterFreq(2), 1000);

            % - save the relevant data in an output array for future use
            signal_out(electrode_i,:,trial_i) = signal_in_flitered;
            test{electrode_i}(trial_i,:) = signal_in_flitered;

        catch % NaN out, if no event time
            signal_out(electrode_i,:,trial_i) = nan(1,length(time_win));
            test{electrode_i}(trial_i,:) = nan(1,length(time_win));
        end

    end
end

%% Test: check alignment

line_colors = magma(n_electrodes); % Generate a colormap

figure; hold on % Open a figure window

for electrode_i = 1:n_electrodes % Loop through the electrodes
    plot(time_win,nanmean(test{electrode_i}),'color',line_colors(electrode_i,:)) % Plot the average LFP
end

xlabel(['Time from ' align_event ' (ms)']) % X-axis label
ylabel(['Voltage (\muV)']) % Y-axis label

%% Curate: clean faulty channels

fault_ch = 22;

% Average channels before and after together

    for trial_i = 1:n_trials
            signal_out(fault_ch,:,trial_i) = nanmean(signal_out(fault_ch+[-1 1],:,trial_i));
            test{fault_ch}(trial_i,:) = nanmean(signal_out(fault_ch+[-1 1],:,trial_i));
    end

%% Analysis: Run laminar toolbox (Westerberg)
analysis = SUITE_LAM(signal_out, 'spc', 0.05, 'times', time_win); % Run laminar toolbox

%% Analysis: Estimate power using matlab bandpower function
clear power_depth*

for electrode_i = 1:n_electrodes % For each electrode
    power_depth_alpha(electrode_i,1) = bandpower(nanmean(squeeze(signal_out(electrode_i,:,:))'),1000,[8 30]); % Get alpha/beta power
    power_depth_gamma(electrode_i,1) = bandpower(nanmean(squeeze(signal_out(electrode_i,:,:))'),1000,[40 120]); % Get gamma power
end

% Normalize power of each contact to maximum power across electrode
power_depth_alpha = power_depth_alpha./max(power_depth_alpha); % for alpha
power_depth_gamma = power_depth_gamma./max(power_depth_gamma); % for gamma

%% Figure: generate a summary figure

f_h = figure('Renderer', 'painters', 'Position', [100 300 1400 450]); hold on; % open figure window

% Power spectral density
ax1 = subplot(1, 3, 1);
P_PSD_BASIC(analysis.PSD_NORM, analysis.PSD_F, f_h, ax1)
xlabel('Frequency (Hz)'); ylabel('Contact'); 
title('Power spectral density')

% Alpha/beta & gamma power x depth
ax2 = subplot(1, 3, 2); hold on
plot(power_depth_alpha, 1:n_electrodes,'LineWidth',2)
plot(power_depth_gamma, 1:n_electrodes,'LineWidth',2)
legend({'alpha','gamma'})
set(gca,'XLim',[0 1],'YLim',[1 n_electrodes],'YDir','reverse')
xlabel('Normalized power'); ylabel('Electrode contact')
title('Power x depth')

% Cross-contact correlation
ax3 = subplot(1, 3, 3);
P_CORRE_BASIC(analysis.CORRE, 1:64, f_h, ax3)
title('Cross-contact correlation')

sgtitle(session_name,'Interpreter','None')  % Whole figure title
