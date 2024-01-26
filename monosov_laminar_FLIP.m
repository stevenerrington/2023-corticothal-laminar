clear all; clc

area = 'vlpfc';
%% Load data
% Load neural data (local field potentials)
clear data_in; load_data = load(['/Users/stevenerrington/Desktop/Test/LemmyKim-08162023-001_MYInfoPavChoice_' area '_lfp.mat']); % Load LFP data
data_in = load_data.([area '_data']); clear load_data
load('/Users/stevenerrington/Desktop/Test/LemmyKim-08162023-001_MYInfoPavChoice_pds.mat'); % Load behavioral data

%% Parameterization
% Set time windows for extraction
time_win = [-500:1000];
align_event = 'timetargeton'; %timefpon, timetargeton, timereward
filterFreq = [1 280];

% Get experiment parameters
n_electrodes = size(data_in,1);
n_trials = size(data_in,2);

%% Restructure data for spectrolaminar analysis
% - for each electrode
for electrode_i = 1:n_electrodes
    fprintf(['   |-- Analysing electrode %i of %i \n'], electrode_i, n_electrodes);

    % - across each trial
    for trial_i = 1:n_trials

        try % if there is an alignment time
            % - get the time to align on
            onset_time = round(beh.(align_event)(trial_i)*1000) ;

            % - get the LFP data for the given
            signal_in = data_in{electrode_i, trial_i} (onset_time+time_win,1)';
            signal_in_flitered = filter_signal(signal_in, filterFreq(1) , filterFreq(2), 1000);

            % - save the relevant data in an output array for future use
            signal_out(electrode_i,:,trial_i) = signal_in_flitered; % nchans x trialtime x ntrials

        catch % NaN out, if no event time
            signal_out(electrode_i,:,trial_i) = nan(1,length(time_win));
        end

    end
end

valid_trials = find(~isnan(beh.(align_event)));

%% Power analysis
[signal_normalized, signal_nonnormalized] = lfp_to_normpower(signal_out(:,:,valid_trials));

[startinglowfreq,endinglowfreq,...
    startinghighfreq,endinghighfreq,...
    goodnessvalue,superficialchannel,deepchannel,...
    highfreqmaxchannel,lowfreqmaxchannel,...
    crossoverchannel] = ...
    FLIPAnalysis(signal_nonnormalized,0:size(signal_nonnormalized,1)-1,1:size(signal_nonnormalized,2),1);



