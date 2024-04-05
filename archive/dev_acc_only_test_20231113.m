clear all; clc

% Directories
if ispc()
    dirs.root = 'C:\Users\Steven\Documents\GitHub\2023-corticothal-laminar';
    dirs.data = 'Y:\Steven\2023-acc-vlpfc-thal\data';
    dirs.data = 'C:\Users\Steven\Desktop\monosov_laminar\acc_only';
else
    dirs.root = '/Users/stevenerrington/Desktop/Projects/2023-corticothal-laminar';
    dirs.data = '/Volumes/Share2/Steven/2023-acc-vlpfc-thal/data';
end

% Add paths
addpath(dirs.root);

data_mat_files = dir([dirs.data '\*acc_lfp.mat']); % list of all .mat files in dir
n_mat = size(data_mat_files,1); % number of files
area_label = 'acc';


%% Load data
for session_i = 1:n_mat
    
    file_in = data_mat_files(session_i).name;
    try
        clear data_in
        % Load neural data (local field potentials)
        fprintf(['Loading ' area_label ' data | ' file_in ' | %i of %i... \n'],session_i,n_mat);
        clear data_in; load_data = load(fullfile(dirs.data,file_in)); % Load LFP data
        data_in = load_data.([area_label '_data']); clear load_data
        
        % Load behavioral data
        clear beh
        load(fullfile(dirs.data,[file_in(1:end-12) '_pds.mat'])); % Load behavioral data
        
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
        clear csd_signal_out lfp_contact_out
        % Restructure data for CSD/spectrolaminar analysis
        % - for each electrode
        for electrode_i = 1:n_electrodes           
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
        freq_in = [3:120]; freq_step = 1; fs = 1000;
        fSteps = freq_in(1):freq_step:freq_in(end);

        
        clear laminar_analysis spectrolaminar_analysis
        laminar_analysis = SUITE_LAM(csd_signal_out, 'spc', 0.05, 'times', time_win); % Run laminar toolbox
        spectrolaminar_analysis = getSpectroLaminar(nanmean(csd_signal_out,3), fs, freq_in, freq_step);
        
        filename = {file_in};
        
        laminar_table(session_i,:) = ...
            table(filename,beh,laminar_analysis,spectrolaminar_analysis,{csd_signal_out});
    catch
        filename_error{session_i} = {file_in};
        
    end
end

for session_i = 1:size(laminar_table,1)
    
    clear power_depth_alpha power_depth_gamma
    
    for electrode_i = 1:n_electrodes % For each electrode
        power_depth_alpha(electrode_i,1) = bandpower(nanmean(squeeze(laminar_table.Var5{session_i,1}(electrode_i,:,:))'),1000,[8 30]); % Get alpha/beta power
        power_depth_gamma(electrode_i,1) = bandpower(nanmean(squeeze(laminar_table.Var5{session_i,1}(electrode_i,:,:))'),1000,[40 120]); % Get gamma power
    end
    
    % Normalize power of each contact to maximum power across electrode
    power_depth_alpha_out(:,session_i) = power_depth_alpha./max(power_depth_alpha); % for alpha
    power_depth_gamma_out(:,session_i) = power_depth_gamma./max(power_depth_gamma); % for gamma
    
end


for session_i = 1:size(laminar_table,1)
f_h = figure('Renderer', 'painters', 'Position', [100 300 800 450]); hold on; % open figure window

% Power spectral density
ax1 = subplot(1, 4, [1 2 3]);
P_PSD_BASIC(laminar_table(session_i,:).laminar_analysis.PSD_NORM,...
    laminar_table(session_i,:).laminar_analysis.PSD_F, f_h, ax1)
xlabel('Frequency (Hz)'); ylabel('Contact'); 
title('Power spectral density')
colormap(parula)

ax2 = subplot(1, 4, 4); hold on
plot(power_depth_alpha_out(:,session_i), 1:64,'LineWidth',2)
plot(power_depth_gamma_out(:,session_i), 1:64,'LineWidth',2)
legend({'alpha','gamma'})
set(gca,'XLim',[0 1],'YLim',[1 64],'YDir','reverse')
xlabel('Normalized power'); ylabel('Electrode contact')
title('Power x depth')

sgtitle(laminar_table(session_i,:).filename,'Interpreter','None')  % Whole figure title


end
