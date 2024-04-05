
figure('Renderer', 'painters', 'Position', [100 100 1000 600]); % Open a figure window

count = 0;

for session_i = 1:3
    for area_i = 1:3
        count = count + 1;
        line_colors = magma(n_electrodes); % Generate a colormap
        
        subplot(3,3,count); hold on
        
        for electrode_i = 1:n_electrodes % Loop through the electrodes
            plot(time_win,nanmean(csd_signal_save{session_i,area_i}(electrode_i,:,:),3),...
                'color',line_colors(electrode_i,:)) % Plot the average LFP
        end
        
        xlabel(['Time from ' align_event ' (ms)']) % X-axis label
        ylabel(['Voltage (\muV)']) % Y-axis label
    end
end


%% Analysis: Clean faulty channels
fault_ch = {[],[22],[]; [],[],[]; [],[],[]};

for session_i = 1:3
    session_name = session_list{session_i};
    fprintf(['Analysing session: ' session_name '   \n']);
    
    for area_i = 1:3
        area_label = area_list{area_i};
        fprintf(['   |- Running analysis on ' area_label ' data \n']);
        
        n_trials = size(csd_signal_save{session_i,area_i},3);
        
        if ~isempty(fault_ch{session_i,area_i})
            % Average channels before and after together
            for trial_i = 1:n_trials
                csd_signal_save{session_i,area_i}(fault_ch{session_i,area_i},:,trial_i)...
                    = nanmean(csd_signal_save{session_i,area_i}(fault_ch{session_i,area_i}+[-1 1],:,trial_i));
            end
        end
        
    end
end

%% Curation: concatenate data across all trials to increase power

csd_signal_data_horz = [];

for session_i = 1:3
    session_name = session_list{session_i};
    fprintf(['Analysing session: ' session_name '   \n']);
    
    for area_i = 1:3
        area_label = area_list{area_i};
        fprintf(['   |- Running analysis on ' area_label ' data \n']);
        
        clear splitA keepIdx
        splitA = num2cell(csd_signal_save{session_i,area_i}, [1 2]); %split A keeping dimension 1 and 2 intact
        csd_signal_data_horz{session_i,area_i} = horzcat(splitA{:});
        
        [~,keepIdx] = find(~isnan(csd_signal_data_horz{session_i,area_i}(1,:)));
        
        csd_signal_data_horz{session_i,area_i} = ...
            csd_signal_data_horz{session_i,area_i}...
            (:,keepIdx);
    end
    
end

%% Analysis: Run laminar toolbox (Westerberg)
clear laminar_analysis
for session_i = 1:3
    session_name = session_list{session_i};
    fprintf(['Analysing session: ' session_name '   \n']);
    
    for area_i = 1:3
        area_label = area_list{area_i};
        fprintf(['   |- Running analysis on ' area_label ' data \n']);
        [laminar_analysis{session_i,area_i}.PSD,...
            laminar_analysis{session_i,area_i}.PSD_NORM,...
            laminar_analysis{session_i,area_i}.PSD_F]...
            = D_PSD_BASIC(csd_signal_data_horz{session_i,area_i});
    end
end

%% Figure: Plot Westerberg PSD

f_h = figure('Renderer', 'painters', 'Position', [100 300 1400 450]); hold on; % open figure window

count = 0;

for session_i = 1:3
    for area_i = 1:3
        area_label = area_list{area_i};
        count = count + 1;
        
        ax1 = subplot(3,3,count); hold on
        P_PSD_BASIC(laminar_analysis{session_i,area_i}.PSD_NORM,...
            laminar_analysis{session_i,area_i}.PSD_F,...
            f_h, ax1)
        title(area_label)

    end
end

%% Analysis: Run spectrolaminar suite (Sajad)
% % Define parameters
% freq_in = [3:120]; freq_step = 1; fs = 1000;
% fSteps = freq_in(1):freq_step:freq_in(end);
% 
% % Determine the array index for gamma and alpha values
% clear gamma_index alpha_index
% gamma_freq = [40:120]; gamma_index = find(ismember(fSteps,gamma_freq));
% alpha_freq = [8:25]; alpha_index = find(ismember(fSteps,alpha_freq));
% 
% session_i = 1;
% area_i = 1;
% 
% spectrolaminar_analysis = getSpectroLaminar(csd_signal_data_horz{session_i,area_i}, fs, freq_in, freq_step);
% plot_psd = H_SMOOTHD1_50um(spectrolaminar_analysis.norm_peak_1);


%% Analysis: get bandpower for alpha/gamma ranges

clear power_depth_*
for session_i = 1:3
    session_name = session_list{session_i};
    fprintf(['Analysing session: ' session_name '   \n']);
    
    for area_i = 1:3
        area_label = area_list{area_i};
        fprintf(['   |- Running analysis on ' area_label ' data \n']);
        for electrode_i = 1:n_electrodes % For each electrode
            power_depth_alpha(electrode_i,1) = bandpower(csd_signal_data_horz{session_i,area_i}(electrode_i,:),1000,[8 30]); % Get alpha/beta power
            power_depth_gamma(electrode_i,1) = bandpower(csd_signal_data_horz{session_i,area_i}(electrode_i,:),1000,[40 120]); % Get gamma power
        end
        % Normalize power of each contact to maximum power across electrode
        power_depth_alpha_out(session_i,area_i,:) = power_depth_alpha./max(power_depth_alpha); % for alpha
        power_depth_gamma_out(session_i,area_i,:) = power_depth_gamma./max(power_depth_gamma); % for gamma
    end
    
end

%% Figure: plot bandpower
figure('Renderer', 'painters', 'Position', [100 100 1000 600]); % Open a figure window

count = 0;

for session_i = 1:3
    for area_i = 1:3
        area_label = area_list{area_i};
        count = count + 1;
        
        subplot(3,3,count); hold on
        plot(squeeze(power_depth_alpha_out(session_i,area_i,:)), 1:n_electrodes,'r','LineWidth',2)
        plot(squeeze(power_depth_gamma_out(session_i,area_i,:)), 1:n_electrodes,'b','LineWidth',2)
        legend({'alpha','gamma'})
        set(gca,'XLim',[0 1],'YLim',[1 n_electrodes],'YDir','reverse')
        xlabel('Normalized power'); ylabel('Electrode contact')
        
        title(area_label)

        
    end
end
