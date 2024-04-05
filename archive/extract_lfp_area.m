clc; clear all; close all

% Parameters & general settings
acc_channels = [1:64];
vlpfc_channels = [65:128];
thalamus_channels = [129:192];

% Define output directory
save_dir = 'C:\Users\Steven\Desktop\monosov_laminar'; % where data will be saved

% Get a list of all files in a directory
lfp_dir = 'Y:\MONKEYDATA\Slayer2\MYsimRecord\AIFP\Thalamus'; % main lfp storage dir
pds_dir = 'Y:\MONKEYDATA\Slayer2\MYsimRecord\PDS\Thalamus'; % main beh storage dir

data_mat_files = dir([lfp_dir '\*.mat']); % list of all .mat files in dir
n_mat = size(data_mat_files,1); % number of files

clear my_sim_map
% Loop through data files in directory.
for file_i = 1:n_mat
    % Get mat file name
    data_file = data_mat_files(file_i).name;
    
    fprintf('Analysing session %i of %i | %s    \n',...
        file_i,n_mat,data_file(1:end-22));
    
    % Load data file
    clear AIFP
    load(fullfile(lfp_dir,data_file))
    
    % ------------------------------------------------------------------
    % 2023-11-06: SE for SA
    % For eye data. I thought it was stored in a separate
    % variable/structure within the datafile, but it appears not. If you
    % run this code again, then you should be able to call the eye data and
    % save it separately.
    
    % Example-----------------------------
    % Add this to the scripts below to extract eyes:
    %  load(fullfile(lfp_dir,data_file))
    %  eye_data = AIFP.XXXX.eye; (I'm not sure what the variable is, as
    %  it's taking so long to load in...)
    %  save(fullfile(save_dir,[data_file(1:end-22) '_eye.mat']),'eye_data', '-v7.3');
    %  ------------------------------------------------------------------
    
    try
        % Break down large LFP file into separate areas for future processing
        clear acc_data vlpfc_data thalamus_data
        acc_data = AIFP.FP.trialbase(acc_channels,:); % approx 5GB
        vlpfc_data = AIFP.FP.trialbase(vlpfc_channels,:); % approx 5GB
        thalamus_data = AIFP.FP.trialbase(thalamus_channels,:); % approx 5GB
        
        % Save individual area files
        save(fullfile(save_dir,[data_file(1:end-22) '_acc_lfp.mat']),'acc_data', '-v7.3');
        save(fullfile(save_dir,[data_file(1:end-22) '_vlpfc_lfp.mat']),'vlpfc_data', '-v7.3');
        save(fullfile(save_dir,[data_file(1:end-22) '_thalamus_lfp.mat']),'thalamus_data', '-v7.3');
        
        % Load behavioral datafile (drawn from first PDS file)
        clear all_beh PDS
        all_beh = dir([pds_dir '\' data_file(1:end-22) '*.mat']); % list of all .mat files in dir
        load(fullfile(pds_dir,all_beh(1).name),'PDS')
        
        % Save filenames and behavior for table call
        beh = PDS;
        acc_file = [data_file(1:end-22) '_acc_lfp.mat'];
        vlpfc_file = [data_file(1:end-22) '_vlpfc_lfp.mat'];
        thalamus_file = [data_file(1:end-22) '_thalamus_lfp.mat'];
        
        % Create data map for current file
        my_sim_map(file_i,:) = table(file_i,{data_file(1:end-22)},{pds_dir},{lfp_dir},...
            beh,{acc_file},{vlpfc_file},{thalamus_file},...
            'VariableNames',{'file_i','session','pds_dir','lfp_dir','beh','acc_file','vlpfc_file','thal_file'});
    catch
        my_sim_map(file_i,:) = table(file_i,{data_file(1:end-22)},{pds_dir},{lfp_dir},...
            NaN,NaN,NaN,NaN,...
            'VariableNames',{'file_i','session','pds_dir','lfp_dir','beh','acc_file','vlpfc_file','thal_file'});
    end
    
    
end

