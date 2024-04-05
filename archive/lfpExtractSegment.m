clear
dataFolder = 'E:\Dropbox\SEF CMAND UNITS & LFP LAYER PAPER\SAJAD - ERRINGTON - SCHALL - DataAnalysis\DATA_LFPS';
fileList = dir([dataFolder '\*.mat']);
fileNameFormat = 'lfpData_probe_X.mat';
segmentType = 'segmented'; % can be 'segmented' (only part of the full session) or 'fullTrace' (the entire session's data - slower)
%% spectrolaminar plot parameters:
Fs = 24414.0625 / 24;
freqRange = 1:120;
normType = 'peak_1';
freqStepSize = 1;
%% extract data:
validLFPfile = 0;
probeIdxList = 1:168;
for probeIdxIdx = 1:numel(probeIdxList)
    probeIdx = probeIdxList(probeIdxIdx);
    clear LFPdata LFPinfo LFP_depthArranged
%     fileName = fileList(fileIdx).name;
%     probeIdx = str2num( fileName( ~ismember( fileName, fileNameFormat ) ) );
    load(['E:\Dropbox\SEF CMAND UNITS & LFP LAYER PAPER\SAJAD - ERRINGTON - SCHALL - DataAnalysis\DATA_LFPS\lfpData_probe_' int2str(probeIdx) '.mat'], 'LFPdata', 'LFPinfo')
%     if exist('LFPData', 'var')
        validLFPfile = validLFPfile + 1;
        if strcmpi(segmentType, 'segmented') == 1
            segmentSamples = round( size( LFPdata.lfp, 2 )/2 + ( size( LFPdata.lfp, 2 ) / 20 )*[-0.5 0.5] ); % 1:size( LFPdata.lfp, 2 ); %
            LFP_depthArranged = LFPdata.lfp( :, segmentSamples(1):segmentSamples(end) );
        else
            LFP_depthArranged = LFPdata.lfp;
            segmentSamples = nan;
        end
        [spectroLaminar.data.(segmentType)]  = getSpectroLaminar( LFP_depthArranged, Fs, freqRange, freqStepSize, normType );
        spectroLaminar.info.probeIdx = probeIdx;
        spectroLaminar.info.monkey = LFPinfo.monkey{1};
        spectroLaminar.info.Area = LFPinfo.Area{1};
        spectroLaminar.info.grid = LFPinfo.grid{1};
        spectroLaminar.info.probeDepth = LFPinfo.probeDepth(1);
        spectroLaminar.info.fileName = LFPinfo.sessionName{1};
        spectroLaminar.info.probeTag = LFPinfo.probeTag(1);
        spectroLaminar.info.Fs = Fs;
        spectroLaminar.info.freqRange = freqRange;
        spectroLaminar.info.normType = normType;
        spectroLaminar.info.freqStepSize = freqStepSize;
        spectroLaminar.info.segmentSamples = segmentSamples;
        save(['E:\Dropbox\SEF CMAND UNITS & LFP LAYER PAPER\SAJAD - ERRINGTON - SCHALL - DataAnalysis\DATA_SPECTROLAMINAR\' segmentType '\spectroLaminar_probe_' int2str(probeIdx) '.mat'], 'spectroLaminar')
        spectroLaminarMaster.info.probeIdx(validLFPfile) = probeIdx;
        spectroLaminarMaster.info.monkey{validLFPfile} = spectroLaminar.info.monkey;
        spectroLaminarMaster.info.Area{validLFPfile} = spectroLaminar.info.Area;
        spectroLaminarMaster.info.grid{validLFPfile} = spectroLaminar.info.grid;
        spectroLaminarMaster.info.probeDepth(validLFPfile) = spectroLaminar.info.probeDepth;
        spectroLaminarMaster.info.fileName{validLFPfile} = spectroLaminar.info.fileName;
        spectroLaminarMaster.info.probeTag(validLFPfile) = spectroLaminar.info.probeTag;
        spectroLaminarMaster.info.info.Fs(validLFPfile) = spectroLaminar.info.Fs;
        spectroLaminarMaster.info.info.freqRange{validLFPfile} = spectroLaminar.info.freqRange;
        spectroLaminarMaster.info.info.normType{validLFPfile} = spectroLaminar.info.normType;
        spectroLaminarMaster.info.info.freqStepSize(validLFPfile) = spectroLaminar.info.freqStepSize;
        spectroLaminarMaster.info.info.segmentSamples{validLFPfile} = spectroLaminar.info.segmentSamples;
        spectroLaminarMaster.data.(segmentType).raw(:,:,validLFPfile) = spectroLaminar.data.(segmentType).raw;
        spectroLaminarMaster.data.(segmentType).(['norm_' spectroLaminar.info.normType])(:,:,validLFPfile) = spectroLaminar.data.(segmentType).(['norm_' spectroLaminar.info.normType]);
%         spectroLaminarMaster.data.fullTrace.raw(:,:,validLFPfile) = spectroLaminar.data.fullTrace.raw;
%         spectroLaminarMaster.data.fullTrace.(['norm_' spectroLaminar.info.normType])(:,:,validLFPfile) = spectroLaminar.data.fullTrace.(['norm_' spectroLaminar.info.normType]);
%     end
end
save(['E:\Dropbox\SEF CMAND UNITS & LFP LAYER PAPER\SAJAD - ERRINGTON - SCHALL - DataAnalysis\DATA_SPECTROLAMINAR\' segmentType '\spectroLaminarMaster.mat'], 'spectroLaminarMaster', '-v7.3')