function [spectroLaminar]  = getSpectroLaminar( LFP_depthArranged, Fs, freqRange, stepSize, normType )

% timeSeries - LFP or EEG
% Fs - sampling frequency
% freqRange - the range of frequencies considered
% stepSize - the width of frequency inclusion in each bin

%%% Example inputs:
% timeSeries = lfp.data.(sprintf('LFP_%d', k))(100000:300000);
% freqSet = [3:90];
% stepSize = 1;
% Fs = lfp.info.samplingFreq;   % Sampling frequency
%%
if nargin < 3
    freqRange = [3:150];
end
if nargin < 4
    stepSize = 1;
end
if nargin < 5
    normType = 'peak_1';
end
%%
for ch = 1:size( LFP_depthArranged, 1)
    LFP = LFP_depthArranged(ch,:);
    T = 1/Fs;                     % Sample time
    L = length(LFP);          % Length of signal
    t = (0:L-1)*T;                % Time vector
    NFFT = 2^nextpow2(L); % Next power of 2 from length of y
    Y = fft(LFP,NFFT)/L;
    f = Fs/2*linspace(0,1,NFFT/2+1);
    absY=2*abs(Y(1:NFFT/2+1));
    %frequency window start
    fSteps = freqRange(1):stepSize:freqRange(end);
    for i=1:numel(fSteps)
        fBIN = [];
        fWinStart(i)= fSteps(i) - stepSize/2;
        fWinEnd(i) = fSteps(i) + stepSize/2;
        fBIN = [find( f>fWinStart(i), 1 )   find( f<fWinEnd(i), 1, 'last')];
        spectroLaminar.raw(ch,i) = mean( absY( fBIN(1):fBIN(2) ) );
    end
end

%% now that we have the power ber frequency bin for each channel, let's normalize the data:
% normalizing such that the depth at which maximum frequency is observed
% becomes 1:
if strcmpi(normType, 'peak_1')
    normSL = nan( size( spectroLaminar.raw ) );
     for i=1:numel(fSteps)
         normSL(:,i) = spectroLaminar.raw(:,i) / max( spectroLaminar.raw(:,i) );
     end
     % CAN ADD ADDITIONAL WAYS TO NORMALIZE
end
spectroLaminar.(['norm_' normType]) = normSL;



