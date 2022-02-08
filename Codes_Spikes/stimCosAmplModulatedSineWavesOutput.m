function SS = stimCosAmplModulatedSineWavesOutput(myScreenInfo,Pars)
% CosAmplModulatedSineWaves makes cosine amplitude modulated sine waves
% (for sound) 
%
% SS = CosAmplModulatedSineWaves(myScreenInfo,Pars) returns an object SS of type ScreenStim
%
% SS = CosAmplModulatedSineWaves(myScreenInfo) uses the default parameters
% 2021-05 EvB

% based on:
%2011-02 MC
% 2013-10 MK (significantly) modified stimWaveOutput for specific needs (optogenetic stimulation) and to allow 2 channels

%%

if nargin < 1
    error('Must at least specify myScreenInfo');
end

if nargin < 2
    Pars = [];
end

%% The parameters and their definition
pp = cell(1,1);
pp{1}  = {'dur',         'Stimulus duration (s*10)',                20,1,600};
pp{2}  = {'ampl',        'Max Amplitude (gain*1000)',               1000,0,1000}; % e.g. if you know speaker:output freq-->dB
pp{3}  = {'SoundFreq',   'Frequency of pulses (kHz*10)',               1,1,260};
pp{4}  = {'ModFreq',     'Amplitude Modulation Frequency (10*Hz)',  50, 10, 120};
pp{5}  = {'Outputtype',  'SoundCard(1), NIBord(2)',                 2, 1, 2};  %maybe this could be soundcard vs. 
pp{6} = {'TriggerType',  'Manual(1), External(2),Immediate(3)',       3,1,3}; %Syncing? Only used with NI board


x = XFile('stimCosAmplModulatedSineWavesOutput',pp);
% x.Write; % call this ONCE: it writes the .x file

%% Parse the parameters

if isempty(Pars)
    Pars = x.ParDefaults;
end

dur      = Pars(1)/10   % s
ampl     = Pars(2)/1000 % Gain of sound 
SoundFreq    = Pars(3)*100 %Hz frequency of sound
ModFreq    = Pars(4)/10 %Cosine Amplitude Modulation in Hz 
Outputtype = Pars(5) % 1 for Soundcard, 2 for NI board
TriggerType = Pars(6) % 1 for Soundcard, 2 for NI board


%% Make the stimulus

rigInfo = RigInfoGet;
SS = ScreenStim; % initialization
SS.Type = 'stimCosAmplModulatedSineWavesOutput';
SS.Parameters = Pars;
devices = PsychPortAudio('GetDevices');
if Outputtype == 1
    % Use soundcard
    %Change values
    fs = rigInfo.WaveInfo.SampleRate;     
    if ~isempty(rigInfo.SoundCardCalibrationFile)
        %load in calibration file and multiply by the gain factor
        ampl = ampl;  % Gain factor
    end
elseif Outputtype==2
    % Use NI AO
    %Change values
    fs = rigInfo.WaveInfo.SampleRate;
    indAO = rigInfo.SoundAOChannel;
    ampl=ampl*10;
end

%% building the sound waveform
t = (0:1/fs:dur)'; %Make for totaldur seconds
SoundWaveForm = sin(2*pi*SoundFreq*t).*(cos(2*pi*ModFreq*t - pi) + 1)/2*ampl;


%% putting things together
if Outputtype==2
    SS.WaveStim.WaveForm = cat(2, repmat(zeros(size(SoundWaveForm)), 1, indAO), SoundWaveForm);
    SS.WaveStim.Waves = cat(2, repmat(zeros(size(SoundWaveForm)), 1, indAO), SoundWaveForm);    
    % figure; plot(t,SS.WaveStim.SoundWaveForm)
    SS.WaveStim.SampleRate = fs;
    switch TriggerType
        case 1
            SS.WaveStim.TriggerType = 'Manual';
        case 2
            SS.WaveStim.TriggerType = 'HwDigital';
        case 3
            SS.WaveStim.TriggerType = 'Immediate';
    end
elseif Outputtype==1
    SS.WaveSoundcard = SoundWaveForm;
end
%% a blank visual stimulus
SS.nTextures = 1;
SS.nFrames = round(myScreenInfo.FrameRate*dur );
if ~mod(SS.nFrames, 2)
    % we want an odd number of frames if the SyncSquare is Flickering.
    SS.nFrames=SS.nFrames+1;
end
SS.Orientations = zeros(1,SS.nFrames);
SS.Amplitudes = zeros(1,SS.nFrames);
SS.nImages = 1;
SS.ImagePointers = Screen('MakeTexture', myScreenInfo.windowPtr, 0, [], 0, 1);
SS.ImageSequence = ones(1,SS.nFrames);
SS.SourceRects = repmat([1; 1; 1; 1],[1 1 SS.nFrames]);
SS.DestRects   = repmat([1; 1; 1; 1],[1 1 SS.nFrames]);

%debug only
% fprintf('stimOptiWaveOutput : nFrames = %d\n', SS.nFrames);
return

%% To test the code 

myScreenInfo = ScreenInfo(RigInfoGet);
SS = stimCosAmplModulatedSineWavesOutput(myScreenInfo); %#ok<UNRCH>
SS.Show(myScreenInfo)
% addpath('C:\Users\Experiment\Documents\GitHub\Linear Track Behav - fUSiRig\Used\StimMapping')
% show( SS, myScreenInfo);


