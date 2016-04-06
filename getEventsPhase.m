function out = getEventsPhase(subj, expt)
out = [];

if strcmp(expt,'tacs_er')
    out.info = subjFileInfo(subj);
elseif strcmp(expt,'tacs_enc_xdiva')
    out.info = subjFileInfo_xdiva(subj);
end

stimData=load([out.info.dataPath out.info.stimulationFileName]);
eegData =load([out.info.dataPath out.info.eegEncodingFileName]);

%% Timing and event markers

stimSamples = stimData(:,9)-stimData(1,9);
stimTime    = stimSamples/1000;  % in seconds
eegSamples  = eegData(:,13)-eegData(1,13);
eegTime     = eegSamples/1000;  % in seconds
eventCodesSamps = eegSamples(eegData(:,12)~=0);
eventCodes  = eegData(eegData(:,12)~=0,12);

%%  get the phase of the samples

stimElec = stimData(:,1)/1e3; % in mA
stimElecAngle = angle(hilbert(stimElec));
stimElecAngle3Quant = quant(stimElecAngle,pi/3);

%% get the phase of each event

FixationIdx         = find(eventCodes==2); % fixations
FirstStimSamps      = eventCodesSamps(FixationIdx+1);
out.TrueAngleStims  = stimElecAngle(FirstStimSamps);
out.QuantAngleStims = stimElecAngle3Quant(FirstStimSamps);

save([out.info.dataPath 'EventsPhase'],'out')
