
function [enc_out,msg]=tACS_EncodingMain_CueResponse(thePath)
% tACS stimulus encoding presentation script.
%
% this scripts presents stimuli to be encoded for a time depending on
% frequency of stimulation. the subject's task is to identify the color
% of a fixation marker super-imposed on the stimuli of interest. stimuli
% , stimulus order and stimulus conditions are pre-determined by the
% tacs_er make list script. thePath indicates the path of the tacs_Eer
% structure.
%
% subjects performance on the cue color identification is stored peridiocally.
%
%------------------------------------------------------------------------%
% Author:       Alex Gonzalez
% Created:      Aug 1th, 2015
% LastUpdate:   Sept 15th, 2015
%------------------------------------------------------------------------%

%% Set-up

% clear all the screens
close all;
sca;

% load the task
fileName = strcat(thePath.subjectPath,'/tacs_er_task.mat');
if exist(fileName,'file')
    load(fileName);
else
    tacs_er =tACS_ER_makelist_RepeatStims(thePath);
end

% For debugging:
% PsychDebugWindowConfiguration

% output structure
enc_out = [];

% Presentation Parameters
PresParams = [];
PresParams.stimFrequency        = 6;
PresParams.stimDurationInCycles = 0.5;
PresParams.stimDurationInSecs   = 1/PresParams.stimFrequency*PresParams.stimDurationInCycles;
PresParams.cueDurationInSecs    = PresParams.stimDurationInSecs;
PresParams.ITI_Range            = [1.5 2]; % variable ITI in secs
PresParams.MaxResponseTime      = 1.5;       % maximum to make perceptual decision
PresParams.PreStimFixColor      = [1 1 1];
PresParams.PreStimFixColorStr   = 'WHITE';
PresParams.CueColor1            = [0.77 0.05 0.2];
PresParams.CueColor1Str         = 'RED';
PresParams.CueColor2            = [0.2 0.1385 1];
PresParams.CueColor2Str         = 'RED';
PresParams.lineWidthPix         = 5;       % Set the line width for our fixation cross
PresParams.Nmasks               = 50;      % number of noise masks
PresParams.nsMaskSize           = [255 255];  % noise mask size (same as stimuli)
PresParams.SaveEvNtrials        = 50;

% determine cue response mapping depending on subject number and active
% Keyboard.
laptopResponseKeys = ['k','l'];
keypadResponseKeys = ['1','2'];
if mod(tacs_er.subjNum,2)
    responseMap = [1,2];
else
    responseMap = [2,1];
end
laptopResponseKeys = laptopResponseKeys(responseMap);
keypadResponseKeys = keypadResponseKeys(responseMap);

PresParams.CueColorsID{1} = PresParams.CueColor1;
PresParams.CueColorsID{2} = PresParams.CueColor2;

enc_out.PresParams  = PresParams;
enc_out.expInfo     = tacs_er;

stimNames = tacs_er.EncStimNames;
nTrials   = numel(stimNames);

%%

% Initialize trial timing structure
TimingInfo = [];
TimingInfo.preStimMaskFlip = cell(nTrials,1);
TimingInfo.stimPresFlip = cell(nTrials,1);
TimingInfo.postStimMaskFlip = cell(nTrials,1);
TimingInfo.trialRT          = zeros(nTrials,1);
TimingInfo.trialKeyPress    = cell(nTrials,1);

try
    
    %---------------------------------------------------------------------%
    % Screen and additional presentation parameters
    %---------------------------------------------------------------------%
    % Get keyboard number
    [activeKeyboardID, laptopKeyboardID, pauseKey, resumeKey] = getKeyboardOr10key;
    % initialize Keyboard Queue
    KbQueueCreate(activeKeyboardID);
    % Start keyboard queue
    KbQueueStart(activeKeyboardID);
    
    if laptopKeyboardID==activeKeyboardID
        PresParams.RespToCue1 = laptopResponseKeys(1);
        PresParams.RespToCue2 = laptopResponseKeys(2);
    else
        PresParams.RespToCue1 = keypadResponseKeys(1);
        PresParams.RespToCue2 = keypadResponseKeys(2);
    end
    
    % initialie window
    [window, windowRect] = initializeScreen;
    screenXpixels = windowRect(3);
    screenYpixels = windowRect(4);
    
    % Get the centre coordinate of the window
    [xCenter, yCenter] = RectCenter(windowRect);
    
    % Get coordinates for fixation cross
    fixCrossCoords = fixCross(xCenter, yCenter,screenXpixels,screenYpixels);
    
    % Query the frame duration
    ifi = Screen('GetFlipInterval', window);
    
    % Get the durations in frames
    % variable pre-stimulus noise mask duration
    ITIsFrames         = randi(round(PresParams.ITI_Range/ifi),nTrials,1);
    
    % fixed stimulus duration
    stimDurFrames      = round(PresParams.stimDurationInSecs/ifi);
    % post-stim max response period duration
    MaxRespFrames       = round(PresParams.MaxResponseTime /ifi);
    
    % creat noise masks
    noiseTextures = cell(PresParams.Nmasks ,1);
    for ii = 1:PresParams.Nmasks
        noiseTextures{ii} = Screen('MakeTexture', window, rand(PresParams.nsMaskSize(1),PresParams.nsMaskSize(2)));
        tstring = sprintf('Loading Noise Masks %g %%', floor(ii/PresParams.Nmasks *100));
        DrawFormattedText(window,tstring,'center','center',255,50);
        Screen('Flip',window);
    end
    
    % pre-make image textures
    imgTextures = cell(nTrials,1);
    for ii = 1:nTrials
        imgTextures{ii}=Screen('MakeTexture', window, tacs_er.Stimuli(stimNames{ii}));
        tstring = sprintf('Loading Stimuli  %g %%',floor(ii/nTrials*100));
        DrawFormattedText(window,tstring,'center','center',255,50);
        Screen('Flip',window);
    end
    
    % Set encoding cue color
    cueColors = zeros(nTrials,3);
    for ii = 1:nTrials
        if tacs_er.EncStimCue(ii)==1
            cueColors(ii,:) = PresParams.CueColorsID{1};
        else
            cueColors(ii,:) = PresParams.CueColorsID{2};
        end
    end
    
    %---------------------------------------------------------------------%
    % Participant Instructions
    %---------------------------------------------------------------------%
    tstring = ['Instructions\n\n' ...
        'You will be presented with a ' PresParams.PreStimFixColorStr ' Fixation  Cross. ' ...
        'Your task is to respond with ' PresParams.RespToCue1 ' for the ' PresParams.CueColor1Str ' Fixation and '...
        'with ' PresParams.RespToCue2 ' for ' PresParams.CueColor2Str ' Fixations. ' ...
        'You will have ' num2str(PresParams.MaxResponseTime) ' seconds to respond '...
        'as quickly and as accurately as possible. \n\n'...
        'Press ''' resumeKey ''' to begin the experiment.'];
    
    
    DrawFormattedText(window,tstring, 'wrapat', 'center', 255, 75, [],[],[],[],[xCenter*0.1,0,screenXpixels*0.8,screenYpixels]);
    Screen('Flip',window);
    
    % resume if Resume Key is pressed
    WaitTillResumeKey(resumeKey,activeKeyboardID)
    
    %%
    %---------------------------------------------------------------------%
    % Trials
    %---------------------------------------------------------------------%
    % Set timing for each flip in a trial
    stimFlipDurSecs = (stimDurFrames - 0.5) * ifi;
    
    % Maximum priority level
    topPriorityLevel = MaxPriority(window);
    Priority(topPriorityLevel);
    
    % iterate through trials
    for tt = 1:nTrials
        
        % empty flip var
        flip     = [];
        
        % Pre-stimulus noise mask (variable ITI); store the first one
        Screen('DrawTexture', window, noiseTextures{randi(PresParams.Nmasks)}, [], [], 0);
        Screen('DrawLines', window, fixCrossCoords,PresParams.lineWidthPix, PresParams.PreStimFixColor, [0 0], 2);
        [flip.VBLTimestamp, flip.StimulusOnsetTime, flip.FlipTimestamp, flip.Missed, flip.Beampos,] ...
            = Screen('Flip', window);
        TimingInfo.preStimMaskFlip{tt}=flip;
        vbl = flip.VBLTimestamp;
        
        for ii=1:(ITIsFrames(tt)-1)
            Screen('DrawTexture', window, noiseTextures{randi(PresParams.Nmasks)}, [], [], 0);
            Screen('DrawLines', window, fixCrossCoords,PresParams.lineWidthPix, PresParams.PreStimFixColor, [0 0], 2);
            vbl = Screen('Flip', window, vbl + 0.5*ifi);
        end
        
        % Checks if the Pause Key has been pressed.
        CheckForPauseKey(pauseKey,resumeKey,activeKeyboardID)
        KbQueueFlush(activeKeyboardID);
        
        % Draw Stimulus for stimFlipDurSecs
        Screen('DrawTexture', window, imgTextures{tt}, [], [], 0);
        Screen('DrawLines', window, fixCrossCoords,PresParams.lineWidthPix, cueColors(tt,:), [0 0], 2);
        [flip.VBLTimestamp, flip.StimulusOnsetTime, flip.FlipTimestamp, flip.Missed, flip.Beampos,] ...
            = Screen('Flip', window, vbl + 0.5*ifi);
        TimingInfo.stimPresFlip{tt}=flip;
        trialTime = GetSecs;
        vbl = flip.VBLTimestamp;
        
        % Draw Post-Stim Noise
        Screen('DrawTexture', window, noiseTextures{randi(PresParams.Nmasks)}, [], [], 0);
        [flip.VBLTimestamp, flip.StimulusOnsetTime, flip.FlipTimestamp, flip.Missed, flip.Beampos,] ...
            = Screen('Flip', window, vbl + stimFlipDurSecs);
        
        TimingInfo.postStimMaskFlip{tt}=flip;
        vbl = flip.VBLTimestamp;
        
        % Re-draw noise mask until response or until max resp time
        for ii = 1:(MaxRespFrames-1)
            [pressed,firstPress] = KbQueueCheck(activeKeyboardID);
            
            if pressed
                TimingInfo.trialKeyPress{tt} = KbName(firstPress);
                TimingInfo.trialRT(tt) = firstPress(find(firstPress,1))-trialTime;
                break
            end
            Screen('DrawTexture', window, noiseTextures{randi(PresParams.Nmasks)}, [], [], 0);
            vbl  = Screen('Flip', window,vbl + 0.5*ifi);
        end
        
        % if no response.
        if ~pressed
            TimingInfo.trialRT(tt) = nan;
        end
        KbQueueFlush(activeKeyboardID);
        
        % save every PresParams.SaveEvNtrials
        if mod(tt,PresParams.SaveEvNtrials)==0
            tempName = sprintf('/tacs_er.s%i.encoding.%s.mat;', thePath.subjNum, datestr(now,'dd.mm.yyyy.HH.MM'));
            save([thePath.subjectPath,tempName],'TimingInfo');
        end
        
        % Discard used image texture
        Screen('Close', imgTextures{tt})
    end
    
    %---------------------------------------------------------------------%
    % End of Experiment. Store data, and Close.
    %---------------------------------------------------------------------%
    
    % store additional outputs
    enc_out.TimingInfo = TimingInfo;
    
    % save
    fileName = 'tacs_er.encoding.mat';
    cnt = 0;
    while 1
        savePath = strcat(thePath.subjectPath,'/',fileName);
        if ~exist(savePath,'file')
            save(savePath,'enc_out')
            break
        else
            cnt = cnt+1;
            warning(strcat(fileName,' already existed.'))
            fileName = strcat('tacs_er.encoding','-',num2str(cnt),'.mat');
            warning(strcat('saving as ', filenName))
        end
    end        
    
    % End of Experiment string
    tstring = ['End of Experiment.\n \n' ...
        'Press ''' resumeKey ''' to exit.'];
    
    DrawFormattedText(window,tstring, 'center', 'center', 255, 40);
    Screen('Flip',window);
    WaitTillResumeKey(resumeKey,activeKeyboardID)

    msg='allGood';
catch msg
    sca
    keyboard
end

% Clear the screen
Priority(0);
sca;
KbQueueStop(activeKeyboardID);
Screen('CloseAll');
ShowCursor;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% auxiliary functions and definitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %-------------------------------------------------------------------------%
% fixCrossCoords
% Set Fixation Cross Coordinates
%-------------------------------------------------------------------------%
function fixCrossCoords = fixCross(xCenter, yCenter,screenXpixels,screenYpixels)

fixCrossXlength = max(0.02*screenXpixels,0.02*screenYpixels); % max of 2% screen dims
fixCrossYlength = fixCrossXlength;

LeftExtent  = xCenter-fixCrossXlength/2;
RightExtent = xCenter+fixCrossXlength/2 ;
BottomExtent = yCenter+fixCrossYlength/2 ;
TopExtent   =  yCenter- fixCrossYlength/2 ;

fixCrossXCoords   = [LeftExtent RightExtent; yCenter yCenter];
fixCrossYCoords   = [xCenter xCenter; BottomExtent TopExtent];

fixCrossCoords       = [fixCrossXCoords fixCrossYCoords];

end

%-------------------------------------------------------------------------%
% WaitTillResumeKey
% Wait until Resume Key is pressed on the keyboard
%-------------------------------------------------------------------------%
function WaitTillResumeKey(resumeKey,activeKeyboardID)

KbQueueFlush(activeKeyboardID);
while 1
    [pressed,firstPress] = KbQueueCheck(activeKeyboardID);
    if pressed
        if strcmp(resumeKey,KbName(firstPress));
            break
        end
    end
    WaitSecs(0.1);
end
KbQueueFlush(activeKeyboardID);
end

%-------------------------------------------------------------------------%
% CheckForPauseKey
% Check if the resume key has been pressed, and pause exection until resume
% key is pressed.
%-------------------------------------------------------------------------%
function CheckForPauseKey(pauseKey,resumeKey,activeKeyboardID)

[pressed,firstPress] = KbQueueCheck(activeKeyboardID);
if pressed
    if strcmp(pauseKey,KbName(firstPress));
        WaitTillResumeKey(resumeKey,activeKeyboardID)
    end
end
end
