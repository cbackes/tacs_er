function [enc_out,msg]=tACS_Encoding_FlashStim_OddBallCueResp(thePath)
% tACS stimulus encoding presentation script.
%
% this scripts presents stimuli to be encoded for a time depending on
% frequency of stimulation. the stimulus will flash at that frequency
% for a given duration. on a percentage of the trials the fixation c
% cross will change colors, which the subjects will need to identify
% it will only happen at a random point during the presentation of the
% stimuli. stimuli, stimulus order and stimulus conditions are pre-determined by the
% tacs_er make list  oddball. thePath indicates the path of the tacs_Eer
% structure.
%
% subjects performance on the cue color identification is stored peridiocally.
%
%------------------------------------------------------------------------%
% Author:       Alex Gonzalez
% Created:      Nov 3th, 2015
% LastUpdate:   Nov 4th, 2015
%------------------------------------------------------------------------%

%% Set-up

% clear all the screens
close all;
sca;

% load the task
fileName = strcat(thePath.subjectPath,'/tacs_er.task.mat');
if exist(fileName,'file')
    load(fileName);
else
    tacs_er =tACS_ER_makelist_Oddball(thePath);
end

% For debugging:
% PsychDebugWindowConfiguration

% Presentation Parameters
PresParams = [];
switch tacs_er.exptType
    case 'behav_v9'
    otherwise
        error('task not supported')
end

PresParams.stimFrequency        = 6;
PresParams.stimDurationInCycles  = 0.5;
PresParams.stimDurationInSecs   = 1/PresParams.stimFrequency*PresParams.stimDurationInCycles;
PresParams.totalStimDuration    = 2.5;
PresParams.nCycles              = PresParams.stimFrequency*PresParams.totalStimDuration;
PresParams.stimPresCycles       = ones(PresParams.nCycles,1); % stimuli will be presented every cycle
PresParams.OddBallnCycles       = 2;
PresParams.OddBallCycleRange    = [2 PresParams.nCycles-PresParams.OddBallnCycles-1]; % range in which the oddball might be presented


PresParams.cueDurationInSecs    = PresParams.stimDurationInSecs;
PresParams.ITI_Range            = [0.5 1]; % variable ITI in secs
PresParams.PostStimTime         = 0.5;     % time after stim (with no fixation).
PresParams.PreStimFixColor      = [0.53 0.1 0.53];
PresParams.PreStimFixColorStr   = 'PURPLE';
PresParams.CueColor1            = [0.77 0.05 0.2];
PresParams.CueColor1Str         = 'RED';
PresParams.CueColor2            = [0.2 0.1385 1];
PresParams.CueColor2Str         = 'BLUE';
PresParams.lineWidthPix         = 5;       % Set the line width for our fixation cross
PresParams.SaveEvNtrials        = 50;
PresParams.PauseEvNtrials       = tacs_er.nEncStim; % after evry block

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

stimNames = tacs_er.EncStimNames;
nTrials   = numel(stimNames);

%%

% Initialize trial timing structure
TimingInfo = [];
TimingInfo.preStimMaskFlip = cell(nTrials,1);
TimingInfo.stimPresFlip = cell(nTrials,PresParams.nCycles);
TimingInfo.postStimMaskFlip = cell(nTrials,1);
TimingInfo.trialRT          = nan(nTrials,1);
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
    postStimFrames     = round(PresParams.PostStimTime/ifi);
    
    % pre-make image textures
    imgTextures = cell(nTrials,1);
    for ii = 1:nTrials
        imgTextures{ii}=Screen('MakeTexture', window, tacs_er.Stimuli(stimNames{ii}));
        pgrStr = sprintf('Loading Stimuli  %g %%',floor(ii/nTrials*100));
        DrawFormattedText(window,pgrStr,'center','center',255,50);
        Screen('Flip',window);
    end
    
    % Set oddball color cue color
    cueColors = nan(nTrials,3);
    for ii = 1:nTrials
        if tacs_er.EncOddBallConds(ii)==1
            cueColors(ii,:) = PresParams.CueColorsID{1};
        elseif tacs_er.EncOddBallConds(ii)==2
            cueColors(ii,:) = PresParams.CueColorsID{2};
        end
    end
    
    % set oddball cycle presentation
    oddBallCycle = randi(PresParams.OddBallCycleRange,nTrials,1);
    PresParams.oddBallCycle =oddBallCycle;
    
    %---------------------------------------------------------------------%
    % Participant Instructions
    %---------------------------------------------------------------------%
    InstStr = ['Instructions\n\n' ...
        'You will be presented with a ' PresParams.PreStimFixColorStr ' Fixation  Cross. ' ...
        'A trial will start with the presentation of an image that will be shown on and off for ' num2str(PresParams.totalStimDuration) ...
        ' seconds. On some trials, the centrally presented fixation cross will change colors for a short duration. ' ...
        'Your task is to respond with ''' PresParams.RespToCue1 ''' for the ' PresParams.CueColor1Str ' Fixation and '...
        'with ''' PresParams.RespToCue2 ''' for ' PresParams.CueColor2Str ' Fixations. \n\n' ...
        'If no questions, \n'...
        'Press ''' resumeKey ''' to begin the experiment.'];
    
    PauseStr = ['Rest Pause. \\ '...
        'To continue press the ''' resumeKey ''' key.' ];
    
    DrawFormattedText(window,InstStr, 'wrapat', 'center', 255, 75, [],[],[],[],[xCenter*0.1,0,screenXpixels*0.8,screenYpixels]);
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
        
        % Pre-stimulus (variable ITI); store the first one
        Screen('DrawLines', window, fixCrossCoords,PresParams.lineWidthPix, PresParams.PreStimFixColor, [0 0], 2);
        [flip.VBLTimestamp, flip.StimulusOnsetTime, flip.FlipTimestamp, flip.Missed, flip.Beampos,] ...
            = Screen('Flip', window);
        TimingInfo.preStimMaskFlip{tt}=flip;
        vbl = flip.VBLTimestamp;
        
        for ii=1:(ITIsFrames(tt)-1)
            Screen('DrawLines', window, fixCrossCoords,PresParams.lineWidthPix, PresParams.PreStimFixColor, [0 0], 2);
            vbl = Screen('Flip', window, vbl + 0.5*ifi);
        end
        
        % Checks if the Pause Key has been pressed.
        CheckForPauseKey(pauseKey,resumeKey,activeKeyboardID)
        KbQueueFlush(activeKeyboardID);
        
        % cycles
        for ii =1:PresParams.nCycles
            % Draw Stimulus for stimFlipDurSecs
            Screen('DrawTexture', window, imgTextures{tt}, [], [], 0);
            if ii==1
                Screen('DrawLines', window, fixCrossCoords,PresParams.lineWidthPix, PresParams.PreStimFixColor, [0 0], 2);
                [flip.VBLTimestamp, flip.StimulusOnsetTime, flip.FlipTimestamp, flip.Missed, flip.Beampos,] ...
                    = Screen('Flip', window, vbl + 0.5*ifi);
                trialTime = GetSecs;
                vbl = flip.VBLTimestamp;
                
                Screen('DrawLines', window, fixCrossCoords,PresParams.lineWidthPix, PresParams.PreStimFixColor, [0 0], 2);
                vbl  = Screen('Flip', window, vbl + stimFlipDurSecs);
            elseif (oddBallCycle(tt)-ii>=0 && oddBallCycle(tt)-ii<PresParams.OddBallnCycles) && tacs_er.EncOddBallTrials(tt)
                Screen('DrawLines', window, fixCrossCoords,PresParams.lineWidthPix, cueColors(tt,:), [0 0], 2);
                [flip.VBLTimestamp, flip.StimulusOnsetTime, flip.FlipTimestamp, flip.Missed, flip.Beampos,] ...
                    = Screen('Flip', window, vbl + stimFlipDurSecs);
                vbl = flip.VBLTimestamp;
                
                Screen('DrawLines', window, fixCrossCoords,PresParams.lineWidthPix, cueColors(tt,:), [0 0], 2);
                vbl  = Screen('Flip', window, vbl + stimFlipDurSecs);
            else
                Screen('DrawLines', window, fixCrossCoords,PresParams.lineWidthPix, PresParams.PreStimFixColor, [0 0], 2);
                [flip.VBLTimestamp, flip.StimulusOnsetTime, flip.FlipTimestamp, flip.Missed, flip.Beampos,] ...
                    = Screen('Flip', window, vbl + stimFlipDurSecs);
                vbl = flip.VBLTimestamp;
                
                Screen('DrawLines', window, fixCrossCoords,PresParams.lineWidthPix, PresParams.PreStimFixColor, [0 0], 2);
                vbl  = Screen('Flip', window, vbl + stimFlipDurSecs);
            end
            TimingInfo.stimPresFlip{tt,ii}=flip;
        end
        [pressed,firstPress] = KbQueueCheck(activeKeyboardID);
        if pressed && tacs_er.EncOddBallTrials(tt)
            TimingInfo.trialKeyPress{tt} = KbName(firstPress);
            TimingInfo.trialRT(tt) = firstPress(find(firstPress,1))-trialTime;
        end
        
        % Draw Post-Stim Blank
        [flip.VBLTimestamp, flip.StimulusOnsetTime, flip.FlipTimestamp, flip.Missed, flip.Beampos,] ...
            = Screen('Flip', window, vbl + stimFlipDurSecs);
        TimingInfo.postStimMaskFlip{tt}=flip;
        vbl = flip.VBLTimestamp;
        for ii = 1:(postStimFrames-1)
            vbl  = Screen('Flip', window,vbl + 0.5*ifi);
        end
        
        % save every PresParams.SaveEvNtrials
        if mod(tt,PresParams.SaveEvNtrials)==0
            tempName = sprintf('/tacs_er.s%i.encoding.%s.mat;', thePath.subjNum, datestr(now,'dd.mm.yyyy.HH.MM'));
            save([thePath.subjectPath,tempName],'TimingInfo');
        end
        
        % Pause and wait for resume every PresParams.PauseEvNtrials
        if mod(tt,PresParams.PauseEvNtrials)==0
            DrawFormattedText(window,PauseStr, 'center', 'center', 255, 75, [],[],[],[],[xCenter*0.1,0,screenXpixels*0.8,screenYpixels]);
            Screen('Flip',window);
            WaitTillResumeKey(resumeKey,activeKeyboardID)
        end
        
        % Discard used image texture
        Screen('Close', imgTextures{tt})
    end
    
    %---------------------------------------------------------------------%
    % End of Experiment. Store data, and Close.
    %---------------------------------------------------------------------%
    
    % store additional outputs
    % output structure
    enc_out = [];
    enc_out.PresParams = PresParams;
    tacs_er.Stimuli = []; % don't re-store stimuli
    enc_out.exptInfo  = tacs_er;
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
    EndStr = ['End of Experiment.\n \n' ...
        'Press ''' resumeKey ''' to exit.'];
    
    DrawFormattedText(window,EndStr, 'center', 'center', 255, 40);
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
