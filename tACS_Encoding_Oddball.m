
function [enc_out,msg]=tACS_Encoding_Oddball(thePath)
% tACS stimulus encoding presentation script.
%
% this scripts presents stimuli to be encoded for a time depending on
% frequency of stimulation. the subject's task is to responde when an
% oddball dot appears in one of 4 locations in the image. the subject will
% respond left or right for the oddball.
% stimulus order and stimulus conditions are pre-determined by the
% tacs_er make list script. thePath indicates the path of the tacs_er
% structure.
%
% subjects performance on the oddball target task is stored periodically
%
%------------------------------------------------------------------------%
% Author:       Alex Gonzalez
% Created:      Oct 12th, 2015
% LastUpdate:   Oct 12th, 2015
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
    tacs_er =tACS_ER_makelist_RS_Oddball(thePath);
end

% For debugging:
% PsychDebugWindowConfiguration

% Presentation Parameters
PresParams = [];
PresParams.stimFrequency        = 5;

PresParams.stimDurationInCycles = 0.5;
PresParams.stimDurationInSecs   = 1/PresParams.stimFrequency*PresParams.stimDurationInCycles;
PresParams.OddBallDurationInSecs= PresParams.stimDurationInSecs/2; % only appear on second half of total duration
PresParams.ITI_Range            = [1 1.5];   % variable ITI in secs
PresParams.MaxResponseTime      = 1.5;       % maximum to make oddball detection
PresParams.PostStimTime         = 0.5;       % time after stim (with no fixation).
PresParams.PreStimFixColor      = [1 1 1];
PresParams.PreStimFixColorStr   = 'WHITE';
PresParams.OddBallColor         = [0.77 0.05 0.2];
PresParams.OddBallColorStr      = 'RED';
PresParams.OddBallSize          = 10;  % pixels
PresParams.lineWidthPix         = 5;       % Set the line width for our fixation cross
PresParams.Nmasks               = 50;      % number of noise masks
PresParams.nsMaskSize           = tacs_er.stimSize;  % noise mask size (same as stimuli)
PresParams.SaveEvNtrials        = 50;
PresParams.PauseEvNtrials       = tacs_er.nEncStim; % after evry block

laptopResponseKeys = ['k','l'];
keypadResponseKeys = ['1','2'];

enc_out.PresParams  = PresParams;

stimNames = tacs_er.EncStimNames;
nTrials   = numel(stimNames);

%%

% Initialize trial timing structure
TimingInfo = [];
TimingInfo.preStimMaskFlip  = cell(nTrials,1);
TimingInfo.stimPresFlip     = cell(nTrials,1);
TimingInfo.oddBallPresFlip  = cell(nTrials,1);
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
        PresParams.RespToLeft = laptopResponseKeys(1);
        PresParams.RespToRight = laptopResponseKeys(2);
    else
        PresParams.RespToLeft = keypadResponseKeys(1);
        PresParams.RespToRight = keypadResponseKeys(2);
    end
    
    % initialie window
    [window, windowRect] = initializeScreen;
    screenXpixels = windowRect(3);
    screenYpixels = windowRect(4);
    
    % Get the centre coordinate of the window
    [xCenter, yCenter] = RectCenter(windowRect);
    
    % Get coordinates for fixation cross
    fixCrossCoords = fixCross(xCenter, yCenter,screenXpixels,screenYpixels);
    
    % Coordinates for oddballs.
    dotDist = PresParams.nsMaskSize(1)*0.15;
    OB_Coords{1}  = [-dotDist+xCenter -dotDist+yCenter]; % upper left
    OB_Coords{2}  = [-dotDist+xCenter dotDist+yCenter]; % bottom left
    OB_Coords{3}  = [dotDist+xCenter -dotDist+yCenter]; % bottom right
    OB_Coords{4}  = [dotDist+xCenter dotDist+yCenter]; % bottom right
    
    % Query the frame duration
    ifi = Screen('GetFlipInterval', window);
    
    % Get the durations in frames
    % variable pre-stimulus noise mask duration
    ITIsFrames         = randi(round(PresParams.ITI_Range/ifi),nTrials,1);
    
    % fixed stimulus duration
    stimDurFrames      = round(PresParams.stimDurationInSecs/ifi);
    oddBallDurFrames   = round(PresParams.OddBallDurationInSecs/ifi);
    postStimFrames     = round(PresParams.PostStimTime/ifi);
    
    % post-stim max response period duration
    MaxRespFrames       = round(PresParams.MaxResponseTime /ifi);
    
    % creat noise masks
    noiseTextures = cell(PresParams.Nmasks ,1);
    for ii = 1:PresParams.Nmasks
        noiseTextures{ii} = Screen('MakeTexture', window, rand(PresParams.nsMaskSize(1),PresParams.nsMaskSize(2)));
        pgrStr = sprintf('Loading Noise Masks %g %%', floor(ii/PresParams.Nmasks *100));
        DrawFormattedText(window,pgrStr,'center','center',255,50);
        Screen('Flip',window);
    end
    
    % pre-make image textures
    imgTextures = cell(nTrials,1);
    for ii = 1:nTrials
        imgTextures{ii}=Screen('MakeTexture', window, tacs_er.Stimuli(stimNames{ii}));
        pgrStr = sprintf('Loading Stimuli  %g %%',floor(ii/nTrials*100));
        DrawFormattedText(window,pgrStr,'center','center',255,50);
        Screen('Flip',window);
    end
    
    %---------------------------------------------------------------------%
    % Participant Instructions
    %---------------------------------------------------------------------%
    InstStr = ['Instructions\n\n' ...
        'You will be presented with a ' PresParams.PreStimFixColorStr ' Fixation  Cross, ' ...
        'on a white noise background. The background will then change into an image, ' ...
        'indicating the begining of a trial. '...
        'On some trials, there will be a ' PresParams.OddBallColorStr ' dot in one of four possible ' ...
        'locations imposed on the image. ' ...
        'Your task is to respond with ''' PresParams.RespToLeft ''' if the dot is to the left of the image,' ...
        'or with ''' PresParams.RespToRight ''' if it is to the right of the image. \n' ...
        'You will have ' num2str(PresParams.MaxResponseTime) ' seconds to respond '...
        'as quickly and as accurately as possible. \n\n'...
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
    stimFlipDurSecs   = (stimDurFrames - 0.5) * ifi;
    oddBallOnsetSecs   = (stimDurFrames-oddBallDurFrames-0.5) * ifi;
    oddBallFlipDurSecs = (oddBallDurFrames-0.5) * ifi;
    
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
        [flip.VBLTimestamp, flip.StimulusOnsetTime, flip.FlipTimestamp, flip.Missed, flip.Beampos,] ...
            = Screen('Flip', window, vbl + 0.5*ifi);
        TimingInfo.stimPresFlip{tt}=flip;
        vbl = flip.VBLTimestamp;
        
        % if an oddball trial...
        if tacs_er.EncOddBallTrials(tt)            
            % Draw Oddball on top of stimulus for stimFlipDurSecs-oddballDursecs
            Screen('DrawTexture', window, imgTextures{tt}, [], [], 0);
            Screen('DrawDots', window, OB_Coords{tacs_er.EncOddBallLocs(tt)},PresParams.OddBallSize, PresParams.OddBallColor, [0 0],0);
            [flip.VBLTimestamp, flip.StimulusOnsetTime, flip.FlipTimestamp, flip.Missed, flip.Beampos,] ...
                = Screen('Flip', window, vbl + oddBallOnsetSecs);
            TimingInfo.oddBallPresFlip{tt}=flip;
            trialTime = GetSecs;
            vbl = flip.VBLTimestamp;
            
            % Draw Post-Stim Noise
            Screen('DrawTexture', window, noiseTextures{randi(PresParams.Nmasks)}, [], [], 0);
            [flip.VBLTimestamp, flip.StimulusOnsetTime, flip.FlipTimestamp, flip.Missed, flip.Beampos,] ...
                = Screen('Flip', window, vbl + oddBallFlipDurSecs);
            
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
        else
            % Draw Post-Stim Noise
            Screen('DrawTexture', window, noiseTextures{randi(PresParams.Nmasks)}, [], [], 0);
            [flip.VBLTimestamp, flip.StimulusOnsetTime, flip.FlipTimestamp, flip.Missed, flip.Beampos,] ...
                = Screen('Flip', window, vbl + stimFlipDurSecs);
            TimingInfo.postStimMaskFlip{tt}=flip;
            vbl = flip.VBLTimestamp;
            for ii = 1:(postStimFrames-1)
                Screen('DrawTexture', window, noiseTextures{randi(PresParams.Nmasks)}, [], [], 0);
                vbl  = Screen('Flip', window,vbl + 0.5*ifi);
            end
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
