
function [ret_out,msg]=tACS_RetrievalMain(thePath)
% tACS experiment image recognition memory presentation script
% This script takes the 'tacs_er' structure that contains the stimuli and
% order that shall be used at retrieval. These images are then loaded and
% saved prior to presentation. The task is recognition memory.
% Subjects can respond to images as 'old', 'new', or 'unsure'. If the
% subjects respond old or new, they are then given a confidence scale to
% indicate their confidence in the decision.
%
%------------------------------------------------------------------------%
% Author:       Alex Gonzalez
% Created:      Aug 20th, 2015
% LastUpdate:   Sept 15, 2015
% TO DO :       (1)EEG markers
%------------------------------------------------------------------------%

% clear all the screens
close all;
sca;

% load the task
fileName = strcat(thePath.subjectPath,'/tacs_er_task.mat');
if exist(fileName,'file')
    load(fileName);
else
    error('no task created, must run encoding task first!')
end

%PsychDebugWindowConfiguration;

% Presentation Parameters
PresParams  = [];
PresParams.stimDurationInSecs   = 3;
PresParams.ITI_Range            = [1.5 2]; % variable ITI in secs
PresParams.MaxResponseTime      = 3;       % maximum to make recognition decision
PresParams.MaxConfDecInSecs     = 5; % max time to make confidene decision
PresParams.SaveEvNtrials        = 20; % save progress every X# of trials.
PresParams.lineWidthPix         = 5;       % Set the line width for our fixation cross
PresParams.dotColor             = [1 1 1];
PresParams.fixCrossColor        = [1 1 1];
PresParams.textColor            = [1 1 1];
PresParams.ConfidenceScale      = 0; % confidence scale flag
PresParams.ConfBarColor         = [0.2 0.1385 1];

% determine numbers for recognition decision
% depending on subject number and active Keyboard.
laptopResponseKeys = ['j','k','l'];
keypadResponseKeys = ['1','2','3'];
RespConds         = {'old','unsure','new'};

if mod(tacs_er.subjNum,2)
    responseMap = [1,2,3];
else
    responseMap = [3,2,1];
end
laptopResponseKeys = laptopResponseKeys(responseMap);
keypadResponseKeys = keypadResponseKeys(responseMap);
RespConds           = RespConds (responseMap);

PresParams.RespConds = RespConds;

% get the IDs of the trials
stimNames = tacs_er.RetStimNames;
nTrials   = numel(stimNames);

%%

% Initialize trial timing structure
TimingInfo = [];
TimingInfo.preStimFixFlip   = cell(nTrials,1);
TimingInfo.stimPresFlip     = cell(nTrials,1);
TimingInfo.trialRT          = nan(nTrials,1);
TimingInfo.trialKeyPress    = cell(nTrials,1);
TimingInfo.CondResp         = cell(nTrials,1);
TimingInfo.Confidence       = nan(nTrials,1);
TimingInfo.ConfResp         = cell(nTrials,1);

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
    
    % get correct mapping to keyboard
    if laptopKeyboardID==activeKeyboardID
        PresParams.RespButtons  = laptopResponseKeys;
    else
        PresParams.RespButtons = keypadResponseKeys;
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
    
    % Get the fixation time for every trial (random ITI)
    % in seconds.
    ITIsSecs        = randi(round(PresParams.ITI_Range/ifi),nTrials,1)*ifi;
    
    % post-stim max confidence response period duration
    MaxConfDescFrames  = round(PresParams.MaxConfDecInSecs /ifi);
    
    % pre-make image textures
    imgTextures = cell(nTrials,1);
    for ii = 1:nTrials
        imgTextures{ii}=Screen('MakeTexture', window, tacs_er.Stimuli(stimNames{ii}));
        loadStr = sprintf('Loading Stimuli  %g %%',floor(ii/nTrials*100));
        DrawFormattedText(window,loadStr,'center','center',255,50);
        Screen('Flip',window);
    end
    
    % Get coordindates for confidence bar
    [ConfidenceBarCoords] = ConfBarParams(xCenter, yCenter,screenXpixels,screenYpixels);
    RightExtent = max(ConfidenceBarCoords(1,:));
    LeftExtent  = min(ConfidenceBarCoords(1,:));
    TopExtent   = max(ConfidenceBarCoords(2,:));
    BottomExtent= min(ConfidenceBarCoords(2,:));
    CenterYpos  = ConfidenceBarCoords(2,1);
    
    %---------------------------------------------------------------------%
    % Participant Instructions
    %---------------------------------------------------------------------%
    if PresParams.ConfidenceScale
        InstString = ['Instructions\n\n' ...
            'You will be presented with images that you might recognized from the previous experiment. '...
            'Your task is to indentify which images were presented before and which ones are new '...
            'by pressing a button. For ' RespConds{1} ' images you will be pressing the '...
            PresParams.RespButtons(1) ' key, for ' RespConds{3} ' images you will press the '...
            PresParams.RespButtons(3) ' key. If you are unsure, you can press the '...
            PresParams.RespButtons(2) ' key. You will have ' num2str(PresParams.MaxResponseTime) ...
            ' seconds to respond, and please do so as quickly and as accurately as possible. '...
            'If you identify the image as old or new, you will also be indicating your '...
            'confidence in that decision by clicking on a scale with the mouse. \n\n'...
            'If no questions, \n'...'
            'Press ''' resumeKey ''' to begin the experiment.'];
    else
        InstString = ['Instructions\n\n' ...
            'You will be presented with images that you might recognized from the previous experiment. '...
            'Your task is to indentify which images were presented before and which ones are new '...
            'by pressing a button. For ' RespConds{1} ' images you will be pressing the '...
            PresParams.RespButtons(1) ' key, for ' RespConds{3} ' images you will press the '...
            PresParams.RespButtons(3) ' key. If you are unsure, you can press the '...
            PresParams.RespButtons(2) ' key. You will have ' num2str(PresParams.MaxResponseTime) ...
            ' seconds to respond, and please do so as quickly and as accurately as possible. '...
            'If you identify the image as old or new, you will also be indicating your '...
            'confidence by pressing ' PresParams.RespButtons(1) ', ' PresParams.RespButtons(2) ' and ' ...
            PresParams.RespButtons(3) ' for low, medium and high confidence, respectively. \n\n'...
            'If there are no questions, \n'...'
            'Press ''' resumeKey ''' to begin the experiment.'];
    end
    
    NoRespText = 'No response recorded, please answer quicker!';
    
    WrongKeyText = ['Please use only the following keys: \n \n' ...
        PresParams.RespButtons(1) ' for ' RespConds{1} '\n' ...
        PresParams.RespButtons(2) ' for ' RespConds{2} '\n' ...
        PresParams.RespButtons(3) ' for ' RespConds{3} '\n\n'...
        'Press ''' resumeKey ''' to continue'];
    
    WrongConfKeyText = ['Please use only the following keys: \n \n' ...
        PresParams.RespButtons(1) ' for low confidence \n' ...
        PresParams.RespButtons(2) ' for mid confidence \n' ...
        PresParams.RespButtons(3) ' for high confidence \n\n'...
        'Press ''' resumeKey ''' to continue'];
    
    DrawFormattedText(window,InstString, 'wrapat', 'center', 255, ...
        75, [],[],[],[],[xCenter*0.1,0,screenXpixels*0.8,screenYpixels]);
    Screen('Flip',window);
    
    % resume if Resume Key is pressed
    WaitTillResumeKey(resumeKey,activeKeyboardID)
    %%
    %---------------------------------------------------------------------%
    % Trials
    %---------------------------------------------------------------------%
    
    % Maximum priority level
    topPriorityLevel = MaxPriority(window);
    Priority(topPriorityLevel);
    
    % iterate through trials
    for tt = 1:nTrials
        
        % empty flip var
        flip     = [];
        
        % Pre-stimulus fixation (variable ITI).
        Screen('DrawLines', window, fixCrossCoords,PresParams.lineWidthPix, PresParams.fixCrossColor, [0 0], 2);
        [flip.VBLTimestamp, flip.StimulusOnsetTime, flip.FlipTimestamp, flip.Missed, flip.Beampos,] ...
            = Screen('Flip', window);
        TimingInfo.preStimFixFlip{tt}=flip;
        vbl = flip.VBLTimestamp;
        
        % Re-draw for last frame of ITI, taking into account the previous
        % presentation
        itiDur = vbl - 1.5*ifi + ITIsSecs(tt);
        Screen('DrawLines', window, fixCrossCoords,PresParams.lineWidthPix, PresParams.fixCrossColor, [0 0], 2);
        vbl = Screen('Flip', window, itiDur);
        
        % Checks if the Pause Key has been pressed.
        CheckForPauseKey(pauseKey,resumeKey,activeKeyboardID)
        KbQueueFlush(activeKeyboardID);
        
        % Draw Stimulus
        Screen('DrawTexture', window, imgTextures{tt}, [], [], 0);
        [flip.VBLTimestamp, flip.StimulusOnsetTime, flip.FlipTimestamp, flip.Missed, flip.Beampos,] ...
            = Screen('Flip', window, vbl + 0.5*ifi);
        TimingInfo.stimPresFlip{tt}=flip;
        trialTime = GetSecs;
        vbl = flip.VBLTimestamp;
        
        % Wait for Response
        [secs,key]=KbQueueWait2(activeKeyboardID,PresParams.stimDurationInSecs-2*ifi);
        if secs<inf && numel(key)==1
            TimingInfo.trialKeyPress{tt} = key;
            TimingInfo.trialRT(tt) = secs-trialTime;
            
            switch key
                case PresParams.RespButtons(1)
                    TimingInfo.CondResp{tt} = RespConds{1};
                case PresParams.RespButtons(2)
                    TimingInfo.CondResp{tt} = RespConds{2};
                case PresParams.RespButtons(3)
                    TimingInfo.CondResp{tt} = RespConds{3};
                otherwise
                    TimingInfo.CondResp{tt} = 'wrongkey';
                    DrawFormattedText(window, WrongKeyText, 'center' , 'center');
                    Screen('Flip', window, vbl + 0.5*ifi);
                    WaitTillResumeKey(resumeKey,activeKeyboardID)
            end
        else
            DrawFormattedText(window,NoRespText, 'center' , 'center');
            Screen('Flip', window, vbl + 0.5*ifi);
            WaitSecs(1);
        end
        
        if strcmp(TimingInfo.CondResp{tt},'old') || strcmp(TimingInfo.CondResp{tt},'new')
            if PresParams.ConfidenceScale
                % If response is old or new, probe confidence.
                confResp = 0; % confidence response flag.
                SetMouse(xCenter,CenterYpos,window);
                for ii = 1:(MaxConfDescFrames)
                    % Draw Confidence Bar
                    Screen('DrawLines', window, ConfidenceBarCoords,PresParams.lineWidthPix, PresParams.ConfBarColor, [0 0], 2);
                    DrawFormattedText(window, [' Confidence for ' TimingInfo.CondResp{tt}], 'center' , TopExtent-0.1*screenYpixels, PresParams.textColor);
                    DrawFormattedText(window, ' 0 ', LeftExtent-20, TopExtent-0.1*screenYpixels, PresParams.textColor);
                    DrawFormattedText(window, '100' , RightExtent-20, TopExtent-0.1*screenYpixels, PresParams.textColor);
                    
                    % Get the current position of the mouse
                    [mx, ~, buttons] = GetMouse(window);
                    
                    % draw dot indicating position of mouse
                    if mx>= RightExtent
                        Screen('DrawDots', window, [RightExtent CenterYpos], 15, PresParams.dotColor, [], 2);
                        mx = RightExtent;
                    elseif mx<= LeftExtent
                        Screen('DrawDots', window, [LeftExtent CenterYpos], 15, PresParams.dotColor, [], 2);
                        mx = LeftExtent;
                    else
                        Screen('DrawDots', window, [mx CenterYpos], 15, PresParams.dotColor, [], 2);
                    end
                    
                    % If there was a click, record. else continue to draw
                    if sum(buttons)
                        confResp = 1;
                        Pct = (mx-LeftExtent)/(RightExtent-LeftExtent);
                        TimingInfo.Confidence(tt) = Pct;
                        
                        buttonPress = sprintf('Response: %.2g',Pct);
                        DrawFormattedText(window,buttonPress,'center', BottomExtent+100, PresParams.textColor);
                        HideCursor();
                        Screen('Flip', window, vbl + 0.5* ifi);
                        SetMouse(xCenter,CenterYpos,window);
                        WaitSecs(0.5);
                        break
                    else
                        vbl  = Screen('Flip', window, vbl +  0.5*ifi);
                    end
                end
                if ~confResp
                    DrawFormattedText(window,NoRespText, 'center' , 'center');
                    Screen('Flip', window, vbl + 0.5*ifi);
                    WaitSecs(1);
                end
            else
                DrawFormattedText(window, [' Confidence for ' TimingInfo.CondResp{tt} '?'], 'center' , TopExtent-0.1*screenYpixels, PresParams.textColor);
                vbl=Screen('Flip', window, vbl + 0.5* ifi);
                % Wait for Response
                [secs,key]=KbQueueWait2(activeKeyboardID,PresParams.MaxConfDecInSecs-2*ifi);
                if secs<inf && numel(key)==1
                    switch key
                        case PresParams.RespButtons(1)
                            TimingInfo.ConfResp{tt} = 'low';
                        case PresParams.RespButtons(2)
                            TimingInfo.ConfResp{tt} = 'mid';
                        case PresParams.RespButtons(3)
                            TimingInfo.ConfResp{tt} = 'high';
                        otherwise
                            TimingInfo.ConfResp{tt} = 'wrongkey';
                            DrawFormattedText(window, WrongConfKeyText, 'center' , 'center');
                            Screen('Flip', window, vbl + 0.5*ifi);
                            WaitTillResumeKey(resumeKey,activeKeyboardID)
                    end
                else
                    DrawFormattedText(window,NoRespText, 'center' , 'center');
                    Screen('Flip', window, vbl + 0.5*ifi);
                    WaitSecs(1);
                end
            end
        end
        % save every PresParams.SaveEvNtrials
        if mod(tt,PresParams.SaveEvNtrials)==0
            tempName = sprintf('/tacs_er.s%i.test.%s.mat', thePath.subjNum, datestr(now,'dd.mm.yyyy.HH.MM'));
            save([thePath.subjectPath,tempName],'TimingInfo');
        end
     
        % Discard used image texture
        Screen('Close', imgTextures{tt})
        
    end
    %---------------------------------------------------------------------%
    % End of Experiment
    %---------------------------------------------------------------------%
    % store additional outputs
    ret_out = [];
    ret_out.PresParams  = PresParams;
    ret_out.expInfo     = tacs_er;
    ret_out.TimingInfo  = TimingInfo;
    
    % save
    fileName = 'tacs_er.test.mat';
    cnt = 0;
    while 1
        savePath = strcat(thePath.subjectPath,'/',fileName);
        if ~exist(savePath,'file')
            save(savePath,'ret_out')
            break
        else
            cnt = cnt+1;
            warning(strcat(fileName,' already existed.'))
            fileName = strcat('tacs_er.test','-',num2str(cnt),'.mat');
            warning(strcat('saving as ', fileName))
        end
    end
    
    InstString = ['End of Experiment.\n \n' ...
        'Press ''' resumeKey ''' to exit.'];
    
    DrawFormattedText(window,InstString, 'center', 'center', 255, 40);
    Screen('Flip',window);
    WaitTillResumeKey(resumeKey,activeKeyboardID)
    KbQueueStop(activeKeyboardID);
    
    msg='allGood';
catch msg
    sca
    ShowCursor
    keyboard
end

% Clear the screen
Priority(0);
sca;
Screen('CloseAll');
ShowCursor;

end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% auxiliary functions and definitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------%
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

%-------------------------------------------------------------------------%
% ConfBarParams
% Confidence Bar Coordinates
%-------------------------------------------------------------------------%
function [ConfidenceBarCoords] = ConfBarParams(xCenter, yCenter,screenXpixels,screenYpixels)

% Here we set the size of our confidence bar
BarLength = 0.5*screenXpixels; % 50% of the width screen
HeightOfBarWhisks = 0.05*screenYpixels; % 5% of the height of screen

LeftExtent  = xCenter-BarLength/2;
MidLeftExtent = xCenter-BarLength/4;
RightExtent = xCenter+BarLength/2 ;
MidRightExtent = xCenter+BarLength/4;
BottomExtent = yCenter+HeightOfBarWhisks/2 ;
TopExtent   =  yCenter- HeightOfBarWhisks/2 ;

HorizontalBarCoords   = [LeftExtent RightExtent; yCenter yCenter];

LeftWhisk             = [LeftExtent LeftExtent ; BottomExtent TopExtent];
RightWhisk            = [RightExtent RightExtent ; BottomExtent TopExtent];
MidLeftWhisk          = [MidLeftExtent MidLeftExtent ; BottomExtent TopExtent];
MidRightWhisk         = [MidRightExtent MidRightExtent ; BottomExtent TopExtent];
CenterWhisk           = [xCenter xCenter ; BottomExtent TopExtent];

ConfidenceBarCoords   = [HorizontalBarCoords LeftWhisk RightWhisk ...
    MidLeftWhisk MidRightWhisk CenterWhisk];
end
