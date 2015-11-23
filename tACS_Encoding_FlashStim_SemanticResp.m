%function [enc_out,msg]=tACS_Encoding_FlashStim_SemanticResp(thePath)
% tACS stimulus encoding presentation script.
%
% this scripts presents stimuli to be encoded for a time depending on
% frequency of stimulation. the stimulus will flash at that frequency
% for a given duration. subjects will be making a shallow encoding decision
% such as scene or face and respond accordingly.
% stimulus order and stimulus conditions are pre-determined by the
% tacs_er make list  oddball. thePath indicates the path of the tacs_er
% structure.
%
% subjects performance on the cue color identification is stored peridiocally.
%
% tACS Marker Codes:
%       1 -> start presentation
%       2 -> pre stimulus fixation
%       3 -> face stimuli
%       4 -> scene stimuli
%       5 -> end of presentation
%------------------------------------------------------------------------%
% Author:       Alex Gonzalez
% Created:      Nov 5th, 2015
% LastUpdate:   Nov 17th, 2015
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
    tacs_er =tACS_ER_makelist_SemanticDesc(thePath);
end

% For debugging:
% PsychDebugWindowConfiguration

% Presentation Parameters
PresParams = [];
switch tacs_er.exptType
    case 'behav_v10'
        PresParams.tACSstim             = 1;
    case 'tacs_enc'
        PresParams.tACSstim             = 1;
    otherwise
        error('task not supported')
end

PresParams.stimFrequency        = 6;
PresParams.stimDurationInCycles  = 0.5;
PresParams.stimDurationInSecs   = 1/PresParams.stimFrequency*PresParams.stimDurationInCycles;
PresParams.totalStimDuration    = 1;
PresParams.nCycles              = PresParams.stimFrequency*PresParams.totalStimDuration;
PresParams.stimPresCycles       = ones(PresParams.nCycles,1); % stimuli will be presented every cycle
PresParams.preStartTime         = 5; % in seconds.

PresParams.ITI_Range            = [0.5 1]; % variable ITI in secs
PresParams.PostStimTime         = 0.5;     % time after stim (with no fixation).
PresParams.PreStimFixColor      = [1 1 1];
PresParams.PreStimFixColorStr   = 'WHITE';
PresParams.lineWidthPix         = 5;      % Set the line width for our fixation cross
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

stimNames = tacs_er.EncStimNames;
nTrials   = numel(stimNames);

% stimulation parameters
if PresParams.tACSstim
    PresParams.StarStimIP           = '10.0.0.42';
    PresParams.StimTemplate         = 'AG-tACS-Fz1mA6Hz-FC12AF34Ret';
    PresParams.preStartTime         = 30; % in seconds.
    
    addpath(genpath('../MatNIC_v2.5/'))
    try
          [ret, status, socket] = MatNICConnect(PresParams.StarStimIP);
          if ret<0
              error(status)
          end
          [ret,LSLOutlet]=MatNICMarkerConnectLSL('alexLSL');          
          if ret<0
              error('could not connect to streaming layer')
          end
          ret = MatNICLoadTemplate(PresParams.StimTemplate,socket);
          if ret<0
              error('could not load template')
          end

          [~,status] = MatNICQueryStatus(socket);
          assert(strcmp(status, 'CODE_STATUS_STIMULATION_READY'),'aborting, stimulation not ready.')
          
          % only works during stimulation
%           [ret,StartImpedances] = MatNICGetImpedance(socket);
%           if ret<0
%               error(' could not obtain impedances');
%           end
%           disp('Initial Impedances');
%           fprintf(StartImpedances);
          
          ret = MatNICStartStimulation(socket);
          if ret<0
              error('could not start stimulation')
          end
    catch msg
        MatNICMarkerCloseLSL(LSLOutlet);
        MatNICAbortStimulation(socket);
        sca
        keyboard;
    end
end    

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
        PresParams.RespToCond1 = laptopResponseKeys(1);
        PresParams.RespToCond2 = laptopResponseKeys(2);
    else
        PresParams.RespToCond1 = keypadResponseKeys(1);
        PresParams.RespToCond2 = keypadResponseKeys(2);
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
            
    %---------------------------------------------------------------------%
    % Participant Instructions
    %---------------------------------------------------------------------%
    InstStr = ['Instructions\n\n' ...
        'A centrally presented ' PresParams.PreStimFixColorStr ' Fixation  Cross will indicate ' ...
        'the start of a new event. Subsequently, you will presented with an image that will be '...
        'shown on and off for ' num2str(PresParams.totalStimDuration) ' seconds. Your task will be ' ...
        'to indicate if the image is a SCENE or a FACE by pressing a button. \n\n'...
        'For FACES press the ''' PresParams.RespToCond1 ''' key. \n'...
        'For SCENES press the ''' PresParams.RespToCond2 ''' key. \n\n'...        
        'If no questions, \n'...
        'Press ''' resumeKey ''' key to begin the experiment.'];
    
    PauseStr = ['Rest Pause. \n\n '...
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
    
    if PresParams.tACSstim
        % start of experiment marker
        sendMarker(1,LSLOutlet)
    end
    
    % Draw blank for a bit    
    Screen('Flip', window);
    WaitSecs(PresParams.preStartTime);    
    
   
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
                
        if PresParams.tACSstim
            % fixation marker
            sendMarker(2,LSLOutlet)
        end
        
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
                [flip.VBLTimestamp, flip.StimulusOnsetTime, flip.FlipTimestamp, flip.Missed, flip.Beampos,] ...
                    = Screen('Flip', window, vbl + 0.5*ifi);
                trialTime = GetSecs;  
            else
                [flip.VBLTimestamp, flip.StimulusOnsetTime, flip.FlipTimestamp, flip.Missed, flip.Beampos,] ...
                    = Screen('Flip', window, vbl + stimFlipDurSecs);                                
            end
            vbl = flip.VBLTimestamp;
            % blank screen for stimFlipDurSecs
            vbl  = Screen('Flip', window, vbl + stimFlipDurSecs);
            TimingInfo.stimPresFlip{tt,ii}=flip;
            
            if PresParams.tACSstim
                % stimulus markers
                if tacs_er.EncStimType(tt)==1
                    sendMarker(3,LSLOutlet)
                else
                    sendMarker(4,LSLOutlet)
                end
            end
        end
        [pressed,firstPress] = KbQueueCheck(activeKeyboardID);
        if pressed
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
    Screen('Flip',window);

    if PresParams.tACSstim
        % end of experiment marker
        sendMarker(5,LSLOutlet)
        MatNICMarkerCloseLSL(LSLOutlet);
        MatNICAbortStimulation(socket);
    end
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
            warning(strcat('saving as ', fileName))
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
    MatNICAbortStimulation(socket);
    MatNICMarkerCloseLSL(LSLOutlet);
    sca
    keyboard
end

% Clear the screen
Priority(0);
sca;
KbQueueStop(activeKeyboardID);
Screen('CloseAll');
ShowCursor;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% auxiliary functions and definitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % %-------------------------------------------------------------------------%
% % fixCrossCoords
% % Set Fixation Cross Coordinates
% %-------------------------------------------------------------------------%
% function fixCrossCoords = fixCross(xCenter, yCenter,screenXpixels,screenYpixels)
% 
% fixCrossXlength = max(0.02*screenXpixels,0.02*screenYpixels); % max of 2% screen dims
% fixCrossYlength = fixCrossXlength;
% 
% LeftExtent  = xCenter-fixCrossXlength/2;
% RightExtent = xCenter+fixCrossXlength/2 ;
% BottomExtent = yCenter+fixCrossYlength/2 ;
% TopExtent   =  yCenter- fixCrossYlength/2 ;
% 
% fixCrossXCoords   = [LeftExtent RightExtent; yCenter yCenter];
% fixCrossYCoords   = [xCenter xCenter; BottomExtent TopExtent];
% 
% fixCrossCoords       = [fixCrossXCoords fixCrossYCoords];
% 
% end

% %-------------------------------------------------------------------------%
% % WaitTillResumeKey
% % Wait until Resume Key is pressed on the keyboard
% %-------------------------------------------------------------------------%
% function WaitTillResumeKey(resumeKey,activeKeyboardID)
% 
% KbQueueFlush(activeKeyboardID);
% while 1
%     [pressed,firstPress] = KbQueueCheck(activeKeyboardID);
%     if pressed
%         if strcmp(resumeKey,KbName(firstPress));
%             break
%         end
%     end
%     WaitSecs(0.1);
% end
% KbQueueFlush(activeKeyboardID);
% end

% %-------------------------------------------------------------------------%
% % CheckForPauseKey
% % Check if the resume key has been pressed, and pause exection until resume
% % key is pressed.
% %-------------------------------------------------------------------------%
% function CheckForPauseKey(pauseKey,resumeKey,activeKeyboardID)
% 
% [pressed,firstPress] = KbQueueCheck(activeKeyboardID);
% if pressed
%     if strcmp(pauseKey,KbName(firstPress));
%         WaitTillResumeKey(resumeKey,activeKeyboardID)
%     end
% end
% end

% %-------------------------------------------------------------------------%
% % sendMarker
% % sends LabStreamLayer markers to NIC on stimulation computer
% %-------------------------------------------------------------------------%
% function sendMarker(marker,LSLOutlet)
%     ret=MatNICMarkerSendLSL(marker,LSLOutlet);
%     if ret<0
%         warning('error sending marker')
%     end
% end