
function [out,msg]=tACS_EncodingMain_CueResponse(tacs_er)
% core script for stimulus presentation on tACS Encoding.

% clear all the screens
close all;
sca;

% output structure
out = [];



% Presentation Parameters
PresParams = [];
PresParams.stimFrequency        = 6;
PresParams.stimDurationInCycles = 0.5;
PresParams.stimDurationInSecs   = 1/PresParams.stimFrequency*PresParams.stimDurationInCycles;
PresParams.cueDurationInSecs    = PresParams.stimDurationInSecs;
%PresParams.noiseFrameInterval   = 1; % define
%PresParams.waitFrameInterval    = 1; % define
PresParams.ITI_Range            = [1 1.5]; % variable ITI in secs
PresParams.MaxResponseTime      = 1;       % maximum to make perceptual decision

% noise mask size -> should be the size of all the stimuli
PresParams.nsMaskSize = [255 255];

% determine cue response mapping depending on subject number.
if mod(tacs_er.subjNum,2)
   PresParams.RespToCue1 = 'K';
   PresParams.RespToCue2 = 'L';
else
   PresParams.RespToCue1 = 'L';
   PresParams.RespToCue2 = 'K';
end

out.PresParams  = PresParams;
out.expInfo     = tacs_er;
%% temporary code
% load images

% Define colors
WHITE   = [1 1 1];
BLACK   = [0 0 0];
GREY    = WHITE/2;
RED     = [1 0 0.8];
BLUE    = [0.8 0 1];
PURPLE  = [1 0 1];

TimingInfo = [];
try    
    
    %---------------------------------------------------------------------%
    % Screen and additional presentation parameters
    %---------------------------------------------------------------------%
    % Get keyboard number
    [activeKeyID, ~, pauseKey, resumeKey] = getKeyboardOr10key;

    % initialie
    [window, windowRect] = initializeScreen;
    
    % pre-make noise masks
    Nmasks = 50;
    for ii = 1:Nmasks
        noiseTextures{ii} = Screen('MakeTexture', window, rand(PresParams.nsMaskSize(1),PresParams.nsMaskSize(2)));
    end
    
    % get stims
    stimNames = tacs_er.EncStimNames;
    nTrials   = numel(stimNames);
       
    % pre-make image textures
    imgTextures = cell(nTrials,1);
    for ii = 1:nTrials
        imgTextures{ii}=Screen('MakeTexture', window, tacs_er.Stimuli(stimNames{ii}));
    end
        
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
    ITIsFrames         = randi(round(PresParams.ITI_Range/ifi),nTrials);
    % fixed stimulus duration
    stimDurFrames      = round(PresParams.stimDurationInSecs/ifi);
    % post-stim max response period duration
    MaxRespFrames       = round(PresParams.MaxResponseTime /ifi);  
    
    % Set the line width for our fixation cross
    lineWidthPix = 4;
    
    % Sync us and get a time stamp
    vbl = Screen('Flip', window);
    
    
    %---------------------------------------------------------------------%
    % Participant Instructions    
    %---------------------------------------------------------------------%
    
    tstring = ['Instructions\n\n' ...
    'You will be presented with a PURPLE Fixation Cross. \n'...
    'The Fixation cross will then turn either RED or BLUE, with a background image of faces and landmarks. \n'...
    'Your task is to respond with ' PresParams.RespToCue1 ' for the RED Fixation and \n'...
    'with ' PresParams.RespToCue2 ' for BLUE Fixations. \n' ...
    'You will have ' num2str(PresParams.MaxResponseTime ) ' second to respond, and please do so \n'...
    'as quickly and as accurately as possible.\n'...    
    'Press R to begin the experiment.'];


    DrawFormattedText(Window,tstring,'center','center',255,50);
    Screen('Flip',Window);

    getKey(resumeKey,activeKeyID);
    
    % iterate through trials
    for tt = 1:nTrials
        
        % Maximum priority level
        topPriorityLevel = MaxPriority(window);
        Priority(topPriorityLevel);
                        
        respFlag  =0;
        Screen('DrawLines', window, fixCrossCoords,lineWidthPix, PURPLE, [0 0], 2);
        vbl  = Screen('Flip', window,vbl + (standardDurFrames - 0.5) * ifi);
        %WaitSecs(3);
        
        % Pre-stimulus noise mask (variable ITI)
        for ii = 1:ITIsFrames(tt)
            Screen('DrawTexture', window, noiseTextures{randi(Nmasks)}, [], [], 0);
            Screen('DrawLines', window, fixCrossCoords,lineWidthPix, PURPLE, [0 0], 2);
            vbl  = Screen('Flip', window,vbl + (PresParams.noiseFrameInterval - 0.5) * ifi);
        end
        
        % Draw Stimulus
        Screen('DrawTexture', window, imgTextures{tt}, [], [], 0);
        vbl  = Screen('Flip', window,vbl + (standardDurFrames-0.5) * ifi);
        x = vbl;
        
        % Keep stimulus on for stimDurationFrames, then draw noise mask
        Screen('DrawTexture', window, noiseTextures{randi(Nmasks)}, [], [], 0);
        Screen('DrawLines', window, fixCrossCoords,lineWidthPix, BLUE, [0 0], 2);
        vbl  = Screen('Flip', window,vbl + (stimDurationFrames - 0.5) * ifi);
        
        stimDur = vbl-x;
        
        % noise mask until response
        for ii = 1:(MaxRespFrames-1)
            % Get the current position of the mouse
            [~,~,buttons] = GetMouse(window,0);
            Screen('DrawTexture', window, noiseTextures{randi(Nmasks)}, [], [], 0);
            Screen('DrawLines', window, fixCrossCoords,lineWidthPix, BLUE, [0 0], 2);
            vbl  = Screen('Flip', window,vbl + (PresParams.noiseFrameInterval - 0.5) * ifi);
            
            if sum(buttons)>0
                respFlag = 1;
                break
            end
        end
        
        % if decision is made in the allowed time, probe confidence
        if respFlag
            % Draw the the confidence bar and wait for postDescIntervalSecs
            Screen('DrawLines', window, ConfidenceBarCoords,lineWidthPix, BLUE, [0 0], 2);
            DrawFormattedText(window, Conditions{ConfRespSide(tt)}{1}, LeftExtent-50, TopExtent-0.15*screenYpixels, WHITE);
            DrawFormattedText(window, Conditions{ConfRespSide(tt)}{2} , RightExtent-50, TopExtent-0.15*screenYpixels, WHITE);
            vbl  = Screen('Flip', window, vbl + (PresParams.waitFrameInterval - 0.5) * ifi);
            WaitSecs(PresParams.postDescIntervalSecs);
            SetMouse(xCenter,yCenter,window);
            for ii = 1:(MaxConfidenceDescDurationFrames-1)
                
                % Draw Confidence Bar
                Screen('DrawLines', window, ConfidenceBarCoords,lineWidthPix, BLUE, [0 0], 2);
                DrawFormattedText(window, Conditions{ConfRespSide(tt)}{1}, LeftExtent-50, TopExtent-0.15*screenYpixels, WHITE);
                DrawFormattedText(window, Conditions{ConfRespSide(tt)}{2} , RightExtent-50, TopExtent-0.15*screenYpixels, WHITE);
            
                % Get the current position of the mouse
                [mx, my, buttons] = GetMouse(window);
                
                % draw dot indicating position of mouse
                if mx>= RightExtent
                    Screen('DrawDots', window, [RightExtent yCenter], 10, WHITE, [], 2);
                    mx = RightExtent;
                elseif mx<= LeftExtent
                    Screen('DrawDots', window, [LeftExtent yCenter], 10, WHITE, [], 2);
                    mx = LeftExtent;
                else
                    Screen('DrawDots', window, [mx yCenter], 10, WHITE, [], 2);
                end
                
                % If there was a click, record. else continue to draw
                if sum(buttons)
                    if mx-xCenter<0
                        Pct = 1- (mx-LeftExtent)/(xCenter-LeftExtent);
                        buttDec = Conditions{ConfRespSide(tt)}{1};                        
                    elseif mx-xCenter>0
                        Pct = 1-(RightExtent-mx)/(RightExtent-xCenter);
                        buttDec = Conditions{ConfRespSide(tt)}{2};
                    else
                        buttDec = 'no Dec';
                        Pct = 0;
                    end
                    buttonPress = sprintf('Button Pressed: %s %.2g',buttDec,Pct);
                    DrawFormattedText(window,buttonPress,'center', BottomExtent+100, WHITE);
                    HideCursor();
                    vbl  = Screen('Flip', window, vbl + (PresParams.waitFrameInterval - 0.5) * ifi);
                    WaitSecs(1);
                    SetMouse(xCenter,yCenter,window);
                    break
                else
                    vbl  = Screen('Flip', window, vbl + (PresParams.waitFrameInterval - 0.5) * ifi);
                end
            end
        end
    end
    msg='allGood';
catch msg
    sca
    keyboard
end

% Clear the screen
Priority(0);
sca;
Screen('CloseAll');

save;
ShowCursor;

end

function ConfidenceBarCoords = barParams(xCenter, yCenter,screenXpixels,screenYpixels)

% Here we set the size of our confidence bar
BarLength = 0.5*screenXpixels; % 50% of the width screen
HeightOfBarWhisks = 0.05*screenYpixels; % 5% of the height of screen

LeftExtent  = xCenter-BarLength/2;
RightExtent = xCenter+BarLength/2 ;
BottomExtent = yCenter+HeightOfBarWhisks/2 ;
TopExtent   =  yCenter- HeightOfBarWhisks/2 ;

HorizontalBarCoords   = [LeftExtent RightExtent; yCenter yCenter];

LeftWhisk             = [LeftExtent LeftExtent ; BottomExtent TopExtent];
RightWhisk            = [RightExtent RightExtent ; BottomExtent TopExtent];
CenterWhisk           = [xCenter xCenter; yCenter+[1 -1]*HeightOfBarWhisks/3];

ConfidenceBarCoords         = [HorizontalBarCoords LeftWhisk RightWhisk CenterWhisk ];

end

function fixCrossCoords = fixCross(xCenter, yCenter,screenXpixels,screenYpixels)

% Here we set the size of fixation cross
fixCrossXlength = 0.02*screenXpixels; % 50% of the width screen
fixCrossYlength = 0.02*screenYpixels; % 50% of the width screen

LeftExtent  = xCenter-fixCrossXlength/2;
RightExtent = xCenter+fixCrossXlength/2 ;
BottomExtent = yCenter+fixCrossYlength/2 ;
TopExtent   =  yCenter- fixCrossYlength/2 ;

fixCrossXCoords   = [LeftExtent RightExtent; yCenter yCenter];
fixCrossYCoords   = [xCenter xCenter; BottomExtent TopExtent];

fixCrossCoords       = [fixCrossXCoords fixCrossYCoords];

end

