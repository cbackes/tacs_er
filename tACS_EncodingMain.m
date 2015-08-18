 
function msg=tACS_EncodingMain()
% core script for stimulus presentation on tACS Encoding.
%
%
%


% Clear the workspace
close all;
sca;

% Set Condition names
Conditions = [];
Conditions{1} = {'Scene','Face'};
Conditions{2} = {'Face','Scene'};

% Presentation Parameters
PresParams = [];
PresParams.stimFrequency = 6;
PresParams.stimDurationInCycles = 0.5;
PresParams.stimDurationInSecs  = 1/PresParams.stimFrequency*PresParams.stimDurationInCycles;
PresParams.noiseFrameInterval = 1;
PresParams.waitFrameInterval  = 1;
PresParams.PreStim_noiseDurationRangeSecs  = [0.5 1.5];    % variable ITI in secs
PresParams.PostStim_MaxNoiseDurationSecs   = 5;            % maximum of two seconds to make perceptual decision
PresParams.MaxConfidenceDescDurationSecs   = 10;
PresParams.postDescIntervalSecs            = 0.2;

% noise mask size -> should be the size of all the stimuli
nsMaskSize = [255 255];

% load images
img = [];
img{1} = imread('./stim/people/Jay_Leno.jpg');
img{2} = imread('./stim/people/Pamela_Anderson.jpg');
img{3} = imread('./stim/landmarks/Panama_Canal.jpg');
img{4} = imread('./stim/landmarks/Space_Needle.jpg');

nTrials = numel(img);
% Get Stimulus conditions; 1=face, 2=scene
SceFacConds = [ones(nTrials/2,1);2*ones(nTrials/2,1)]; % insert better code here :)

% Shuffle Trials
stimOrd = randperm(nTrials);
stims   = img(stimOrd);
SceFacConds =SceFacConds(stimOrd); 

nFaces      = sum(SceFacConds==1);
FaceStimId  = find(SceFacConds);
nScences    = sum(SceFacConds==2);
SceneStimId  = find(SceFacConds);

% Counterbalanced Left/Right Decisions with Stimulus Condition
temp1 = perms(FaceStimId);
temp2 = perms(SceneStimId);
LeftStimsIDs = [temp1(1:nFaces/2)'; temp2(1:nScences/2)'];

% 1 for Left-Scene/Right-Face; 2 for Left-Face/Right-Scene as determined by
% Conditions.
ConfRespSide = 2*ones(nTrials,1);
ConfRespSide(LeftStimsIDs)=1;

% Define colors
WHITE   = [1 1 1];
BLACK   = [0 0 0];
GREY    = WHITE/2;
RED     = [1 0 0];
BLUE    = [0 0 1];

try
    
    [window, windowRect] = initializeScreen;
    
    % pre-make noise masks
    Nmasks = 20;
    for ii = 1:Nmasks
        noiseTextures{ii} = Screen('MakeTexture', window, rand(nsMaskSize(1),nsMaskSize(2)));
    end
    
    imgTextures = [];
    for ii = 1:nTrials
        imgTextures{ii}=Screen('MakeTexture', window, stims{ii});
    end
    
    % Here we set the initial position of the mouse to be in the centre of the
    % screen
    %SetMouse(xCenter, yCenter, window);
    
    screenXpixels = windowRect(3);
    screenYpixels = windowRect(4);
    
    % Get the centre coordinate of the window
    [xCenter, yCenter] = RectCenter(windowRect);
    
    % Get coordinates for fixation cross
    fixCrossCoords = fixCross(xCenter, yCenter,screenXpixels,screenYpixels);
    
    % Get coordindates for confidence bar
    ConfidenceBarCoords = barParams(xCenter, yCenter,screenXpixels,screenYpixels);
    RightExtent = max(ConfidenceBarCoords(1,:));
    LeftExtent  = min(ConfidenceBarCoords(1,:));
    TopExtent   = max(ConfidenceBarCoords(2,:));
    BottomExtent= min(ConfidenceBarCoords(2,:));
    
    % Query the frame duration
    ifi = Screen('GetFlipInterval', window);
    
    % Get the durations in frames
    standardDurFrames                   = 1;
    stimDurationFrames                  = round(PresParams.stimDurationInSecs/ifi);
    PreStim_noiseDurationFrames         = randi(round(PresParams.PreStim_noiseDurationRangeSecs/ifi),nTrials);
    PostStim_MaxNoiseDurationFrames     = round(PresParams.PostStim_MaxNoiseDurationSecs/ifi);
    MaxConfidenceDescDurationFrames     = round(PresParams.MaxConfidenceDescDurationSecs/ifi);
    
    % Set the line width for our fixation cross
    lineWidthPix = 4;
    
    % Sync us and get a time stamp
    vbl = Screen('Flip', window);
    
    % iterate through trials
    for tt = 1:nTrials
        
        % Maximum priority level
        topPriorityLevel = MaxPriority(window);
        Priority(topPriorityLevel);
        
        respFlag  =0;
        Screen('DrawLines', window, fixCrossCoords,lineWidthPix, RED, [0 0], 2);
        vbl  = Screen('Flip', window,vbl + (standardDurFrames - 0.5) * ifi);
        WaitSecs(3);
        
        % Pre-stimulus noise mask (variable ITI)
        for ii = 1:PreStim_noiseDurationFrames(tt)
            Screen('DrawTexture', window, noiseTextures{randi(Nmasks)}, [], [], 0);
            Screen('DrawLines', window, fixCrossCoords,lineWidthPix, RED, [0 0], 2);
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
        for ii = 1:(PostStim_MaxNoiseDurationFrames-1)
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

