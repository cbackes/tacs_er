
function debugResumeCommand

% Get keyboard number
[activeKeyboardID, ~, pauseKey, resumeKey] = getKeyboardOr10key;
% initialize Keyboard Queue
KbQueueCreate(activeKeyboardID);
% Start keyboard queue
KbQueueStart(activeKeyboardID);

[secs,key]=KbQueueWait2(activeKeyboardID);

e
%WaitTillResumeKey(resumeKey,activeKeyboardID)

end
%%

% Wait until Resume Key is pressed on the keyboard
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

% Check if the resume key has been pressed, and pause exection until resume
% key is pressed.
function CheckForPauseKey(pauseKey,resumeKey,activeKeyboardID)

[pressed,firstPress] = KbQueueCheck(activeKeyboardID);
if pressed
    if strcmp(pauseKey,KbName(firstPress));
        WaitTillResumeKey(resumeKey,activeKeyboardID)
    end
end
end