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