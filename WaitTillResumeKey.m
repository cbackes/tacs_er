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
