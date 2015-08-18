
function [keys RT k] = qkeys(startTime,dur,deviceNumber,term)
% If dur==-1 or if term==1, terminates after first keypress
% use term=1 to self-terminate trials with a max length of dur

myStart = GetSecs;

KbQueueCreate(deviceNumber);
KbQueueStart();

if ~exist('term','var')
    term = 0;
end
    
if (dur == -1)||(term==1)
    while 1
        [k.pressed, k.firstPress, k.firstRelease, k.lastPress, k.lastRelease]=...
            KbQueueCheck();
        if k.pressed
            break
        end
        if term
            if (GetSecs-startTime)>dur
                break
            end
        end
        WaitSecs(0.001);
    end
else
    WaitSecs('UntilTime',startTime+dur);
end

KbQueueStop();

if (dur ~= -1)&&(term==0)
    [k.pressed, k.firstPress, k.firstRelease, k.lastPress, k.lastRelease]=...
        KbQueueCheck();
end

if k.pressed == 0
    keys = 'noanswer';
    RT = 0;
else
    keys = KbName(k.firstPress);
    f = find(k.firstPress);
    RT = k.firstPress(f)-myStart;
end