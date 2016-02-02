% test images for xdiva

nXDivaPreludeFrames = 10;
stimDurationInFrames = 5;
nXDivaFrames          = 80;
stimFrequency       = 6;
stimSize            = [300 300];

ZeroPhaseBlackFrameIDs       = (nXDivaPreludeFrames+1):...
    stimDurationInFrames*2:nXDivaFrames;
ZeroPhaseBlackFrameIDs(stimFrequency+1:end)=[];

ZeroPhaseWhiteFrameIDs       = ZeroPhaseBlackFrameIDs+stimDurationInFrames;

images = zeros([stimSize,1,3], 'uint8');
images(:,:,1,1) = 1*ones(stimSize);
images(:,:,1,2) = 255*ones(stimSize);

fileName = 'test_phase_xdiva_';
for tt = 1:6           
        
    imageSequence = zeros(nXDivaFrames,1,'uint32');    
    imageSequence(1)=2;
    imageSequence(ZeroPhaseBlackFrameIDs+tt*2-2)=1;
    imageSequence(ZeroPhaseWhiteFrameIDs+tt*2-2)=2;
    
    save(['~/Desktop/' fileName num2str(tt) '.mat'],'images','imageSequence')
end

% 1 frame precision aka. +-18 degrees.

