
function [tacs_er] = tACS_ER_makelist_SemanticDesc_xDivaStim(thePath)
% make lists for tacs encoding & retrieval (tacs_er) experiment
% This design:
% 1) One encoding block with flashing stimuli.
% 2) Responses will be made to the stimulus itself, no cue conditions.
%
% Important variables
% subjNum       -> subject  ID
% thePath       -> path for all stimulus and data directories
%
% nEncStim      -> # of UNIQUE encoding stimuli
% nEncTrials    -> # of TOTAL trials (should be multiple of nEncStim)
% nEncCond      -> # of conditions at encoding (that are critical for retrieval)
% nEncBlocks    -> # of encoding blocks
%
% nRetTrials    -> # of TOTAL retrieval trials
% nFoilTrials   -> # of NEW trials (foils)
% nRetConds     -> # of retrieval conditions
% nRetBlocks    -> # of retrieval blocks
% nRetCondTrials-> VECTOR of # trials per Retrieal Condition
% maxNumConsecutiveOld -> MAX # of consecutive old trials
%
% RandStream    -> random stream seed
%
% EncStimType       -> VECTOR trial of indicating type (1  faces, 2 scene)
% EncStimTypeIDs    -> {face,scene}
% EncCondCodeIDs    -> {'Face0','Face90','Face180','Scn0','Scn90','Scn180'}
% EncCondTrialCode  -> VECTOR trial encoding indicating EncCondCodeID
% EncStimNames      -> VECTOR trial stimulus ID
% EncStimUniqueIDs  -> VECTOR trial with unique ID for each stim
% EncBlockID        -> VECTOR trial block ID
%
% RetCondIDs        -> {'OldFace','OldScene','NewFace','NewScene'};
% RetCondTrialCode  -> VECTOR trial 1:4, indicating RetCondID
% RetStimNames      -> VECTOR trial of retrieaval IDs
% RetBlockID        -> VECTOR of trial block ID
%
% Stimuli           -> key, matrix map of stimuli's names to image matrix.

%------------------------------------------------------------------------%
% Author:       Alex Gonzalez
% Created:      Jan, 2016
% LastUpdate:   April 6, 2016
%------------------------------------------------------------------------%

%% Set-up

% parameter list
nUniqueFaceStimuli = 150;
nUniqueSceneStimuli = 150;
nEncStim    = nUniqueSceneStimuli+ nUniqueFaceStimuli;
nEncPhases = 5;                 % constrained by xdiva design.
nEncConds  = 2*nEncPhases;     % faces/scenes x phase (5 different phases)
EncPhases  = 0:(360/nEncPhases):359;

switch thePath.exptType
    case {'tacs_enc_xdiva'}
        stimSize        = [600 600];
        nEncBlocks      = 1;
    otherwise
        error('incorrect expt!')
end

nEncTrials = nEncBlocks*nEncStim;
nRetTrials = nEncStim*1.5;   % encoding trials + retrieval trials
nFoilTrials =nRetTrials-nEncStim;
nRetConds  = 4;     % old and new conditons per stim type (Face/scene)
nRetBlocks = nEncBlocks;    % blocks
maxNumConsecutiveOld = 8;   % maximum number of old trials in a row

% % reset stream (to avoid duplicate lists)
% s = RandStream.create('mt19937ar','seed',sum(100*clock));
% if strfind(version,'R2014b')>0
%     RandStream.setGlobalStream(s)
% end

% changed on 4/6/16 to make it work on different matlab versions. 
s = RandStream.create('mt19937ar','seed',thePath.subjNum);
RandStream.setGlobalStream(s)

%% load data names
% scenes
cd(fullfile(thePath.stim,'notfamous_scenes/cropped'));
temp = dir('*.jpg');
count = 1;
ScenesMat = cell(numel(temp),1);
for n = 1:size(temp,1)
    ScenesMat{n} = zeros(stimSize(1),stimSize(2),'uint8');
    try
        x=imread(temp(n).name); %check if readable
        ScenesMat{n} = imresize(x(:,:,1),2); % only resizing case for xdiva
        SceneNames(count) = temp(n);
        count = count + 1;
    catch
        fprintf(['\n' temp(n).name ' unreadable\n']);
    end
end
%get the names in a cell array and shuffle
SceneNames = {SceneNames.name}';
[SceneNames,index] = Shuffle(SceneNames);
ScenesMat = ScenesMat(index);


% faces
cd(fullfile(thePath.stim,'notfamous_people/cropped'));
temp = dir('*.jpg');
count = 1;
FacesMat = cell(numel(temp),1);
for n = 1:size(temp,1)
    try
        x=imread(temp(n).name);  %check if readable
        FacesMat{n} = imresize(x(:,:,1),2); % only resizing case for xdiva
        FaceNames(count) = temp(n);
        count = count + 1;
    catch
        fprintf(['\n' temp(n).name ' unreadable\n']);
    end
end
%get the names in a cell array and shuffle
FaceNames = {FaceNames.name}';
[FaceNames,index] = Shuffle(FaceNames);
FacesMat = FacesMat(index);

StimObj = containers.Map( [FaceNames; SceneNames], [FacesMat; ScenesMat]);

%% Encoding
% Equal numbers of faces and scenes: #N encoding trials/2 per type
% Each block will have one of each of stimuli, for a total of # of blocks x
% # of unique stimuli.
%
% Oddball can be in 4 locations
% Each stimuli will be associated with one, and then allocated across blocks
% such that the same stimuli appears only once with the associted oddball
%
% EncStimType -> 1 for faces, 2 for scenes
% EncStimCue  -> 1 or 2
% EncStimNames -> name of the stimuli as it apperas on the database
% EncBlockID   -> block ID

% order of scenes and faces per block
EncStimTypeIDs = {'Face','Scene'};
EncStimType = zeros(nEncTrials,1);
EncBlockID  = zeros(nEncTrials,1);
EncStimNames = cell(nEncTrials,1);
nEncTrialsBlock = nEncTrials/nEncBlocks;
nFacesEncBlock  = nEncTrialsBlock/2;
nScenesEncBlock = nFacesEncBlock;

% Get Encoding Stimuli
temp = Shuffle(FaceNames);
EncFaces = temp(1:nFacesEncBlock);
temp = Shuffle(SceneNames);
EncScenes = temp(1:nScenesEncBlock);
EncStimUniqueIDs = zeros(nEncTrials,2); EncStimUniqueIDs(:,1)=1:nEncTrials;

% assign face/scence trials
for bb = 1:nEncBlocks
    TrialBlockID = ((bb-1)*nEncTrialsBlock+1):bb*nEncTrialsBlock;
    FaceTrials  = datasample(TrialBlockID,nFacesEncBlock,'replace',false);
    SceneTrials = setdiff(TrialBlockID,FaceTrials);
    EncBlockID(TrialBlockID) = bb;
    EncStimType(FaceTrials)  = 1;
    EncStimType(SceneTrials) = 2;
    
    % assign individial stimuli to trials
    FaceBlockUniqueID        = randperm(nFacesEncBlock,nFacesEncBlock);
    ScnBlockUniqueID         = randperm(nScenesEncBlock,nScenesEncBlock);
    EncStimUniqueIDs(FaceTrials,2) = FaceBlockUniqueID;
    EncStimUniqueIDs(SceneTrials,2) = ScnBlockUniqueID+nFacesEncBlock;
    
    EncStimNames(FaceTrials)  = EncFaces(FaceBlockUniqueID);
    EncStimNames(SceneTrials) = EncScenes(ScnBlockUniqueID);
    
end

%% assign encoding conditions
% max 12 condtions for a 2 x 6 design
%
% ** Note, not counterbalancing by stimuli at the moment, that would have
% to be done across subjects. **
nEncCondTrials      = nEncTrials/nEncConds;
nEncCondBlockTrials = nEncCondTrials/nEncBlocks;
EncCondTrialCode    = zeros(nEncTrials,1);
faceConds = 1:(nEncConds/2);
sceneConds = (nEncConds/2+1):nEncConds;
EncCondCodeIDs = cell(nEncConds,1);
EncCondCodeIDs(faceConds) = strcat('Face',cellfun(@num2str,num2cell(EncPhases'),'UniformOutput',false));
EncCondCodeIDs(sceneConds) = strcat('Scn',cellfun(@num2str,num2cell(EncPhases'),'UniformOutput',false));

if mod(nEncCondBlockTrials,1)~=0
    error('uneven encoding condition trials')
end

for bb = 1:nEncBlocks
    TrialBlockID = ((bb-1)*nEncTrialsBlock+1):bb*nEncTrialsBlock;
    
    % for face trials
    nFaceConds=numel(faceConds);
    availableFaceTrials = intersect(find(EncStimType==1),TrialBlockID);
    for ii = 1:nFaceConds
        condiTrials=datasample(availableFaceTrials,nEncCondBlockTrials,'replace',false);
        EncCondTrialCode(condiTrials) = faceConds(ii);
        availableFaceTrials = setdiff(availableFaceTrials,condiTrials);
    end
    % for scenes trials
    nSceneConds=numel(sceneConds);
    availableSceneTrials = intersect(find(EncStimType==2),TrialBlockID);
    for ii = 1:nSceneConds
        condiTrials=datasample(availableSceneTrials,nEncCondBlockTrials,'replace',false);
        EncCondTrialCode(condiTrials) = sceneConds(ii);
        availableSceneTrials = setdiff(availableSceneTrials,condiTrials);
    end
end

% check for equal number of conditions throughout
assert(sum(histc(EncCondTrialCode,1:nEncConds)==nEncCondTrials)==nEncConds,'unequal number of trials produced')

%% assign retrieval conditions
% Condtions
% 1 -> old face
% 2 -> old scene
% 3 -> new face
% 4 -> new scene

RetCondIDs        = {'OldFace','OldScene','NewFace','NewScene'};
nRetTrialsBlock   = nRetTrials/nRetBlocks;
nRetTrialTypesBlock = [nUniqueFaceStimuli nUniqueSceneStimuli nFoilTrials/2 nFoilTrials/2]'/nRetBlocks;
nRetCondTrials      = nRetTrialTypesBlock*nRetBlocks;

% assign condition such that each is equally likely per block.
% also retry the sampling such that there are not more than maxNumConsecutiveOld
% i.e., no that many consecutive old trials.
counter = 1;
while true
    RetCondTrialCode = zeros(nRetTrials,1);
    RetBlockID = zeros(nRetTrials,1);
    for rr = 1:nRetBlocks
        TrialBlockID = ((rr-1)*nRetTrialsBlock+1):rr*nRetTrialsBlock;
        availableTrials = TrialBlockID;
        for cc = 1:nRetConds
            condTrials  = datasample(availableTrials,nRetTrialTypesBlock(cc),'replace',false);
            RetCondTrialCode(condTrials)=cc;
            availableTrials  = setdiff(availableTrials,condTrials);
        end
        RetBlockID(TrialBlockID)=rr;
    end
    
    % quick way to the number of consecutive old trials
    temp=conv(ones(maxNumConsecutiveOld,1),double((RetCondTrialCode<=2)));
    if sum(temp>=maxNumConsecutiveOld)==0
        break
    end
    counter=counter+1;
    if counter >10000
        error('could not arrive at a stable sequence in 10000 attempts')
    end
end
% make sure that each condition is equally represented
assert(sum(histc(RetCondTrialCode,1:nRetConds)==nRetCondTrials)==nRetConds,'unequal number of trials produced')
%% assign stimuli based on retrieval condition

RetStimNames = cell(nRetTrials,1);
temp         = setdiff(FaceNames,EncFaces);
FoilFaces    = Shuffle(temp(1:nRetCondTrials(3)));
temp         = setdiff(SceneNames,EncScenes);
FoilScenes   = Shuffle(temp(1:nRetCondTrials(4)));

RetStimNames(RetCondTrialCode==1) = Shuffle(EncFaces); % old faces
RetStimNames(RetCondTrialCode==2) = Shuffle(EncScenes); % old scenes
RetStimNames(RetCondTrialCode==3) = Shuffle(FoilFaces); % new faces
RetStimNames(RetCondTrialCode==4) = Shuffle(FoilScenes); % new scenes

%% store necessary variables
tacs_er = [];
tacs_er.subjNum     = thePath.subjNum;
tacs_er.exptType    = thePath.exptType;
tacs_er.thePath     = thePath;

% parameterts for making stimuli list
tacs_er.nEncStim = nEncStim;
tacs_er.nEncTrials = nEncTrials;
tacs_er.nEncConds  = nEncConds;
tacs_er.nEncPhases = nEncPhases;
tacs_er.EncPhases  = EncPhases;
tacs_er.nEncBlocks = nEncBlocks;

tacs_er.stimSize    = stimSize;
tacs_er.nRetTrials = nRetTrials;
tacs_er.nFoilTrials = nFoilTrials;
tacs_er.nRetConds  = nRetConds;
tacs_er.nRetBlocks = nRetBlocks;
tacs_er.maxNumConsecutiveOld = maxNumConsecutiveOld;

tacs_er.RandStream = s;

% encoding
tacs_er.EncStimType     = EncStimType;
tacs_er.EncStimTypeIDs  = EncStimTypeIDs;
tacs_er.EncCondTrialCode= EncCondTrialCode;
tacs_er.EncCondCodeIDs  = EncCondCodeIDs;
tacs_er.EncStimNames    = EncStimNames;
tacs_er.EncStimUniqueIDs=EncStimUniqueIDs;
tacs_er.EncBlockID      = EncBlockID;
x = tacs_er.EncCondTrialCode;
condNums = x;
condNums(x>nEncPhases)=x(x>nEncPhases)-nEncPhases;
tacs_er.EncTrialPhaseConds  = condNums;
tacs_er.EncTrialPhase   = tacs_er.EncPhases(condNums);

% retrieval
tacs_er.RetCondIDs      = RetCondIDs;
tacs_er.RetCondTrialCode= RetCondTrialCode;
tacs_er.RetStimNames    = RetStimNames;
tacs_er.RetBlockID      = RetBlockID;
tacs_er.nRetCondTrials  = nRetCondTrials;

% sanity check: verify that old labels match stimuli
assert(sum(ismember(tacs_er.RetStimNames,tacs_er.EncStimNames) == ...
    (tacs_er.RetCondTrialCode<=2))==tacs_er.nRetTrials,'Labels do not match stimuli')

% encodign conditions at retrieval
[~,i]=ismember(tacs_er.RetStimNames,tacs_er.EncStimNames);
tacs_er.EncCondAtRet = zeros(tacs_er.nRetTrials,1);
tacs_er.EncCondAtRet(i>0) = tacs_er.EncCondTrialCode(i(i>0));

% Stimulus Object
tacs_er.Stimuli = StimObj;
cd(thePath.main)

%% Stimulus Presentation Parameters
PresParams                      = [];
PresParams.VideoFrameRate       = 60;
PresParams.nXDivaFrames         = 90;
PresParams.nXDivaPreludeFrames  = 10; % leave first 10 spots blank
PresParams.stimFrequency        = 6;
PresParams.stimDurationInCycles  = 0.5;
PresParams.totalTrialDuration   = 60; % in seconds;
PresParams.stimDurationInSecs   = 1/PresParams.stimFrequency*PresParams.stimDurationInCycles;
PresParams.nCycles              = PresParams.stimFrequency*PresParams.totalTrialDuration;
PresParams.stimPresCycles       = ones(PresParams.nCycles,1); % stimuli will be presented every cycle
PresParams.DegreesPerFrame      = 360/(PresParams.VideoFrameRate/PresParams.stimFrequency);
PresParams.VideoFrameDurSecs    = 1/PresParams.VideoFrameRate;
PresParams.stimDurationInFrames = PresParams.VideoFrameRate*PresParams.stimDurationInCycles/PresParams.stimFrequency;

%% Save individual images for xDiva presentation and prepare sequence


PresParams.FixCrossFrameIDs           = 1; % on the first frame
PresParams.FixCrossMinNFrames         = 20;% corresponding to 333ms or 2 cycles of 6Hz
PresParams.InstructionsLength         = 2*60*PresParams.VideoFrameRate;
% bug here noticed on the mod function. it should be thePath.subjectNum! 
PresParams.InstructionSet             = mod(thePath.subjectPath,2)==0 +1;
tacs_er.EncFaceRespID                 = mod(thePath.subjectPath,2)==0 +1;
tacs_er.EncSceneRespID                = mod(thePath.subjectPath,2)==1 +1;

PresParams.tACSPreTaskWarmUpDur       = 3*60*PresParams.VideoFrameRate;     
tacs_er.PresParams = PresParams;
% Save into subjects path
fileName = 'tacs_er_xdiva.task.mat';
if ~exist([thePath.subjectPath '/' fileName],'file')
    save([thePath.subjectPath '/' fileName],'tacs_er')
    mkdir([thePath.subjectPath '/xdiva/'])
    overwriteFlag=1;
else
    warning('tacs task for this subject already created')
    while 1
        s = input('Overwrite? (Y/N)','s');
        if strcmp(s,'Y')
            overwriteFlag=1;
            save([thePath.subjectPath '/' fileName],'tacs_er')
            break
        elseif strcmp(s,'N')
            break
        end
    end
end

if overwriteFlag
    ZeroPhaseImgFrameIDs       = (PresParams.FixCrossMinNFrames+1):...
        PresParams.stimDurationInFrames*2:PresParams.nXDivaFrames;
    ZeroPhaseImgFrameIDs(PresParams.stimFrequency+1:end)=[];
    
    ZeroPhaseBlankFrameIDs       = ZeroPhaseImgFrameIDs+PresParams.stimDurationInFrames;
    
    images = zeros([stimSize,1,3], 'uint8');
    images(:,:,1,2) = 128*ones(stimSize,'uint8');
    images(:,:,1,3) = fixCrossImage(stimSize(1),round(stimSize(1)/20));
    
    
    % Save individual stims into subjects path
    fileName = 'tacs_enc_xdiva_';
    for tt = 1:nEncTrials
        
        images(:,:,1,1) = StimObj(EncStimNames{tt});
        
        condNum = tacs_er.EncTrialPhaseConds(tt);
        % n phases go from 1 to nEncPhases, in jumps of 2 frames=72deg for 6Hz
        % stimulation
        imageSequence = zeros(PresParams.nXDivaFrames,1,'uint32');
        imageSequence(1)=3;
        imageSequence(ZeroPhaseImgFrameIDs+(condNum-1)*2)=1;
        imageSequence(ZeroPhaseBlankFrameIDs+(condNum-1)*2)=2;
        
        save([thePath.subjectPath '/xdiva/' fileName num2str(tt) '.mat'],'images','imageSequence')
    end
    
    % save instructions (changes by subj number)
    if PresParams.InstructionSet==1
        img  = rgb2gray(imread([thePath.stim '/EncInstructions_xDiva1.png']));
    else
        img  = rgb2gray(imread([thePath.stim '/EncInstructions_xDiva2.png']));
    end
    images = zeros([size(img),1,1], 'uint8');
    images(:,:,1,1) = img;        
    
    % 2 mins worth of instructions.
    imageSequence = zeros(PresParams.InstructionsLength,1,'uint32');
    imageSequence(1)=1;    
    save([thePath.subjectPath '/xdiva/Instructions_1.mat'],'images','imageSequence')
    
    % save pre-task warm up sequence    
    images = zeros([200,200,1,1], 'uint8');
    images(:,:,1,1) = 128*ones([200 200],'uint8');
    
    imageSequence = zeros(PresParams.tACSPreTaskWarmUpDur,1,'uint32');
    imageSequence(1)=1;    
    save([thePath.subjectPath '/xdiva/tACSPreTask_1.mat'],'images','imageSequence')
    
end



end

function im=fixCrossImage(imgSize,fixCrossSize)
% all in pixels
% only supports square imgages (for now).
center = round(imgSize(1)/2);
horizontalMarkIDs = [(center-fixCrossSize):(center+fixCrossSize)];
verticalMarkIDs   = horizontalMarkIDs;

im = 128*ones(imgSize(1),imgSize(1),'uint8');
im([-1:1]+center,horizontalMarkIDs) = 255;
im(verticalMarkIDs,[-1:1]+center)   = 255;

return

end

