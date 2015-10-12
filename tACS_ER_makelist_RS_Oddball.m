
function [tacs_er] = tACS_ER_makelist_RS_Oddball(thePath)
% make lists for tacs encoding & retrieval (tacs_er) experiment
% This design:
% 1) all stimuli presented once per encoding block. Designed for 4 blocks.
% each stim will appear with an oddball at least once (throw-out trials)
% 2) responses are made to the location of the odd ball
% oddballs are counter balanced for block and location in the screen
%
% Important variables
% subjNum       -> subject  ID
% thePath       -> path for all stimulus and data directories
%
% nEncStim      -> # of UNIQUE encoding stimuli
% nEncTrials    -> # of TOTAL trials (should be multiple of nEncStim)
% nEncCond      -> # of conditions at encoding (that are critical for retrieval)
% nEncBlocks    -> # of encoding blocks
% nEncOddBalls  -> # of encoding oddballs
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
% EncOddBallTrials  -> VECTOR trial (BOOL) for trials that have an oddball
% EncOddBallLocs    -> VECTOR trial (1-4) possible locations for oddballs
% (NaNs for trials not having an oddball).

% 
% RetCondIDs        -> {'OldFace','OldScene','NewFace','NewScene'};
% RetCondTrialCode  -> VECTOR trial 1:4, indicating RetCondID
% RetStimNames      -> VECTOR trial of retrieaval IDs
% RetBlockID        -> VECTOR of trial block ID
%
% Stimuli           -> key, matrix map of stimuli's names to image matrix.

%------------------------------------------------------------------------%
% Author:       Alex Gonzalez
% Created:      Oct 8th, 2015
% LastUpdate:   Oct 8th, 2015
%------------------------------------------------------------------------%

%% Set-up

% parameter list
nUniqueFaceStimuli = 120;
nUniqueSceneStimuli = 120;
nEncStim    = nUniqueSceneStimuli+ nUniqueFaceStimuli;
nEncPhases = 6;
nEncConds  = 2*nEncPhases;     % faces/scenes x phase (6 different phases)
EncPhases  = 0:(360/nEncPhases):359;
nOddBalls  = nEncStim; % one per unique stim

switch thePath.exptType
    case 'behav_v7'
        stimSize        = [300 300];
        nEncBlocks      = 4;
        nOddBallLocs    = 4;
end

nEncTrials = nEncBlocks*nEncStim;
nRetTrials = 360;   % encoding trials + retrieval trials
nFoilTrials =nRetTrials-nEncStim;
nRetConds  = 4;     % old and new conditons per stim type (Face/scene)
nRetBlocks = nEncBlocks;    % blocks
maxNumConsecutiveOld = 8;   % maximum number of old trials in a row

% reset stream (to avoid duplicate lists)
s = RandStream.create('mt19937ar','seed',sum(100*clock));
if strfind(version,'R2014b')>0
    RandStream.setGlobalStream(s)
end

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
        ScenesMat{n} = x(:,:,1);
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
        FacesMat{n} = x(:,:,1);
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
nOddBallsBlock  = nOddBalls/nEncBlocks;
nOddBallsBlockLocs = nOddBallsBlock/nOddBallLocs; % # OB per location per block (60 in current design)
OddBallTrials   = false(nEncTrials,1);
OddBallLocs     = nan(nEncTrials,1);
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
    
    %select stims for odd balls and keep track of the ones by block
    if bb==1
        % all are available in the first block
        availableFaceIDs  = EncStimUniqueIDs(FaceTrials);
        availableSceneIDs = EncStimUniqueIDs(SceneTrials);
    end
    oddballFaceIDs = datasample(availableFaceIDs,nOddBallsBlock/2,'replace',false);
    oddballScnIDS  = datasample(availableSceneIDs,nOddBallsBlock/2,'replace',false);
    
    OddBallTrials(TrialBlockID) = ismember(EncStimUniqueIDs(TrialBlockID,2),oddballFaceIDs) | ...
        ismember(EncStimUniqueIDs(TrialBlockID,2),oddballScnIDS);
    
    availableFaceIDs    = setdiff(availableFaceIDs,oddballFaceIDs);
    availableSceneIDs   = setdiff(availableSceneIDs,oddballScnIDS);
    
    % assign locations for oddballs
    x = find(OddBallTrials(TrialBlockID));
    for ll = 1:nOddBallLocs
        y = datasample(x,nOddBallsBlockLocs,'replace',false);
        OddBallLocs(TrialBlockID(y)) = ll;
        x = setdiff(x,y);
    end
end

assert(sum(histc(OddBallLocs,1:4)==nOddBallsBlock)==nOddBallLocs,'uneven oddball!')
assert(sum(OddBallTrials)==nOddBalls,'oddball numbers do not match!!')
assert(numel(unique(EncStimNames(OddBallTrials)))==nEncStim,'not unique stimuli for each oddball!')

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
tacs_er.nEncOddBalls = nOddBalls;

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
tacs_er.EncOddBallTrials = OddBallTrials;
tacs_er.EncOddBallLocs   = OddBallLocs;

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

% Save into subjects path
fileName = 'tacs_er.task.mat';
if ~exist([thePath.subjectPath '/' fileName],'file')
    save([thePath.subjectPath '/' fileName],'tacs_er')
else
    if ~(thePath.subjNum==0)
        warning('tacs task for this subject already created')
        while 1
            s = input('Overwrite? (Y/N)','s');
            if strcmp(s,'Y')
                save([thePath.subjectPath '/' fileName],'tacs_er')
                break
            elseif strcmp(s,'N')
                break
            end
        end
    end
end

end


