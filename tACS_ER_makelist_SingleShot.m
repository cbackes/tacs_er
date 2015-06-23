
function [tacs_er] = tACS_ER_makelist_SingleShot(thePath)

% make lists for tacs encoding & retrieval (tacs_er) experiment

% lists parameters
nEncTrials = 180;   % 180 trials at encoding
nEncConds  = 6;        % 2 x 3 design
nEncBlocks = 2;        % experiment has two blocks of nTrials/2

nRetTrials = nEncTrials*2;   % encoding trials + retrieval trials
nRetConds  = 4;              % old and new conditons per stim type (Face/scene)
nRetBlocks = nEncBlocks;    % blocks
maxNumConsecutiveOld = 5;   % maximum number of old trials in a row
% reset stream (to avoid duplicate lists)
s = RandStream.create('mt19937ar','seed',sum(100*clock));
if strfind(version,'R2014b')>0
    RandStream.setGlobalStream(s)
else
    RandStream.setDefaultStream(s);
end

%% load data names
cd(fullfile(thePath.stim,'landmarks'));
temp = dir('*.jpg');
count = 1;
for n = 1:size(temp,1)
    try
        imread(temp(n).name); %check if readable
        SceneNames(count) = temp(n);
        count = count + 1;
    catch
        fprintf(['\n' temp(n).name ' unreadable\n']);
    end
end

%get the names in a cell array and shuffle
SceneNames = {SceneNames.name}';
SceneNames = Shuffle(SceneNames);

cd(fullfile(thePath.stim,'people'));
temp = dir('*.jpg');
count = 1;
for n = 1:size(temp,1)
    try
        imread(temp(n).name);  %check if readable
        FaceNames(count) = temp(n);
        count = count + 1;
    catch
        fprintf(['\n' temp(n).name ' unreadable\n']);
    end
end
%get the names in a cell array and shuffle
FaceNames = {FaceNames.name}';
FaceNames = Shuffle(FaceNames);

%% Encoding
% 1st counter balance:
% get face and landmarks equal splits in each block
% codes 1 for faces, 2 for scenes
EncStimTypeIDs = {'Face','Scene'};

EncStimType = zeros(nEncTrials,1);
nEncTrialsBlock = nEncTrials/nEncBlocks;
nFacesEncBlock  = nEncTrialsBlock/2;
nScenesEncBlock = nFacesEncBlock;
for rr = 1:nEncBlocks
    TrialBlockID = ((rr-1)*nEncTrialsBlock+1):rr*nEncTrialsBlock;
    
    availableFaceTrials  = datasample(TrialBlockID,nFacesEncBlock,'replace',false);
    SceneTrials = setdiff(TrialBlockID,availableFaceTrials);
    EncStimType(availableFaceTrials)=1;
    EncStimType(SceneTrials)=2;
    % choose nTrialsBlock/2 faces
    
end
% assign actual stimuli to the list. note that since the order was shuffle
% obtaining the first N trials will result in a random order of stimuli
EncStimNames = cell(nEncTrials,1);
EncFaces     = 1:nFacesEncBlock*nEncBlocks;
EncScenes    = 1:nScenesEncBlock*nEncBlocks;
EncStimNames(EncStimType==1) = FaceNames(EncFaces);
EncStimNames(EncStimType==2) = SceneNames(EncScenes);

%% assign encoding conditions
% six conditions 2 x 3 design
% 1) faces at 0deg
% 2) faces at 90deg
% 3) faces at 180deg
% 4) scenes at 0deg
% 5) scenes at 90deg
% 6) scenes at 180deg
%
% ** Note, we are not counterbalancing these by block
nEncCondTrials      = nEncTrials/nEncConds;
EncCondTrialCode    = zeros(nEncTrials,1);
EncCondCodeIDs      = {'Face0','Face90','Face180','Scn0','Scn90','Scn180'};

% for face trials
faceConds = 1:3; 
nFaceConds=numel(faceConds);
availableFaceTrials = find(EncStimType==1);
for ii = 1:nFaceConds
    condiTrials=datasample(availableFaceTrials,nEncCondTrials,'replace',false);
    EncCondTrialCode(condiTrials) = faceConds(ii);
    availableFaceTrials = setdiff(availableFaceTrials,condiTrials);
end
% for scenes trials
sceneConds = 4:6; 
nSceneConds=numel(sceneConds);
availableSceneTrials = find(EncStimType==2);
for ii = 1:nSceneConds
    condiTrials=datasample(availableSceneTrials,nEncCondTrials,'replace',false);
    EncCondTrialCode(condiTrials) = sceneConds(ii);
    availableSceneTrials = setdiff(availableSceneTrials,condiTrials);
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
nRetTrialTypesBlock = nRetTrialsBlock/nRetConds;
nRetCondTrials   = nRetTrials/nRetConds;
% assign condition such that each is equally likely per block.
% also retry the sampling such that there are not more than maxNumConsecutiveOld
% i.e., no that many consecutive old trials.

counter = 1;
while true
    RetCondTrialCode = zeros(nRetTrials,1);
    for rr = 1:nRetBlocks
        TrialBlockID = ((rr-1)*nRetTrialsBlock+1):rr*nRetTrialsBlock;
        availableTrials = TrialBlockID;
        for cc = 1:nRetConds
            condTrials  = datasample(availableTrials,nRetTrialTypesBlock,'replace',false);
            RetCondTrialCode(condTrials)=cc;
            availableTrials  = setdiff(availableTrials,condTrials);
        end        
    end
    
    % quick way to the number of consecutive old trials
    temp=conv(ones(maxNumConsecutiveOld,1),double((RetCondTrialCode<=2)));
    if sum(temp>=maxNumConsecutiveOld)==0
        break
    end
    counter=counter+1;
    if counter >1000
        error('could not arrive at a stable sequence in 500 attempts')
    end
end
% make sure that each condition is equally represented
assert(sum(histc(RetCondTrialCode,1:nRetConds)==nRetCondTrials)==nRetConds,'unequal number of trials produced')
%% assign stimuli based on retrieval condition and block
% for each block, make sure that trials early at encoding are as likely as
% trials late at encoding.

% e.g. the first retrieval block should contain an equal number of early
% face and scene trials, and late face and scene trials. in this case early
% and late simply refer to the Encoding Block. 
% ** Code only works for two encoding conditions at the moment.
% ** encoding phase condition not counterbalanced across blocks, though it
% should havve very minimal effects.

RetStimNames = cell(nRetTrials,1);
EarlyEncBlockIDs = 1:nEncTrialsBlock;
LateEncBlockIDs  = (nEncTrialsBlock+1):nEncTrials;

% split encoding stimuli by block
EarlyFaces = Shuffle(find(EncStimType(EarlyEncBlockIDs)==1));
LateFaces = Shuffle(nEncTrialsBlock+find(EncStimType(LateEncBlockIDs)==1));
EarlyScenes = Shuffle(find(EncStimType(EarlyEncBlockIDs)==2));
LateScenes  = Shuffle(nEncTrialsBlock+find(EncStimType(LateEncBlockIDs)==2));

% Trials at retrieval split by block
FirstRetBlockTrials = 1:nRetTrials/2;
SecondRetBlockTrials = nRetTrials/2+1 : nRetTrials;

% assignment of old faces
temp1 = 1:floor(numel(EarlyFaces)/2);   
temp2 = 1:ceil(numel(LateFaces)/2);
temp3 = setdiff(1:numel(EarlyFaces),temp1); 
temp4 = setdiff(1:numel(LateFaces),temp2); 
FirstBlockIDs = Shuffle([EarlyFaces(temp1);LateFaces(temp2)]);
SecondBlockIDs = Shuffle([EarlyFaces(temp3);LateFaces(temp4)]);
RetStimNames(find(RetCondTrialCode(FirstRetBlockTrials)==1)) = EncStimNames(FirstBlockIDs);
RetStimNames(nRetTrialsBlock+find(RetCondTrialCode(SecondRetBlockTrials)==1)) = EncStimNames(SecondBlockIDs);

% assigment of old scenes
temp1 = 1:floor(numel(EarlyScenes)/2); 
temp2 = 1:ceil(numel(LateScenes)/2);
temp3 = setdiff(1:numel(EarlyScenes),temp1); 
temp4 = setdiff(1:numel(LateScenes),temp2); 
FirstBlockIDs = Shuffle([EarlyScenes(temp1);LateScenes(temp2)]);
SecondBlockIDs = Shuffle([EarlyScenes(temp3);LateScenes(temp4)]);
RetStimNames(find(RetCondTrialCode(FirstRetBlockTrials)==2)) = EncStimNames(FirstBlockIDs);
RetStimNames(nRetTrialsBlock+find(RetCondTrialCode(SecondRetBlockTrials)==2)) = EncStimNames(SecondBlockIDs);

% assign new faces
temp = FaceNames(~ismember(FaceNames,EncStimNames));
RetStimNames(RetCondTrialCode==3)=Shuffle(temp(1:nRetCondTrials));

% assign new scenes
temp = SceneNames(~ismember(SceneNames,EncStimNames));
RetStimNames(RetCondTrialCode==4)=Shuffle(temp(1:nRetCondTrials));


%% store necessary variables
tacs_er = [];

% parameterts for making stimuli list 
tacs_er.nEncTrials = nEncTrials;
tacs_er.nEncConds  = nEncConds;
tacs_er.nEncBlocks = nEncBlocks;

tacs_er.nRetTrials = nRetTrials;
tacs_er.nRetConds  = nRetConds;
tacs_er.nRetBlocks = nRetBlocks;    
tacs_er.maxNumConsecutiveOld = maxNumConsecutiveOld;

tacs_er.RandStream = s;

% encoding 
tacs_er.EncStimType = EncStimType; 
tacs_er.EncStimTypeIDs = EncStimTypeIDs;
tacs_er.EncCondTrialCode = EncCondTrialCode;
tacs_er.EncCondCodeIDs = EncCondCodeIDs;
tacs_er.EncStimNames  = EncStimNames;

% retrieval
tacs_er.RetCondIDs = RetCondIDs;
tacs_er.RetCondTrialCode = RetCondTrialCode;
tacs_er.RetStimNames = RetStimNames;

% sanity check: verify that old labels match stimuli
assert(sum(ismember(tacs_er.RetStimNames,tacs_er.EncStimNames) == ...
    (tacs_er.RetCondTrialCode<=2))==tacs_er.nRetTrials,'Labels do not match stimuli')

% encodign conditions at retrieval
[~,i]=ismember(tacs_er.RetStimNames,tacs_er.EncStimNames);
tacs_er.EncCondAtRet = zeros(tacs_er.nRetTrials,1);
tacs_er.EncCondAtRet(i>0) = tacs_er.EncCondTrialCode(i(i>0));

%%
dfwef
cd(thePath.main)

end


