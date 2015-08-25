function thePath = tACS_ER_path(subjNum,exptType)
% tACS_ER_path
% tACS encoding & retrieval task path structure
% returns path information
% subjNum   -> subject number (integer)
% exptType  -> experiment type
%   'behav'     --> behavioral only
%   'eeg'       --> eeg at encoding and retrieval (no tacs)
%   'eeg_enc'    --> eeg at encoding only (no tacs)
%   'tacs_enc'   --> tacs at encoding
%
% subject 0 reserved for debugging/testing

%------------------------------------------------------------------------%
% Author:       Alex Gonzalez (from similar lab copies)
% Created:      May 25, 2015
% LastUpdate:   Aug 20, 2015
%------------------------------------------------------------------------%

exptOptions = {'behav','eeg','eeg_enc','tacs_enc'};

if ~any(strcmp(exptOptions,exptType))
    error('Experiment type not available; please see help tACS_ER_path')
end

basepath = pwd;
cd(basepath);
addpath(basepath);

thePath = [];
thePath.subjNum = subjNum;
thePath.exptType = exptType;
thePath.main    = basepath;
thePath.stim    = fullfile(pwd,'stim');
if ~exist(thePath.stim,'dir')
    error('At the wrong directory, necessary files not found.')
end
thePath.scripts = fullfile(pwd,'scripts');
thePath.data    = fullfile(pwd,'data');

subjectPath = strcat(thePath.data,'/',exptType,'/s',num2str(subjNum));

if exist(subjectPath,'dir')
    if ~(subjNum==0)
        warning('Subject directory already present.')
    end
else
    mkdir(subjectPath);
end
thePath.subjectPath = subjectPath;

addpath(genpath(thePath.stim));
addpath(genpath(thePath.scripts));
addpath(genpath(thePath.subjectPath));

return


