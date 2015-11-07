function thePath = tACS_ER_path(subjNum,exptType)
% tACS_ER_path
% tACS encoding & retrieval task path structure
% returns path information
% subjNum   -> subject number (integer)
% exptType  -> experiment type
%   'behav'     --> behavioral original
%   'behav_v3'  --> non famous stimuli, cue response
%   'behav_v4'  --> passive viewing.
%   'behav_v5'  --> 4 presentations
%   'behav_v6'  --> v4 with @ 5Hz
%   'behav_v7'  --> 4 presentations, left/right oddball in one of
%   presentations
%   'behav_v8'  --> 4 presentations with face/scene decision on stimuli
%   'behav_v9'  --> multiple presentations whithin a trial, cue color oddball response
%   'behav_v10' --> v9, with semantic decision
%   'eeg'       --> eeg at encoding and retrieval (no tacs)
%   'eeg_enc'   --> eeg at encoding only (no tacs)
%   'tacs_enc'  --> tacs at encoding
%
% subject 0 reserved for debugging/testing

%------------------------------------------------------------------------%
% Author:       Alex Gonzalez (from similar lab copies)
% Created:      May 25, 2015
% LastUpdate:   Nov 6s, 2015
%------------------------------------------------------------------------%

exptOptions = {'behav','behav_v3','behav_v4','behav_v5','behav_v6',...
    'behav_v7','behav_v8','behav_v9','behav_v10','eeg','eeg_enc','tacs_enc'};

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


