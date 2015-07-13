function [thePath] = tACS_ER_path()

% tACS_ER -> tACS encoding & retrieval task
% return path info
% Alex Gonzalez
% May 2015
%

% set path
basepath = pwd;
cd(basepath);
addpath(genpath(pwd));

thePath.data = fullfile(pwd,'data');
thePath.stim = fullfile(pwd,'stim');
thePath.scripts = fullfile(pwd,'scripts');
thePath.main = basepath;

return


