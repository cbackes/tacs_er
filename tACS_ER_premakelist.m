function [] = tACS_ER_premakelist(slist)

% generate lists for tACS_ER experiment in advance
% Alex Gonzalez
% May 2015

thePath = tACS_ER_path;
addpath(thePath.scripts);

for sx = 1:length(slist)   
    tacs_er = tACS_ER_makelist(thePath,1,1); 
    save(fullfile(thePath.data,sprintf('tACS_ER.init.%s.mat',slist{sx})),'tacs_er');        
end