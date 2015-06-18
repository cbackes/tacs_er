% script to crop and recolor face images

opt = 2;
if opt==1
    stimDir = '~/Google Drive/Research/tACS/tACS_ER_task/stim/people/uncropped/';
    newDir  = '~/Google Drive/Research/tACS/tACS_ER_task/stim/people/additional/';
    mkdir(newDir)
elseif opt==2
    stimDir = '~/Google Drive/Research/tACS/tACS_ER_task/stim/landmarks/uncropped/';
    newDir  = '~/Google Drive/Research/tACS/tACS_ER_task/stim/landmarks/additional/';
    mkdir(newDir)
end

% load directory contents and names of the files
list = dir(stimDir);
fileNames = {list.name}';
nFiles = numel(fileNames);

% resize pixel dimensions
width = 225;
height = 225;

for ff = 1:nFiles
    fileName = [stimDir fileNames{ff}];
    
    fe = exist(fileName,'file');
    if fe>0 && fe ~=7
        try
            IM = imread(fileName);
            % convert to greyscale
            if ~ismatrix(IM)
                IM = rgb2gray(IM);
            end
            IM2 = imresize(IM,[width height]);
            
            fileName = [newDir fileNames{ff}];
            imwrite(IM2,fileName,'jpeg')
        catch ME
        end
    end
end
