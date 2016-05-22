

dataPath = '~/Google Drive/Research/tACS/tACS_ER_task/data/tacs_enc_xdiva/';

nSubjs      = 32;
nTrials     = 300;
colNames    = {'TrialNum','Resp1','Resp2','NoAn','RTs'};

enc_out =[];
enc_out.colNames = colNames;
enc_out.nTrials = nTrials;
for ss = 1:nSubjs
    dataMat     = nan(nTrials,numel(colNames));
    try
        load([dataPath 'xdivafiles/s' num2str(ss) '_xdiva/Exp_MATL_PD1010_3_Cz/RTSeg_s002.mat'],'TimeLine');
        ValidTrials = [TimeLine.cndNmb]==3;
        temp=TimeLine(ValidTrials);
        dataMat(:,strcmp(colNames,'TrialNum'))   = [temp.trlNmb];
        dataMat(:,strcmp(colNames,'Resp1'))      = strcmp({temp.respString},'1');
        dataMat(:,strcmp(colNames,'Resp2'))      = strcmp({temp.respString},'2');
        dataMat(:,strcmp(colNames,'NoAn'))       = strcmp({temp.respString},'Mis');
        dataMat(:,strcmp(colNames,'RTs'))        = [temp.respTimeSec];
        dataMat(strcmp({temp.respString},'Mis'),strcmp(colNames,'RTs')) = nan;
    catch msg
        keyboard
    end
    enc_out.dataMat = dataMat;
    save([dataPath '/s' num2str(ss) '/encData_xdiva.mat'],'enc_out')
end

