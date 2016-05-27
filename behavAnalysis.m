function behav_out = behavAnalysis(expt,nSubjs)
% behavioral analysis function for tACS_ER task
% analyzes both encoding and retrieval data for tacs_enc
%
% thePath input determines the experimental folder to load data from
%
%------------------------------------------------------------------------%
% Author:         Alex Gonzalez
% Created:        Dec 1, 2015
% LastUpdate:     April 21, 2016
%------------------------------------------------------------------------%
%
%
% encoding task:
% Perceptual Discrimination
% 1) Face hit rate
% 2) Scene hit rate
% 3) Discrimination rate (accuracy) (stats across subjects)
% 4) RTs of for Faces and Scenes; (stats across subjects)
%
% retrieval task:
% 1) hit rate, fa rate, cr rate, miss rate
% 2) hit rate, fa rate, cr rate, miss rate for Faces and Scenes
% independently
% 3) RTs
% 4) Confidence
%

dataPath = ['~/Google Drive/Research/tACS/tACS_ER_task/data/' expt '/'];
behav_out     = [];
behav_out.encSubj = cell(nSubjs,1);
behav_out.retSubj = cell(nSubjs,1);

% encoding stats pre-allocation
behav_out.encSummary = [];
behav_out.encSummary.goodSubj           = nan(nSubjs,1);
behav_out.encSummary.meanAcc            = nan(nSubjs,1);
behav_out.encSummary.FaceHR             = nan(nSubjs,1);
behav_out.encSummary.SceneHR            = nan(nSubjs,1);
behav_out.encSummary.medianRTs          = nan(nSubjs,1);
behav_out.encSummary.medianFaceRTs      = nan(nSubjs,1);
behav_out.encSummary.medianSceneRTs     = nan(nSubjs,1);
behav_out.encSummary.FaceScene_RTs_TVals = nan(nSubjs,1);
        
        
% retrieval stats pre-allocation
behav_out.retSummary =[];
behav_out.retSummary.dPrime             = nan(nSubjs,1);
behav_out.retSummary.dPrime_C           = nan(nSubjs,1);
behav_out.retSummary.Face_dPrime        = nan(nSubjs,1);
behav_out.retSummary.Face_dPrime_C      = nan(nSubjs,1);
behav_out.retSummary.Scene_dPrime       = nan(nSubjs,1);
behav_out.retSummary.Scene_dPrime_C     = nan(nSubjs,1);

% reaction time summaries
behav_out.retSummary.medianHit_RTs     = nan(nSubjs,1);
behav_out.retSummary.medianCRs_RTs     = nan(nSubjs,1);
behav_out.retSummary.medianFA_RTs      = nan(nSubjs,1);
behav_out.retSummary.medianMiss_RTs    = nan(nSubjs,1);
behav_out.retSummary.HitCRs_RTs_TVals   = nan(nSubjs,1);


% confidence summaries
behav_out.retSummary.dPrimeConf             = nan(nSubjs,3);
behav_out.retSummary.dPrimeConf_C           = nan(nSubjs,3);
behav_out.retSummary.meanAccuracyByConf     = nan(nSubjs,3);
behav_out.retSummary.nRespByConf            = nan(nSubjs,3);
behav_out.retSummary.nH_nMiss_nFA_nCRs      = nan(nSubjs,3,4);
behav_out.retSummary.medianHit_RTsConf     = nan(nSubjs,3);
behav_out.retSummary.medianCRs_RTsConf         = nan(nSubjs,3);
behav_out.retSummary.medianFA_RTsConf          = nan(nSubjs,3);
behav_out.retSummary.medianMiss_RTsConf        = nan(nSubjs,3);

% trials
behav_out.EncRet = [];
behav_out.EncRet.EncHitTrialIDs         = cell(nSubjs,1);
behav_out.EncRet.EncMissTrialIDs        = cell(nSubjs,1);
behav_out.EncRet.EncHitRTs              = cell(nSubjs,1);
behav_out.EncRet.EncMissRTs             = cell(nSubjs,1);

behav_out.Trials =[]; nEncTrials = 300; nRetTrials = 450;
behav_out.Trials.EncTrialPhaseCond      = nan(nSubjs,nEncTrials);
behav_out.Trials.EncTrialPhase          = nan(nSubjs,nEncTrials);
behav_out.Trials.RetStimIDAtEnc         = nan(nSubjs,nEncTrials);
behav_out.Trials.StimTypeEnc            = nan(nSubjs,nEncTrials); %
behav_out.Trials.EncCondAtRet           = nan(nSubjs,nRetTrials);
behav_out.Trials.EncStimIDAtRet         = nan(nSubjs,nRetTrials);


for ss = 1:nSubjs
    try
        if strcmp(expt,'tacs_enc_xdiva')        
            load([dataPath 's' num2str(ss) '/encData_xdiva.mat'])
            load([dataPath 's' num2str(ss) '/tacs_er_xdiva.task.mat'])
            behav_out.encSubj{ss}               = EncAnalysisBySubj_xdiva(enc_out,tacs_er);       
        else
            load([dataPath 's' num2str(ss) '/tacs_er.encoding.mat'])
             behav_out.encSubj{ss}              = EncAnalysisBySubj(enc_out);       
        end
        behav_out.encSummary.goodSubj(ss)           = behav_out.encSubj{ss}.goodSubj;
        behav_out.encSummary.meanAcc(ss)            = behav_out.encSubj{ss}.meanAccuracy;
        behav_out.encSummary.FaceHR(ss)             = behav_out.encSubj{ss}.FaceHitRate;
        behav_out.encSummary.SceneHR(ss)            = behav_out.encSubj{ss}.SceneHitRate;
        behav_out.encSummary.medianRTs(ss)          = behav_out.encSubj{ss}.RTsStats.medianRT;
        behav_out.encSummary.medianFaceRTs(ss)      = behav_out.encSubj{ss}.FaceRTsStats.medianRT;
        behav_out.encSummary.medianSceneRTs(ss)     = behav_out.encSubj{ss}.SceneRTsStats.medianRT;
        behav_out.encSummary.FaceScene_RTs_TVals(ss) ...
            = behav_out.encSubj{ss}.RTsStats_FvsSc.tstat;
        behav_out.encSummary.HRByPhase(ss,:,:)          = behav_out.encSubj{ss}.HRbyPhase;        
        behav_out.encSummary.medianRTsByPhase(ss,:,:)   = behav_out.encSubj{ss}.medianRTsByPhase;
        
    catch msg   
        fprintf('issues processing subject %i \n',ss)
    end
    try
        
        if strcmp(expt,'tacs_enc_xdiva')
            load([dataPath 's' num2str(ss) '/tacs_enc_xdiva.test.mat'])
            behav_out.Trials.EncTrialPhaseCond(ss,:)        = ret_out.expInfo.EncTrialPhaseConds;
            behav_out.Trials.EncTrialPhase(ss,:)            = ret_out.expInfo.EncTrialPhase;
            behav_out.Trials.EncCondAtRet(ss,:)             = ret_out.expInfo.EncCondAtRet;
            
        elseif strcmp(expt,'tacs_enc')
            load([dataPath 's' num2str(ss) '/tacs_er.test.mat'])
            load([dataPath 's' num2str(ss) '/EventsPhase.mat']);
            behav_out.Trials.EncTrialPhase(ss,:) = out.TrueAngleStims;
        end
        behav_out.retSubj{ss}                           = RetAnalysisBySubj(ret_out);
        
        % IDs of test items during encoding.
        behav_out.Trials.EncStimIDAtRet(ss,:)           = behav_out.retSubj{ss}.EncStimIDAtRet;
        [s,i]                                           = sort(behav_out.Trials.EncStimIDAtRet(ss,:));
        behav_out.Trials.RetStimIDAtEnc(ss,:)           = i(s>0);
        behav_out.Trials.StimTypeEnc(ss,:)              = ret_out.expInfo.EncStimType;
        
        % DPrimes and accuracy summary measures
        behav_out.retSummary.dPrime(ss)                 = behav_out.retSubj{ss}.dPrime;
        behav_out.retSummary.dPrime_C(ss)               = behav_out.retSubj{ss}.dPrime_C;
        behav_out.retSummary.Face_dPrime(ss)            = behav_out.retSubj{ss}.Face_dPrime;
        behav_out.retSummary.Face_dPrime_C(ss)          = behav_out.retSubj{ss}.Face_dPrime_C;
        behav_out.retSummary.Scene_dPrime(ss)           = behav_out.retSubj{ss}.Scene_dPrime;
        behav_out.retSummary.Scene_dPrime_C(ss)         = behav_out.retSubj{ss}.Scene_dPrime_C;
        behav_out.retSummary.meanAccuracyByConf(ss,:)   = behav_out.retSubj{ss}.meanAccuracyByConf;
        behav_out.retSummary.dPrimeConf(ss,:)           = behav_out.retSubj{ss}.dPrimeConf;
        behav_out.retSummary.dPrimeConf_C(ss,:)         = behav_out.retSubj{ss}.dPrimeConf_C;
        
        behav_out.retSummary.Face_dPrimeConf(ss,:)      = behav_out.retSubj{ss}.Face_dPrimeConf;
        behav_out.retSummary.Face_dPrimeConf_C(ss,:)    = behav_out.retSubj{ss}.Face_dPrimeConf_C;
        behav_out.retSummary.Scene_dPrimeConf(ss,:)     = behav_out.retSubj{ss}.Scene_dPrimeConf;
        behav_out.retSummary.Scene_dPrimeConf_C(ss,:)   = behav_out.retSubj{ss}.Scene_dPrimeConf_C;
        
        % numbers of trials per conditions
        behav_out.retSummary.nH_nMiss_nFA_nCRs(ss,:,:)  = behav_out.retSubj{ss}.nH_nMiss_nFA_nCRs;
        behav_out.retSummary.nRespByConf(ss,:)          = behav_out.retSubj{ss}.nRespByConf;
        
        % reaction times        
        behav_out.retSummary.medianHit_RTs(ss)          = behav_out.retSubj{ss}.Hit_RTsStats.medianRT;
        behav_out.retSummary.medianCRs_RTs(ss)          = behav_out.retSubj{ss}.CRs_RTsStats.medianRT;
        behav_out.retSummary.medianFA_RTs(ss)           = behav_out.retSubj{ss}.FA_RTsStats.medianRT;
        behav_out.retSummary.medianMiss_RTs(ss)         = behav_out.retSubj{ss}.Misses_RTsStats.medianRT;
        behav_out.retSummary.HitCRs_RTs_TVals(ss)       = behav_out.retSubj{ss}.HitCRs_RTs.tstat;
        
        behav_out.retSummary.medianFaceScnHit_RTs(ss,:) = [behav_out.retSubj{ss}.Face_Hit_RTsStats.medianRT, behav_out.retSubj{ss}.Scene_Hit_RTsStats.medianRT];
        behav_out.retSummary.medianFaceScnCRs_RTs(ss,:) = [behav_out.retSubj{ss}.Face_CRs_RTsStats.medianRT, behav_out.retSubj{ss}.Scene_CRs_RTsStats.medianRT];
        behav_out.retSummary.medianFaceScnFA_RTs(ss,:)  = [behav_out.retSubj{ss}.Face_FA_RTsStats.medianRT,  behav_out.retSubj{ss}.Scene_FA_RTsStats.medianRT];
        behav_out.retSummary.medianFaceScnMiss_RTs(ss,:)= [behav_out.retSubj{ss}.Face_Misses_RTsStats.medianRT,behav_out.retSubj{ss}.Scene_Misses_RTsStats.medianRT];
        
        behav_out.retSummary.medianHit_RTsConf(ss,:)    = [behav_out.retSubj{ss}.Hit_RTsStatsConf.medianRT];
        behav_out.retSummary.medianCRs_RTsConf(ss,:)    = [behav_out.retSubj{ss}.CRs_RTsStatsConf.medianRT];
        behav_out.retSummary.medianFA_RTsConf(ss,:)     = [behav_out.retSubj{ss}.FA_RTsStatsConf.medianRT];
        behav_out.retSummary.medianMiss_RTsConf(ss,:)   = [behav_out.retSubj{ss}.Misses_RTsStatsConf.medianRT];
        
        % trials
        behav_out.EncRet.EncHitTrialIDs{ss}         = behav_out.retSubj{ss}.EncHitTrialIDs;
        behav_out.EncRet.EncMissTrialIDs{ss}        = behav_out.retSubj{ss}.EncMissTrialIDs;
        behav_out.EncRet.EncHCTrials{ss}            = behav_out.retSubj{ss}.EncHCTrials;
        behav_out.EncRet.EncMCTrials{ss}            = behav_out.retSubj{ss}.EncMCTrials;
        behav_out.EncRet.EncLCTrials{ss}            = behav_out.retSubj{ss}.EncLCTrials;
        behav_out.EncRet.EncHitRTs{ss}              = behav_out.retSubj{ss}.EncHitRTs;
        behav_out.EncRet.EncMissRTs{ss}             = behav_out.retSubj{ss}.EncMissRTs;
        
    catch msg
        disp(msg)
    end
    
end
save([dataPath '/Summary/BehavSummary'],'behav_out')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliary functions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Encoding Analyses by subject
function enc=EncAnalysisBySubj(enc_out)
enc =[];
enc.encCond = enc_out.exptInfo.EncStimType; % face/scene trials
enc.RTs     = enc_out.TimingInfo.trialRT;
keyPress = enc_out.TimingInfo.trialKeyPress;

cond1Key=enc_out.PresParams.RespToCond1;
cond2Key=enc_out.PresParams.RespToCond2;

% get simple responses information
binOut = tabulateBinaryResponses(keyPress,enc.encCond,cond1Key,cond2Key);

enc.CorrectResp     = binOut.Correct;
enc.InCorrectResp   = binOut.InCorrect;
enc.meanAccuracy    = binOut.mACC;
enc.FaceHitRate     = binOut.Cond1HR;
enc.SceneHitRate    = binOut.Cond2HR;
enc.FaceCorrectResp = binOut.Cond1Correct;
enc.SceneCorrectResp= binOut.Cond2Correct;
enc.otherEncInfo    = binOut;

% RT info (for correct)
enc.RTsStats        = analyzeRTs(enc.RTs(enc.CorrectResp));
enc.FaceRTsStats    = analyzeRTs(enc.RTs(enc.FaceCorrectResp));
enc.SceneRTsStats   = analyzeRTs(enc.RTs(enc.SceneCorrectResp));
[~,a1,~,a2] = ttest2(enc.RTs(enc.FaceCorrectResp),enc.RTs(enc.SceneCorrectResp));
enc.RTsStats_FvsSc.p=a1;
enc.RTsStats_FvsSc.tstat=a2.tstat;


end

function enc = EncAnalysisBySubj_xdiva(enc_out,tacs_er)
enc =[];
enc.goodSubj = mean(isnan(enc_out.dataMat(:,5)))<0.5; % liberal threshold for responses.
enc.encCond  = tacs_er.EncStimType; % face/scene trials


% Get Response IDs based on Subjects Responses. Categorization task should be
% close to perfect performance.
enc.StimIDResp1 = mode(enc.encCond(enc_out.dataMat(:,strcmp(enc_out.colNames,'Resp1'))==1));
enc.StimIDResp2 = mode(enc.encCond(enc_out.dataMat(:,strcmp(enc_out.colNames,'Resp2'))==1));

enc.encResp  = enc_out.dataMat(:,strcmp(enc_out.colNames,'Resp1'))*enc.StimIDResp1 + ...
                enc_out.dataMat(:,strcmp(enc_out.colNames,'Resp2'))*enc.StimIDResp2;


if ~(isnan(enc.StimIDResp1) || isnan(enc.StimIDResp2))
    enc.encRespTab = array2table(crosstab([nan ;enc.encCond],[0;enc.encResp]),'VariableNames', ...
    {'NoResp','FaceResp','SceneResp'},'rowNames',{'FaceStim','SceneStim'});
else
     enc.encRespTab = array2table(nan(2,3),'VariableNames', ...
    {'NoResp','FaceResp','SceneResp'},'rowNames',{'FaceStim','SceneStim'});
end
    
deltaTime = (tacs_er.PresParams.VideoFrameDurSecs*2)*(tacs_er.EncTrialPhaseConds-1);
enc.RTs     = enc_out.dataMat(:,strcmp(enc_out.colNames,'RTs'))-deltaTime;

% get simple responses information
binOut = tabulateBinaryResponses2(enc.encResp,enc.encCond);

enc.CorrectResp     = binOut.Correct;
enc.InCorrectResp   = binOut.InCorrect;
enc.meanAccuracyAll = binOut.mACC;
enc.meanAccuracy    = mean(binOut.Correct(~binOut.noResp));
enc.FaceHitRate     = enc.encRespTab.FaceResp(1)/sum(enc.encRespTab.FaceResp);
enc.SceneHitRate    = enc.encRespTab.SceneResp(2)/sum(enc.encRespTab.SceneResp);
enc.FaceCorrectResp = binOut.Cond1Correct;
enc.SceneCorrectResp= binOut.Cond2Correct;
enc.otherEncInfo    = binOut;

enc.HRbyPhase           = zeros(2,tacs_er.nEncPhases); %  faces, scenes

for ii = 1:tacs_er.nEncPhases    
    trials = tacs_er.EncTrialPhaseConds==ii & ~(enc.encResp==0);        
    temp = tabulateBinaryResponses2(enc.encResp(trials),enc.encCond(trials)); 
    enc.HRbyPhase(1,ii) = temp.mACC;
    enc.HRbyPhase(2,ii) = temp.Cond1HR;
    enc.HRbyPhase(3,ii) = temp.Cond2HR;
end
% RT info (for correct)
enc.RTsStats        = analyzeRTs(enc.RTs(enc.CorrectResp));
enc.FaceRTsStats    = analyzeRTs(enc.RTs(enc.FaceCorrectResp));
enc.SceneRTsStats   = analyzeRTs(enc.RTs(enc.SceneCorrectResp));
[~,a1,~,a2] = ttest2(enc.RTs(enc.FaceCorrectResp),enc.RTs(enc.SceneCorrectResp));
enc.RTsStats_FvsSc.p=a1;
enc.RTsStats_FvsSc.tstat=a2.tstat;

enc.medianRTsByPhase    = zeros(3,tacs_er.nEncPhases); %  faces, scenes
for ii = 1:tacs_er.nEncPhases        
    trials = tacs_er.EncTrialPhaseConds==ii & ~(enc.encResp==0);  
    enc.medianRTsByPhase(1,ii) = nanmedian(enc.RTs(trials));
    for jj = 1:2
        trials = tacs_er.EncTrialPhaseConds==ii & ~(enc.encResp==0) & (enc.encCond==jj);  
        enc.medianRTsByPhase(jj+1,ii) = nanmedian(enc.RTs(trials));
    end
end

%enc.PhaseRTAnova            = [];
% [~,a,b]      = anova1(enc.RTs(~(enc.encResp==0)), ...
%     tacs_er.EncTrialPhaseConds(~(enc.encResp==0)));
% enc.PhaseRTAnova.anovaTable = a;
% enc.PhaseRTAnova.groupIno   = b;
% enc.PhaseRTAnova.FStat  = a{2,5};
% enc.PhaseRTAnova.Pval   = a{2,6};


end

% Retrieval Analyses
function ret=RetAnalysisBySubj(ret_out)

ret = [];
ret.retCond     = ret_out.expInfo.RetCondTrialCode;
% 1 old face, 2 old scene, 3 new face, 4 new scene
ret.OldTrials   = ret.retCond==1 | ret.retCond==2;
ret.NewTrials   = ret.retCond==3 | ret.retCond==4;
ret.FaceTrials  = ret.retCond==1 | ret.retCond==3;
ret.SceneTrials = ret.retCond==2 | ret.retCond==4;
[~,ret.EncStimIDAtRet]=ismember(ret_out.expInfo.RetStimNames,ret_out.expInfo.EncStimNames);


ret.RTs         = ret_out.TimingInfo.trialRT;

keyPress        = ret_out.TimingInfo.trialKeyPress;

cond1Key=ret_out.PresParams.RespButtons(1);
cond2Key=ret_out.PresParams.RespButtons(3);

% get overall response information
binOut = tabulateBinaryResponses(keyPress,1*ret.OldTrials+2*ret.NewTrials,cond1Key,cond2Key);

ret.CorrectResp     = binOut.Correct;
ret.InCorrectResp   = binOut.InCorrect;
ret.meanAccuracy    = binOut.mACC;

ret.HitRate         = binOut.Cond1HR;
ret.CorrRejRate     = binOut.Cond2HR;
ret.FARate          = binOut.Cond2FAR;
ret.MissRate        = binOut.Cond1FAR;

ret.Hits            = binOut.Cond1Correct;
ret.CRs             = binOut.Cond2Correct;
ret.Misses          = binOut.Cond1InCorrect;
ret.FA              = binOut.Cond2InCorrect;

[ret.dPrime, ret.dPrime_C] = calc_dPrime(sum(ret.Hits),sum(ret.Misses),sum(ret.FA),sum(ret.CRs));

% RT analyses
ret.Hit_RTsStats     = analyzeRTs(ret.RTs(ret.Hits));
ret.CRs_RTsStats     = analyzeRTs(ret.RTs(ret.CRs));
ret.Misses_RTsStats  = analyzeRTs(ret.RTs(ret.Misses));
ret.FA_RTsStats      = analyzeRTs(ret.RTs(ret.FA));

[~,a1,~,a2]          = ttest2(ret.RTs(ret.Hits),ret.RTs(ret.CRs));
ret.HitCRs_RTs.p     = a1;
ret.HitCRs_RTs.tstat =a2.tstat;

% separated by faces & scenes
types = {'Face','Scene'};
for tt = 1:numel(types)
    type = types{tt};
    trials = ret.([type 'Trials']);
    
    binOut = tabulateBinaryResponses(keyPress(trials),1*ret.OldTrials(trials)+2*ret.NewTrials(trials),cond1Key,cond2Key);
    
    ret.([type '_meanACC'])      = binOut.mACC;
    ret.([type '_HitRate'])      = binOut.Cond1HR;
    ret.([type '_CorrRejRate'])  = binOut.Cond2HR;
    ret.([type '_FARate'])       = binOut.Cond2FAR;
    ret.([type '_MissRate'])     = binOut.Cond1FAR;
    
    [ret.([type '_dPrime']), ret.([type '_dPrime_C'])]= calc_dPrime(sum(binOut.Cond1Correct),....
        sum(binOut.Cond1InCorrect),sum(binOut.Cond2InCorrect),sum(binOut.Cond2Correct));
    
    % RT analyses
    ret.([type '_Hit_RTsStats'])     = analyzeRTs(ret.RTs(ret.Hits & trials));
    ret.([type '_CRs_RTsStats'])     = analyzeRTs(ret.RTs(ret.CRs & trials));
    ret.([type '_Misses_RTsStats'])  = analyzeRTs(ret.RTs(ret.Misses & trials));
    ret.([type '_FA_RTsStats'])      = analyzeRTs(ret.RTs(ret.FA & trials));
end

% Confidence
ret.ConfidenceResp = ret_out.TimingInfo.ConfResp;
ret.Confidence     = 3*strcmp(ret.ConfidenceResp,'high')+ ...
    2*strcmp(ret.ConfidenceResp,'mid')+1*strcmp(ret.ConfidenceResp,'low');
ret.nRespByConf    = histc(ret.Confidence,1:3);

% Tabulate Hits, FA, CRs and Misses as a function of confidence
for cc = 1:3
    % get overall response information
    trials = ret.Confidence==cc;
    
    binOut = tabulateBinaryResponses(keyPress(trials),...
        1*ret.OldTrials(trials)+2*ret.NewTrials(trials),...
        cond1Key,cond2Key);
    
    ret.meanAccuracyByConf(cc)    = binOut.mACC;
    
    ret.HitRateByConf(cc)         = binOut.Cond1HR;
    ret.CorrRejRateByConf(cc)     = binOut.Cond2HR;
    ret.FARateByConf(cc)          = binOut.Cond2FAR;
    ret.MissRateByConf(cc)        = binOut.Cond1FAR;
    
    x1           = binOut.Cond1Correct; % hits
    x2           = binOut.Cond1InCorrect; % misses
    x3           = binOut.Cond2InCorrect; % # FA
    x4           = binOut.Cond2Correct; % # CRs

    % #s of hits/misses/fa/CRs
    ret.nH_nMiss_nFA_nCRs(cc,:)   = sum([x1 x2 x3 x4]);
    
    % compute dprime
    ret.dPrimeConf(cc)      = binOut.dPrime;
    ret.dPrimeConf_C(cc)    = binOut.dPrimeC;   
    
    % get RT stats per confidence.
    rts = ret.RTs(trials);
    ret.Hit_RTsStatsConf(cc)     = analyzeRTs(rts(x1));
    ret.CRs_RTsStatsConf(cc)     = analyzeRTs(rts(x2));
    ret.Misses_RTsStatsConf(cc)  = analyzeRTs(rts(x3));
    ret.FA_RTsStatsConf(cc)      = analyzeRTs(rts(x4));
    
    % compute dPrime for faces and scenes
    types = {'Face','Scene'};
    for tt = 1:numel(types)
        type = types{tt};
        trials2 = ret.([type 'Trials']) & trials;
        binOut = tabulateBinaryResponses(keyPress(trials2),...
        1*ret.OldTrials(trials2)+2*ret.NewTrials(trials2),...
        cond1Key,cond2Key);
        
        ret.([type '_dPrimeConf'])(cc)   = binOut.dPrime;
        ret.([type '_dPrimeConf_C'])(cc) = binOut.dPrimeC;
    end
end

ret.EncHitTrialIDs     = ret.EncStimIDAtRet(ret.Hits);
ret.EncHCHitTrialIDs   = ret.EncStimIDAtRet(ret.Hits & ret.Confidence==3);
ret.EncMCHitTrialIDs   = ret.EncStimIDAtRet(ret.Hits & ret.Confidence==2);
ret.EncLCHitTrialIDs   = ret.EncStimIDAtRet(ret.Hits & ret.Confidence==1);
ret.EncMissTrialIDs    = ret.EncStimIDAtRet(ret.Misses );
ret.EncHCMissTrialIDs  = ret.EncStimIDAtRet(ret.Misses & ret.Confidence==3);
ret.EncMCMissTrialIDs  = ret.EncStimIDAtRet(ret.Misses & ret.Confidence==2);
ret.EncLCMissTrialIDs  = ret.EncStimIDAtRet(ret.Misses & ret.Confidence==1);
ret.EncHCTrials        = ret.EncStimIDAtRet(ret.OldTrials & ret.Confidence==3);
ret.EncMCTrials        = ret.EncStimIDAtRet(ret.OldTrials & ret.Confidence==2);
ret.EncLCTrials        = ret.EncStimIDAtRet(ret.OldTrials & ret.Confidence==1);
ret.EncHitRTs          = ret.RTs(ret.Hits);
ret.EncMissRTs         = ret.RTs(ret.Misses);


end

% Other functions.
function out=tabulateBinaryResponses(responses,correctAns,key1,key2)

out = [];
fields = {'multipleResp','noResp','Cond1Correct','Cond2Correct',...
    'Cond1InCorrect','Cond2InCorrect','Cond1HR','Cond2HR','Cond1FAR',...
    'Cond2FAR','Correct','InCorrect','mACC','dPrime','dPrimeC'};
for ii = 1:numel(fields)
    out.(fields{ii}) = [];
end

% handle special cases
% for more than one answer, use the first button response
if any(cellfun(@iscell,responses))
    ids = find(cellfun(@iscell,responses));
    for ii = 1:numel(ids)
        responses(ids(ii))=responses{ids(ii)}(1);
    end
    out.multipleResp = ids;
end

% no response trials
if any(cellfun(@isempty,responses))
    out.noResp=find(cellfun(@isempty,responses));
end

try
    % tabulate correct/incorrect by cond
    out.Cond1Correct   = strcmp(responses,key1)&(correctAns==1); % hits
    out.Cond1InCorrect = strcmp(responses,key2)&(correctAns==1); % misses
    out.Cond2Correct   =  strcmp(responses,key2)&(correctAns==2); % FA
    out.Cond2InCorrect = strcmp(responses,key1)&(correctAns==2); % CRs
    
    out.Cond1HR  = sum(out.Cond1Correct)/sum(correctAns==1);  % correct response to cond 1
    out.Cond2HR  = sum(out.Cond2Correct)/sum(correctAns==2);   % correct response to cond 2
    out.Cond1FAR = sum(out.Cond1InCorrect)/(sum(out.Cond1InCorrect)+sum(out.Cond1Correct)); % responded as condition 1 (when it was condition 2)
    out.Cond2FAR = sum(out.Cond2InCorrect)/(sum(out.Cond2InCorrect)+sum(out.Cond2Correct)); % responded as condition 2 (when it was condition 1)
    
    out.Correct   = out.Cond1Correct | out.Cond2Correct;
    out.InCorrect = out.Cond1InCorrect | out.Cond2InCorrect;
    out.mACC = mean(out.Correct);
    
    x1           = out.Cond1Correct; % hits
    x1s          = sum(x1);
    x2           = out.Cond1InCorrect; % misses
    x2s          = sum(x2);
    x3           = out.Cond2InCorrect; % # FA
    x3s          = sum(x3); % # FA
    x4           = out.Cond2Correct; % # CRs
    x4s          = sum(x4);
    
    % compute dprime
    [out.dPrime, out.dPrimeC] = calc_dPrime(x1s,x2s,x3s,x4s);
    
catch msg
    keyboard
end

end

% version when resp and ans match.
function out=tabulateBinaryResponses2(resp,stim)

out = [];
fields = {'multipleResp','noResp','Cond1Correct','Cond2Correct',...
    'Cond1InCorrect','Cond2InCorrect','Cond1HR','Cond2HR','Cond1FAR',...
    'Cond2FAR','Correct','InCorrect','mACC','dPrime','dPrimeC'};
for ii = 1:numel(fields)
    out.(fields{ii}) = [];
end

% no response trials
out.noResp = (resp==0);

try
    % tabulate correct/incorrect by cond
    out.Cond1Correct   = (resp==1)&(stim==1); % hits
    out.Cond1InCorrect = (resp==2)&(stim==1); % misses
    out.Cond2Correct   = (resp==2)&(stim==2); % FA
    out.Cond2InCorrect = (resp==1)&(stim==2); % CRs
    
    out.Cond1HR  = sum(out.Cond1Correct)/sum(stim==1);  % correct response to cond 1
    out.Cond2HR  = sum(out.Cond2Correct)/sum(stim==2);   % correct response to cond 2
    out.Cond1FAR = sum(out.Cond1InCorrect)/(sum(out.Cond1InCorrect)+sum(out.Cond1Correct)); % responded as condition 1 (when it was condition 2)
    out.Cond2FAR = sum(out.Cond2InCorrect)/(sum(out.Cond2InCorrect)+sum(out.Cond2Correct)); % responded as condition 2 (when it was condition 1)
    
    out.Correct   = out.Cond1Correct | out.Cond2Correct;
    out.InCorrect = out.Cond1InCorrect | out.Cond2InCorrect;
    out.mACC = mean(out.Correct);
    
    x1           = out.Cond1Correct; % hits
    x1s          = sum(x1);
    x2           = out.Cond1InCorrect; % misses
    x2s          = sum(x2);
    x3           = out.Cond2InCorrect; % # FA
    x3s          = sum(x3); % # FA
    x4           = out.Cond2Correct; % # CRs
    x4s          = sum(x4);
    
    % compute dprime
    [out.dPrime, out.dPrimeC] = calc_dPrime(x1s,x2s,x3s,x4s);
    
catch msg
    keyboard
end

end
function out=analyzeRTs(RTs)
out.meanRT = nanmean(RTs);
out.sdRT = nanstd(RTs);
out.medianRT = nanmedian(RTs);
out.kurtosisRT = kurtosis(RTs);
end

