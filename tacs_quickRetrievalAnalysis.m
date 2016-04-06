
expt = 'tacs_enc_xdiva';
mainPath = ['~/Google Drive/Research/tACS/tACS_ER_task/data/' expt '/'];
filename = 'tacs_enc_xdiva.test.mat';

RespTypes = {'old','new','unsure'};
subjs = [1:17];
nSubjs = numel(subjs);
dPFaces = zeros(nSubjs,1);
dPScns  =  zeros(nSubjs,1);

BehavTableMeans = zeros(4,3,nSubjs);
BehavTableNums =  zeros(4,3,nSubjs);

for kk = 1:nSubjs
    load([mainPath 's' num2str(subjs(kk)),'/',filename])
    
    for ii = 1:4
        for jj = 1:numel(RespTypes)
            BehavTableNums(ii,jj,kk) = sum(strcmp(ret_out.TimingInfo.CondResp(ret_out.expInfo.RetCondTrialCode==ii),RespTypes{jj}));
            BehavTableMeans(ii,jj,kk)=mean(strcmp(ret_out.TimingInfo.CondResp(ret_out.expInfo.RetCondTrialCode==ii),RespTypes{jj}));
        end
    end
    
    nFaceHits  = BehavTableNums(1,1,kk);
    nFaceMiss = BehavTableNums(1,2,kk);
    nFaceFA    = BehavTableNums(3,1,kk);
    nFaceCR    = BehavTableNums(3,2,kk);
    dPFaces(kk) = calc_dPrime(nFaceHits,nFaceMiss,nFaceFA,nFaceCR);
    
    
    nScnHits  = BehavTableNums(2,1,kk);
    nScnMiss = BehavTableNums(2,2,kk);
    nScnFA    = BehavTableNums(4,1,kk);
    nScnCR    = BehavTableNums(4,2,kk);
    dPScns(kk) = calc_dPrime(nScnHits,nScnMiss,nScnFA,nScnCR);
    
    
end
nUnsures=sum(squeeze(BehavTableNums(:,3,:)));