
function out = subjFileInfo_xdiva(subj)

out = [];
out.dataPath                = ['~/Google Drive/Research/tACS/tACS_ER_task/data/tacs_enc_xdiva/s' num2str(subj) '/' ];
out.infoStimFileName        = [];
out.stimulationFileName     = [];
out.eegEncodingFileName     = [];
out.eegRetrievalInfoFileName= [];
out.eegRetrievalFileName    = [];

switch  subj
    case 1        
        out.infoStimFileName        = '20160215095702_s1.info';
        out.stimulationFileName     = '20160215095701_s1.stim';
        out.eegEncodingFileName     = '20160215095702_s1.easy';
        
        out.eegRetrievalInfoFileName= '';
        out.eegRetrievalFileName    = '';
    case 2        
        out.infoStimFileName        = '20160220103208_s2.info';
        out.stimulationFileName     = '20160220103207_s2.stim';
        out.eegEncodingFileName     = '20160220103208_s2.easy';

        out.eegRetrievalInfoFileName= '20160220110256_s2.info';
        out.eegRetrievalFileName    = '20160220110256_s2.easy';
    case 3        
        out.infoStimFileName        = '20160220121937_s3.info';
        out.stimulationFileName     = '20160220121936_s3.stim';
        out.eegEncodingFileName     = '20160220121937_s3.easy';

        out.eegRetrievalInfoFileName= '20160220123951_s3.info';
        out.eegRetrievalFileName    = '20160220123951_s3.easy';
    case 4
        out.infoStimFileName        = '20160222092623_s4.info';
        out.stimulationFileName     = '20160222092623_s4.stim';
        out.eegEncodingFileName     = '20160222092623_s4.easy';
        
        out.eegRetrievalInfoFileName= '20160222094629_s4.info';
        out.eegRetrievalFileName    = '20160222094629_s4.easy';
    case 5
        out.infoStimFileName        = '20160225092800_s5.info';
        out.stimulationFileName     = '20160225092759_s5.stim';
        out.eegEncodingFileName     = '20160225092800_s5.easy';
        
        out.eegRetrievalInfoFileName= '20160225094731_s5.info';
        out.eegRetrievalFileName    = '20160225094731_s5.easy';
    case 6
        out.infoStimFileName        = '20160225125237_s6.info';
        out.stimulationFileName     = '20160225125236_s6.stim';
        out.eegEncodingFileName     = '20160225125237_s6.easy';
        
        out.eegRetrievalInfoFileName= '20160225131151_s6.info';
        out.eegRetrievalFileName    = '20160225131151_s6.easy';
    case 7
        out.infoStimFileName        = '20160226110221_s7.info';
        out.stimulationFileName     = '20160226110220_s7.stim';
        out.eegEncodingFileName     = '20160226110221_s7.easy';
        
        out.eegRetrievalInfoFileName= '20160226112153_s7.info';
        out.eegRetrievalFileName    = '20160226112153_s7.easy';
    case 8
        out.infoFileName            = '20160228113509_s8.info';
        out.stimulationFileName     = '20160228113508_s8.stim';
        out.eegEncodingFileName     = '20160228113509_s8.easy';
        
        out.eegRetrievalInfoFileName= '20160228115418_s8.info';
        out.eegRetrievalFileName    = '20160228115418_s8.easy';
    case 9
        out.infoStimFileName        = '20160228130559_s9.info';
        out.stimulationFileName     = '20160228130558_s9.stim';
        out.eegEncodingFileName     = '20160228130559_s9.easy';
        
        out.eegRetrievalInfoFileName= '20160228133034_s9.info';
        out.eegRetrievalFileName    = '20160228133034_s9.easy';
    case 10
        out.infoStimFileName        = '20160301111434_s10.info';
        out.stimulationFileName     = '20160301111433_s10.stim';
        out.eegEncodingFileName     = '20160301111434_s10.easy';

        out.eegRetrievalInfoFileName= '20160301113218_s10.info';
        out.eegRetrievalFileName    = '20160301113218_s10.easy';
    case 11
        out.infoStimFileName        = '20160301130748_s11.info';
        out.stimulationFileName     = '20160301130747_s11.stim';
        out.eegEncodingFileName     = '20160301130748_s11.easy';

        out.eegRetrievalInfoFileName= '20160301132838_s11.info';
        out.eegRetrievalFileName    = '20160301132838_s11.easy';
    case 12
        out.infoFileName            = '20160303111420_s12.info';
        out.stimulationFileName     = '20160303111420_s12.stim';
        out.eegEncodingFileName     = '20160303111420_s12.easy';
        
        out.eegRetrievalInfoFileName= '20160303113304_s12.info';
        out.eegRetrievalFileName    = '20160303113304_s12.easy';
    case 13
        out.infoStimFileName        = '20160305135814_s13.info';
        out.stimulationFileName     = '20160305135813_s13.stim';
        out.eegEncodingFileName     = '20160305135814_s13.easy';
        
        out.eegRetrievalInfoFileName= '20160305141652_s13.info';
        out.eegRetrievalFileName    = '20160305141652_s13.easy';
    case 14
        out.infoStimFileName        = '20160313102344_s14.info';
        out.stimulationFileName     = '20160313102343_s14.stim';
        out.eegEncodingFileName     = '20160313102344_s14.easy';

        out.eegRetrievalInfoFileName= '20160313104242_s14.info';
        out.eegRetrievalFileName    = '20160313104242_s14.easy';
    case 15
        out.infoStimFileName        = '20160313121124_s15.info';
        out.stimulationFileName     = '20160313121123_s15.stim';
        out.eegEncodingFileName     = '20160313121124_s15.easy';

        out.eegRetrievalInfoFileName= '20160313123140_s15.info';
        out.eegRetrievalFileName    = '20160313123140_s15.easy';
    case 16
        out.infoStimFileName        = '20160316093009_s16.info';
        out.stimulationFileName     = '20160316093008_s16.stim';
        out.eegEncodingFileName     = '20160316093009_s16.easy';

        out.eegRetrievalInfoFileName= '20160316094719_s16.info';
        out.eegRetrievalFileName    = '20160316094719_s16.easy';
    case 17
        out.infoStimFileName        = '20160316110653_s17.info';
        out.stimulationFileName     = '20160316110652_s17.stim';
        out.eegEncodingFileName     = '20160316110653_s17.easy';

        out.eegRetrievalInfoFileName= '20160316112606_s17.info';
        out.eegRetrievalFileName    = '20160316112606_s17.easy';
    case 18
        out.infoStimFileName        = '20160317110642_s18.info';
        out.stimulationFileName     = '20160317110641_s18.stim';
        out.eegEncodingFileName     = '20160317110642_s18.easy';

        out.eegRetrievalInfoFileName= '20160317112526_s18.info';
        out.eegRetrievalFileName    = '20160317112526_s18.easy';    
    case 19
        out.infoStimFileName        = '20160321093700_s19.info';
        out.stimulationFileName     = '20160321093659_s19.stim';
        out.eegEncodingFileName     = '20160321093700_s19.easy';

        out.eegRetrievalInfoFileName= '20160321095943_s19.info';
        out.eegRetrievalFileName    = '20160321095943_s19.easy';    
    case 20
        out.infoStimFileName        = '20160321111104_s20.info';
        out.stimulationFileName     = '20160321111104_s20.stim';
        out.eegEncodingFileName     = '20160321111104_s20.easy';

        out.eegRetrievalInfoFileName= '20160321112949_s20.info';
        out.eegRetrievalFileName    = '20160321112949_s20.easy';    
    otherwise
        error('subject does not exist')
end

