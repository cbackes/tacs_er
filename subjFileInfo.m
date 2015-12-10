function out = subjFileInfo(subj)

out = [];
out.dataPath                = ['~/Google Drive/Research/tACS/tACS_ER_task/data/tacs_enc/s' num2str(subj) '/' ];
out.infoStimFileName        = [];
out.stimulationFileName     = [];
out.eegEncodingFileName     = [];
out.eegRetrievalInfoFileName= [];
out.eegRetrievalFileName    = [];

switch  subj
    case 0
        warning([' although data was collected for this subject, markers did not '...
            'work. '])
    case 1
        warning('no subject')
    case 2
        warning('self-paced retrieval task. changed from subject 4 on')
        out.infoStimFileName        = '20151121103036_s2.info';
        out.stimulationFileName     = '20151121103035_s2.stim';
        out.eegEncodingFileName     = '20151121103036_s2.easy';
        out.eegRetrievalInfoFileName= '20151121104805_s2.info';
        out.eegRetrievalFileName    = '20151121104805_s2.easy';
    case 3
        warning('self-paced retrieval task. changed from subject 4 on')
        out.infoStimFileName        = '20151121132947_s3.info';
        out.stimulationFileName     = '20151121132946_s3.stim';
        out.eegEncodingFileName     = '20151121132947_s3.easy';
        out.eegRetrievalInfoFileName= '20151121134627_s3.info';
        out.eegRetrievalFileName    = '20151121134627_s3.easy';
    case 4
        out.infoStimFileName        = '20151124132824_s4.info';
        out.stimulationFileName     = '20151124132824_s4.stim';
        out.eegEncodingFileName     = '20151124132824_s4.easy';
        out.eegRetrievalInfoFileName= '20151124134534_s4.info';
        out.eegRetrievalFileName    = '20151124134534_s4.easy';
    case 5
        out.infoStimFileName        = '20151124152044_s5.info';
        out.stimulationFileName     = '20151124152043_s5.stim';
        out.eegEncodingFileName     = '20151124152044_s5.easy';
        out.eegRetrievalInfoFileName= '20151124153808_s5.info';
        out.eegRetrievalFileName    = '20151124153808_s5.easy';
    case 6
        out.infoStimFileName        = '20151128103023_s6.info';
        out.stimulationFileName     = '20151128103023_s6.stim';
        out.eegEncodingFileName     = '20151128103023_s6.easy';
        out.eegRetrievalInfoFileName= '20151128104644_s6.info';
        out.eegRetrievalFileName    = '20151128104644_s6.easy';
    case 7
        out.infoStimFileName        = '20151202092812_s7.info';
        out.stimulationFileName     = '20151202092812_s7.stim';
        out.eegEncodingFileName     = '20151202092812_s7.easy';
        out.eegRetrievalInfoFileName= '20151202094551_s7.info';
        out.eegRetrievalFileName    = '20151202094551_s7.easy';
    case 8
        out.infoFileName            = '20151202112009_s8.info';
        out.stimulationFileName     = '20151202112009_s8.stim';
        out.eegEncodingFileName     = '20151202112009_s8.easy';
        out.eegRetrievalInfoFileName= '20151202113931_s8.info';
        out.eegRetrievalFileName    = '20151202113931_s8.easy';
    case 9
        out.infoStimFileName        = '20151202132023_s9.info';
        out.stimulationFileName     = '20151202132023_s9.stim';
        out.eegEncodingFileName     = '20151202132023_s9.easy';
        out.eegRetrievalInfoFileName= '20151202133758_s9.info';
        out.eegRetrievalFileName    = '20151202133758_s9.easy';
    case 10
        out.infoStimFileName        = '20151203142200_s10.info';
        out.stimulationFileName     = '20151203142200_s10.stim';
        out.eegEncodingFileName     = '20151203142200_s10.easy';
        out.eegRetrievalInfoFileName= '20151203144126_s10.info';
        out.eegRetrievalFileName    = '20151203144126_s10.easy';
    case 11
        out.infoStimFileName        = '20151206122419_s11.info';
        out.stimulationFileName     = '20151206122419_s11.stim';
        out.eegEncodingFileName     = '20151206122419_s11.easy';
        out.eegRetrievalInfoFileName= '20151206123950_s11.info';
        out.eegRetrievalFileName    = '20151206123950_s11.easy';
    otherwise
        error('subject does not exist')
end

