% ensure that the time and date of execution are not stored in the provenance information
global ft_default
ft_default.trackcallinfo = 'no';

% allow the user to specify the path where additional data is present, e.g. the channel layout or anatomy files
%if exist('path', 'var')
%    addpath(path)
%end

%%
%-------------------------------------------------------------------------------
% The subject list file should be an ascii file with 5 columns and each row corresponding to a single subject
% COLUMNS:
% 1:   Subject ID
% 2:   0 or 1. flag signifying if this subject has good RESTING STATE data from ALL 3 scans
% 3:   0 or 1. flag signifying if this subject has good MOTOR TASK data from BOTH scans
% 4:   0 or 1. flag signifying if this subject has good WORKING MEMORY data from BOTH scans
% 5:   0 or 1. flag signifying if this subject has good STORY MATH data from BOTH scans
if ~exist( 'subjlistfile','var')
    error('subjlistfile should be specified')
end
subjlistmat=load(subjlistfile);
if size(subjlistmat,2)~=5,
    error('subjlistfile should have 5 columns');
end

%%
%-------------------------------------------------------------------------------
if ~exist( 'avgspecdir','var')
    error('avgspecdir should be specified. This is the directory where the files with the averaged spectra will be saved for all subjects have been saved. i.e. /HCP/scratch/meg/frankenstein/workspace/hcp/findpeaks/analysis/');
end
if ~exist( 'savealphabetadir','var')
    error('savealphabetadir should be specified. This is the directory where the files with peak alpha and beta frequencis for each subject will be saved . i.e. /HCP/scratch/meg/frankenstein/workspace/hcp/findpeaks/genericpeaks/');
end
% Files with the averaged spectra are saved with the naming
%[subjID,'_ABpeakfile.mat'];

%{
if runSite==1,
    avgspecdir=['/mnt/hps/slurm/michalareasg/pipeline/TestDetectFreqPeaks/analysis/'];
    savealphabetadir=['/mnt/hps/slurm/michalareasg/pipeline/TestDetectFreqPeaks/genericpeaks'];
    
elseif runSite==2,
    avgspecdir=['/HCP/scratch/meg/frankenstein/workspace/hcp/findpeaks/analysis/'];
    savealphabetadir=['/HCP/scratch/meg/frankenstein/workspace/hcp/findpeaks/genericpeaks'];
end;
%}


subjLogInfo=subjlistmat;
indWithIn=find((subjLogInfo(:,2)>0)|(subjLogInfo(:,3)>0)|(subjLogInfo(:,4)>0)|(subjLogInfo(:,5)>0));
numList=unique(subjLogInfo(:,1));

subjList=cellfun(@(x) num2str(x) , num2cell(numList),'Uniformoutput',false);

Nsubj=length(indWithIn);

for iSubj=1:Nsubj
    %subjID=subjList{indWithIn(iSubj)};
    
    subjID=num2str(subjLogInfo(indWithIn(iSubj),1));
    
    hasRestin=subjLogInfo(indWithIn(iSubj),2);
    hasMotort=subjLogInfo(indWithIn(iSubj),3);
    hasWrkmem=subjLogInfo(indWithIn(iSubj),4);
    hasStoryM=subjLogInfo(indWithIn(iSubj),5);
    
    
    inspectra_Restin=[];
    inspectra_Motort=[];
    inspectra_Wrkmem=[];
    inspectra_StoryM=[];

    if hasMotort,
        spefile=[avgspecdir,'casespectra_',subjID,'_Motort.mat'];
        load(spefile,'casespectra','allfreqs');
        inspectra_Motort=casespectra;
    end
    if hasWrkmem,
        spefile=[avgspecdir,'casespectra_',subjID,'_Wrkmem.mat'];
        load(spefile,'casespectra','allfreqs');
        inspectra_Wrkmem=casespectra;
    end
    if hasStoryM,
        spefile=[avgspecdir,'casespectra_',subjID,'_StoryM.mat'];
        load(spefile,'casespectra','allfreqs');
        inspectra_StoryM=casespectra;
    end
    if hasRestin,
        spefile=[avgspecdir,'casespectra_',subjID,'_Restin.mat'];
        load(spefile,'casespectra','allfreqs');
        inspectra_Restin=casespectra;
    end
    
    
    infreqs=allfreqs;
    
    
    [parH]=hcp_avgspecmanualalphabetagui(subjID,inspectra_Restin,inspectra_Motort,inspectra_Wrkmem,inspectra_StoryM,infreqs,savealphabetadir);
    
    
    is6in=0;
    while ~is6in
    inNumber=input('Input number 6 to continue:');
    if inNumber==6,
       is6in=1; 
    end
    end
    
    
end