%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2011-2014 by the Human Connectome Project, WU-Minn Consortium (1U54MH091657)
%
% This file is part of megconnectome.
%
% megconnectome is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% megconnectome is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with megconnectome.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setup the execution environment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opengl software;

% ensure that the time and date of execution are not stored in the provenance information
global ft_default
ft_default.trackcallinfo = 'no';

% allow the user to specify the path where additional data is present, e.g. the channel layout or anatomy files
if exist('path', 'var')
    addpath(path)
end

if ~exist('filename', 'var')
    error('filename should be specified')
end

% the filename is assumed to be something like
% 'rawdatadir/Phase1MEG/Subjects/SUBJECTID/Experiments/SUBJECTID_MEG/Scans/1-Rnoise_MNN_V1/Resources/4D/c,rfDC'

if ~exist('logfilename', 'var')   % This parameter is used in the case that an additional Logfile is required to slipt the trials
    logfilename='';
end

tok = tokenize(filename, '/');

if ~exist('subjectid', 'var')
    subjectid = tok{end-7};
end

if ~exist('experimentid', 'var')
    experimentid = tok{end-5};
end

if ~exist('scanid', 'var')
    scanid = tok{end-3};
end

if ~exist('pipelinedatadir', 'var')
    pipelinedatadir = hcp_pathdef;
end

% print the matlab and megconnectome version to screen for provenance
ver('megconnectome')

% print the value of all local variables to screen for provenance
w = whos;
w = {w.name};
w = setdiff(w, {'w', 'ans'});
for i=1:length(w)
    fprintf(hcp_printstruct(w{i}, eval(w{i})));
end

% change to the location of the processed data (input and output)
cd(pipelinedatadir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% execute the pipeline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hcp_check_pipelineoutput('baddata', 'subject', subjectid, 'experiment', experimentid, 'scan', scanid);
hcp_check_pipelineoutput('icaclass', 'subject', subjectid, 'experiment', experimentid, 'scan', scanid);

badchansuffix='baddata_badchannels'; % mnemonic for file that contains bad channel and segment info
badsegmsuffix='baddata_badsegments'; % mnemonic for file that contains bad channel and segment info
icainfosuffix='icaclass'; % mnemonic to save results


tok2         = tokenize(scanid, '-');
scanmnem     = tok2{2};

%------- The following is just for the cases where the suffix "_Run1 or
%Run2" has been added to the scanid in order to differentiate between 2
%different runs of the same paradigm. i.e. The way Robert has saved data in
%his database for subject CP10168.
indRunStr=regexp(scanmnem,'_Run');
if ~isempty(indRunStr),
    scanmnem=scanmnem(1:indRunStr(1)-1);
end


resultprefix = sprintf('%s_%s', experimentid, scanid);


% the location of the dataset, i.e. the c,rfDC file with full path
fname = filename;
logfname=logfilename;
%======================================

%--- IS the following used anywhere?????
%     if (strcmp(subjectid ,'CP10128')|strcmp(subjectid ,'CP10129'))
%         cfg.lineFreq        = [50 100];
%     else
%         cfg.lineFreq        = [60 120];
%     end
%----------------------------------------------



%==========================================================
%==========================================================
%----- Extract precleaned trialinfo  for all data groups

trialFunName=['trialfun_',scanmnem];
%trialFunName=['testtrialfun_',scanmnem];
trialDefFuncFile=['alltrialdefparams_',scanmnem];
trialDefFuncHandle=str2func(trialDefFuncFile);
%------------------------------------------
savesuffix_general='tmegpreproc';
savesuffix_trialinfo='trialinfo';
%--- Extract trials definition
%hcp_read_ascii(trialDefFile); % loads  allTrlCfgs
allTrlCfgs=trialDefFuncHandle();
Ncasegroups=length(allTrlCfgs);
outputTrialInfoFile        = [resultprefix,'_',savesuffix_general,'_',savesuffix_trialinfo];
%outputTrialInfoSummaryFile = [resultprefix,'_',savesuffix_general,'_raw',savesuffix_trialinfo,'_QC.txt'];
trlInfo=[];
trlInfo.lockNames=[];
trlInfo.lockTrl=[];

montage         = hcp_exgmontage(subjectid, experimentid, scanid);

for iCaseGroup=1:Ncasegroups
    iCaseGroup
    
    trlCfg                 = allTrlCfgs{iCaseGroup};
    trlCfg.datafile        = fname;
    trlCfg.trialfun        = trialFunName;
    
    trlInfo.lockNames{iCaseGroup}=trlCfg.trialdef.lockMnem;
    
    if ~isempty(regexp(scanid,'Motor'))
        trlCfg.trialdef.montage=montage;
    end
    if ~isempty(logfname)
        trlCfg.expResultsFile=logfname;
    end
    %----------------------------------------
    
    %------------------------------------------
    eval(['[trl,trlInfoColDescr,trialSummary,scanStartSamp,scanEndSamp ]=',trialFunName,'(trlCfg);']);
    %cfgDefTr = ft_definetrial(trlCfg);
    %trl=cfgDefTr.trl;
    
    %     if size(trl,2)>3,
    %         trialinfo=trl(:,4:end);
    %     else
    %         trialinfo=[];
    %     end
    trlInfo.lockTrl{iCaseGroup}=trl;
    trlInfo.trlColDescr{iCaseGroup}=trlInfoColDescr;
end
%==========================================================
% Save the summary ascii file
%hcp_write_ascii(outputTrialInfoSummaryFile,'trialSummary'); % Trial Summary should be the same no matter what the input arguments are. This is because the trial definition function creates all information about the trial definition. This is what the summary contains. The trl field only contains the trials defined by the input arguments.
%==========================================================
%----- Extract and clean data groups

Ncasegroups=length(trlInfo.lockTrl);

cleanTrlInfo=trlInfo;
for iCaseGroup=1:Ncasegroups
    iCaseGroup
    
    datagroupid=trlInfo.lockNames{iCaseGroup};
   
    cfg  =  [];
    cfg.badchanfile     = [resultprefix,'_',badchansuffix,'.txt'];
    cfg.badsegmfile     = [resultprefix,'_',badsegmsuffix,'.txt'];
    cfg.icainfofile     = [resultprefix,'_',icainfosuffix];
    %cfg.caseprefix      = resultprefix;
    cfg.datafile        = filename;
    cfg.badsegmode      = allTrlCfgs{iCaseGroup}.badsegmode; %tmpCfg.badsegmode='remfull'; % This is used by tmeg_preproc to know if trials should be fully removed, or bad segments should be replaced by nans(i.e. 'BU' case)
    cfg.montage         = montage; %hcp_exgmontage(subjectid, experimentid, scanid);
    cfg.outputfile      = [resultprefix,'_',savesuffix_general,'_',datagroupid];
    
    
    if strcmp(subjectid,'CP10128') || strcmp(subjectid,'CP10129')
        cfg.lineFreq        = [50 100];
    else
        cfg.lineFreq        = [60 120];
    end
    
    cfg.trl=trlInfo.lockTrl{iCaseGroup};
    
    cleanTrl=hcp_extract_allfromrun(cfg);
    cleanTrlInfo.lockTrl{iCaseGroup}=cleanTrl;
    cleanTrlInfo.trlColDescr{iCaseGroup}=[trlInfo.trlColDescr{iCaseGroup}; [num2str(size(trlInfo.trlColDescr{iCaseGroup},1)+1),'. hastrialNANs']];
    
    %hcp_check_pipelineoutput('tmegpreproc', 'subject', subjectid, 'experiment', experimentid, 'scan', scanid,'datagroupid',datagroupid);
end

trlInfo=cleanTrlInfo;
%save(cleantrialinfofile,'trlInfo');


hcp_write_matlab(outputTrialInfoFile,'trlInfo');
%hcp_write_ascii(outputTrialInfoSummaryFile,'trialSummary'); % Trial Summary should be the same no matter what the input arguments are. This is because the trial definition function creates all information about the trial definition. This is what the summary contains. The trl field only contains the trials defined by the input arguments.

hcp_check_pipelineoutput('tmegpreproc', 'subject', subjectid, 'experiment', experimentid, 'scan', scanid,'datagroup',trlInfo.lockNames);
