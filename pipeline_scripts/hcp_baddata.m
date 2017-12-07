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

tok2=tokenize(scanid,'-');
scanmnem=tok2{2};

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
pipelineprefix='baddata';
resultprefix = sprintf('%s_%s_%s', experimentid, scanid, pipelineprefix);

% dsname = fileparts(filename);

% Read in the raw data in one piece
cfg = [];
cfg.dataset = filename;
dataraw=ft_preprocessing(cfg);

dataraw = ft_selectdata(dataraw,'channel','MEG');

Nsamples=size(dataraw.trial{1},2);

isTask=strcmp(scanmnem,'Wrkmem')|strcmp(scanmnem,'Motort')|strcmp(scanmnem,'StoryM');
isResting=strcmp(scanmnem,'Restin');
%=========================================================
if isTask,
    trialFunName=['trialfun_',scanmnem];
    trialDefFuncFile=['alltrialdefparams_',scanmnem];
    trialDefFuncHandle=str2func(trialDefFuncFile);
    %------------------------------------------
    savesuffix_trialinfo='trialinfo';
    %--- Extract trials definition
    allTrlCfgs=trialDefFuncHandle();
    Ncasegroups=length(allTrlCfgs);
    outputTrialInfoSummaryFile = [resultprefix,'_raw',savesuffix_trialinfo,'_QC.txt'];
    montage         = hcp_exgmontage(subjectid, experimentid, scanid);
    iCaseGroup=1;
    %------------------------------------------
    trlCfg                 = allTrlCfgs{iCaseGroup};
    trlCfg.datafile        = filename;
    trlCfg.trialfun        = trialFunName;
    if ~isempty(regexp(scanid,'Motor'))
        trlCfg.trialdef.montage=montage;
    end
    %------------------------------------------
    try
        eval(['[trl,trlInfoColDescr,trialSummary,scanStartSamp,scanEndSamp,warninfo]=',trialFunName,'(trlCfg);']);
        % Save the summary ascii file
        hcp_write_ascii(outputTrialInfoSummaryFile,'trialSummary'); % Trial Summary should be the same no matter what the input arguments are. This is because the trial definition function creates all information about the trial definition. This is what the summary contains. The trl field only contains the trials defined by the input arguments.
        ErrorFile=['ERROR_',outputTrialInfoSummaryFile];
        delete([ErrorFile,'*']); % Remove any error files present from previous runs
        WarningFile=['WARNING_',outputTrialInfoSummaryFile];
        if isempty(warninfo)
            delete([WarningFile,'*']); % Remove any warning files present from previous runs
        else
            hcp_write_ascii(WarningFile,'warninfo');
        end
        
    catch
        ErrorMessage={'Something went wrong with trialdefinition function. hcp_baddata continued without any information from these functions i.e. manual bad segments at the start and end'};
        ErrorFile=['ERROR_',outputTrialInfoSummaryFile];
        delete([ErrorFile,'*']);
        WarningFile=['WARNING_',outputTrialInfoSummaryFile];
        delete([WarningFile,'*']); % Remove any warning files present from previous runs
        delete([outputTrialInfoSummaryFile,'*']);
        hcp_write_ascii(ErrorFile,'ErrorMessage'); % Trial Summary should be the same no matter what the input arguments are. This is because the trial definition function creates all information about the trial definition. This is what the summary contains. The trl field only contains the trials defined by the input arguments.
        scanStartSamp=nan;
        scanEndSamp=nan;
    end
elseif isResting,
    idealDurSecs=5*60;  % 5 minutes Restting state
    actualDurSecs=(Nsamples./dataraw.fsample);
    durDif=actualDurSecs-idealDurSecs;
    if durDif<0,
        scanStartSamp=nan;
        scanEndSamp=nan;
        error('The resting state scan is shorter than 5 minutes');
    elseif durDif>0,
        frontPadSamps=floor((2/3)*durDif*dataraw.fsample);
        idealDurSamps=floor(idealDurSecs*dataraw.fsample);
        backPadSamps=Nsamples-frontPadSamps-idealDurSamps;
        if backPadSamps<0,
            error('There is a problem with splitting the excess time in 2/3 front and 1/3 back pad')
        end
        scanStartSamp=frontPadSamps;
        scanEndSamp=frontPadSamps+idealDurSamps;
    else
        scanStartSamp=1;
        scanEndSamp=Nsamples;
    end
    
    
    
end
%=========================================================
%=========================================================
%=========================================================
if ~exist([resultprefix '_manual_badsegments.txt'],'file');
    %--------------
    manual_badsegments=[];
    
    if isTask|isResting
        if scanStartSamp>1
            manual_badsegments=[manual_badsegments
                [1 scanStartSamp] ];
        end
        if scanEndSamp<Nsamples
            manual_badsegments=[manual_badsegments
                [scanEndSamp Nsamples] ];
        end
        
    end
    hcp_write_ascii(sprintf('%s_manual_badsegments.txt', resultprefix), 'manual_badsegments');
    %--------------
else
    %--------------
    manual_badsegments=[];
    
    hcp_read_ascii([resultprefix '_manual_badsegments.txt']);
    if isTask|isResting
        if scanStartSamp>1
            manual_badsegments=[manual_badsegments
                [1 scanStartSamp] ];
        end
        if scanEndSamp<Nsamples
            manual_badsegments=[manual_badsegments
                [scanEndSamp Nsamples] ];
        end
        
        
        tmpbadsample=zeros(1,Nsamples);
        for iLine=1:size(manual_badsegments,1)
            tmpbadsample(manual_badsegments(iLine,1):manual_badsegments(iLine,2))=1;
        end
        tmpflank  = diff([0 tmpbadsample 0]);
        tmpmanual_badsegments=[];
        tmpmanual_badsegments(:,1) = find(tmpflank== 1);
        tmpmanual_badsegments(:,2) = find(tmpflank==-1) - 1;
        manual_badsegments=tmpmanual_badsegments;
        
        hcp_write_ascii(sprintf('%s_manual_badsegments.txt', resultprefix), 'manual_badsegments');
    end
    %--------------
end
%=========================================================
%=========================================================
%=========================================================


if ~exist([resultprefix '_manual_badchannels.txt'],'file');
    manual_badchannels=[];
    hcp_write_ascii(sprintf('%s_manual_badchannels.txt', resultprefix), 'manual_badchannels');
else
    hcp_read_ascii([resultprefix '_manual_badchannels.txt']);
end




%plotlayout='4D248.mat';
corrthreshold=0.4;
zthreshold=20;
stdratiothreshold=0.5;



options_correl = {'resultprefix', resultprefix, 'corrthreshold',      corrthreshold};
options_zscore = {'resultprefix', resultprefix, 'zthreshold',         zthreshold, 'skipped_intervals', manual_badsegments};
options_std    = {'resultprefix', resultprefix, 'stdratiothreshold',  stdratiothreshold};


tic;
[badchannel_neighcorr, dataclean1] = hcp_qc_neighcorrel(dataraw, options_correl);
baddata_time.corel = toc;


tic;
[badsegment_zscore,    dataclean2] = hcp_qc_zscore(dataclean1, options_zscore);
baddata_time.zscore = toc;
% we don't need this any more
clear dataclean1

tic;
[badchannel_neighstd ]= hcp_qc_neighstdratio(dataclean2, options_std);
baddata_time.std = toc;
% we don't need this any more
clear dataclean2


bandpass = [1 150]; % band pass frequency
if strcmp(subjectid,'CP10128') || strcmp(subjectid,'CP10129')
    % these were recorded in Glasgow
    bandstop = [49 51 ; 99 101]; % band stop frequency
else
    % these were recorded in St Louis
    bandstop = [59 61 ; 119 121]; % band stop frequency
end


channels          = dataraw.label'; % Here should be the channels of data to be analysed
knownbadchannels  = unique(cat(2, badchannel_neighcorr, badchannel_neighstd, manual_badchannels)); % previously known bad channels
if(~isempty(knownbadchannels)) knownbadchannels=strcat('-',knownbadchannels); end
%     resamplefs  = 500;
%     options_ICA = {'resultprefix', resultprefix, 'resamplefs', resamplefs, 'channels', channels,'knownbadchannels', knownbadchannels, 'skipped_intervals', manual_badsegments, 'bandpass', bandpass, 'bandstop', bandstop,'modality',modality};
options_ICA = {'resultprefix', resultprefix, 'channels', channels, 'knownbadchannels', knownbadchannels, 'skipped_intervals', manual_badsegments, 'bandpass', bandpass, 'bandstop', bandstop};

tic;
[badsegment_ica, badchannel_ica] = hcp_ICA_qualitycheck_pipeline(dataraw, options_ICA);

disp('after ICA qc')
baddata_time.ica=toc;

% we don't need this any more
clear dataraw

badchannel_neighcorr
badchannel_neighstd
badchannel_ica

% combine the bad channels and compute their correspondence
disp('bad ch concatenation')
badchannel_all = cat(2, badchannel_neighcorr, badchannel_neighstd, badchannel_ica, manual_badchannels);
disp('1')
badchannel_all = unique(badchannel_all);
if isempty(badchannel_all)% This is because in Matlab 2013a unique of an empty matrix return a 0x1 empty matrix
    badchannel_ica=[{'A75'} {'A246'}]; % Francesco: we shoud define a list of known brocken channels
    badchannel_all=[{'A75'} {'A246'}];
end
% disp('2')
% badchannel_correspondence = 0;
% disp('3')
% for j=1:length(badchannel_all)
%     if ismember(badchannel_all{j}, badchannel_neighcorr) && ismember(badchannel_all{j}, badchannel_neighstd) && ismember(badchannel_all{j}, badchannel_ica)
%         badchannel_correspondence = badchannel_correspondence + 1;
%     end
% end % for all bad channels
% disp('4')
% if isempty(badchannel_all)
%     badchannel_correspondence = 1; % all methods agree that there is not any bad channel
% else
%     badchannel_correspondence = badchannel_correspondence./length(badchannel_all); % compute the ratio
% end
disp('5')
% collect the results in a structure
badchannel = [];
badchannel.neigh_corr      = badchannel_neighcorr;
badchannel.neigh_stdratio  = badchannel_neighstd;
badchannel.ica             = badchannel_ica;
badchannel.manual          = manual_badchannels;
badchannel.all             = badchannel_all;
% badchannel.correspondence  = badchannel_correspondence;

% write it to an ascii file, this works just like a normal save command
hcp_write_ascii(sprintf('%s_badchannels.txt', resultprefix), 'badchannel');

disp('bad segments concatenation')
% collect the results in a structure
maxsmp1 = max(badsegment_zscore(:));
maxsmp2 = max(badsegment_ica(:));
maxsmp3 = max(manual_badsegments(:));
maxsmp  = max([maxsmp1, maxsmp2, maxsmp3]);
badsample = zeros(1,maxsmp);
if ~isempty(badsegment_zscore)
    for j=1:size(badsegment_zscore,1)
        badsample(badsegment_zscore(j,1):badsegment_zscore(j,2)) = 1;
    end
end
if ~isempty(badsegment_ica)
    for j=1:size(badsegment_ica,1)
        badsample(badsegment_ica(j,1):badsegment_ica(j,2)) = 1;
    end
end
if ~isempty(manual_badsegments)
    for j=1:size(manual_badsegments,1)
        badsample(manual_badsegments(j,1):manual_badsegments(j,2)) = 1;
    end
end

if ~isempty(badsample)
    flank  = diff([0 badsample 0]);
    badsegment_all=[];
    badsegment_all(:,1) = find(flank== 1);
    badsegment_all(:,2) = find(flank==-1) - 1;
    
    badsegment = [];
    badsegment.zscore    = badsegment_zscore;
    badsegment.ica       = badsegment_ica;
    badsegment.manual    = manual_badsegments;
    badsegment.all       = badsegment_all;
else
    badsegment = [];
    badsegment.zscore    = [];
    badsegment.ica       = [];
    badsegment.manual    = [];
    badsegment.all       = [];
end

if(~isempty(badsegment.ica) && ~isempty(badsegment.manual))
flagica=zeros(1,badsegment.manual(end,2));
for i=1:size(badsegment.ica,1)
    flagica(badsegment.ica(i,1):badsegment.ica(i,2))=1;
end

for i=1:size(badsegment.manual,1)
flagica(badsegment.manual(i,1):badsegment.manual(i,2))=0;
end
diffflag=diff(flagica);
junk=find(diffflag==1);
badsegmentica=zeros(size(junk,2),2);
badsegmentica(:,1)=junk'+1;
junk=find(diffflag==-1);
badsegmentica(:,2)=junk';
badsegment.ica=badsegmentica;
end


% write it to an ascii file, this works just like a normal save command
hcp_write_ascii(sprintf('%s_badsegments.txt', resultprefix), 'badsegment');
% hcp_write_ascii(sprintf('%s_times.txt', resultprefix), 'baddata_time');

% ensure that the expected pipeline output is present
hcp_check_pipelineoutput('baddata', 'subject', subjectid, 'experiment', experimentid, 'scan', scanid);

clear badsegment badsegment_all badchannel badchannel_neighcorr badsegment_zscore badchannel_neighstd badsegment_ica badchannel_ica channels 
