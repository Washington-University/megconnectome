function [badsegment_zscore,    dataclean] = hcp_qc_zscore  (datain, options_zscore)

% This function identifies segments of bad data. Id does so in three steps. Step
% A: use the z-vale of the data one-sample difference. A threshold is defined by
% the used and all the points above this threshold are identified. Then they are
% clustered so that two consequent artifact points are no more than 500msec apart.
% By this procedure squid jumps are identified. Also muscle artifacts can be
% identified by this step. Step B: The data is searched for cliiping(channels with
% flat signal) by fieldtrip's cliiping artifact detection. Step C: his step is
% using fieltrip's muscle artifact detection. This procedure is optimized (data
% filtering) for muscle artifacts detection. Artifacts identified here might
% coincide with artifacts in stepA. At the end all bad segments from the three
% steps are concatanated in a singe sequence. Also the data is split in pseudo
% trials of 1 sec duration and the ones containg artifacts are excluded and the
% remaining clean dataset is passed at the output for subsequent use.
%
% INPUTS
%     datain : This data structure should contain the
%     entire dataset in one continous trial.
%
%     options_correl: This is a row cell which key-value pairs in which the
%     parameters for the neighbour correlation are set. These paramters
%     are:
%
%               resultprefix:   This is a char containing the subject and
%                               experiment ID and which is used to name the plots created
%                               by this function
%
%               zthreshold:  The threshold of z-value of the one-sample
%               data time derivative. The same   threshold is used in step A
%               and step C.
%
% OUPUTS:
%               badsegment_zscore: This the bad segments definition matric.
%               First column is the start sample and the second column is
%               the end sample.
%
%               dataclean  : This is a data structure with pseudo trials of 1 sec duration.
%                     From this data structure the pseudo trials that coincide with artifacts are removed.

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

if length(datain.trial)>1
    error(['input dataset should have the entire data set in a single trial']);
    return;
end

zthreshold     = ft_getopt(options_zscore, 'zthreshold');        % skipped intervals previously identified
skipped_intervals = ft_getopt(options_zscore, 'skipped_intervals') ;
%---Distance between successive artifacts for cluster definition

Fsample=datain.fsample;
clusterStep=0.5; % msec
artifPad=0.1;
clusterStepSamps=floor(clusterStep*Fsample); % msec
artifPadSamps=floor(artifPad*Fsample); % msec

%data1=ft_selectdata(datain,'channel',{'MEG','MEGREF','-AF*','-M*'});
data1=datain;
Nchannels=length(data1.label);

%=============================================================
% if skipped intervals have ben provided by the user then put nans there so
% they are not taking into consideration in the z-score computations
if ~isempty(skipped_intervals)
    Nskipped=size(skipped_intervals,1);
    for iSk=1:Nskipped,
        data1.trial{1}(:,skipped_intervals(iSk,1):skipped_intervals(iSk,2))=nan;
    end
end
%=============================================================

%=============================================================
%  Z value artifact identification. It also counts how many channels have
%  benn over the threshold. It also uses a time window to define clusters
%  of astrtifacts.
%===================================================
diffData=abs(diff(data1.trial{1},1,2));
meanDiff=nanmean(diffData,2);
stdDiff=nanstd(diffData,0,2);

Ntimes=size(diffData,2);
diffData=(diffData-repmat(meanDiff,1,Ntimes))./repmat(stdDiff,1,Ntimes);

zData=diffData;
clear diffData meanDiff stdDiff Ntimes;

%{
zData=diffData;
for iChan=1:Nchannels,
  % iChan
  zData(iChan,:)=(diffData(iChan,:)-meanDiff(iChan))./stdDiff(iChan);
end
%}

[i1,j1] = find(zData>zthreshold);
%==========================================================
uniqChInd=unique(i1);
uniqTimeInd=unique(j1);
chanInfo=zeros(length(uniqChInd),2);
timeInfo=zeros(length(uniqTimeInd),2);
for iCh=1:length(uniqChInd)
    chanInfo(iCh,1)=uniqChInd(iCh);
    chanInfo(iCh,2)=sum(i1==uniqChInd(iCh));
end
for iTime=1:length(uniqTimeInd)
    timeInfo(iTime,1)=uniqTimeInd(iTime);
    timeInfo(iTime,2)=sum(j1==uniqTimeInd(iTime));
end
%------------------------------
%------------------------------
if ~isempty(timeInfo)
    clusterInfo=[];
    countClust=1;
    clustChanInd=[];
    clustChanInd{countClust}=[i1(j1==timeInfo(1,1))]';
    clusterInfo(countClust,:)=timeInfo(1,1)*[1 1];
    for iTime=1:size(timeInfo,1)-1,
        tmpChanInd=i1(j1==timeInfo(iTime+1,1));
        if timeInfo(iTime+1,1)-timeInfo(iTime,1)<=clusterStepSamps
            clusterInfo(countClust,2)=timeInfo(iTime+1,1);
            clustChanInd{countClust}=unique([clustChanInd{countClust}   tmpChanInd']);
        else
            countClust=countClust+1;
            clustChanInd{countClust}=tmpChanInd';
            clusterInfo(countClust,:)=timeInfo(iTime+1,1)*[1 1];
        end
    end
    %------------------------------
    sinClustInfo=[];
    countSing=1;
    for iClust=1:size(clusterInfo,1),
        if length(clustChanInd{iClust})==1,
            sinClustInfo(countSing,:)=[iClust clustChanInd{iClust}];
            countSing=countSing+1;
        end
        
    end
    dualClustInfo=[];
    countDual=1;
    for iClust=1:size(clusterInfo,1),
        if length(clustChanInd{iClust})==2,
            dualClustInfo(countDual,:)=[iClust clustChanInd{iClust}];
            countDual=countDual+1;
        end
    end
    
    
else
    
    clusterInfo=[];
    
end

NtotSamps=size(zData,2);
if ~isempty(clusterInfo)
    for iClust=1:size(clusterInfo,1)
        if (clusterInfo(iClust,1)-artifPadSamps) >0
            clusterInfo(iClust,1)=clusterInfo(iClust,1)-artifPadSamps;
        end
        if (clusterInfo(iClust,2)+artifPadSamps) < NtotSamps
            clusterInfo(iClust,2)=clusterInfo(iClust,2)+artifPadSamps;
        end
    end
end
clear zData;
%================================================================
%--- Perform Filetrip artifact detection for clipping and muscle
cfg=[];
cfg.continuous='no';
cfg.trl=[data1.sampleinfo 0];% This is done to overcome issue with ft_artifact_zvalue.m
%--- Settings for Clipping
cfg.artfctdef.clip.pretim   = 0.01; % pre-artifact rejection-interval in seconds
cfg.artfctdef.clip.psttim   = 0.01; % post-artifact rejection-interval in seconds
cfg.artfctdef.clip.thresh   = 0.05; % minimum duration in seconds of a datasegment with consecutive identical samples to be considered as 'clipped'
%--- Settings for Muscle
cfg.artfctdef.muscle.cutoff      = zthreshold;
cfg.artfctdef.muscle.trlpadding  = 0;
cfg.artfctdef.muscle.fltpadding  = 0;
cfg.artfctdef.muscle.artpadding  = 0.01;
cfg.artfctdef.muscle.bpfilter    = 'yes';
cfg.artfctdef.muscle.bpfreq      = [110 140];
cfg.artfctdef.muscle.bpfiltord   = 4;
cfg.artfctdef.muscle.bpfilttype  = 'but';
cfg.artfctdef.muscle.hilbert     = 'yes';
cfg.artfctdef.muscle.boxcar      = 0.2;
%---------------------------------------------------

% Downsample to avoid huge memory consumption by hilbert transform
if data1.fsample>1200
    % Downsample to avoid huge memory consumption by hilbert transform
    defDSratio=3;
    downcfg=[];
    downcfg.resamplefs = data1.fsample./defDSratio;
    downcfg.detrend ='no';
    datadown=ft_resampledata(downcfg,data1);
    
    cfg.trl=[1 length(datadown.time{1}) 0];
    
    %---------------------------------------------------
    cfg.artfctdef.clip.channel=datadown.label;
    [cfgClip,artClip]=ft_artifact_clip(cfg,datadown);
    if ~isempty(artClip) % Recover
        artClip=floor(artClip);
        artClip=[(defDSratio*artClip(:,1))-(defDSratio-1) defDSratio*artClip(:,2)];
    end
    %----------------------------------------------
    disp('FINDING MUSCLE ARTIFACTS');
    cfg.artfctdef.muscle.channel=datadown.label;
    [cfgMuscle,artMuscle]=ft_artifact_muscle(cfg,datadown);
    disp('DONE WITH FINDING MUSCLE ARTIFACTS');
    if ~isempty(artMuscle) % Recover
        artMuscle=[(defDSratio*artMuscle(:,1))-(defDSratio-1) defDSratio*artMuscle(:,2)];
    end
    clear datadown;
else
    %---------------------------------------------------
    cfg.artfctdef.clip.channel=data1.label;
    [cfgClip,artClip]=ft_artifact_clip(cfg,data1);
    artClip=floor(artClip);
    %----------------------------------------------
    
    cfg.artfctdef.muscle.channel=data1.label;
    [cfgMuscle,artMuscle]=ft_artifact_muscle(cfg,data1);
end
%========================================================
% combine the bad segments in a structure
badsegall=[clusterInfo; artClip; artMuscle];

if ~isempty(badsegall)
    [i1,j1]=find(badsegall(:,1)<1);
    if ~isempty(i1),
        badsegall(i1,1)=1;
    end
    [i1,j1]=find(badsegall(:,2)>data1.sampleinfo(2));
    if ~isempty(i1)
        badsegall(i1,2)=data1.sampleinfo(2);
    end
    maxsmp = max(badsegall(:));
    badsmp = zeros(1,maxsmp);
    for i=1:size(badsegall,1)
        badsmp(badsegall(i,1):badsegall(i,2)) = 1;
    end
    badsmp = diff([0 badsmp 0]);
    onset  = find(badsmp== 1);
    offset = find(badsmp==-1) - 1;
    badsegment_zscore = [onset(:) offset(:)];
else
    badsegment_zscore = [];
end

% try
%   badsegall=[clusterInfo; artClip; artMuscle];
%   if ~isempty(badsegall)
%     [i1,j1]=find(badsegall(:,1)<1);
%     if ~isempty(i1),
%       badsegall(i1,1)=1;
%     end
%     [i1,j1]=find(badsegall(:,2)>NtotSamps);
%     if ~isempty(i1),
%       badsegall(i1,1)=NtotSamps;
%     end
%
%     badsegall=sortrows(badsegall,1);
%
%     dummyBin=zeros(1,badsegall(end,2));
%     for iBad=1:size(badsegall,1)
%       dummyBin(badsegall(iBad,1):badsegall(iBad,2))=1;
%     end
%     dummyDif=[0 diff(dummyBin)];
%     artiStart=find(dummyDif>0);
%     artiEnd=[find(dummyDif<0)  badsegall(end,2)];
%     if length(artiStart)==length(artiEnd)
%       badsegment_zscore=[artiStart' artiEnd'];
%     elseif length(artiStart)-length(artiEnd)==1
%       badsegment_zscore=[artiStart(1:end-1)' artiEnd'];
%     elseif length(artiStart)-length(artiEnd)==-1
%       badsegment_zscore=[[1 artiStart]' artiEnd'];
%     end
%
%   else
%     badsegment_zscore=[];
%   end
%
% catch
%   errorfile=['errorcase_',resultprefix];
%   fh=fopen(errorfile,'w+');
%   fprintf('%s',resultprefix);
%   fclose(fh);
% end % try

%========================================================
% -- SPLIT DATA IN PSEUDO TRIALS OF ONE SECOND
%=============================================
defTrDur=1; %sec
defTrSamps=floor(defTrDur.*Fsample);
nSamples=size(data1.trial{1},2);

Ntrials=floor(nSamples./defTrSamps);

trlStart=defTrSamps.*((0:Ntrials-1)')+1;
trlEnd=defTrSamps.*((1:Ntrials)');
trlDur=zeros(Ntrials,1);

trl=[trlStart trlEnd trlDur];


cfg=[];
cfg.trl=trl;
data1=ft_redefinetrial(cfg,data1);

cfg=[];
cfg.demean='yes';
data1=ft_preprocessing(cfg,data1);
%===================================================

%========================================================
%---- Reject Artifact intervals
artCfg=[];
artCfg.artfctdef.zvalue.artifact = [skipped_intervals;  badsegment_zscore]; % remove maually set skipped intervals and identified bad segments
dataclean = ft_rejectartifact(artCfg, data1); clear data1;
%=======================================
clearvars -except badsegment_zscore dataclean



