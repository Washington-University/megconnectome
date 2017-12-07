function[outStatus]=hcp_srcavgdics_contrasts(inCfg)

% This function loads time frequency data and computes its trial average  as well
% as the average of its planar gradient
% outdatafile is the file where the averaged data will be save
% outinfofile is the ZIP file where any plotted figures will be saved

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
outStatus=-1;
global saveExtraDir;

%================================================
% Get inputs
if isfield(inCfg,'savedir')
    saveExtraDir=inCfg.savedir;
else
    saveExtraDir=[];
end

subjectid=inCfg.subjectid;
experimentid=inCfg.experimentid;
multiscanid=inCfg.multiscanid; % cell with the scan id of  both files. They shoudl be in scannning order
contrastlist=inCfg.contrastlist;
anatomydir = inCfg.anatomydir;
gridtype=inCfg.gridtype; %Can be '3D'(using 8mm) or '2D'%
bandinfo=inCfg.bandinfo;
beamflambda =inCfg.beamflambda;    % char i.e. '100%' with the regularisation lambda for the beamformer
%------------------------------------------
% Arrange input frequency bands
freqBands=reshape([bandinfo{:,2}],2,size(bandinfo,1))';
freqBandCents=mean(freqBands,2);
freqBandCents=round(freqBandCents*2)./2;
freqBandSmoothF=(freqBands(:,2)-freqBands(:,1));
freqBandNames=bandinfo(:,1);

%================================================
% Get basic information about input contrasts
Nfiles=length(contrastlist); % 1 (Motor) or 2(WM and SM) files
Ncontr=length(contrastlist{1}); % Check how many contrasts from the first file's contrasts


contrnames=[];
datagroups=[];
for iC=1:Ncontr,
    contrnames{iC}=contrastlist{1}{iC}.mnemprint;
    datagroups=unique([datagroups,contrastlist{1}{iC}.lockmode]);
end

Ngroups=length(datagroups);

%=============================================================
% Sort contralist list so that contrasts with the same datagroup of the
% first condition are in sequence. This is to avoid reloading data files
% all the time
sortcontrlist=[];
countcontr=1;
for iGroup=1:Ngroups
    for iC=1:Ncontr,
        if strcmp(contrastlist{1}{iC}.lockmode{1},datagroups{iGroup}),
            sortcontrlist{1}{countcontr}=contrastlist{1}{iC};
            if Nfiles==2
                sortcontrlist{2}{countcontr}=contrastlist{2}{iC};
            end
            countcontr=countcontr+1;
        end
    end
end
sortcontrnames=[];
sortgroupnames=[];
isCompCntr=[]; % It says if it is a comparison between conditions
isDifGroups=[]; % It says if there are 2 different groups in each contrast
for iC=1:Ncontr,
    sortcontrnames{iC}=sortcontrlist{1}{iC}.mnemprint;
    sortgroupnames{iC}=sortcontrlist{1}{iC}.lockmode;
    isCompCntr(iC,1)=length(sortgroupnames{iC})-1;
    isDifGroups(iC,1)=length(unique(sortgroupnames{iC}))-1;
end
if Nfiles==2
    is2Files=1;
else
    is2Files=0;
end
%=============================================================
%=============================================================
%=============================================================
%% LOOP OVER CONTRASTS
gridAllLF=[]; % Initialize variable where the source space with LeadFields will be placed
prevGroups={'' ''};
for iC=1:Ncontr
    
    %=====================================================================
    % Srcavglcmv is only computing single condition contrasts. Check and
    % if comparison between 2 different conditions the skip to next
    if isCompCntr(iC),
        disp(['Srcavglvmc pipeline is only computing single condition contrasts , NOT between different conditions. Skipping contrast: ', sortcontrnames{iC}]);
        continue;
    end
    %--------------------------------------------------------------------
    
    tmpTok=tokenize(multiscanid{1},'-');
    scanmnem=tmpTok{2};
    %------- The following is just for the cases where the suffix "_Run1 or
    %Run2" has been added to the scanid in order to differentiate between 2
    %different runs of the same paradigm. i.e. The way Robert has saved data in
    %his database for subject CP10168.
    indRunStr=regexp(scanmnem,'_Run');
    if ~isempty(indRunStr),
        scanmnem=scanmnem(1:indRunStr(1)-1);
    end
    
    
    
    
    %---------------------------------------------------------
    % Check if the current data group is the same as the one loaded for
    % the preceeding processed contrast
    prevMatchIndx=[];
    tmpIndx=find(strcmp(sortgroupnames{iC},prevGroups));
    if ~isempty(tmpIndx)
        prevMatchIndx{1}=tmpIndx(1);
    else
        prevMatchIndx{1}=[];
    end
    %---------------------------------------------------------
    
    % if the same don't load anything else load data group for current contrast
    if ~isempty(prevMatchIndx{1})
        eval(['tmpdata1A=alldata',num2str(prevMatchIndx{1}),'A;']);
        eval(['tmpdata1B=alldata',num2str(prevMatchIndx{1}),'B;']);
    else
        tmpdata1A=[];
        tmpdata1B=[];
    end
    
    %---If contrast from new datagroup then load data ----------------
    if isempty(tmpdata1A)
        genprefix = sprintf('%s_%s', experimentid, multiscanid{1}); tmpdatafile=[genprefix,'_tmegpreproc_',sortgroupnames{iC}{1}];
        hcp_read_matlab(tmpdatafile,'data');
        
        %=================================================================
        %{
        % ADDHOCK REMOVING SUPINE BALANCING - TEMPORARY - REMOVE BEFORE COMMITING
        tmpgrad1 = ft_datatype_sens(ft_apply_montage(data.grad, data.grad.balance.invcomp, 'inverse', 'yes', 'keepunused', 'yes'));
        tmpgrad2 = ft_datatype_sens(ft_apply_montage(tmpgrad1, tmpgrad1.balance.pca, 'inverse', 'yes', 'keepunused', 'yes'));
        tmpgrad3 = ft_datatype_sens(ft_apply_montage(tmpgrad2, tmpgrad2.balance.Supine, 'inverse', 'yes', 'keepunused', 'yes'));
        tmpgrad2 = ft_datatype_sens(ft_apply_montage(tmpgrad3, data.grad.balance.pca,'balancename','pca', 'keepunused', 'yes'));
        tmpgrad1 = ft_datatype_sens(ft_apply_montage(tmpgrad2, data.grad.balance.invcomp, 'balancename','invcomp', 'keepunused', 'yes'));
        tmpgrad1.balance.invcomp=tmpgrad1.balance.invcomp2;
        tmpgrad1.balance.pca=tmpgrad1.balance.pca2;
        tmpgrad1.balance.previous{1}='pca';
        tmpgrad1.balance.current='invcomp';
        tmpgrad1.balance=rmfield(tmpgrad1.balance,{'invcomp1','invcomp2','pca1','pca2'});
        tmpgrad1.balance.pca=rmfield(tmpgrad1.balance.pca,{'chantypeorg','chantypenew','chanunitorg','chanunitnew'});
        tmpgrad1.balance.invcomp=rmfield(tmpgrad1.balance.invcomp,{'chantypeorg','chantypenew','chanunitorg','chanunitnew'});
        data.grad=tmpgrad1;
        %}
        %=================================================================
        alldata1A=data;
        Fsample=data.fsample;
        if isempty(gridAllLF), % If this is the first processed contrast then the grid variable is empty and the grid should be loaded and processed for inverse solution.
            %--------------------------
            % At the moment the grad is taken from the first scan and it is assumed that the
            % subject has not moved significantly in the other one. Maybe a check should be made and if the difference is above a threshold then a filter should be derived for each file separately;
            
            grad=ft_convert_units(data.grad,'cm');
            fileGrid2D=[anatomydir,'/',experimentid,'_anatomy_sourcemodel_2d.mat'];
            fileGrid3D=[anatomydir,'/',experimentid,'_anatomy_sourcemodel_3d8mm.mat'];
            fileGrid3DStandard=['standard_sourcemodel3d8mm']; % Add in the compilation of the executable
            %fileAnat=[anatomydir,'/',experimentid,'_anatomy_anatomical.nii'];
            fileAnatTemplMcGill='mni_icbm152_t1_tal_nlin_sym_09a.nii';
            fileVol=[anatomydir,'/',experimentid,'_anatomy_headmodel.mat'];
            
            mriTempl=ft_read_mri(fileAnatTemplMcGill);
            
            %mriSubj=ft_read_mri(fileAnat);
            if strcmp(gridtype,'2D')
                hcp_read_matlab(fileGrid2D,'sourcemodel2d');
                grid=ft_convert_units(sourcemodel2d,'cm');
            elseif strcmp(gridtype,'3D')
                hcp_read_matlab(fileGrid3D,'sourcemodel3d');
                grid=ft_convert_units(sourcemodel3d,'cm');
                hcp_read_matlab(fileGrid3DStandard,'sourcemodel');
                gridStand=sourcemodel;% already in cm
            end
            
            hcp_read_matlab(fileVol,'headmodel');
            vol=ft_convert_units(headmodel,'cm');
            
            %----- Create Leadfields ---------------------------------------
            allChansMEG=ft_channelselection({'MEG'},grad.label);
            cfg=[];
            cfg.grid=grid; % Grid for Individual's Brain in MEG sensor space
            cfg.vol=vol;
            cfg.grad=grad;
            cfg.reducerank      = 2; %(default = 3 for EEG, 2 for MEG)
            cfg.normalize       = 'yes' ; %Normalise Leadfield: 'yes' for beamformer , 'no' for MNE
            cfg.normalizeparam  = 1; %depth normalization parameter (default = 0.5)
            cfg.feedback='no';
            cfg.channel=allChansMEG;
            if strcmp(gridtype,'2D') %JUST TO MAKE SURE ENTIRE CORTICAL SHEET IN SOURCE SPACE
                cfg.inwardshift=-1;
            end
            gridAllLF= ft_prepare_leadfield(cfg);
            gridAllLF.label=allChansMEG;
            
            
            
        end
        clear data;
        
        % ------ Load data from the second session of the same task ------
        
        if is2Files
            genprefix = sprintf('%s_%s', experimentid, multiscanid{2}); tmpdatafile=[genprefix,'_tmegpreproc_',sortgroupnames{iC}{1}];
            hcp_read_matlab(tmpdatafile,'data');
            
            %=================================================================
            %{
            % ADDHOCK REMOVING SUPINE BALANCING - TEMPORARY - REMOVE BEFORE COMMITING
            tmpgrad1 = ft_datatype_sens(ft_apply_montage(data.grad, data.grad.balance.invcomp, 'inverse', 'yes', 'keepunused', 'yes'));
            tmpgrad2 = ft_datatype_sens(ft_apply_montage(tmpgrad1, tmpgrad1.balance.pca, 'inverse', 'yes', 'keepunused', 'yes'));
            tmpgrad3 = ft_datatype_sens(ft_apply_montage(tmpgrad2, tmpgrad2.balance.Supine, 'inverse', 'yes', 'keepunused', 'yes'));
            tmpgrad2 = ft_datatype_sens(ft_apply_montage(tmpgrad3, data.grad.balance.pca,'balancename','pca', 'keepunused', 'yes'));
            tmpgrad1 = ft_datatype_sens(ft_apply_montage(tmpgrad2, data.grad.balance.invcomp, 'balancename','invcomp', 'keepunused', 'yes'));
            tmpgrad1.balance.invcomp=tmpgrad1.balance.invcomp2;
            tmpgrad1.balance.pca=tmpgrad1.balance.pca2;
            tmpgrad1.balance.previous{1}='pca';
            tmpgrad1.balance.current='invcomp';
            tmpgrad1.balance=rmfield(tmpgrad1.balance,{'invcomp1','invcomp2','pca1','pca2'});
            tmpgrad1.balance.pca=rmfield(tmpgrad1.balance.pca,{'chantypeorg','chantypenew','chanunitorg','chanunitnew'});
            tmpgrad1.balance.invcomp=rmfield(tmpgrad1.balance.invcomp,{'chantypeorg','chantypenew','chanunitorg','chanunitnew'});
            data.grad=tmpgrad1;
            %}
            %=================================================================
            
            
            alldata1B=data; clear data;
            
            
        else
            alldata1B=[];
        end
        
        
        
    else
        alldata1A=tmpdata1A;clear tmpdata1A;
        alldata1B=tmpdata1B;clear tmpdata1B;
    end
    %====================================================================
    %====================================================================
    % Up to this point data for the data group of the current contrast
    % has been either loaded ot retained from previous contrast.  This data
    % contains all possible conditions defined in it. The trials of the
    % current contrast will be extracted from this data.
    %====================================================================
    %====================================================================
    
    
    
    %====================================================================
    %====================================================================
    % Now Start processing the contrast at hand
    
    %------------------------------------------------------
    % Get data from constrast conditions and merge from different files
    sel=sortcontrlist{1}{iC}.selection{1};
    datacntr1A=ft_selectdata(alldata1A,'rpt',sel);
    if is2Files,
        sel=sortcontrlist{2}{iC}.selection{1};
        datacntr1B=ft_selectdata(alldata1B,'rpt',sel);
        cfg=[]; %
        [indChA,indChB]=match_str(datacntr1A.label,datacntr1B.label); datacntr1A=ft_selectdata(datacntr1A,'channel',indChA); datacntr1B=ft_selectdata(datacntr1B,'channel',indChB);
        datacntr1=ft_appenddata(cfg,datacntr1A,datacntr1B); clear datacntr1A datacntr2A
    else
        datacntr1=datacntr1A;
    end
    %------------------------------------------------------
    %=========================================================
    % Get Spatial Filter Type;
    % 0 : filter computed from covariance matrix of the average over trials
    % 1 : filter computed from covariance matrix from all data points
    invFiltFlag=0;
    
    invfiltertype=sortcontrlist{1}{iC}.invfiltertype;
    if isempty(invfiltertype)
        invFiltFlag=1;
    else
        if strcmp(invfiltertype,'avg')
            invFiltFlag=0;
        elseif strcmp(invfiltertype,'all')
            invFiltFlag=1;
        else
            error(['invfiltertype can be :  avg or all. Change it in contrast ']);
        end
    end
    
    %------------------------------------------------------
    %====================================================
    %====================================================
    
    
    %===================================================================
    % Get timing and baseline settings for current contrast
    %-------------------------------------------
    % --- Check Time Settings
    timeperiods=[];
    isTimeWin=[];
    timepoints=[];
    Ntimepoints=[];
    Ntimeperiods=[];
    
    iCond=1;
    timeperiods{iCond}=sortcontrlist{1}{iC}.timeperiods{iCond};
    if ~isempty(timeperiods{iCond})
        Ntimeperiods=size(timeperiods{iCond},1);
        if size(timeperiods{iCond},2)==2,  % If timeperiods has 2 columns then it means that the processing is to be preformed in periods rather than in time points
            isTimeWin{iCond}=1;
            timepoints{iCond}=mean(timeperiods{iCond},2)';
        else                               % Processing to be performed at time points
            isTimeWin{iCond}=0;
            timepoints{iCond}=timeperiods{iCond}';
        end
        Ntimepoints=length(timepoints{iCond});
        
    else
        isTimeWin{iCond}=0;
        timepoints{iCond}=[]; % Leave empty as timeperiods to signify that the original time axis must be kept
        %Ntpts=length(avg1.time);
    end
    
    %-----------------------------------------------
    % Get baseline Settings
    baseline{1}=sortcontrlist{1}{iC}.baseline{1};
    hasBaseline=0;
    if ~isempty(baseline{1}),
        hasBaseline=1;
    end
    %==================================================================
    
    %==================================================================
    
    datacntrtrim1=datacntr1; clear datacntr1;
    Ntrials1=length(datacntrtrim1.trial);
    
    %==================================================================
    %Check if the datasets have fixed or variable size variables (this is important for the computation of the covariavce matrix)
    isFixedLen1=[];
    
    npointsMat1=cellfun(@(x) size(x,2),datacntrtrim1.trial); %,'UniformOutput','false');
    if length(unique(npointsMat1))==1,
        isFixedLen1=1;
    else
        isFixedLen1=0;
    end
    
    if (~isFixedLen1)&(~isempty(timeperiods{1}))
        error('The trials for the specified timeperiods do not have fixed time lengths. This is not yet supported. WHat is supported is variable length trials with the entire trials used. For this case, do not set the timeperiods field in the contrast.');
    end
    
    %=========================================================
    %---------------------------------
    % If is a time window then just sample the data so no more than 200
    % points are processed for speed.
    newtimepoints=[];
    if ~isFixedLen1
        %  Here it is assumed that each trial has variable length bit each
        %  trial is starting from the same point (i.e. Story blocks in story/math task)
        
        %Check how many points you get with 25 msec interval between timepoints
        iniTimeStep=0.025;
        maxAllowedTimes=20000;
        
        tmpTimeStep=iniTimeStep;
        isMoreThanMax=1;
        while isMoreThanMax,
            uniqidealTimes=[];
            idealTimesPerTrl=[];
            NtotalIdealTimes=0;
            for iTrl=1:Ntrials1,
                curTrlTime=datacntrtrim1.time{:};
                totPeriodRange=[curTrlTime(1) curTrlTime(end)];
                if hasBaseline
                    totPeriodRange=[[min(baseline{1}(2),min(totPeriodRange))] [max(baseline{1}(2),max(totPeriodRange))]];
                end
                if (totPeriodRange(1)<0)&(totPeriodRange(end)>=0)
                    tmpIdealTimes=[fliplr([0:-tmpTimeStep:totPeriodRange(1)]) [tmpTimeStep:tmpTimeStep:totPeriodRange(end)]];
                elseif (totPeriodRange(1)>=0)&(totPeriodRange(end)>=0)
                    tmpStartTime=round(totPeriodRange(1)./tmpTimeStep)*tmpTimeStep;
                    tmpIdealTimes=[tmpStartTime:tmpTimeStep:totPeriodRange(end)];
                elseif (totPeriodRange(1)<0)&(totPeriodRange(end)<=0)
                    tmpEndTime=round(totPeriodRange(end)./tmpTimeStep)*tmpTimeStep;
                    tmpIdealTimes=[fliplr([tmpEndTime:-tmpTimeStep:totPeriodRange(1)])];
                end
                idealTimesPerTrl{iTrl}=tmpIdealTimes;
                uniqidealTimes=unique([uniqidealTimes tmpIdealTimes]);
                NtotalIdealTimes=NtotalIdealTimes+length(tmpIdealTimes);
            end
            if NtotalIdealTimes > maxAllowedTimes,
                overNumRatio=NtotalIdealTimes./maxAllowedTimes;
                tmpTimeStep=(ceil(ceil(1000*tmpTimeStep*overNumRatio)./5)*5)./1000;
            else
                isMoreThanMax=0;
            end
            
        end
        newtimepoints{1}=uniqidealTimes;
    else
        
        if isTimeWin{1}
            if ~isFixedLen1
                error('This is time windowed case but the data has no fixed length. No supported yet');
            end
            iniTimeStep=0.025;
            maxAllowedTimes=20000;
            
            tmpTimeStep=iniTimeStep;
            isMoreThanMax=1;
            while isMoreThanMax,
                uniqidealTimes=[];
                idealTimes=[];
                NtotalIdealTimes=0;
                
                curTrlTime=datacntrtrim1.time{1};
                totPeriodRange=[curTrlTime(1)  curTrlTime(end)];
                if hasBaseline
                    totPeriodRange=[[min(baseline{1}(2),min(totPeriodRange))] [max(baseline{1}(2),max(totPeriodRange))]];
                end
                if (totPeriodRange(1)<0)&(totPeriodRange(end)>=0)
                    tmpTotalTimes=[fliplr([0:-tmpTimeStep:totPeriodRange(1)]) [tmpTimeStep:tmpTimeStep:totPeriodRange(end)]];
                elseif (totPeriodRange(1)>=0)&(totPeriodRange(end)>=0)
                    tmpTotalTimes=[totPeriodRange(1):tmpTimeStep:totPeriodRange(end)];
                elseif (totPeriodRange(1)<0)&(totPeriodRange(end)<0)
                    tmpTotalTimes=[fliplr([totPeriodRange(end):-tmpTimeStep:totPeriodRange(1)])];
                end
                
                for iPer=1:Ntimeperiods
                    [indIn]=find((tmpTotalTimes>=timeperiods{1}(iPer,1))&(tmpTotalTimes<timeperiods{1}(iPer,2)));
                    tmpIdealTimes=tmpTotalTimes(indIn);
                    
                    idealTimes=[idealTimes tmpIdealTimes];
                    
                end
                idealTimes=unique(idealTimes);
                NtotalIdealTimes=Ntrials1*length(idealTimes);
                if NtotalIdealTimes > maxAllowedTimes,
                    overNumRatio=NtotalIdealTimes./maxAllowedTimes;
                    tmpTimeStep=(ceil(ceil(1000*tmpTimeStep*overNumRatio)./5)*5)./1000;
                else
                    isMoreThanMax=0;
                end
            end
            newtimepoints{1}=idealTimes;
        else
            if isempty(timepoints{1})
                iniTimeStep=0.025;
                maxAllowedTimes=20000;
                
                tmpTimeStep=iniTimeStep;
                isMoreThanMax=1;
                while isMoreThanMax,
                    uniqidealTimes=[];
                    idealTimes=[];
                    NtotalIdealTimes=0;
                    
                    curTrlTime=datacntrtrim1.time{1};
                    totPeriodRange=[curTrlTime(1)  curTrlTime(end)];
                    if hasBaseline
                        totPeriodRange=[[min(baseline{1}(2),min(totPeriodRange))] [max(baseline{1}(2),max(totPeriodRange))]];
                    end
                    if (totPeriodRange(1)<0)&(totPeriodRange(end)>=0)
                        tmpTotalTimes=[fliplr([0:-tmpTimeStep:totPeriodRange(1)]) [tmpTimeStep:tmpTimeStep:totPeriodRange(end)]];
                    elseif (totPeriodRange(1)>=0)&(totPeriodRange(end)>=0)
                        tmpTotalTimes=[totPeriodRange(1):tmpTimeStep:totPeriodRange(end)];
                    elseif (totPeriodRange(1)<0)&(totPeriodRange(end)<0)
                        tmpTotalTimes=[fliplr([totPeriodRange(end):-tmpTimeStep:totPeriodRange(1)])];
                    end
                    idealTimes=tmpTotalTimes;
                    idealTimes=unique(idealTimes);
                    NtotalIdealTimes=Ntrials1*length(idealTimes);
                    if NtotalIdealTimes > maxAllowedTimes,
                        overNumRatio=NtotalIdealTimes./maxAllowedTimes;
                        tmpTimeStep=(ceil(ceil(1000*tmpTimeStep*overNumRatio)./5)*5)./1000;
                    else
                        isMoreThanMax=0;
                    end
                end
                newtimepoints{1}=idealTimes;
                
                
            else
                % Here it is assumes that the baseline period is INCLUDED in
                % the timepoints
                
                newtimepoints{1}=timepoints{1};
                NtotalIdealTimes=Ntrials1*Ntimepoints;
            end
        end
    end
    
    %=========================================================
    %=========================================================
    % Set output type
    outtype='pow';
    curCM=[];
    if isfield(sortcontrlist{1}{iC},'connemetric'),
        if strcmp(sortcontrlist{1}{iC}.connemetric,'emgcoh')
            outtype='emgcoh';
        end
    end
    
    emgRefChan=[];
    if strcmp(outtype,'emgcoh')
        emgRefChan=['EMG_',sortcontrlist{1}{iC}.mnemtrl{1}];
    end
    
    refchan =  emgRefChan;
    %=========================================================
    %=========================================================
    % Perform time frequency analysis perf Frequency and do source
    % localization
    %------------------------------------------------------
    Nbands=size(freqBands,1);
    
    %membCase='LH';
    %timeperiods{1}=[-0.6:0.015:0.6];
    
    
    %--- LOOP over Frequencies
    for iBand=1:Nbands,
        
        
        curWin=freqBandSmoothF(iBand)./2;
        curFreq=freqBandCents(iBand);
        curBand=freqBands(iBand,:);
        curBandName=freqBandNames{iBand};
        
        %{
        K=3;
        tw=1./curWin;
        fw=curWin;
        
        K = 2*tw*fw-1;
        fw=(K+1)./(2*tw);
        tw=(K+1)./(2*fw);
        %}
        %=========================================================
        % Define freq analysis setting and do freq analysis
        if curFreq<30
            K=1;
        else
            K=3;
        end
        tw=min(3./curWin,0.5);
        fw=floor((K+1)./(2*tw));
        
        %------------------------------------
        % Check padding
        trlTimeRangeMax=max(cellfun(@(x) diff(x([1 end])), datacntrtrim1.time ));
        if trlTimeRangeMax >8
            tfreqPad=8*(ceil(trlTimeRangeMax./8));
        else
            tfreqPad=8;
        end
        
        %--------------------------------------
        % Perform Time Freq
        cfg=[];
        cfg.output='fourier';
        cfg.pad=tfreqPad;
        cfg.method='mtmconvol';
        cfg.foi   = curFreq;  %vector 1 x numfoi, frequencies of interest
        cfg.taper      = 'hanning';      %'dpss', 'hanning' or many others, see WINDOW (default = 'dpss')
        cfg.tapsmofrq  = fw;
        cfg.t_ftimwin  = tw;   %vector 1 x numfoi, length of time window (in seconds)
        cfg.toi        = newtimepoints{1};% timeperiods{1};  %vector 1 x numtoi, the times on which the analysis windows
        cfg.keeptrials='yes';
        tmpfreq=ft_freqanalysis(cfg,datacntrtrim1);
        
        channelsMEGEMG=tmpfreq.label;
        channelsMEG=ft_channelselection('MEG',tmpfreq.label);
        
        
        %=================================================================
        % Compute Total ross spectral density
        if (~isFixedLen1)
            %---------------------------------------------------------------
            % If trial dont have a fixed length compute cross spectral
            % density from all the data points
            tmpTotFour=reshape(squeeze(tmpfreq.fourierspctrm),length(tmpfreq.label),Ntrials1*length(newtimepoints{1}));
            [indNan1,indNan2]=find(isnan(tmpTotFour));
            tmpTotFour(:,unique(indNan2))=[];
            tmpTotSpec=tmpTotFour*ctranspose(tmpTotFour)./Ntrials1*length(newtimepoints{1}); clear tmpTotFour;
            tmpspec=ft_checkdata(tmpfreq,'cmbrepresentation','fullfast'); clear tmpfreq;
            tmpspecavg=ft_selectdata(tmpspec,'avgovertime','yes');
            tmpspecavg.crossspctrm=tmpTotSpec; clear tmpTotSpec;
            
            specavg4src=tmpspecavg;
            specavg4src=rmfield(specavg4src,'time');
            specavg4src.dimord='chan_chan_freq';
            %---------------------------------------------------------------
            % If Baseline is given compute cross spec for baseline. Here is
            % is assumed the the baseline is at the beginning of the trial
            % so all trials will contain the baseline
            if hasBaseline
                indActFirst=find(tmpspecavg.time>baseline{1}(2),1,'first');
                indActLast=length(tmpspecavg.time);
                indBaseFirst=find(tmpspecavg.time>=baseline{1}(1),1,'first');
                indBaseLast=find(tmpspecavg.time<=baseline{1}(2),1,'last');
                timeActFirst=tmpspecavg.time(indActFirst);
                timeActLast=tmpspecavg.time(indActLast);
                timeBaseFirst=tmpspecavg.time(indBaseFirst);
                timeBaseLast=tmpspecavg.time(indBaseLast);
            else
                indActFirst=1;
                indActLast=length(tmpspecavg.time);
                indBaseFirst=nan;
                indBaseLast=nan;
                timeActFirst=tmpspecavg.time(indActFirst);
                timeActLast=tmpspecavg.time(indActLast);
                timeBaseFirst=nan;
                timeBaseLast=nan;
            end
            if hasBaseline
                specavg4src_base=ft_selectdata(tmpspec,'toilim',[timeBaseFirst timeBaseLast],'avgovertime','yes');
                specavg4src_base=rmfield(specavg4src_base,'time');
                specavg4src_base.dimord='chan_chan_freq';
            else
                specavg4src_base=nan;
            end
            
            
        else
            %------------------------------------------------
            % If trial have fixed length then just average for each time
            % point in the trial first
            tmpspec=ft_checkdata(tmpfreq,'cmbrepresentation','fullfast'); clear tmpfreq;
            
            %------------------------------------------------
            % Compute total cross spectral density
            tmpspecavg=ft_selectdata(tmpspec,'avgovertime','yes');
            specavg4src=tmpspecavg;
            specavg4src=rmfield(specavg4src,'time');
            specavg4src.dimord='chan_chan_freq';
            %---------------------------------------------------------
            % Compute cross spectral density for baseline
            if hasBaseline
                indActFirst=find(tmpspec.time>baseline{1}(2),1,'first');
                indActLast=length(tmpspec.time);
                indBaseFirst=find(tmpspec.time>=baseline{1}(1),1,'first');
                indBaseLast=find(tmpspec.time<=baseline{1}(2),1,'last');
                timeActFirst=tmpspec.time(indActFirst);
                timeActLast=tmpspec.time(indActLast);
                timeBaseFirst=tmpspec.time(indBaseFirst);
                timeBaseLast=tmpspec.time(indBaseLast);
            else
                indActFirst=1;
                indActLast=length(tmpspec.time);
                indBaseFirst=nan;
                indBaseLast=nan;
                timeActFirst=tmpspec.time(indActFirst);
                timeActLast=tmpspec.time(indActLast);
                timeBaseFirst=nan;
                timeBaseLast=nan;
            end
            if hasBaseline
                specavg4src_base=ft_selectdata(tmpspec,'toilim',[timeBaseFirst timeBaseLast],'avgovertime','yes');
                specavg4src_base=rmfield(specavg4src_base,'time');
                specavg4src_base.dimord='chan_chan_freq';
            else
                specavg4src_base=nan;
            end
            
        end;
        
        %-----------------------------------------
        % Do source analysis
        srcCfg=[];
        %srcCfg.latency=timeperiods{1}([1 2]);
        srcCfg.realfilter='yes';
        srcCfg.vol=vol;
        srcCfg.grad=grad;
        srcCfg.method='dics';
        srcCfg.frequency = curFreq;
        srcCfg.fixedori='yes';
        srcCfg.feedback='text';
        srcCfg.lambda=beamflambda;
        srcCfg.projectnoise='yes';
        srcCfg.keepfilter='yes';
        srcCfg.keepleadfield='yes';
        srcCfg.keepmom='no';
        %srcCfg.refchan =  ['EMG_LH'];
        
        [indA,indB]=match_str(channelsMEG,gridAllLF.label);
        
        inIndices=find(gridAllLF.inside);
        NsrcIn=length(inIndices);
        Nsrc=size(gridAllLF.pos,1);
        
        gridCase=gridAllLF;
        gridCase.leadfield(inIndices)=cellfun(@(x,y) x(y,:),gridAllLF.leadfield(inIndices),repmat({indB},1,NsrcIn),'UniformOutput',false);
        
        srcCfg.grid=gridCase;
        srcCfg.grid.label=gridCase.label(indB);
        disp('...beamforming...');
        
        %[indA,indB]=match_str(specavg4src.label,channelsMEG);
        
        if strcmp(outtype,'emgcoh')
            srcCfg.refchan=refchan;
        end
        sourceLocAll=ft_sourceanalysis(srcCfg,specavg4src);
        clear specavg4src;
        srcCfg.filter=sourceLocAll.avg.filter; % any subsequent localization is performed with computed filters
        
        %--------------------------------------------------------------------
        % Assign output source structure
        sourceLocOut=sourceLocAll;
        sourceLocOut.freq=curFreq;
        sourceLocOut.band=curBand;
        sourceLocOut.bandname=curBandName;
        
        if (~isFixedLen1),
            totalMetric=sourceLocAll.avg.pow;
        else
            
            if (isTimeWin{1})
                totalMetric=nan(Nsrc,Ntimeperiods);
                
                for iPer=1:Ntimeperiods
                    perspecavg=ft_selectdata(tmpspec,'toilim',timeperiods{1}(iPer,:),'avgovertime','yes');
                    perspecavg=rmfield(perspecavg,'time');
                    perspecavg.dimord='chan_chan_freq';
                    tmpSourceLoc=ft_sourceanalysis(srcCfg,perspecavg);clear perspecavg;
                    totalMetric(:,iPer)=tmpSourceLoc.avg.pow;
                end
                sourceLocOut.time=mean(timeperiods{1},2);
                sourceLocOut.period=timeperiods{1};
                
            else
                if ~isempty(timepoints{1})
                    
                    NelecTot=length(channelsMEGEMG);
                    dummyInd=1:NelecTot;
                    [indA,indB]=match_str(tmpspec.label,channelsMEG);
                    dummyInd(indA)=[];
                    indMEG=indA;
                    indEMG=dummyInd;
                    channelsEMG=tmpspec.label(indEMG);
                    
                    
                    NelecMEG=length(indMEG);
                    NelecEMG=length(indEMG);
                    
                    totFilter=reshape([sourceLocAll.avg.filter{inIndices}],NelecMEG,NsrcIn)';
                    totFilterT=totFilter';
                    
                    
                    Ntimes=length(tmpspec.time);
                    
                    %---- Loop over time and project cross spec density matrix ---------
                    
                    totalMetric=nan(Nsrc  , Ntimeperiods);
                    
                    tic;
                    for iTime=1:Ntimes,
                        iTime
                        tmpSens_X=squeeze(tmpspec.crsspctrm(:,:,:,iTime));
                        tmpSensM_X=tmpSens_X(indMEG,indMEG);
                        if strcmp(outtype,'pow')
                            tmpSrc_X=totFilter*tmpSensM_X*(totFilterT);  clear tmpSensM_X;
                            tmpSrc_A=real(diag(tmpSrc_X))';              clear tmpSrc_X;
                            totalMetric(inIndices,iTime)=tmpSrc_A;
                            
                        elseif strcmp(outtype,'emgcoh')
                            
                            tmpSensE_X=tmpSens_X(indEMG,indEMG);
                            tmpSensEM_X=tmpSens_X(indEMG,indMEG);
                            %tmpSensME_X=tmpSens_X(indMEG,indEMG); clear tmpSens_X;
                            
                            
                            tmpSrc_X=totFilter*tmpSensM_X*(totFilterT);  clear tmpSensM_X;
                            tmpSrc_A=real(diag(tmpSrc_X))';              clear tmpSrc_X;
                            tmpSensE_A=real(diag(tmpSensE_X))';          clear tmpSensE_X;
                            tmpSrcEM_X=tmpSensEM_X*totFilterT;           clear tmpSensEM_X;
                            %tmpSrcME_X=totFilter*tmpSensME_X;
                            
                            refchanInd = match_str(channelsEMG,refchan);
                            
                            
                            totalMetric(inIndices,iTime)=(abs(tmpSrcEM_X(refchanInd,:)))./(sqrt(tmpSrc_A).*sqrt(tmpSensE_A(refchanInd)));
                            
                            clear tmpSrcEM_X tmpSensE_A ;
                            
                            
                        end
                    end
                    
                    sourceLocOut.time=mean(timeperiods{1},2);
                    sourceLocOut.period=timeperiods{1};
                else
                    sourceLocOut.period='all';
                    if strcmp(outtype,'pow')
                        totalMetric=sourceLocOut.avg.pow;
                    elseif strcmp(outtype,'emgcoh')
                        totalMetric=sourceLocOut.avg.coh;
                        sourceLocOut.avg=rmfield(sourceLocOut.avg,'coh');
                    end
                    
                    
                end
            end
            
        end
        clear sourceLocAll;
        if strcmp(outtype,'pow')
            sourceLocOut.avg.pow=totalMetric;
        elseif strcmp(outtype,'emgcoh')
            sourceLocOut.avg=rmfield(sourceLocOut.avg,'pow');
            sourceLocOut.avg.emgcoh=totalMetric;
        end
        
        
        %--------------------------------------------------------------------
        %Assign baseline source structure if baseline is given
        if hasBaseline
            if strcmp(outtype,'emgcoh')
                srcCfg.refchan=refchan;
                sourceLocBase=ft_sourceanalysis(srcCfg,specavg4src_base);
                sourceLocBase.avg.emgcoh=sourceLocBase.avg.coh;
                sourceLocBase.avg=rmfield(sourceLocBase.avg,{'coh','pow'});
            else
                sourceLocBase=ft_sourceanalysis(srcCfg,specavg4src_base);
                
            end
            
            sourceLocBase.time=mean(baseline{1});
            sourceLocBase.period=baseline{1};
            sourceLocBase.freq=curFreq;
            sourceLocBase.band=curBand;
            sourceLocBase.bandname=curBandName;
            
            
        else
            sourceLocBase=[];
        end
        clear tmpspec tmpspecavg specavg4src specavg4src_base;
        %==========================================
        % Now the sourceLocOut should contain the source Localisation results for the desired time points/ periods
        % and if baseline provided sourceLocBase should have the baseline power
        % ============================================================
        %% PLOT SOME FIGURES
        % If 3D plot on tweaked spm8's T1.nii
        % If 2D plot on the subject's surface in Head Space (Shouldn't we use the freesurfer template or another template of cortex?)
        %
        
        %----------------------------------
        % Need to select a time period to plot
        source2Plot=sourceLocOut;
        if ~isFixedLen1 % This case power has been computed from the entire data
            time2plot=[nan nan];
            source2Plot=sourceLocOut;
        else
            if isTimeWin{1},
                if Ntimeperiods>1 % If more than one time periods find the first one with time >0
                    indFirst=find(sourceLocOut.time>0,1,'first');
                    if isempty(indFirst), % if not higher than 0 then take the first
                        indFirst=1;
                    end
                    %source2Plot=sourceLocOut;
                    if strcmp(outtype,'emgcoh')
                        source2Plot.avg.pow=source2Plot.avg.emgcoh;
                    end
                    source2Plot.time=sourceLocOut.time(indFirst);
                    source2Plot.period=sourceLocOut.period(indFirst,:);
                    source2Plot.avg.pow=sourceLocOut.avg.pow(:,indFirst);
                    time2plot=source2Plot.period;
                else
                    time2plot=sourceLocOut.period;
                    source2Plot=sourceLocOut;
                    if strcmp(outtype,'emgcoh')
                        source2Plot.avg.pow=source2Plot.avg.emgcoh;
                    end
                    
                end
                
            else
                if ~isempty(timepoints{1})
                    % in Time Points. Assuming time has a part >0 relative to reference.
                    %   Plotting what happens in 350 msec
                    indFirst=find(sourceLocOut.time>0,1,'first');
                    if isempty(indFirst), % if not higher than 0 then take the first
                        indFirst=1;
                    end
                    timeStart=sourceLocOut.time(indFirst);
                    timeEnd=timeStart+0.350;
                    indEnd=find(sourceLocOut.time>timeEnd,1,'last');
                    
                    %source2Plot=sourceLocOut;
                    if strcmp(outtype,'emgcoh')
                        source2Plot.avg.pow=source2Plot.avg.emgcoh;
                    end
                    source2Plot.time=mean([timeStart timeEnd]);
                    source2Plot.period=[timeStart timeEnd];
                    source2Plot.avg.pow=nanmean(source2Plot.avg.pow(:,indFirst:indEnd),2);
                    time2plot=source2Plot.period;
                else
                    time2plot=[nan nan];
                    if strcmp(outtype,'emgcoh')
                        source2Plot.avg.pow=source2Plot.avg.emgcoh;
                    end
                end
            end
        end
        
        
        %---------------------------------------
        % Just set the baseline source
        source2PlotBase=sourceLocBase;
        if hasBaseline
            if strcmp(outtype,'emgcoh')
                source2PlotBase.avg.pow=source2PlotBase.avg.emgcoh;
            end
            [tmps1,tmps2]=size(source2PlotBase.avg.pow);
            if tmps2>tmps1,
                source2PlotBase.avg.pow=source2PlotBase.avg.pow';
            end
        end
        
        %---------------------------------------
        % Do the plot
        if strcmp(gridtype,'3D')
            source2Plot.pos=gridStand.pos;
            plotsource3D(source2Plot,mriTempl,invFiltFlag,experimentid, scanmnem,sortcontrlist{1}{iC}, curBandName, time2plot,source2PlotBase);
        elseif strcmp(gridtype,'2D')
            plotsource2D(source2Plot,invFiltFlag,experimentid, scanmnem,sortcontrlist{1}{iC},curBandName,time2plot,source2PlotBase);
        end
        
        % ============================================================
        %% SAVE results
        
        %--- Save the Source results
        source=rmfield(sourceLocOut,{'cfg','freq'}); % save some space
        %------------------
        if strcmp(gridtype,'2D')
            source.brainstructure=sourcemodel2d.brainstructure;
            source.brainstructurelabel=sourcemodel2d.brainstructurelabel; %{'CORTEX_LEFT' ;  'CORTEX_RIGHT'};
        end
        %------------------
        if isfield(source,'time')
            source.dimord='pos_time';
            sourcetype='dtseries';
        else
            source.dimord='pos';
            sourcetype='dscalar';
        end
        
        if strcmp(outtype,'emgcoh')
            source.emgcoh=source.avg.emgcoh;
            source.avg=rmfield(sourceLocOut.avg,'emgcoh'); % save some space
            outparameter='emgcoh';
            
        else
            %==========================================
            % Apply Weighting to the data so that the colorscale works in Workbench
            % The scaling is (10^21)^2 the square stands for power
            scaleWeight=(10^21)^2;
            source.power=scaleWeight.*source.avg.pow;
            %==========================================
            source.avg=rmfield(sourceLocOut.avg,'pow'); % save some space
            outparameter='power';
        end
        saveFnameData=[saveExtraDir,experimentid,'_',scanmnem,'_',sortcontrlist{1}{iC}.mnemprint,'_[FB-',curBandName,']'];
        
        
        %{
        % --- Save the baseline
        if ~isempty(sourceLocBase)
            tmpsource=rmfield(sourceLocBase,'cfg'); % save some space
            if isfield(tmpsource,'time')
                tmpsource.dimord='pos_time';
                tmpsourcetype='dtseries';
            else
                tmpsource.dimord='pos';
                tmpsourcetype='dscalar';
            end
            
            if strcmp(outtype,'emgcoh')
                tmpsource.baselineemgcoh=tmpsource.avg.emgcoh;
                tmpsource.avg=rmfield(sourceLocBase.avg,'emgcoh'); % save some space
                outparameter='baselineemgcoh';
                source.baselineemgcoh=tmpsource.baselineemgcoh;
                
            else
                tmpsource.baselinepower=tmpsource.avg.pow;
                tmpsource.avg=rmfield(sourceLocBase.avg,'pow'); % save some space
                outparameter='baselinepower';
                source.baselinepower=tmpsource.baselinepower;
            end
            
            source.baselinetime=sourceLocBase.time;
            source.baselineperiod=sourceLocBase.period;
            
            
            
            % save it as a cifti file
            hcp_write_cifti([saveFnameData '.' outparameter],tmpsource, 'parameter',outparameter, 'type', tmpsourcetype,'precision','double');
            
            clear tmpsource;
            
        end
        %}
        
        
        %Save in Mat
        %hcp_write_matlab(saveFnameData,'source');
        clear sourceLocOut source2Plot sourceLocOutBase source2PlotBase;
        % save it as a cifti file
        hcp_write_cifti([saveFnameData,'.',outparameter],source, 'parameter', outparameter, 'type', sourcetype,'precision','double');  clear source;
        
        
    end
    clear datacntrtrim1;
    
    % ============================================================
    %% VERY IMPORTANT - update previoius data group variable for next iteration
    
    prevGroups=sortgroupnames{iC};
    
    
    
end



end
%--- END OF MAIN --------------------------
%=====================================================
%=====================================================
%=====================================================
function[]=plotsource3D(source2Plot,mri2interp,invFiltType,experimentid, scanmnem,curContrast,curBandName, timeofplot,source2PlotBase)
global saveExtraDir;


source2Plot=ft_convert_units(source2Plot,'mm');
dumSrc=source2Plot;clear source2Plot;
%------------------------------------------------------------------
% If baseline provided then subtract from condition source power
if ~isempty(source2PlotBase)
    Ntimepoints=size(dumSrc.avg.pow,2);
    dumSrc.avg.pow=dumSrc.avg.pow-repmat(source2PlotBase.avg.pow,1,Ntimepoints);
    basetimeofplot=source2PlotBase.period;
else
    basetimeofplot=[nan nan];
end

%----------------------------------------------------------------
% Interpolate with MRI


intCfg=[];
intCfg.parameter = 'avg.pow';
[dumSrcInterp] = ft_sourceinterpolate(intCfg,dumSrc, mri2interp);

%----------------------------------------------------------------
% Define the plot color limits
inIndices=find(dumSrc.inside);
[maxabsVal,maxabsInd]=max(abs(dumSrc.avg.pow(inIndices)));
maxabsIndTot=inIndices(maxabsInd);

[minVal,minInd]=min((dumSrc.avg.pow(inIndices)));
[maxVal,maxInd]=max((dumSrc.avg.pow(inIndices)));
if sign(maxVal)~=sign(minVal)
    subclim=0.75*maxabsVal*[-1 1];
    %subAlim=0.05*maxabsVal*[-1 1];
else
    meanVal=mean([minVal maxVal]);
    demMinVal=minVal-meanVal;
    demMaxVal=maxVal-meanVal;
    subclim=meanVal+0.75.*[demMinVal demMaxVal];
    %subAlim=meanVal+0.05.*[demMinVal demMaxVal];
end

%subclim=0.75*maxabsVal*[-1 1];


%---------------------------------------------------------
% Plotting settings
plotCfg=[];
plotCfg.nslices       = 25;
plotCfg.method        = 'slice';
plotCfg.funcolormap   ='jet';
plotCfg.funcolorlim   = subclim;
%plotCfg.opacitylim   = subAlim;
%plotCfg.maskparameter  = 'avg.pow';
plotCfg.funparameter  = 'avg.pow';

ft_sourceplot(plotCfg,dumSrcInterp);
h1=gcf;

clear dumSrc dumSrcInterp;


set(h1,'papertype','A4');
set(h1,'paperunits','centimeters')
papersize=get(h1,'PaperSize');
paperposition = [1 1 papersize-1 ];
set(h1,'papersize',papersize);
set(h1,'paperposition',paperposition);
set(h1,'position',[10 10 15*papersize]);

fax(1)=gca;
set(fax(1),'position',[0.25 0.1 0.5 0.5]);

Toph1=axes();
set(Toph1,'position',[0.1 0.8 0.8 0.19]);axis off
dispStringTop1=sprintf('%s',[' experimentid:  ',regexprep(experimentid,'_','\\_')]);
dispStringTop1=sprintf('%s\n%s',dispStringTop1,[' scan:  ',regexprep(scanmnem,'_','\\_')]);
dispStringTop1=sprintf('%s\n%s',dispStringTop1,[' contrast:  ',regexprep(curContrast.mnemprint,'_','\\_')]);
dispStringTop1=sprintf('%s\n%s',dispStringTop1,[' Freq Band:  ',curBandName]);
%if invFiltType>0,
%    dispStringTop1=sprintf('%s\n%s',dispStringTop1,[' Filter type:  ',curContrast.invfiltertype]);
%    dispStringTop1=sprintf('%s\n%s',dispStringTop1,[' Contrast Operation:  ',curContrast.operation]);
%end
if ~isempty(timeofplot)
    dispStringTop1=sprintf('%s\n%s',dispStringTop1,[' Averaged Power in period:  ',num2str(timeofplot(1)),' to ',num2str(timeofplot(2))]);
else
    dispStringTop1=sprintf('%s\n%s',dispStringTop1,[' Averaged Power in period:  all ']);
end
%if invFiltType==0,
dispStringTop1=sprintf('%s\n%s',dispStringTop1,[' Averaged Baseline period:  ',num2str(basetimeofplot(1)),' to ',num2str(basetimeofplot(2))]);
dispStringTop1=sprintf('%s\n%s',dispStringTop1,[' Baseline operation:  ',num2str(curContrast.baselinetype)]);
%end
htTop1=text(0,0.2,dispStringTop1);
set(htTop1,'Fontsize',10)

saveFnameImage=[saveExtraDir,experimentid,'_',scanmnem,'_',curContrast.mnemprint,'_[FB-',curBandName,']_plot.png'];
disp('Saving the figure');
hcp_write_figure(saveFnameImage, h1);
disp('Done Saving');
close(h1);

end

%=====================================================
%=====================================================
%=====================================================
%=====================================================
function[]=plotsource2D(source2Plot,invFiltType,experimentid, scanmnem,curContrast,curBandName,timeofplot,source2PlotBase)
global saveExtraDir;

dumSrc=source2Plot;clear source2Plot;
%------------------------------------------------------------------
% If baseline provided then subtract from condition source power
if ~isempty(source2PlotBase)
    Ntimepoints=size(dumSrc.avg.pow,2);
    dumSrc.avg.pow=dumSrc.avg.pow-repmat(source2PlotBase.avg.pow,1,Ntimepoints);
    basetimeofplot=source2PlotBase.period;
else
    basetimeofplot=[nan nan];
end


%--------------------------------------------------------
% Define plot color limits
inIndices=find(dumSrc.inside);
[maxabsVal,maxabsInd]=max(abs(dumSrc.avg.pow(inIndices)));
maxabsIndTot=inIndices(maxabsInd);

[minVal,minInd]=min((dumSrc.avg.pow(inIndices)));
[maxVal,maxInd]=max((dumSrc.avg.pow(inIndices)));
if sign(maxVal)~=sign(minVal)
    subclim=0.75*maxabsVal*[-1 1];
else
    meanVal=mean([minVal maxVal]);
    demMinVal=minVal-meanVal;
    demMaxVal=maxVal-meanVal;
    subclim=meanVal+0.75.*[demMinVal demMaxVal];
end

%subclim=0.75*maxabsVal*[-1 1];


cfg = [];
cfg.method         = 'surface';
cfg.funparameter   = 'avg.pow';
cfg.funcolorlim=subclim;
ft_sourceplot(cfg,dumSrc);
figsurf=gcf;
axsurf=gca;

clear dumSrc;

h1=figure;
set(h1,'papertype','A4');
set(h1,'paperunits','centimeters')
papersize=get(h1,'PaperSize');
paperposition = [1 1 papersize-1 ];
set(h1,'papersize',papersize);
set(h1,'paperposition',paperposition);
set(h1,'position',[10 10 15*papersize]);

fxa1=axes('Position',[0.1 0.3 0.4 0.3]);axis off;
fxa2=axes('Position',[0.5 0.3 0.4 0.3]);axis off;
fxa3=axes('Position',[0.1 0 0.4 0.3]);axis off;
fxa4=axes('Position',[0.5 0 0.4 0.3]);axis off;


axes(fxa1);
curchild=get(axsurf,'children');
copyobj(curchild([3 2 1]),fxa1);

view(-90,0);
set(gca,'CLim',subclim);set(gca,'ALimMode','manual');set(gca,'ALim',[0 1]);
hcol1=colorbar('West');
set(hcol1,'Position',[[0.05 0.3 0.05 0.28]]);

axes(fxa2);
copyobj(curchild(1:2),fxa2);
view(-90,90);
set(gca,'CLim',subclim);set(gca,'ALimMode','manual');set(gca,'ALim',[0 1]);

axes(fxa3);
curchild=get(axsurf,'children');
copyobj(curchild(1:2),fxa3);
view(180,0);
set(gca,'CLim',subclim);set(gca,'ALimMode','manual');set(gca,'ALim',[0 1]);

axes(fxa4);
copyobj(curchild(1:2),fxa4);
view(0,0);
set(gca,'CLim',subclim);set(gca,'ALimMode','manual');set(gca,'ALim',[0 1]);

close(figsurf);


%===============================================================================
Toph1=axes();
set(Toph1,'position',[0.1 0.8 0.8 0.19]);axis off
dispStringTop1=sprintf('%s',[' experimentid:  ',regexprep(experimentid,'_','\\_')]);
dispStringTop1=sprintf('%s\n%s',dispStringTop1,[' scan:  ',regexprep(scanmnem,'_','\\_')]);
dispStringTop1=sprintf('%s\n%s',dispStringTop1,[' contrast:  ',regexprep(curContrast.mnemprint,'_','\\_')]);
dispStringTop1=sprintf('%s\n%s',dispStringTop1,[' Freq Band:  ',curBandName]);
%if invFiltType>0,
%    dispStringTop1=sprintf('%s\n%s',dispStringTop1,[' Filter type:  ',curContrast.invfiltertype]);
%    dispStringTop1=sprintf('%s\n%s',dispStringTop1,[' Contrast Operation:  ',curContrast.operation]);
%end
if ~isempty(timeofplot)
    dispStringTop1=sprintf('%s\n%s',dispStringTop1,[' Averaged Power in period:  ',num2str(timeofplot(1)),' to ',num2str(timeofplot(2))]);
else
    dispStringTop1=sprintf('%s\n%s',dispStringTop1,[' Averaged Power in period:  all ']);
end
%if invFiltType==0,
dispStringTop1=sprintf('%s\n%s',dispStringTop1,[' Averaged Baseline period:  ',num2str(basetimeofplot(1)),' to ',num2str(basetimeofplot(2))]);
dispStringTop1=sprintf('%s\n%s',dispStringTop1,[' Baseline operation:  ',num2str(curContrast.baselinetype)]);

%end
htTop1=text(0,0.2,dispStringTop1);
set(htTop1,'Fontsize',10)

saveFnameImage=[saveExtraDir,experimentid,'_',scanmnem,'_',curContrast.mnemprint,'_[FB-',curBandName,']_plot.png'];
disp('Saving the figure');
hcp_write_figure(saveFnameImage, h1);
disp('Done Saving');
close(h1);

end

