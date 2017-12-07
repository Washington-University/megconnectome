function[outStatus]=hcp_srcavglcmv_contrasts(inCfg)

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
global saveExtraDir; % This variable is externally set by the user or by default to the pipelinedir and contains the directory where the pipeline results will be saved

%================================================
% Get inputs
if isfield(inCfg,'savedir')
    saveExtraDir=inCfg.savedir;
else
    saveExtraDir=[];
end

subjectid=inCfg.subjectid;
experimentid=inCfg.experimentid; % Experiment ID
multiscanid=inCfg.multiscanid;   % Cell with the scan id of both files. They shoudl be in scannning order
contrastlist=inCfg.contrastlist; % Cell with list of contrasts
anatomydir = inCfg.anatomydir;   % Directory where anatomy file are source space models are located.
gridtype=inCfg.gridtype;         % Can be '3D'(using 8mm) or '2D'%

%===================================================================
% Set basic filtering settings for sensor level preprocessing

% -----------------
% Settings for filtering when localization is gonna be based on the avg of
% trials - ERFs are clearer when data is filtered above 40 Hz
flowcfg_avg=[];
flowcfg_avg.hpfilter='yes';
flowcfg_avg.hpfreq=[1];
flowcfg_avg.hpfiltord=4;

fhighcfg_avg=[];
fhighcfg_avg.lpfilter='yes';
fhighcfg_avg.lpfreq=[40];

% -----------------
% Settings for filtering when localization is gonna be based on all the trial data
% Keeping Frequencies up to 100 Hz
flowcfg_all=[];
flowcfg_all.hpfilter='yes';
flowcfg_all.hpfreq=[0.5];
flowcfg_all.hpfiltord=4;

fhighcfg_all=[];
fhighcfg_all.lpfilter='yes';
fhighcfg_all.lpfreq=[100];
%=====================================================================
%================================================
% Get basic information about input contrasts
Nfiles=length(contrastlist);     % Number of files to be read for a given task. Should be 2 for each Task Type
Ncontr=length(contrastlist{1});  % Check how many contrasts from the first file's contrasts

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
    
    
    %=============================================================
    % Load the time data for the corresponding group and from both
    % files
    
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
            fileGrid2D=[anatomydir,'/',experimentid,'_anatomy_sourcemodel_2d.mat'];     % Subject's 2D cortical sheet source model
            fileGrid3D=[anatomydir,'/',experimentid,'_anatomy_sourcemodel_3d8mm.mat'];  % Subject's 3D volumetric source model with 8mm resolution
            fileGrid3DStandard=['standard_sourcemodel3d8mm'];                           % template 3D grid
            %fileAnat=[anatomydir,'/',experimentid,'_anatomy_anatomical.nii'];           % Subject's MRI
            fileAnatTemplMcGill='mni_icbm152_t1_tal_nlin_sym_09a.nii';                  % Template T1 MRI from McGill
            fileVol=[anatomydir,'/',experimentid,'_anatomy_headmodel.mat'];             % Subject's brain volume
            
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
            cfg.normalize       = 'yes' ; %Normalise Leadfield: 'yes' for beamformer
            cfg.normalizeparam  = 1;      %depth normalization parameter (default = 0.5).
            cfg.feedback='no';
            cfg.channel=allChansMEG;
            if strcmp(gridtype,'2D') %JUST TO MAKE SURE ENTIRE CORTICAL SHEET IN SOURCE SPACE
                cfg.inwardshift=-1;
            end
            gridAllLF= ft_prepare_leadfield(cfg);
            gridAllLF.label=allChansMEG;
            
            
            
        end
        clear data;
        
        
        %----- As the data was just loaded apply filtering at sensor level data
        %alldata1A=ft_preprocessing(flowcfg,alldata1A);
        %alldata1A=ft_preprocessing(fhighcfg,alldata1A);
        
        
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
            %----- As the data was just loaded apply filtering at sensor level data
            %alldata1B=ft_preprocessing(flowcfg,alldata1B);
            %alldata1B=ft_preprocessing(fhighcfg,alldata1B);
        else
            alldata1B=[];
        end
        
        
        
    else % if the datagroup of the current contrast is the same as the previous, then just keep using te same data.
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
    
    %=========================================================
    %-------------------------------------------------------
    % Apply filters according to if the source localization will be based
    % on the ERF or on the entire data
    if invFiltFlag==0,  % avg
        flowcfg=flowcfg_avg;
        fhighcfg=fhighcfg_avg;
    elseif invFiltFlag==1, % all
        flowcfg=flowcfg_all;
        fhighcfg=fhighcfg_all;
    end
    datacntr1=ft_preprocessing(flowcfg,datacntr1);
    datacntr1=ft_preprocessing(fhighcfg,datacntr1);
    
    
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
    
    
    %===================================================================
    % Select data portion that spans from the earliest to the latest
    % timepoints or timeperiod boundaries asked( and baseline boundaries if provided).
    % Here is is assumed that timepoints and time periods are in ascending
    % order.
    
    Ntrials1=length(datacntr1.trial);
    datacntrtrim1=[];
    
    iCond=1;
    eval(['tmpdat=datacntr',num2str(iCond),';']);
    eval(['tmpNtrials=Ntrials',num2str(iCond),';']);
    
    earliestTimeBase=nan;
    latestTimeBase=nan;
    if hasBaseline
        earliestTimeBase=baseline{iCond}(1,1);
        latestTimeBase=baseline{iCond}(1,2);
    end
    if isTimeWin{iCond}==1,
        
        earliestTimeAct=timeperiods{iCond}(1,1);
        latestTimeAct=timeperiods{iCond}(Ntimeperiods,2);
        
        earliestTime=min([earliestTimeBase earliestTimeAct latestTimeBase latestTimeAct]);
        latestTime=max([earliestTimeBase earliestTimeAct latestTimeBase latestTimeAct]);
        
        tmpdattrim=ft_selectdata(tmpdat,'toilim',[earliestTime latestTime]);
        
    else
        if ~isempty(timepoints{iCond})
            earliestTimeAct=timepoints{iCond}(1);
            latestTimeAct=timepoints{iCond}(Ntimepoints);
            
            earliestTime=min([earliestTimeBase earliestTimeAct latestTimeBase latestTimeAct]);
            latestTime=max([earliestTimeBase earliestTimeAct latestTimeBase latestTimeAct]);
            
            tmpdattrim=ft_selectdata(tmpdat,'toilim',[earliestTime latestTime]);
            
        else
            tmpdattrim=tmpdat;
        end
        
    end
    eval(['datacntrtrim',num2str(iCond),'=tmpdattrim;']); clear tmpdat; clear tmpdattrim;
    clear datacntr1;
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
    % Do a quick check if the data has variable length and filter type is
    % avg . This aint possible
    if (~isFixedLen1)&(invFiltFlag==0)
        error('Variable length trials but inverse type filter is avg. This cannot be processed');
    end
    %=========================================================
    
    
    
    
    %=========================================================
    %=========================================================
    % Compute total covariance Matrix
    % If invFiltFlag =0 (avg) then compute covariance from averaged trials
    % If invFiltFlag =1 (all) then compute covariance from NON-averaged trials
    datap=datacntrtrim1;
    
    %--------------------------------------------------------
    %First check if there are any nans in the data - This check is mostly thought for the Story Math task where trials span entire blocks with variable length and bad segmetns have been replaced with nans
    tmpdat=datap;
    for iTr=1:Ntrials1,
        [tmpIndx1,tmpIndx2]=find(isnan(tmpdat.trial{iTr}));
        if ~isempty(tmpIndx2)
            tmpdat.trial{iTr}(:,unique(tmpIndx2))=[];% Remove nans so they do not affect the computation
            tmpdat.time{iTr}(unique(tmpIndx2))=[];
        end
    end
    
    tmptotdat=[tmpdat.trial{:}];
    if isempty(tmptotdat)
        error('All your data is nan - Check in the early steps');
    end
    clear tmptotdat datap;
    %--------------------------------------------------------
    % Then demean and apply baseline if available
    cfg=[];
    cfg.demean='yes';
    %cfg.channel='MEG';
    if hasBaseline,
        cfg.baselinewindow=baseline{1};
    end
    tmpdat=ft_preprocessing(cfg,tmpdat);
    %--------------------------------------------------------
    
    %--------------------------------------------------------
    % Compute covariance
    if invFiltFlag==0, % filter to be computed from trial avg covariance
        
        cfg=[];
        cfg.covariance='no';
        cfg.keeptrials='no';
        cfg.removemean='no';
        cfg.vartrllength=2;
        datavg1=ft_timelockanalysis(cfg,tmpdat);
        
        cfg=[];
        cfg.covariance='yes';
        datavg1=ft_timelockanalysis(cfg,datavg1);
        datavgall=datavg1;
        
    elseif invFiltFlag==1, % filter to be computed from all trial covariance
        cfg=[];
        cfg.covariance='yes';
        cfg.keeptrials='no';
        cfg.removemean='no';
        cfg.vartrllength=2;
        datavg1=ft_timelockanalysis(cfg,tmpdat);
        datavgall=datavg1;
    end
    clear datavg1;
    
    if ~isFixedLen1
        tmpconcatdat=[tmpdat.trial{:}];
        [ind1,ind2]=find(isnan(tmpconcatdat));
        ind2un=unique(ind2);
        tmpconcatdat(:,ind2)=[];
        tmpNconcat=size(tmpconcatdat,2);
        tmpCov=(tmpconcatdat*transpose(tmpconcatdat))./tmpNconcat;
        datavgall.cov=tmpCov;
        clear tmpconcatdat ind1 ind2 ind2un tmpNconcat tmpCov;
    end
    
    
    %--------------------------------------------------------
    
    %=========================================================
    %=========================================================
    %=========================================================
    % Now Compute Inverse Solution
    
    
    %----- Set initial source localization settings
    srcCfg=[];
    srcCfg.vol=vol;
    srcCfg.method='lcmv';
    srcCfg.fixedori='yes';
    srcCfg.feedback='text';
    srcCfg.lambda='100%';
    srcCfg.projectnoise='no';
    srcCfg.keepfilter='yes';
    srcCfg.keepleadfield='yes';
    srcCfg.keepmom='no';
    
    
    %--- Select only channels of current data set in already computed leadfields -----------
    [indA,indB]=match_str(datavgall.label,gridAllLF.label);
    Nsrc=length(gridAllLF.inside);
    gridCase=gridAllLF;
    gridCase.leadfield(gridCase.inside)=cellfun(@(x,y) x(y,:),gridAllLF.leadfield(gridAllLF.inside),repmat({indB},1,Nsrc),'UniformOutput',false);
    
    %----- Add grid to source localization settings
    srcCfg.grid=gridCase;
    
    %----- Compute Inverse Solution from Covariance matrix ------------
    
    %  this provides the spatial filters that will project the data for the
    %  time points and time periods of interest
    datavgall.grad=grad;
    disp('beamforming starting...');
    datavgall.grad=grad;
    sourceLocAll=ft_sourceanalysis(srcCfg,datavgall);
    disp('...beamforming ended');
    
    
    %===================================================================
    %===================================================================
    % Here the Inverse Solution has been computed and data can now be
    % projected through the Spatial Filters
    
    NsrcTotal=size(sourceLocAll.pos,1);
    NsrcIn=length(sourceLocAll.inside);
    if (~isFixedLen1) % if data has no fixed length then only the total power from the entire data will be provided
        sourceLocOut=rmfield(sourceLocAll,'time');
        sourceLocOut.period='all';
    else
        
        sourceLocOut=sourceLocAll;
        if isTimeWin{1}  % If output is in time periods compute covariance for each and project
            sourceLocOut.avg.pow=nan(NsrcTotal,Ntimeperiods);
            for iPeriod=1:Ntimeperiods
                
                cfg=[];
                cfg.covariance='yes';
                cfg.covariancewindow=timeperiods{1}(iPeriod,:);
                cfg.keeptrials='no';
                cfg.removemean='no';
                cfg.vartrllength=2;
                if invFiltFlag==0,
                    tmpDatAvg=ft_timelockanalysis(cfg,datavgall);
                else
                    tmpDatAvg=ft_timelockanalysis(cfg,tmpdat);
                end
                srcCfg.keepmom='no';
                srcCfg.filter=sourceLocAll.avg.filter;
                tmpSourceLoc=ft_sourceanalysis(srcCfg,datavgall);
                sourceLocOut.avg.pow(sourceLocOut.inside,iPeriod)=tmpSourceLoc.avg.pow(tmpSourceLoc.inside);
                sourceLocOut.time=mean(timeperiods{1},2);
                sourceLocOut.period=timeperiods{1};
                clear tmpSourceLoc tmpDatAvg;
            end
            
        else % if output in time points then project all initial data and
            % extract time points while also smoothing with a window equal
            % to the average inter-point distance
            if ~isempty(timepoints{1})
                
                srcCfg.keepmom='yes';
                srcCfg.filter=sourceLocAll.avg.filter;
                tmpSourceLoc=ft_sourceanalysis(srcCfg,datavgall);
                Ndattimes=length(datavgall.time);
                tmpPowAll=nan(NsrcIn,Ndattimes);
                tmpPowAll(tmpSourceLoc.inside,:)=reshape([tmpSourceLoc.avg.mom{tmpSourceLoc.inside}].^2,Ndattimes,NsrcIn)';
                
                
                
                dattimes=datavgall.time;
                smoothHalfTime=mean(diff(timepoints{iCond})./2);
                smoothHalfWin=floor(smoothHalfTime*Fsample); % TODO: When songle points are selected a smoothing window of 30msec is applies around them to avoid outlier effects. The length of this window is hard coded. This can be made dynamic. FIX!!!
                tmpIndMat=[];
                for iTime=1:Ntimepoints,
                    tmpIndMat(iTime)=nearest(dattimes,timepoints{iCond}(iTime));
                end
                tmpIndMat=unique(tmpIndMat);
                NtimesOut=length(tmpIndMat);
                
                % Downsample to get the desired points
                sourceLocOut.time=tmpSourceLoc.time(tmpIndMat);
                tmpPowDown=nan(NsrcIn,NtimesOut);
                for iTime=1:NtimesOut,
                    tmpPowDown(:,iTime)=mean(tmpPowAll(:,(max(1,(tmpIndMat(iTime)-smoothHalfWin)):min(size(tmpPowAll,2),(tmpIndMat(iTime)+smoothHalfWin)))),2);
                end
                sourceLocOut.avg.pow=nan(NsrcTotal,NtimesOut);
                sourceLocOut.avg.pow(sourceLocOut.inside,:)=tmpPowDown;
                clear tmpPowDown tmpPowAll dattimes;
            else
                sourceLocOut=rmfield(sourceLocOut,'time');
            end
            
            
        end
        
        
        
    end
    
    %%
    %------------------------------------------------------------------------
    % If Baseline has been provided then compute the power in the baseline
    % period
    if hasBaseline
        sourceLocBase=sourceLocAll;
        sourceLocBase.avg.pow=nan(NsrcTotal,1);
        sourceLocBase.time=mean(baseline{1});
        sourceLocBase.period=baseline{1};
        cfg=[];
        cfg.covariance='yes';
        cfg.covariancewindow=baseline{1};
        cfg.keeptrials='no';
        cfg.removemean='no';
        cfg.vartrllength=2;
        if invFiltFlag==0,
            tmpDatAvg=ft_timelockanalysis(cfg,datavgall);
        else
            tmpDatAvg=ft_timelockanalysis(cfg,tmpdat);
        end
        srcCfg.keepmom='no';
        tmpSourceLoc=ft_sourceanalysis(srcCfg,datavgall);
        sourceLocBase.avg.pow(sourceLocOut.inside)=tmpSourceLoc.avg.pow(tmpSourceLoc.inside);
        clear tmpSourceLoc tmpDatAvg;
    else
        sourceLocBase=[];
    end
    %----------------------------------
    % Do some cleaning after source localization
    clear sourceLocAll;
    clear tmpdat;
    clear datavgall;
    %------------------------------------------------------------------------
    
    
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
                source2Plot=sourceLocOut;
                source2Plot.time=sourceLocOut.time(indFirst);
                source2Plot.period=sourceLocOut.period(indFirst,:);
                source2Plot.avg.pow=sourceLocOut.avg.pow(:,indFirst);
                time2plot=source2Plot.period;
            else
                time2plot=sourceLocOut.period;
                source2Plot=sourceLocOut;
            end
            
        else % in Time Points. Assuming time has a part >0 relative to reference.
            %   Plotting what happens in 350 msec
            if ~isempty(timepoints{1})
                indFirst=find(sourceLocOut.time>0,1,'first');
                if isempty(indFirst), % if not higher than 0 then take the first
                    indFirst=1;
                end
                timeStart=sourceLocOut.time(indFirst);
                timeEnd=timeStart+0.350;
                indEnd=find(sourceLocOut.time>timeEnd,1,'last');
                
                source2Plot=sourceLocOut;
                source2Plot.time=mean([timeStart timeEnd]);
                source2Plot.period=[timeStart timeEnd];
                source2Plot.avg.pow=mean(sourceLocOut.avg.pow(:,indFirst:indEnd),2);
                time2plot=source2Plot.period;
            else
                source2Plot=sourceLocOut;
                time2plot=[nan nan];
            end
            
        end
    end
    
    %---------------------------------------
    % Just set the baseline source
    source2PlotBase=sourceLocBase;
    %---------------------------------------
    % Do the plot
    if strcmp(gridtype,'3D')
        source2Plot.pos=gridStand.pos;
        plotsource3D(source2Plot,mriTempl,invFiltFlag,experimentid, scanmnem,sortcontrlist{1}{iC},time2plot,source2PlotBase);
    elseif strcmp(gridtype,'2D')
        plotsource2D(source2Plot,invFiltFlag,experimentid, scanmnem,sortcontrlist{1}{iC},time2plot,source2PlotBase);
    end
    
    
    % ============================================================
    %% SAVE results
    
    %--- Save the Source results
    source=rmfield(sourceLocOut,'cfg'); % save some space
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
    
    
    %==========================================
    % Apply Weighting to the data so that the colorscale works in Workbench
    % The scaling is (10^21)^2 the square stands for power
    scaleWeight=(10^21)^2;
    source.power=scaleWeight.*source.avg.pow;
    %==========================================
    
    source.avg=rmfield(sourceLocOut.avg,'pow'); % save some space
    saveFnameData=[saveExtraDir,experimentid,'_',scanmnem,'_',sortcontrlist{1}{iC}.mnemprint];
    
    %{
    % --- Save the baseline
    if ~isempty(sourceLocBase)
        tmpsource=rmfield(sourceLocBase,'cfg'); % save some space
        tmpsource.baselinepower=sourceLocBase.avg.pow;
        tmpsource.avg=rmfield(sourceLocBase.avg,'pow'); % save some space
        source.baselinepower=sourceLocBase.avg.pow;
        
        if isfield(tmpsource,'time')
            tmpsource.dimord='pos_time';
            source.baselinetime=sourceLocBase.time;
            source.baselineperiod=sourceLocBase.period;
            tmpsourcetype='dtseries';
        else
            tmpsource.dimord='pos';
            tmpsourcetype='dscalar';
        end
        
        % save it as a cifti file
        hcp_write_cifti([saveFnameData 'baselinepower'],tmpsource, 'parameter','baselinepower', 'type', tmpsourcetype,'precision','double');
        
        clear tmpsource;
        
    end
    %}
    %save as mat file
    hcp_write_matlab(saveFnameData,'source'); clear sourceLocOut source2Plot sourceLocBase source2PlotBase;
    % save it as a cifti file
    hcp_write_cifti([saveFnameData '.power'], source, 'parameter', 'power', 'type', sourcetype,'precision','double'); clear source;
    
    
    
    % ============================================================
    %% VERY IMPORTANT - update previoius data group variable for next iteration
    
    prevGroups=sortgroupnames{iC};
    
    
    % ============================================================
    % end of processing for current contrast
    
    
end % End of loop over contrasts



end % End of main function


%--- END OF MAIN --------------------------
%=====================================================
%=====================================================
%=====================================================
%=====================================================
%=====================================================
%=====================================================
%=====================================================
%=====================================================
%=====================================================
function[]=plotsource3D(source2Plot,mri2interp,invFiltType,experimentid, scanmnem,curContrast,timeofplot,source2PlotBase)
global saveExtraDir;

source2Plot=ft_convert_units(source2Plot,'mm');
dumSrc=source2Plot; clear source2Plot;
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
[maxabsVal,maxabsInd]=max(abs(dumSrc.avg.pow(dumSrc.inside)));
maxabsIndTot=dumSrc.inside(maxabsInd);

[minVal,minInd]=min((dumSrc.avg.pow(dumSrc.inside)));
[maxVal,maxInd]=max((dumSrc.avg.pow(dumSrc.inside)));
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


%---------------------------------------------------------
% Plotting settings
plotCfg=[];
plotCfg.nslices       = 25;
plotCfg.method        = 'slice';
plotCfg.funcolormap   ='jet';
plotCfg.funcolorlim   = subclim;
plotCfg.funparameter  = 'avg.pow';

ft_sourceplot(plotCfg,dumSrcInterp);
h1=gcf;

clear dumSrc dumSrcInterp;
%---- Put figure in correct format and add text
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
if invFiltType>0,
    dispStringTop1=sprintf('%s\n%s',dispStringTop1,[' Filter type:  ',curContrast.invfiltertype]);
    dispStringTop1=sprintf('%s\n%s',dispStringTop1,[' Contrast Operation:  ',curContrast.operation]);
end
if ~isempty(timeofplot)
    dispStringTop1=sprintf('%s\n%s',dispStringTop1,[' Averaged Power in period:  ',num2str(timeofplot(1)),' to ',num2str(timeofplot(2))]);
else
    dispStringTop1=sprintf('%s\n%s',dispStringTop1,[' Averaged Power in period:  all ']);
end
if invFiltType==0,
    dispStringTop1=sprintf('%s\n%s',dispStringTop1,[' Averaged Baseline period:  ',num2str(basetimeofplot(1)),' to ',num2str(basetimeofplot(2))]);
    dispStringTop1=sprintf('%s\n%s',dispStringTop1,[' Baseline operation:  ',num2str(curContrast.baselinetype)]);
end
htTop1=text(0,0.2,dispStringTop1);
set(htTop1,'Fontsize',10)

%------ Save Figure --------------
saveFnameImage=[saveExtraDir,experimentid,'_',scanmnem,'_',curContrast.mnemprint,'_plot.png'];
disp('Saving the figure');
hcp_write_figure(saveFnameImage, h1);
disp('Done Saving');
close(h1);

end

%=====================================================
%=====================================================
%=====================================================
%=====================================================
function[]=plotsource2D(source2Plot,invFiltType,experimentid, scanmnem,curContrast,timeofplot,source2PlotBase)
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
[maxabsVal,maxabsInd]=max(abs(dumSrc.avg.pow(dumSrc.inside)));
maxabsIndTot=dumSrc.inside(maxabsInd);

[minVal,minInd]=min((dumSrc.avg.pow(dumSrc.inside)));
[maxVal,maxInd]=max((dumSrc.avg.pow(dumSrc.inside)));
if sign(maxVal)~=sign(minVal)
    subclim=0.75*maxabsVal*[-1 1];
else
    meanVal=mean([minVal maxVal]);
    demMinVal=minVal-meanVal;
    demMaxVal=maxVal-meanVal;
    subclim=meanVal+0.75.*[demMinVal demMaxVal];
end


%------------------------------------
% Plot figure;
cfg = [];
cfg.method         = 'surface';
cfg.funparameter   = 'avg.pow';
cfg.funcolorlim=subclim;
ft_sourceplot(cfg,dumSrc);
figsurf=gcf;
axsurf=gca;

clear dumSrc;
%------------------------------
% Set figure in correct format
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


%-----------------------------------------------------------
Toph1=axes();
set(Toph1,'position',[0.1 0.8 0.8 0.19]);axis off
dispStringTop1=sprintf('%s',[' experimentid:  ',regexprep(experimentid,'_','\\_')]);
dispStringTop1=sprintf('%s\n%s',dispStringTop1,[' scan:  ',regexprep(scanmnem,'_','\\_')]);
dispStringTop1=sprintf('%s\n%s',dispStringTop1,[' contrast:  ',regexprep(curContrast.mnemprint,'_','\\_')]);
if invFiltType>0,
    dispStringTop1=sprintf('%s\n%s',dispStringTop1,[' Filter type:  ',curContrast.invfiltertype]);
    dispStringTop1=sprintf('%s\n%s',dispStringTop1,[' Contrast Operation:  ',curContrast.operation]);
end
if ~isempty(timeofplot)
    dispStringTop1=sprintf('%s\n%s',dispStringTop1,[' Averaged Power in period:  ',num2str(timeofplot(1)),' to ',num2str(timeofplot(2))]);
else
    dispStringTop1=sprintf('%s\n%s',dispStringTop1,[' Averaged Power in period:  all ']);
end
if invFiltType==0,
    dispStringTop1=sprintf('%s\n%s',dispStringTop1,[' Averaged Baseline period:  ',num2str(basetimeofplot(1)),' to ',num2str(basetimeofplot(2))]);
    dispStringTop1=sprintf('%s\n%s',dispStringTop1,[' Baseline operation:  ',num2str(curContrast.baselinetype)]);
end
htTop1=text(0,0.2,dispStringTop1);
set(htTop1,'Fontsize',10)



%-----------------------------
%Save Figure

saveFnameImage=[saveExtraDir,experimentid,'_',scanmnem,'_',curContrast.mnemprint,'_plot.png'];
disp('Saving the figure');
hcp_write_figure(saveFnameImage, h1);
disp('Done Saving');
close(h1);

end

