function [outStatus] = hcp_tfavg_contrastsonthefly(inCfg)

% This function loads time frequency data and computes its trial average as well
% as the average of its planar gradient
%
% outdatafile is the file where the averaged data will be saved
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

%================================================
% Get inputs
subjectid=inCfg.subjectid;
experimentid=inCfg.experimentid;
multiscanid=inCfg.multiscanid; % cell with the scan id of  both files. They shoudl be in scannning order
bandinfo=inCfg.bandinfo;
contrastlist=inCfg.contrastlist;
avgmode=inCfg.avgmode;

bands=bandinfo(:,1);
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
% Now sort so that in each group if there is emgcoh all the emgcoh
% contrasts will be grouped together.
sortgroupsFirst=[];
for iGC=1:Ncontr,
    sortgroupsFirst{iGC}=sortgroupnames{iGC}{1};
end
newsortcontrnames=sortcontrnames;
newsortgroupnames=sortgroupnames;
newsortcontrlist=cell(1,2);
for iGr=1:Ngroups,
    [ind1,ind2]=match_str(sortgroupsFirst,datagroups{iGr});
    curList1=sortcontrlist{1}(ind1);
    if is2Files,
        curList2=sortcontrlist{2}(ind1);
    end
    newcurList1=[];
    nonemgList1=[];
    newcurList2=[];
    nonemgList2=[];
    countNew=1;
    for indGC=1:length(curList1)
        if strcmp(curList1{indGC}.connemetric,'emgcoh')
            newcurList1{countNew}=curList1{indGC};
            if is2Files,
                newcurList2{countNew}=curList2{indGC};
            end
            countNew=countNew+1;
        else
            nonemgList1{end+1}=curList1{indGC};
            if is2Files,
                nonemgList2{end+1}=curList2{indGC};
            end
        end
    end
    newsortcontrlist{1}=[newsortcontrlist{1} newcurList1 nonemgList1];
    if is2Files,
        newsortcontrlist{2}=[newsortcontrlist{2} newcurList2 nonemgList2];
    end
end

newsortcontrnames=[];
newsortgroupnames=[];
for iC=1:Ncontr,
    newsortcontrnames{iC}=newsortcontrlist{1}{iC}.mnemprint;
    newsortgroupnames{iC}=newsortcontrlist{1}{iC}.lockmode;
    newisCompCntr(iC,1)=length(newsortgroupnames{iC})-1;
    newisDifGroups(iC,1)=length(unique(newsortgroupnames{iC}))-1;
end
sortcontrlist=newsortcontrlist;
sortcontrnames=newsortcontrnames;
sortgroupnames=newsortgroupnames;
isCompCntr=newisCompCntr;
isDifGroups=newisDifGroups;
%=============================================================
%=============================================================
%% LOOP OVER CONTRASTS
prevGroups={'' ''};
curCM=[];
prevCM=[];
for iC=1:Ncontr
    %curcont=
    
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
    
    curCM=sortcontrlist{1}{iC}.connemetric;
    
    saveFnameData=[experimentid,'_',scanmnem,'_',sortcontrlist{1}{iC}.mnemprint,'_[MODE-',avgmode,']'];
    saveFnameImage=[experimentid,'_',scanmnem,'_',sortcontrlist{1}{iC}.mnemprint,'_[MODE-',avgmode,']_plot.png'];
    %==================================
    % Load the time frequency data for the corresponding group and from both
    % files
    
    NcurGroups=length(sortgroupnames{iC});
    prevMatchIndx=cell(1,NcurGroups);
    for iGr1=1:NcurGroups,
        tmpIndx=find(strcmp(sortgroupnames{iC},prevGroups));
        if ~isempty(tmpIndx)
            prevMatchIndx{iGr1}=tmpIndx(1);
        else
            prevMatchIndx{iGr1}=[];
        end
    end
    if ~isempty(prevMatchIndx{1})
        eval(['tmpdata1A=alldata',num2str(prevMatchIndx{1}),'A;']);
        eval(['tmpdata1B=alldata',num2str(prevMatchIndx{1}),'B;']);
    else
        tmpdata1A=[];
        tmpdata1B=[];
    end
    if NcurGroups>1
        if ~isempty(prevMatchIndx{2})
            eval(['tmpdata2A=alldata',num2str(prevMatchIndx{2}),'A;']);
            eval(['tmpdata2B=alldata',num2str(prevMatchIndx{2}),'B;']);
        else
            tmpdata2A=[];
            tmpdata2B=[];
        end
    end
    
    
    if isempty(tmpdata1A)
        
        genprefix = sprintf('%s_%s', experimentid, multiscanid{1}); tmpdatafile=[genprefix,'_tmegpreproc_',sortgroupnames{iC}{1}];
        load(tmpdatafile,'data');
        grad=data.grad;
        %=====================
        if data.fsample>2000
            newFs=data.fsample/4;
        elseif data.fsample>1000
            newFs=data.fsample/2;
        else
            newFs=data.fsample;
        end
        cfg = [];
        cfg.detrend = 'no';
        cfg.resamplefs = newFs;
        data = ft_resampledata(cfg, data); % clear data;
        
        %=========================
        alldata1A=data; clear data;
        
        if is2Files
            genprefix = sprintf('%s_%s', experimentid, multiscanid{2}); tmpdatafile=[genprefix,'_tmegpreproc_',sortgroupnames{iC}{1}];
            load(tmpdatafile,'data');
            %=====================
            if data.fsample>2000
                newFs=data.fsample/4;
            elseif data.fsample>1000
                newFs=data.fsample/2;
            else
                newFs=data.fsample;
            end
            cfg = [];
            cfg.detrend = 'no';
            cfg.resamplefs = newFs;
            data = ft_resampledata(cfg, data); % clear data;
            
            %=========================
            alldata1B=data;clear data;
        else
            alldata1B=[];
        end
    else
        alldata1A=tmpdata1A;clear tmpdata1A;
        alldata1B=tmpdata1B;clear tmpdata1B;
    end
    
    if isCompCntr(iC)
        if isDifGroups(iC)
            
            if isempty(tmpdata2A)
                genprefix = sprintf('%s_%s', experimentid, multiscanid{1}); tmpdatafile=[genprefix,'_tmegpreproc_',sortgroupnames{iC}{2}];
                load(tmpdatafile,'data');
                %=====================
                if data.fsample>2000
                    newFs=data.fsample/4;
                elseif data.fsample>1000
                    newFs=data.fsample/2;
                else
                    newFs=data.fsample;
                end
                cfg = [];
                cfg.detrend = 'no';
                cfg.resamplefs = newFs;
                data = ft_resampledata(cfg, data); % clear data;
                
                %=========================
                alldata2A=data; clear data;
                if is2Files
                    genprefix = sprintf('%s_%s', experimentid, multiscanid{2}); tmpdatafile=[genprefix,'_tmegpreproc_',sortgroupnames{iC}{2}];
                    load(tmpdatafile,'data');
                    %=====================
                    if data.fsample>2000
                        newFs=data.fsample/4;
                    elseif data.fsample>1000
                        newFs=data.fsample/2;
                    else
                        newFs=data.fsample;
                    end
                    cfg = [];
                    cfg.detrend = 'no';
                    cfg.resamplefs = newFs;
                    data = ft_resampledata(cfg, data); % clear data;
                    
                    %=========================
                    alldata2B=data;clear data;
                else
                    alldata2B=[];
                end
                
                
            else
                alldata2A=tmpdata2A;clear tmpdata2A;
                alldata2B=tmpdata2B;clear tmpdata2B;
            end
            
        else
            alldata2A=alldata1A; clear tmpdata2A;
            alldata2B=alldata1B; clear tmpdata2B;
        end
    else
        alldata2A=[];
        alldata2B=[];
    end
    
    %------------------------------------------------------
    % Get data from constrast conditions and merge from different files
    sel=sortcontrlist{1}{iC}.selection{1};
    if ~isempty(sel)
        datacntr1A=ft_selectdata(alldata1A,'rpt',sel);
    else
        datacntr1A=[];
    end
    
    if is2Files,
        sel=sortcontrlist{2}{iC}.selection{1};
        
        if ~isempty(sel)
            datacntr1B=ft_selectdata(alldata1B,'rpt',sel);
            if ~isempty(datacntr1A)
                cfg=[]; %
                [indChA,indChB]=match_str(datacntr1A.label,datacntr1B.label); datacntr1A=ft_selectdata(datacntr1A,'channel',indChA); datacntr1B=ft_selectdata(datacntr1B,'channel',indChB);
                datacntr1=ft_appenddata(cfg,datacntr1A,datacntr1B); clear datacntr1A datacntr2A
            else
                datacntr1=datacntr1B; clear datacntr1B;
            end
        else
            datacntr1=datacntr1A; clear datacntr1A;
        end
    else
        datacntr1=datacntr1A; clear datacntr1A;
    end
    
    if isempty(datacntr1)
        errorFname=['NOTRIALSFOUND_',experimentid,'_',scanmnem,'_',sortcontrlist{1}{iC}.mnemprint,'_[MODE-',avgmode,'].txt'];
        errorTime=clock;
        errorCase=[experimentid,'_',scanmnem,'_',sortcontrlist{1}{iC}.mnemprint,'_[MODE-',avgmode,']'];
        errorType=['No trials were found for condition 1'];
        hcp_write_ascii(errorFname,'errorTime','errorCase','errorType');
        warning(['No trial found for condition 1 in ',experimentid,'  - ',sortcontrlist{1}{iC}.mnemprint]) ;
        continue
    end
    datacntr1.grad=grad;
    %---------------------
    
    if strcmp(curCM,'emgcoh')
        emgRefChan=['EMG_',sortcontrlist{1}{iC}.mnemtrl{1}];
    end
    
    %--------------------------
    datacntr1_MEG=ft_selectdata(datacntr1,'channel',{'MEG'});
    if strcmp(curCM,'emgcoh')
        datacntr1_EMG=ft_selectdata(datacntr1,'channel',emgRefChan);
    else
        datacntr1_EMG=[];
    end
    
    if strcmp(avgmode,'planar'),
        if ~exist('neighbours','var')
            cfg=[];
            cfg.method='distance';
            cfg.neighbourdist = 0.037;     %  (3.7 cm is the distance to 1st order neighbs)
            neighbours = ft_prepare_neighbours(cfg, datacntr1);
        end
        
        cfg = [];
        cfg.feedback='no';
        cfg.planarmethod = 'sincos';
        cfg.neighbours= neighbours;
        %tmpData.grad=grad;
        datacntr1_MEG = ft_megplanar(cfg, datacntr1_MEG);
        
    end
    clear datacntr1;
    %--------------------------
    
    if isCompCntr(iC)
        sel=sortcontrlist{1}{iC}.selection{2};
        if ~isempty(sel)
            datacntr2A=ft_selectdata(alldata2A,'rpt',sel);
        else
            datacntr2A=[];
        end
        if is2Files,
            sel=sortcontrlist{2}{iC}.selection{2};
            if ~isempty(sel)
                datacntr2B=ft_selectdata(alldata2B,'rpt',sel);
                if ~isempty(datacntr2A)
                    cfg=[]; %
                    [indChA,indChB]=match_str(datacntr2A.label,datacntr2B.label); datacntr2A=ft_selectdata(datacntr2A,'channel',indChA); datacntr2B=ft_selectdata(datacntr2B,'channel',indChB);
                    datacntr2=ft_appenddata(cfg,datacntr2A,datacntr2B); clear datacntr2A datacntr2B
                else
                    datacntr2=datacntr2B; clear datacntr2B;
                end
            else
                datacntr2=datacntr2A; clear datacntr2A;
            end
        else
            datacntr2=datacntr2A; clear datacntr2A;
        end
        
        
        if isempty(datacntr2)
            errorFname=['NOTRIALSFOUND_',experimentid,'_',scanmnem,'_',sortcontrlist{1}{iC}.mnemprint,'_[MODE-',avgmode,'].txt'];
            errorTime=clock;
            errorCase=[experimentid,'_',scanmnem,'_',sortcontrlist{1}{iC}.mnemprint,'_[MODE-',avgmode,']'];
            errorType=['No trials were found for condition 2'];
            hcp_write_ascii(errorFname,'errorTime','errorCase','errorType');
            warning(['No trial found for condition 2 in ',experimentid,'  - ',sortcontrlist{1}{iC}.mnemprint]) ;
            continue
        end
        
        datacntr2.grad=grad;
        %--------------------------
        datacntr2_MEG=ft_selectdata(datacntr2,'channel',{'MEG'});
        if strcmp(curCM,'emgcoh')
            datacntr2_EMG=ft_selectdata(datacntr2,'channel',emgRefChan);
        else
            datacntr2_EMG=[];
        end
        
        if strcmp(avgmode,'planar'),
            cfg = [];
            cfg.feedback='no';
            cfg.planarmethod = 'sincos';
            cfg.neighbours= neighbours;
            %tmpData.grad=grad;
            datacntr2_MEG = ft_megplanar(cfg, datacntr2_MEG);
            
        end
        %--------------------------
        clear datacntr2;
    else
        datacntr2_MEG=[];
        datacntr2_EMG=[];
    end
    %------------------------------------------------------
    % Compute average for each condition
    %==========================================================================
    
    %  -------------- lower frequencies ----------------- %
    foi1 = [1:1:30];      % frequencies from 2.5 to 29 Hz in 1 Hz steps
    %0.5*ones(length(foi1),1)';  % time window: 500 ms
    %  -------------- higher frequencies ----------------- %
    foi2 = [31:2:100];           % frequencies from 30 to 100 Hz in 2.5 Hz steps
    %tw2=0.5*ones(length(foi2),1)';  % time window: 200 ms
    %  -------------- everything together ----------------- %
    foi=[foi1 foi2];
    %tw=[tw1 tw2];
    tw=5./foi;
    tw(tw>0.5)=0.5;
    
    
    
    %==========================================================================
    freqCfg              = [];
    freqCfg.pad          = 8;
    freqCfg.toiStep      = 0.025;
    freqCfg.method       = 'mtmconvol';
    freqCfg.taper        = 'hanning';
    freqCfg.foi          = foi;                          % analysis 2 to 30 Hz in steps of 2 Hz
    freqCfg.t_ftimwin    = tw;                           % length of time window = 0.5 sec
    freqCfg.keeptrials   = 'yes';                          % important keeptrials
    freqCfg.feedback     = 'yes';
    %==========================================================================
    
    timeperiods1=sortcontrlist{1}{iC}.timeperiods{1};
    if isempty(timeperiods1)
        timeperiods1=datacntr1_MEG.time{1};
    end
    if isCompCntr(iC)
        timeperiods2=sortcontrlist{1}{iC}.timeperiods{2};
        if isempty(timeperiods2)
            timeperiods2=datacntr2_MEG.time{1};
        end
    end
    
    
    
    
    if strcmp(curCM,'emgcoh')
        %emgRefChan=['EMG_',sortcontrlist{1}{iC}.mnemtrl{1}];
        
        if strcmp(avgmode,'planar'),
            
            
            freqCfg.output       = 'fourier';
            freqCfg.toi=timeperiods1;
            freqCfg.keeptrials='yes';
            freqCfg.feedback='text';
            
            freqfour_MEG=ft_freqanalysis(freqCfg,datacntr1_MEG);clear datacntr1_MEG;
            freqfour_EMG=ft_freqanalysis(freqCfg,datacntr1_EMG);clear datacntr1_EMG;
            
            combCfg=[];
            combCfg.combinemethod='svd';
            freqfour_MEG=ft_combineplanar(combCfg, freqfour_MEG);
            %freqplanar.grad=grad;
            
            avg1=ft_freqdescriptives([],freqfour_MEG);
            
            cfg=[];cfg.parameter={'fourierspctrm','cumtapcnt'};
            freqfour_MEG=ft_appendfreq(cfg,freqfour_MEG,freqfour_EMG);clear freqfour_EMG;
            
            
            
            
            cfg            = [];
            cfg.method     = 'coh';
            cfg.channelcmb = {'MEG' emgRefChan};
            freqemgcoh   = ft_connectivityanalysis(cfg, freqfour_MEG);clear freqfour_MEG;
            
        else
            
            datacntr1=ft_appenddata([],datacntr1_MEG,datacntr1_EMG);clear datacntr1_MEG datacntr1_EMG
            freqCfg.output       = 'powandcsd';
            freqCfg.toi=timeperiods1;
            freqCfg.keeptrials='no';
            freqCfg.feedback='text';
            freqCfg.channelcmb = ft_channelcombination({'MEG' emgRefChan},datacntr1.label);
            freqpowcsd=ft_freqanalysis(freqCfg,datacntr1);
            tmpcontrdata1=ft_selectdata(freqpowcsd,'channel','MEG');
            avg1=ft_freqdescriptives([],tmpcontrdata1); clear tmpcontrdata1;
            
            cfg            = [];
            cfg.method     = 'coh';
            cfg.channelcmb = ft_channelcombination({'MEG' emgRefChan},datacntr1.label);
            freqemgcoh   = ft_connectivityanalysis(cfg, freqpowcsd);clear freqpowcsd;
            
        end
        
        
        [ind1,ind2]=match_str(avg1.label,freqemgcoh.labelcmb(:,1));
        
        avg1.powspctrm=freqemgcoh.cohspctrm(ind2,:,:); clear freqemgcoh;
        
    else
        
        if strcmp(avgmode,'planar'),
            
            
            freqCfg.output       = 'fourier';
            freqCfg.toi=timeperiods1;
            freqCfg.keeptrials='yes';
            freqCfg.feedback='text';
            
            freqfour_MEG=ft_freqanalysis(freqCfg,datacntr1_MEG);clear datacntr1_MEG;
            if ~isempty(datacntr1_EMG)
                freqfour_EMG=ft_freqanalysis(freqCfg,datacntr1_EMG);clear datacntr1_EMG;
            else
                freqfour_EMG=[];
            end
            combCfg=[];
            combCfg.combinemethod='svd';
            freqfour_MEG=ft_combineplanar(combCfg, freqfour_MEG);
            %freqplanar.grad=grad;
            
            if ~isempty(freqfour_EMG)
                cfg=[];cfg.parameter={'fourierspctrm','cumtapcnt'};
                freqfour_MEG=ft_appendfreq(cfg,freqfour_MEG,freqfour_EMG);
            end
            
            avg1=ft_freqdescriptives([],freqfour_MEG); clear  freqFour_MEG freqFour_EMG;
            
        else
            if ~isempty(datacntr1_EMG)
                datacntr1=ft_appenddata([],datacntr1_MEG,datacntr1_EMG);clear datacntr1_MEG datacntr1_EMG
            else
                datacntr1=datacntr1_MEG;clear datacntr1_MEG datacntr1_EMG
            end
            freqCfg.output       = 'pow';
            freqCfg.toi=timeperiods1;
            freqCfg.keeptrials='no';
            freqCfg.feedback='text';
            avg1=ft_freqanalysis(freqCfg,datacntr1);clear datacntr1;
            
        end
        
        
        
        
        
        if isCompCntr(iC)
            if strcmp(avgmode,'planar'),
                
                
                freqCfg.output       = 'fourier';
                freqCfg.toi=timeperiods2;
                freqCfg.keeptrials='yes';
                freqCfg.feedback='text';
                
                freqfour_MEG=ft_freqanalysis(freqCfg,datacntr2_MEG);clear datacntr2_MEG;
                if ~isempty(datacntr2_EMG)
                    freqfour_EMG=ft_freqanalysis(freqCfg,datacntr2_EMG);clear datacntr2_EMG;
                else
                    freqfour_EMG=[];
                end
                
                combCfg=[];
                combCfg.combinemethod='svd';
                freqfour_MEG=ft_combineplanar(combCfg, freqfour_MEG);
                %freqplanar.grad=grad;
                
                if ~isempty(freqfour_EMG)
                    cfg=[];cfg.parameter={'fourierspctrm','cumtapcnt'};
                    freqfour_MEG=ft_appendfreq(cfg,freqfour_MEG,freqfour_EMG);
                end
                
                avg2=ft_freqdescriptives([],freqfour_MEG); clear  freqFour_MEG freqFour_EMG;
                
            else
                if ~isempty(datacntr2_EMG)
                    datacntr2=ft_appenddata([],datacntr2_MEG,datacntr2_EMG);clear datacntr2_MEG datacntr2_EMG
                else
                    datacntr2=datacntr2_MEG;clear datacntr2_MEG datacntr2_EMG
                end
                freqCfg.output       = 'pow';
                freqCfg.toi=timeperiods2;
                freqCfg.keeptrials='no';
                freqCfg.feedback='text';
                avg2=ft_freqanalysis(freqCfg,datacntr2);clear datacntr2;
                
            end
            
        end
    end
    
    %------------------------------------------------------
    % Apply contrast comparison conditions
    % --- Check Time
    operation=sortcontrlist{1}{iC}.operation;
    if isCompCntr(iC)
        if strcmp(operation,'diff')
            avg1.powspctrm=avg1.powspctrm-avg2.powspctrm;
        elseif strcmp(operation,'relch')
            avg1.powspctrm=(avg1.powspctrm-avg2.powspctrm)./avg1.powspctrm;
        else
            error('For comparisons between conditions the supported operations are diff and relch')
        end
    end
    %--------------------------------------------------------
    % PLOT
    
    if ~isCompCntr(iC)
        if ~strcmp(curCM,'emgcoh')
            baseline=sortcontrlist{1}{iC}.baseline{1};
            if isempty(baseline)
                baseline=[-0.5 0]; % se a default baseline
            end
            cfg=[];
            cfg.baseline      = baseline;
            cfg.baselinetype  = 'relchange';
            avg2plot=ft_freqbaseline(cfg,avg1);
        else
            avg2plot=avg1;
            
        end
        
    else
        avg2plot=avg1;
        
    end
    
    data=avg1;
    hcp_write_matlab(saveFnameData,'data');
    
    plotsensmainareas(sortcontrlist{1}{iC},avg2plot,bandinfo,experimentid,scanmnem,saveFnameImage,avgmode);
    
    prevGroups=sortgroupnames{iC};
    
    
end



end
%--- END OF MAIN --------------------------
%=====================================================
%=====================================================
%=====================================================
%%
function[allfreqpow]=loadsingallfreqpow(experimentid, scanid,groupname,bands,inavgmode)
allfreqpow=[];
for iBand=1:length(bands)
    genprefix = sprintf('%s_%s', experimentid, scanid);
    if strcmp(inavgmode,'mag')
        freqfile=[genprefix,'_tfsens_',groupname,'_',bands{iBand}];
        hcp_read_matlab(freqfile,'freq');
    elseif strcmp(inavgmode,'planar')
        freqfile=[genprefix,'_tfsens_',groupname,'_',bands{iBand},'_planar'];
        hcp_read_matlab(freqfile,'freqplanar');
        freq=freqplanar; clear freqplanar;
    end
    cfg=[];cfg.keeptrials='yes';
    freqpow=ft_freqdescriptives(cfg,freq);
    if iBand==1
        allfreqpow=freqpow;
    else
        cfg=[]; cfg.parameter={'powspctrm'};%,'cumtapcnt'};
        allfreqpow=ft_appendfreq(cfg,allfreqpow,freqpow);
    end
end
end
%%
function[allfreqfour]=loadsingallfreqfour(experimentid, scanid,groupname,bands,inavgmode)
allfreqfour=[];
allcumtapcnt=[];
for iBand=1:length(bands)
    genprefix = sprintf('%s_%s', experimentid, scanid);
    if strcmp(inavgmode,'mag')
        freqfile=[genprefix,'_tfsens_',groupname,'_',bands{iBand}];
        hcp_read_matlab(freqfile,'freq');
    elseif strcmp(inavgmode,'planar')
        freqfile=[genprefix,'_tfsens_',groupname,'_',bands{iBand},'_planar'];
        hcp_read_matlab(freqfile,'freqplanar');
        freq=freqplanar; clear freqplanar;
    end
    
    
    if iBand==1
        allfreqfour=freq;
        allcumtapcnt=freq.cumtapcnt;
    else
        allcumtapcnt=[allcumtapcnt freq.cumtapcnt];
        cfg=[];cfg.parameter={'fourierspctrm'};
        allfreqfour=ft_appendfreq(cfg,allfreqfour,freq);  clear freq;
        
    end
    allfreqfour.cumtapcnt=allcumtapcnt;
end
end

%%
function[]=plotsensmainareas(incontrast,inavg,bandinfo,experimentid,scanmnem,saveFileName,avgmode)

channels=[];
channels.occ = {'A136', 'A137', 'A138', 'A163', 'A164', 'A165', 'A166', 'A185', 'A186', 'A187', 'A203', 'A204', 'A205', 'A219', 'A220', 'A221', 'A238', 'A239'};
channels.occ_L = {'A134', 'A161', 'A182', 'A183', 'A199', 'A200', 'A201', 'A216', 'A217', 'A218', 'A235'};
channels.occ_R = {'A107', 'A139', 'A140', 'A166', 'A167', 'A168', 'A188', 'A189', 'A190', 'A205', 'A206', 'A207', 'A222', 'A223'};
channels.temppf_L = {'A96', 'A97', 'A98', 'A127', 'A128', 'A129', 'A130', 'A155', 'A156', 'A157', 'A179', 'A196'};
channels.temppf_R = {'A112', 'A113', 'A144', 'A145', 'A146', 'A147', 'A148', 'A171', 'A172', 'A173', 'A174', 'A193', 'A209', 'A210', 'A211', 'A227', 'A247'};
channels.front = {'A37', 'A38', 'A60', 'A61', 'A62', 'A63', 'A87', 'A88', 'A89', 'A90', 'A91', 'A92', 'A117', 'A118', 'A119', 'A120', 'A123', 'A124', 'A149', 'A150', 'A151'};
channels.motor_L = {'A7', 'A9', 'A10', 'A11', 'A22', 'A23', 'A24', 'A25', 'A26', 'A41', 'A42', 'A43', 'A44', 'A45', 'A66', 'A67', 'A68', 'A69', 'A70'};
channels.motor_R ={'A13', 'A14', 'A15', 'A16', 'A17', 'A30', 'A31', 'A32', 'A33', 'A34', 'A53', 'A54', 'A55', 'A56', 'A57', 'A80', 'A81', 'A82', 'A83', 'A84'};
channels.pariet= {'A49', 'A73', 'A74', 'A75', 'A76', 'A77', 'A102', 'A103', 'A104', 'A105', 'A106', 'A107', 'A134', 'A135', 'A139', 'A162'};




%==========================================================================
%==========================================================================
% Make figure outline
h1=figure();
set(h1,'papertype','A3');
set(h1,'paperunits','centimeters')
papersize=get(h1,'PaperSize');
paperposition = [1 1 papersize-1 ];


bigFonts1=floor(min(papersize)*0.25);
bigFonts2=floor(min(papersize)*0.5);
bigFonts3=floor(min(papersize)*1);
bigFonts4=floor(min(papersize)*2);

set(h1,'papersize',papersize);
set(h1,'paperposition',paperposition);
set(h1,'position',[10 10 15*papersize]);
%---------------------------------------------
Aha(1,1)=axes();
set(Aha(1,1),'position',[0.1 0.625 0.2 0.1]);
Aha(1,2)=axes();
set(Aha(1,2),'position',[0.375 0.625 0.2 0.1]);
Aha(1,3)=axes();
set(Aha(1,3),'position',[0.65 0.625 0.2 0.1]);
%----
Aha(2,1)=axes();
set(Aha(2,1),'position',[0.1 0.475 0.2 0.1]);
Aha(2,2)=axes();
set(Aha(2,2),'position',[0.375 0.475 0.2 0.1]);
Aha(2,3)=axes();
set(Aha(2,3),'position',[0.65 0.475 0.2 0.1]);
%----
Aha(3,1)=axes();
set(Aha(3,1),'position',[0.1 0.325 0.2 0.1]);
Aha(3,2)=axes();
set(Aha(3,2),'position',[0.375 0.325 0.2 0.1]);
Aha(3,3)=axes();
set(Aha(3,3),'position',[0.65 0.325 0.2 0.1]);
%=======================================================
%---------------------------------------------
Bha(1)=axes();
set(Bha(1),'position',[0.05 0.15 0.15 0.10],'YTick',[],'XTick',[]);
Bha(2)=axes();
set(Bha(2),'position',[0.225 0.15 0.15 0.10],'YTick',[],'XTick',[]);
Bha(3)=axes();
set(Bha(3),'position',[0.4 0.15 0.15 0.10],'YTick',[],'XTick',[]);
Bha(4)=axes();
set(Bha(4),'position',[0.575 0.15 0.15 0.10],'YTick',[],'XTick',[]);
%----

Bha(5)=axes();
set(Bha(5),'position',[0.05 0.01 0.15 0.10],'YTick',[],'XTick',[]);
Bha(6)=axes();
set(Bha(6),'position',[0.225 0.01 0.15 0.10],'YTick',[],'XTick',[]);
Bha(7)=axes();
set(Bha(7),'position',[0.4 0.01 0.15 0.10],'YTick',[],'XTick',[]);
Bha(8)=axes();
set(Bha(8),'position',[0.575 0.01 0.15 0.10],'YTick',[],'XTick',[]);
%----
%=======================================================
%=======================================================
%=======================================================


%plot
%----
axes(Aha(1,1));
tmpavg=ft_selectdata(inavg,'channel',channels.temppf_L,'avgoverchan','yes');colorbar off;
minVal=min(tmpavg.powspctrm(:)); maxVal=max(tmpavg.powspctrm(:)); if sign(minVal)~=sign(maxVal), topozlim=max(abs([minVal maxVal]))*[-1 1]; else  topozlim=[minVal maxVal];end
plotCfg=[]; plotCfg.zlim=topozlim; plotCfg.renderer ='painters';plotCfg.layout='4D248.lay'; plotCfg.ylim=[4 80];ft_singleplotTFR(plotCfg,  tmpavg);colorbar off;
title('TFR avg of Left Temp-PF sens','fontsize',bigFonts2);
%----
axes(Aha(1,2));
tmpavg=ft_selectdata(inavg,'channel',channels.front,'avgoverchan','yes');colorbar off;
minVal=min(tmpavg.powspctrm(:)); maxVal=max(tmpavg.powspctrm(:)); if sign(minVal)~=sign(maxVal), topozlim=max(abs([minVal maxVal]))*[-1 1]; else  topozlim=[minVal maxVal];end
plotCfg=[]; plotCfg.zlim=topozlim; plotCfg.renderer ='painters';plotCfg.layout='4D248.lay'; plotCfg.ylim=[4 80];ft_singleplotTFR(plotCfg,  tmpavg);colorbar off;
title('TFR avg of Frontal sens','fontsize',bigFonts2);
%----
axes(Aha(1,3));
tmpavg=ft_selectdata(inavg,'channel',channels.temppf_R,'avgoverchan','yes');colorbar off;
minVal=min(tmpavg.powspctrm(:)); maxVal=max(tmpavg.powspctrm(:)); if sign(minVal)~=sign(maxVal), topozlim=max(abs([minVal maxVal]))*[-1 1]; else  topozlim=[minVal maxVal];end
plotCfg=[]; plotCfg.zlim=topozlim; plotCfg.renderer ='painters';plotCfg.layout='4D248.lay'; plotCfg.ylim=[4 80];ft_singleplotTFR(plotCfg,  tmpavg);colorbar off;
title('TFR avg of Right Temp-PF sens','fontsize',bigFonts2);
%---------------------------------------------------
axes(Aha(2,1));
tmpavg=ft_selectdata(inavg,'channel',channels.motor_L,'avgoverchan','yes');colorbar off;
minVal=min(tmpavg.powspctrm(:)); maxVal=max(tmpavg.powspctrm(:)); if sign(minVal)~=sign(maxVal), topozlim=max(abs([minVal maxVal]))*[-1 1]; else  topozlim=[minVal maxVal];end
plotCfg=[]; plotCfg.zlim=topozlim; plotCfg.renderer ='painters';plotCfg.layout='4D248.lay'; plotCfg.ylim=[4 80];ft_singleplotTFR(plotCfg,  tmpavg);colorbar off;
title('TFR avg of Left Motor sens','fontsize',bigFonts2);
%----
axes(Aha(2,2));
tmpavg=ft_selectdata(inavg,'channel',channels.pariet,'avgoverchan','yes');colorbar off;
minVal=min(tmpavg.powspctrm(:)); maxVal=max(tmpavg.powspctrm(:)); if sign(minVal)~=sign(maxVal), topozlim=max(abs([minVal maxVal]))*[-1 1]; else  topozlim=[minVal maxVal];end
plotCfg=[]; plotCfg.zlim=topozlim; plotCfg.renderer ='painters';plotCfg.layout='4D248.lay'; plotCfg.ylim=[4 80];ft_singleplotTFR(plotCfg,  tmpavg);colorbar off;
title('TFR avg of Parietal sens','fontsize',bigFonts2);
%----
axes(Aha(2,3));
tmpavg=ft_selectdata(inavg,'channel',channels.motor_R,'avgoverchan','yes');colorbar off;
minVal=min(tmpavg.powspctrm(:)); maxVal=max(tmpavg.powspctrm(:)); if sign(minVal)~=sign(maxVal), topozlim=max(abs([minVal maxVal]))*[-1 1]; else  topozlim=[minVal maxVal];end
plotCfg=[]; plotCfg.zlim=topozlim; plotCfg.renderer ='painters';plotCfg.layout='4D248.lay'; plotCfg.ylim=[4 80];ft_singleplotTFR(plotCfg,  tmpavg);colorbar off;
title('TFR avg of Right Motor sens','fontsize',bigFonts2);
%---------------------------------------------------
axes(Aha(3,1));
tmpavg=ft_selectdata(inavg,'channel',channels.occ_L,'avgoverchan','yes');colorbar off;
minVal=min(tmpavg.powspctrm(:)); maxVal=max(tmpavg.powspctrm(:)); if sign(minVal)~=sign(maxVal), topozlim=max(abs([minVal maxVal]))*[-1 1]; else  topozlim=[minVal maxVal];end
plotCfg=[]; plotCfg.zlim=topozlim; plotCfg.renderer ='painters';plotCfg.layout='4D248.lay'; plotCfg.ylim=[4 80];ft_singleplotTFR(plotCfg,  tmpavg);colorbar off;
title('TFR avg of Left Occ. sens','fontsize',bigFonts2);
%----
axes(Aha(3,2));
tmpavg=ft_selectdata(inavg,'channel',channels.occ,'avgoverchan','yes');colorbar off;
minVal=min(tmpavg.powspctrm(:)); maxVal=max(tmpavg.powspctrm(:)); if sign(minVal)~=sign(maxVal), topozlim=max(abs([minVal maxVal]))*[-1 1]; else  topozlim=[minVal maxVal];end
plotCfg=[]; plotCfg.zlim=topozlim; plotCfg.renderer ='painters';plotCfg.layout='4D248.lay'; plotCfg.ylim=[4 80];ft_singleplotTFR(plotCfg,  tmpavg);colorbar off;
title('TFR avg of Med. Occ. sens','fontsize',bigFonts2);
%----
axes(Aha(3,3));
tmpavg=ft_selectdata(inavg,'channel',channels.occ_R,'avgoverchan','yes');colorbar off;
minVal=min(tmpavg.powspctrm(:)); maxVal=max(tmpavg.powspctrm(:)); if sign(minVal)~=sign(maxVal), topozlim=max(abs([minVal maxVal]))*[-1 1]; else  topozlim=[minVal maxVal];end
plotCfg=[]; plotCfg.zlim=topozlim; plotCfg.renderer ='painters';plotCfg.layout='4D248.lay'; plotCfg.ylim=[4 80];ft_singleplotTFR(plotCfg,  tmpavg);colorbar off;
title('TFR avg of Right Occ. sens','fontsize',bigFonts2);
%----
%---------------------------------------------------


bandStr=[];
for iBand=1:size(bandinfo,1)
    iBand
    axes(Bha(iBand));
    tmpfreqlim=bandinfo{iBand,2};
    tmpavg=ft_selectdata(inavg,'foilim',tmpfreqlim);colorbar off;
    tmpavg=ft_selectdata(tmpavg,'avgovertime','yes');
    %------
    minVal=min(tmpavg.powspctrm(:));
    maxVal=max(tmpavg.powspctrm(:));
    if sign(minVal)~=sign(maxVal),
        topozlim=max(abs([minVal maxVal]))*[-1 1];
    else
        topozlim=[minVal maxVal];
    end
    %------
    plotCfg=[];
    plotCfg.layout='4D248.mat';
    plotCfg.interpmethod='linear';
    plotCfg.marker='on';
    plotCfg.markersymbol='.';
    plotCfg.gridscale=50;
    plotCfg.shading='flat';
    plotCfg.comment='no';
    
    plotCfg.zlim=topozlim;
    
    ft_topoplotTFR(plotCfg,  tmpavg);
    colorbar off;
    title(bandinfo{iBand,1},'fontsize',bigFonts2);
    bandStr = [bandStr, sprintf('%s : %s to %s Hz\n',bandinfo{iBand,1},num2str(bandinfo{iBand,2}(1)),num2str(bandinfo{iBand,2}(2))) ];
end
%{
dispStr=sprintf('Bands: \n%s',bandStr);

txtBoxBand = uicontrol('style','listbox');
set(txtBoxBand,'units','normalized')
set(txtBoxBand,'Position',[0.75 0.01 0.2 0.25]);
set(txtBoxBand,'String',dispStr)
set(txtBoxBand,'Fontsize',6)
set(txtBoxBand,'HorizontalAlignment','left')

%-------------------------------------------------------------------------
print(h1,'-dpng','-r600','testTFmulti1.png');



Bha(8)=axes();
set(Bha(8),'position',[0.575 0.01 0.15 0.10],'YTick',[],'XTick',[]);
%}
dispStr=sprintf('Bands: \n%s',bandStr);
txtH=axes();axis off;
set(txtH,'position',[0.75 0.01 0.24 0.1],'YTick',[],'XTick',[]);
xlim=get(txtH,'Xlim');ylim=get(txtH,'Ylim');
ht1=text(xlim(1),ylim(2),dispStr);
set(ht1,'Fontsize',bigFonts2)

Toph1=axes();
set(Toph1,'position',[0.1 0.95 0.8 0.05]);axis off
dispStringTop1=sprintf('%s\n%s\n%s',[' experimentid:  ',regexprep(experimentid,'_','\\_'),'  scan:  ',regexprep(scanmnem,'_','\\_')],[' contrast:  ',regexprep(incontrast.mnemprint,'_','\\_')],[' avgmode:  ',regexprep(avgmode,'_','\\_')]);
htTop1=text(0,0.2,dispStringTop1);
set(htTop1,'Fontsize',bigFonts2)
%{
Toph2=axes();
set(Toph2,'position',[0.1 0.775 0.4 0.15]); axis off;
htTop2=text(0,0,[' contrast:',incontrast.mnemprint]);
%}
%print(h1,'-dpng',[saveFileName]);
hcp_write_figure(saveFileName, h1,'format','png');
close(h1);
%==========================================================================
% ATTENTION !!! THE FOLLOWING PLOT ASSUMES THAT YOU HAVE A VERY LARGE VNC
% SCREEN . OTHERWISE VNC DIES.
%{
[PATHSTR,NAME,EXT] = fileparts(saveFileName);

saveFileName2=[PATHSTR,NAME,'_multiTFR',EXT];
h2=figure;
%set(0,'Units','centimeters');
%screenSize=get(0,'screensize');
%if screenSize(3)>120,
%    set(h2,'papertype','A2');%A2
%else
    set(h2,'papertype','A3');%A3
%end
set(h2,'paperunits','centimeters')
%set(h2,'papersize',[59.4000 42]);%A2
papersize=get(h2,'PaperSize');
paperposition = [1 1 papersize-1 ];
set(h2,'papersize',papersize);
set(h2,'paperposition',paperposition);
set(h2,'position',[10 10 15*papersize]);

cfg=[];
cfg.layout='4D248.mat';
cfg.showoutline      = 'yes';
ft_multiplotTFR(cfg,inavg);


set(h2,'paperorientation','portrait');

hcp_write_figure(saveFileName2, h2,'format','png');
close(h2);
%}
%==========================================================================
end
