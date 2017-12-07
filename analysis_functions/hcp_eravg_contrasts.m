function [outStatus] = hcp_eravg_contrasts(inCfg)

% This function loads time frequency data and computes its trial average  as well
% as the average of its planar gradient
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
contrastlist=inCfg.contrastlist;


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
prevGroups={'' ''};
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
    
    
    %==================================
    % Load the time data for the corresponding group and from both
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
        hcp_read_matlab(tmpdatafile,'data');
        alldata1A=data;
        
        if iC==1,
            grad=data.grad;
            cfg=[];
            cfg.method='distance';
            cfg.neighbourdist = 0.037;     %  (3.7 cm is the distance to 1st order neighbs)
            neighbours = ft_prepare_neighbours(cfg, data);
        end
        clear data;
        
        if is2Files
            genprefix = sprintf('%s_%s', experimentid, multiscanid{2}); tmpdatafile=[genprefix,'_tmegpreproc_',sortgroupnames{iC}{1}];
            hcp_read_matlab(tmpdatafile,'data');
            alldata1B=data; clear data;
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
                hcp_read_matlab(tmpdatafile,'data');
                alldata2A=data; clear data;
                if is2Files
                    genprefix = sprintf('%s_%s', experimentid, multiscanid{2}); tmpdatafile=[genprefix,'_tmegpreproc_',sortgroupnames{iC}{2}];
                    hcp_read_matlab(tmpdatafile,'data');
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
    if isempty(sel),
        datacntr1A=[];
    else
        datacntr1A=ft_selectdata(alldata1A,'rpt',sel);
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
    
    %--------
    if isempty(datacntr1)
        errorFname=['NOTRIALSFOUND_',experimentid,'_',scanmnem,'_',sortcontrlist{1}{iC}.mnemprint,'_[MODE-',avgmode,'].txt'];
        errorTime=clock;
        errorCase=[experimentid,'_',scanmnem,'_',sortcontrlist{1}{iC}.mnemprint,'_[MODE-',avgmode,']'];
        errorType=['No trials were found for condition 1'];
        hcp_write_ascii(errorFname,'errorTime','errorCase','errorType');
        warning(['No trial found for condition 1 in ',experimentid,'  - ',sortcontrlist{1}{iC}.mnemprint]) ;
        continue
    end
    %--------
    
    if isCompCntr(iC)
        sel=sortcontrlist{1}{iC}.selection{2};
        if isempty(sel),
            datacntr2A=[];
        else
            datacntr2A=ft_selectdata(alldata2A,'rpt',sel);
        end;
        if is2Files,
            sel=sortcontrlist{2}{iC}.selection{2};
            if ~isempty(sel),
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
        
        
        %--------
        if isempty(datacntr2)
            errorFname=['NOTRIALSFOUND_',experimentid,'_',scanmnem,'_',sortcontrlist{1}{iC}.mnemprint,'_[MODE-',avgmode,'].txt'];
            errorTime=clock;
            errorCase=[experimentid,'_',scanmnem,'_',sortcontrlist{1}{iC}.mnemprint,'_[MODE-',avgmode,']'];
            errorType=['No trials were found for condition 2'];
            hcp_write_ascii(errorFname,'errorTime','errorCase','errorType');
            warning(['No trial found for condition 2 in ',experimentid,'  - ',sortcontrlist{1}{iC}.mnemprint]) ;
            continue
        end
        %--------
        
    end
    %------------------------------------------------------
    %====================================================
    %====================================================
    %====================================================
    % DO basic filtering
    
    flowcfg=[];
    flowcfg.hpfilter='yes';
    flowcfg.hpfreq=[1];
    flowcfg.hpfiltord=4;
    
    fhighcfg=[];
    fhighcfg.lpfilter='yes';
    fhighcfg.lpfreq=[40];
    
    
    datacntr1=ft_preprocessing(flowcfg,datacntr1);
    datacntr1=ft_preprocessing(fhighcfg,datacntr1);
    if isCompCntr(iC)
        datacntr2=ft_preprocessing(flowcfg,datacntr2);
        datacntr2=ft_preprocessing(fhighcfg,datacntr2);
    end
    
    %=========================================================
    % Demean
    baseline1=sortcontrlist{1}{iC}.baseline{1};
    dmcfg=[];
    dmcfg.demean='yes';
    dmcfg.baselinewindow=baseline1;
    datacntr1=ft_preprocessing(dmcfg,datacntr1);
    if isCompCntr(iC)
        baseline2=sortcontrlist{1}{iC}.baseline{2};
        dmcfg=[];
        dmcfg.demean='yes';
        dmcfg.baselinewindow=baseline2;
        datacntr2=ft_preprocessing(dmcfg,datacntr2);
    end
    %=========================================================
    % Average
    avg1=ft_timelockanalysis([],datacntr1); clear datacntr1;
    if isCompCntr(iC)
        avg2=ft_timelockanalysis([],datacntr2); clear datacntr2;
    end
    
    
    
    cfg = [];
    cfg.feedback='no';
    %cfg.planarmethod = 'sincos';
    cfg.neighbours= neighbours;
    avg1.grad=grad;
    tmpplanar = ft_megplanar(cfg, avg1);
    
    combCfg=[];
    %combCfg.combinemethod='svd';
    avg1planar=ft_combineplanar(combCfg, tmpplanar);
    if isCompCntr(iC)
        cfg = [];
        cfg.feedback='no';
        %cfg.planarmethod = 'sincos';
        cfg.neighbours= neighbours;
        avg2.grad=grad;
        tmpplanar = ft_megplanar(cfg, avg2);
        
        combCfg=[];
        %combCfg.combinemethod='svd';
        avg2planar=ft_combineplanar(combCfg, tmpplanar);
    end
    
    %====================================================
    %====================================================
    %====================================================
    
    %------------------------------------------------------
    % Apply contrast comparison conditions
    % --- Check Time
    timeperiods=sortcontrlist{1}{iC}.timeperiods{1};
    if ~isempty(timeperiods)
        if size(timeperiods,2)==2,
            isTimeWin=1;
            timepoints=mean(timeperiods,2)';
        else
            isTimeWin=0;
            timepoints=timeperiods';
        end
        Ntpts=size(timeperiods,1);
    else
        isTimeWin=0;
        timepoints=[]; % Leave empty as timeperiods to signify that the original time axis must be kept
        Ntpts=length(avg1.time);
    end
    
    if (isTimeWin)&(Ntpts==1)
        isSingleWin=1;
    elseif (isTimeWin)&(Ntpts>1)
        isSingleWin=0;
        error('eravg for a time structure of multiple time windows is not yet supported. Put just one time window that sets the time limits or leave empty to get all ');
    else
        isSingleWin=0;
    end
    
    %---------------
    
    if isSingleWin
        avg1=ft_selectdata(avg1,'toilim',timeperiods);
        avg1planar=ft_selectdata(avg1planar,'toilim',timeperiods);
        if isCompCntr(iC)
            avg2=ft_selectdata(avg2,'toilim',timeperiods);
            avg2planar=ft_selectdata(avg2planar,'toilim',timeperiods);
        end
    elseif isempty(timeperiods) % original
        disp('Just using the original time axis');
        
    else
        %indtimeord=match_str(tokenize(avg1.dimord,'_'),'time');
        Ndattimes=length(avg1.time);
        tmpIndMat=[];
        for iTime=1:Ntpts,
            tmpIndMat(iTime)=nearest( avg1.time,timepoints(iTime));
        end
        tmpIndMat=unique(tmpIndMat);
        avg1.time=avg1.time(tmpIndMat);
        avg1.avg=avg1.avg(:,tmpIndMat);
        avg1planar.time=avg1planar.time(tmpIndMat);
        avg1planar.avg=avg1planar.avg(:,tmpIndMat);
        
        if isCompCntr(iC)
            avg2.time=avg2.time(tmpIndMat);
            avg2.avg=avg2.avg(:,tmpIndMat);
            avg2planar.time=avg2planar.time(tmpIndMat);
            avg2planar.avg=avg2planar.avg(:,tmpIndMat);
            
        end
        
    end
    %------------------------------------------------------
    %====================================================================
    operation=sortcontrlist{1}{iC}.operation;
    if isCompCntr(iC)
        if strcmp(operation,'diff')
            avg1.avg=avg1.avg-avg2.avg;
            avg1planar.avg=avg1planar.avg-avg2planar.avg;
        else
            error('For eravg comparisons between conditions the supported operation is diff')
        end
    end
    
    %--------------------------------------------------------
    
    % PLOT
    
    avg2plot=avg1;
    avgmode='mag';
    saveFnameData=[experimentid,'_',scanmnem,'_',sortcontrlist{1}{iC}.mnemprint,'_[MODE-',avgmode,']'];
    saveFnameImage=[experimentid,'_',scanmnem,'_',sortcontrlist{1}{iC}.mnemprint,'_[MODE-',avgmode,']_plot.png'];
    plotsensmainareas(sortcontrlist{1}{iC},avg2plot,experimentid,scanmnem,saveFnameImage,avgmode);
    data=avg1;
    hcp_write_matlab(saveFnameData,'data'); clear data;
    
    avg2plot=avg1planar;
    avgmode='planar';
    saveFnameData=[experimentid,'_',scanmnem,'_',sortcontrlist{1}{iC}.mnemprint,'_[MODE-',avgmode,']'];
    saveFnameImage=[experimentid,'_',scanmnem,'_',sortcontrlist{1}{iC}.mnemprint,'_[MODE-',avgmode,']_plot.png'];
    plotsensmainareas(sortcontrlist{1}{iC},avg2plot,experimentid,scanmnem,saveFnameImage,avgmode);
    data=avg1planar;
    hcp_write_matlab(saveFnameData,'data'); clear data;
    
    
    prevGroups=sortgroupnames{iC};
    
    
    
end



end
%--- END OF MAIN --------------------------
%=====================================================
%=====================================================
%=====================================================

%%
function[]=plotsensmainareas(incontrast,inavg,experimentid,scanmnem,saveFileName,avgmode)

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
set(Bha(2),'position',[0.3 0.15 0.15 0.10],'YTick',[],'XTick',[]);
Bha(3)=axes();
set(Bha(3),'position',[0.55 0.15 0.15 0.10],'YTick',[],'XTick',[]);
Bha(4)=axes();
set(Bha(4),'position',[0.8 0.15 0.15 0.10],'YTick',[],'XTick',[]);
%----

Bha(5)=axes();
set(Bha(5),'position',[0.05 0.01 0.15 0.10],'YTick',[],'XTick',[]);
Bha(6)=axes();
set(Bha(6),'position',[0.3 0.01 0.15 0.10],'YTick',[],'XTick',[]);
Bha(7)=axes();
set(Bha(7),'position',[0.55 0.01 0.15 0.10],'YTick',[],'XTick',[]);
Bha(8)=axes();
set(Bha(8),'position',[0.8 0.01 0.15 0.10],'YTick',[],'XTick',[]);
%----
%=======================================================
%=======================================================
%=======================================================
%plot
%----
axes(Aha(1,1));
tmpavg=ft_selectdata(inavg,'channel',channels.temppf_L,'avgoverchan','yes');colorbar off;
plotCfg=[]; plotCfg.layout='4D248.mat'; ft_singleplotER(plotCfg,  tmpavg);colorbar off;
title('ER avg of Left Temp-PF sens','fontsize',bigFonts2);
%----
axes(Aha(1,2));
tmpavg=ft_selectdata(inavg,'channel',channels.front,'avgoverchan','yes');colorbar off;
plotCfg=[]; plotCfg.layout='4D248.mat'; ft_singleplotER(plotCfg,  tmpavg);colorbar off;
title('ER avg of Frontal sens','fontsize',bigFonts2);
%----
axes(Aha(1,3));
tmpavg=ft_selectdata(inavg,'channel',channels.temppf_R,'avgoverchan','yes');colorbar off;
plotCfg=[]; plotCfg.layout='4D248.mat'; ft_singleplotER(plotCfg,  tmpavg);colorbar off;
title('ER avg of Right Temp-PF sens','fontsize',bigFonts2);
%---------------------------------------------------
axes(Aha(2,1));
tmpavg=ft_selectdata(inavg,'channel',channels.motor_L,'avgoverchan','yes');colorbar off;
plotCfg=[]; plotCfg.layout='4D248.mat'; ft_singleplotER(plotCfg,  tmpavg);colorbar off;
title('ER avg of Left Motor sens','fontsize',bigFonts2);
%----
axes(Aha(2,2));
tmpavg=ft_selectdata(inavg,'channel',channels.pariet,'avgoverchan','yes');colorbar off;
plotCfg=[]; plotCfg.layout='4D248.mat'; ft_singleplotER(plotCfg,  tmpavg);colorbar off;
title('ER avg of Parietal sens','fontsize',bigFonts2);
%----
axes(Aha(2,3));
tmpavg=ft_selectdata(inavg,'channel',channels.motor_R,'avgoverchan','yes');colorbar off;
plotCfg=[]; plotCfg.layout='4D248.mat'; ft_singleplotER(plotCfg,  tmpavg);colorbar off;
title('ER avg of Right Motor sens','fontsize',bigFonts2);
%---------------------------------------------------
axes(Aha(3,1));
tmpavg=ft_selectdata(inavg,'channel',channels.occ_L,'avgoverchan','yes');colorbar off;
plotCfg=[]; plotCfg.layout='4D248.mat'; ft_singleplotER(plotCfg,  tmpavg);colorbar off;
title('ER avg of Left Occ. sens','fontsize',bigFonts2);
%----
axes(Aha(3,2));
tmpavg=ft_selectdata(inavg,'channel',channels.occ,'avgoverchan','yes');colorbar off;
plotCfg=[]; plotCfg.layout='4D248.mat'; ft_singleplotER(plotCfg,  tmpavg);colorbar off;
title('ER avg of Med. Occ. sens','fontsize',bigFonts2);
%----
axes(Aha(3,3));
tmpavg=ft_selectdata(inavg,'channel',channels.occ_R,'avgoverchan','yes');colorbar off;
plotCfg=[]; plotCfg.layout='4D248.mat'; ft_singleplotER(plotCfg,  tmpavg);colorbar off;
title('ER avg of Right Occ. sens','fontsize',bigFonts2);
%----
%---------------------------------------------------
%---------------------------------------------------
%---------------------------------------------------
timLim=inavg.time([1 end]);

tmpTopotimes=[timLim(1):(diff(timLim)./8):timLim(end)];
tmpStartTimes =tmpTopotimes(1:end-1);
tmpEndTimes   =tmpTopotimes(2:end);

format short g
for iPer=1:length(tmpStartTimes)
    iPer
    axes(Bha(iPer));
    tmptoilim=[tmpStartTimes(iPer) tmpEndTimes(iPer)];
    tmpavg=ft_selectdata(inavg,'toilim',tmptoilim);colorbar off;
    %tmpavg=ft_selectdata(tmpavg,'avgovertime','yes');
    %minVal=min(tmpavg.avg(:)); maxVal=max(tmpavg.avg(:)); if sign(minVal)~=sign(maxVal), topozlim=max(abs([minVal maxVal]))*[-1 1]; else  topozlim=[minVal maxVal];end
    
    plotCfg=[];
    plotCfg.layout='4D248.mat';
    plotCfg.interpmethod='linear';
    plotCfg.marker='on';
    plotCfg.markersymbol='.';
    plotCfg.gridscale=50;
    plotCfg.shading='flat';
    plotCfg.comment='no';
    
    %plotCfg.zlim=topozlim;
    ft_topoplotER(plotCfg,  tmpavg);
    colorbar off;
    titlStr=sprintf('%2.3f to %2.3f',tmptoilim(1),tmptoilim(2));
    title(titlStr,'fontsize',bigFonts2);
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
hcp_write_figure('testTFmulti1.png', h1, 'resolution', 600);



Bha(8)=axes();
set(Bha(8),'position',[0.575 0.01 0.15 0.10],'YTick',[],'XTick',[]);
%}
%{
dispStr=sprintf('Bands: \n%s',bandStr);
txtH=axes();axis off;
set(txtH,'position',[0.75 0.01 0.24 0.1],'YTick',[],'XTick',[]);
xlim=get(txtH,'Xlim');ylim=get(txtH,'Ylim');
ht1=text(xlim(1),ylim(2),dispStr);
set(ht1,'Fontsize',20)
%}
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
disp('Saving the figure');
hcp_write_figure(saveFileName, h1,'format','png');
disp('Done Saving');
close(h1);
%==========================================================================
%==========================================================================
end
