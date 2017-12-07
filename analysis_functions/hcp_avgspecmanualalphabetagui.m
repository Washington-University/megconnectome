function [parH]=hcp_avgspecmanualalphabetagui(subjID,inspectra_Restin,inspectra_Motort,inspectra_Wrkmem,inspectra_StoryM,infreqs,outSaveDir)


hasRestin=~isempty(inspectra_Restin);
hasMotort=~isempty(inspectra_Motort);
hasWrkmem=~isempty(inspectra_Wrkmem);
hasStoryM=~isempty(inspectra_StoryM);

curspectra_Restin=inspectra_Restin;
curspectra_Motort=inspectra_Motort;
curspectra_Wrkmem=inspectra_Wrkmem;
curspectra_StoryM=inspectra_StoryM;

%curspectra=inspectra;
curfreqs=infreqs;
Nfreqs=length(curfreqs);

%%
range_alpha=[6 14];
range_beta=[15 30];
Xlim_default=[1 50];

peakinfo_alpha=[];
peakinfo_beta=[];



outSaveFile=[outSaveDir,'/',subjID,'_ABpeakfile.mat'];

%%
%========================================================================
% CREATE GUI LAYOUT
parH=[];
parH.Fig =figure('position',[ 13          15        1481         683]);
%-- Main Panels for plotting
parH.panel_plot = uipanel( 'Parent', parH.Fig,'position',[0 0.2 1 0.8]);
parH.panel_fsave = uipanel( 'Parent', parH.Fig,'position',[0 0 1 0.199]);
%-------------------------------------------------
Nsubplots_W=7;
Nsubplots_H=2;
subpan_W=1/Nsubplots_W;
subpan_H=1/Nsubplots_H;

subpan_Wpos=0:subpan_W:1-subpan_W;
subpan_Hpos=1-subpan_H:-subpan_H:0;
%-------------------------------------------------
%===============================
subTitles{1,1}='MOT Post.Raw';
subTitles{1,2}='MOT Post.Relative';
subTitles{1,3}='WM Post.Raw';
subTitles{1,4}='WM Post.Relative';
subTitles{1,5}='STM Post.Raw';
subTitles{1,6}='STM Post.Relative';
subTitles{1,7}='RESTIN Post.Raw';

subTitles=[];
subTitles{2,1}='MOT Ant.Raw';
subTitles{2,2}='MOT Ant.Relative';
subTitles{2,3}='WM Ant.Raw';
subTitles{2,4}='WM Ant.Relative';
subTitles{2,5}='STM Ant.Raw';
subTitles{2,6}='STM Ant.Relative';
subTitles{2,7}='RESTIN Ant.Raw';


%===============================

parH.subplotplanel=[];
for iH=1:Nsubplots_H,
    for iW=1:Nsubplots_W,
        parH.subplotpanel(iH,iW)=uipanel( 'Parent', parH.panel_plot,'position',[subpan_Wpos(iW) subpan_Hpos(iH) subpan_W subpan_H]);
    end
end


parH.subplotaxes=[];
parH.subplottitles=[];
parH.subplotselpeakpos=[];
parH.subplotXhairpos=[];
parH.subplotselpeak=[];
parH.subplotXhair=[];
for iH=1:Nsubplots_H,
    for iW=1:Nsubplots_W,
        parH.subplotaxes(iH,iW)=axes( 'Parent', parH.subplotpanel(iH,iW),'position',[0.1 0.2 0.9 0.7],'ButtonDownFcn',{@getclickposonaxes});
        hold on;
        parH.subplottitles(iH,iW)=title(subTitles{iH,iW});
        parH.subplotselpeakpos{iH,iW}=[0,0];
        parH.subplotXhairpos{iH,iW}=[0,0];
        parH.subplotselpeak(iH,iW)=plot(parH.subplotselpeakpos{iH,iW}(1),parH.subplotselpeakpos{iH,iW}(2),'m+','markersize',30);
        parH.subplotXhair(iH,iW)=plot(parH.subplotXhairpos{iH,iW}(1),parH.subplotXhairpos{iH,iW}(2),'k+','markersize',20);
        
    end
end


%-------------------------------------------------
% Alpha Get Button
Ngetbutton_W=7;
Ngetbutton_H=2;
getbutton_W=0.25;
getbutton_H=0.1;

getbutton_Wpos=0.25; %0+getbutton_W:subpan_W:1-subpan_W+getbutton_W;
getbutton_Hpos=0.01; %01-subpan_H:-subpan_H:0;
%-------------------------------------------------
parH.getbutton_A=[];
for iH=1:Nsubplots_H,
    for iW=1:Nsubplots_W,
        parH.getbutton_A(iH,iW) = uicontrol('Style','pushbutton','Parent',parH.subplotpanel(iH,iW),'Units','normalized','Position',[getbutton_Wpos getbutton_Hpos getbutton_W getbutton_H],'String','Aget','BackgroundColor',[205,92,92]./255,'Callback',{@getA_fromXhair});
    end
end
%-------------------------------------------------
% Beta Get Button
getbutton_Wpos=0.6; %0+getbutton_W:subpan_W:1-subpan_W+getbutton_W;
%-------------------------------------------------
parH.getbutton_B=[];
for iH=1:Nsubplots_H,
    for iW=1:Nsubplots_W,
        parH.getbutton_B(iH,iW) = uicontrol('Style','pushbutton','Parent',parH.subplotpanel(iH,iW),'Units','normalized','Position',[getbutton_Wpos getbutton_Hpos getbutton_W getbutton_H],'String','Bget','BackgroundColor',[100,149,237]./255,'Callback',{@getB_fromXhair});
    end
end


%============================================================
Neditbox_W=7;
Neditbox_H=2;
editbox_W=1/Neditbox_W;
editbox_H=0.6/Neditbox_H;

editbox_Wpos=0:editbox_W:1-editbox_W;
editbox_Hpos=1-editbox_H:-editbox_H:0.4;

parH.editfreq_alpha_panel=[];
for iH=1, %:Neditbox_H,
    for iW=1:Neditbox_W,
        parH.editfreq_alpha_panel(iH,iW)=uipanel( 'Parent', parH.panel_fsave,'position',[editbox_Wpos(iW) editbox_Hpos(iH) editbox_W editbox_H]);
    end
end

parH.editfreq_beta_panel=[];
for iH=2, %:Neditbox_H,
    for iW=1:Neditbox_W,
        parH.editfreq_beta_panel(iH-1,iW)=uipanel( 'Parent', parH.panel_fsave,'position',[editbox_Wpos(iW) editbox_Hpos(iH) editbox_W editbox_H]);
    end
end
%% ====================================
editboxinlet_W=0.6;
editboxinlet_H=0.9;
editboxinlet_Wpos=0.2;
editboxinlet_Hpos=0.05;


parH.editfreq_alpha_inlet=[];
for iH=1, %:Neditbox_H,
    for iW=1:Neditbox_W,
        parH.editfreq_alpha_inlet(iH,iW)=uicontrol('Style','edit','Parent', parH.editfreq_alpha_panel(iH,iW) ,'Units','normalized' ,'position',[editboxinlet_Wpos editboxinlet_Hpos editboxinlet_W editboxinlet_H],'String',num2str(nan));
    end
end


parH.editfreq_beta_inlet=[];
for iH=2, %:Neditbox_H,
    for iW=1:Neditbox_W,
        parH.editfreq_beta_inlet(iH-1,iW)=uicontrol('Style','edit','Parent', parH.editfreq_beta_panel(iH-1,iW) ,'Units','normalized' ,'position',[editboxinlet_Wpos editboxinlet_Hpos editboxinlet_W editboxinlet_H],'String',num2str(nan));
    end
end
%%
%==============================================================================
% Make the button that updates stored list of freqs.
parH.subjIDtext=uicontrol('Style','edit','Parent', parH.panel_fsave ,'Units','normalized' ,'position',[0 0.1 0.09 0.2],'String',subjID);
parH.computepeakbutton = uicontrol('Style','pushbutton','Parent',parH.panel_fsave,'Units','normalized','Position',[0.1 0.1 0.1 0.2],'String','Compute','Callback',{@resetpeakfreqs});
parH.freqrangebutton = uicontrol('Style','pushbutton','Parent',parH.panel_fsave,'Units','normalized','Position',[0.33 0.1 0.1 0.2],'String','Set Range','Callback',{@setallplotfreqrange});
parH.editfreqrangeLow=uicontrol('Style','edit','Parent', parH.panel_fsave ,'Units','normalized' ,'position',[0.23 0.1 0.05 0.2],'String',num2str(Xlim_default(1)));
parH.editfreqrangeHigh=uicontrol('Style','edit','Parent', parH.panel_fsave ,'Units','normalized' ,'position',[0.28 0.1 0.05 0.2],'String',num2str(Xlim_default(2)));
parH.updatebutton = uicontrol('Style','pushbutton','Parent',parH.panel_fsave,'Units','normalized','Position',[0.45 0.1 0.1 0.2],'String','Update','Callback',{@getfromeditboxalphabetapeaks});
parH.exitbutton = uicontrol('Style','pushbutton','Parent',parH.panel_fsave,'Units','normalized','Position',[0.9 0.1 0.1 0.2],'String','EXIT','Callback',{@exitaftersave2file});


parH.textfreqFinalA=uicontrol('Style','text','Parent', parH.panel_fsave ,'Units','normalized' ,'position',[0.63 0.1 0.05 0.2],'String','final A:');
parH.textfreqFinalB=uicontrol('Style','text','Parent', parH.panel_fsave ,'Units','normalized' ,'position',[0.73 0.1 0.05 0.2],'String','final B:');

parH.editfreqFinalA=uicontrol('Style','edit','Parent', parH.panel_fsave ,'Units','normalized' ,'position',[0.68 0.1 0.05 0.2],'String',' ');
parH.editfreqFinalB=uicontrol('Style','edit','Parent', parH.panel_fsave ,'Units','normalized' ,'position',[0.78 0.1 0.05 0.2],'String',' ');


avgspectracell=getavgspectra();
plotavgspectra(avgspectracell);
curavgpeaksraw=computeallpeaks(avgspectracell);
plotavgrawpeaks(curavgpeaksraw);
[peakinfo_alpha,peakinfo_beta]=assignalphabetafromrawpeaks(avgspectracell,curavgpeaksraw);
plotalphabetapeaks(peakinfo_alpha,peakinfo_beta);
putineditboxalphabetapeaks(peakinfo_alpha,peakinfo_beta);


% Initialize the final alpha and beta values with the median values:
alphapeak=[];
betapeak=[];
alphapeak=nanmedian(peakinfo_alpha(3,:));
betapeak=nanmedian(peakinfo_beta(3,:));
set(parH.editfreqFinalA,'String',num2str(alphapeak));
set(parH.editfreqFinalB,'String',num2str(betapeak));

disp('DONE');


%% CREATE FREQ Edit boxes
%{
parH.

parH.freqband_edit_lower = uicontrol('Style','edit','Parent',parH.panel_ctrl_L,'Units','normalized','Position',[0.01 0.3 0.1 0.2],'String',num2str(freqRange(1)));
parH.freqband_edit_upper = uicontrol('Style','edit','Parent',parH.panel_ctrl_L,'Units','normalized','Position',[0.121 0.3 0.1 0.2],'String',num2str(freqRange(2)));
parH.button_freqband_set = uicontrol('Style','pushbutton','Parent',parH.panel_ctrl_L,'Units','normalized','Position',[0.23 0.3 0.2 0.2],'String','set band');
parH.freqband_text_title = uicontrol('Style','text','Parent',parH.panel_ctrl_L,'Units','normalized','Position',[0.01 0.5 0.42 0.2],'String',['Current:  ',num2str(freqRange(1)),' to ',num2str(freqRange(2)),' Hz']);
%======================================================================
%}
%#########################################################################################
%#########################################################################################
%#########################################################################################
%#########################################################################################
%#########################################################################################
%#########################################################################################
%#########################################################################################
%#########################################################################################
%#########################################################################################
%#########################################################################################
%#########################################################################################
%#########################################################################################
%#########################################################################################
%  END OF MAIN FUNCTION - START NESTED FUNCTIONS

    function[avgspectracell]=getavgspectra(hObject,~)
        
        avgspectracell=[];
        dumNanMat=nan(1,Nfreqs);
        for iD1=1:2,
            for iD2=1:7,
                avgspectracell{iD1,iD2}=dumNanMat;
            end
        end
        
        
        if hasMotort,
            % Motort- Raw Post - Avg Motor and Front
            iCase=2; jCase=1;
            curspectra=curspectra_Motort;
            indIn=find((curspectra(:,4)==2)&(ismember(curspectra(:,2),[1]))&(ismember(curspectra(:,3),[2 3 4]))); avgspectracell{iCase,jCase}=mean(curspectra(indIn,5:end),1);
            %------------------
            % Motort- Raw Post - Occ
            iCase=1; jCase=1;
            curspectra=curspectra_Motort;
            indIn=find((curspectra(:,4)==2)&(ismember(curspectra(:,2),[2]))&(ismember(curspectra(:,3),[1]))); avgspectracell{iCase,jCase}=mean(curspectra(indIn,5:end),1);
            %------------------
            % Motort- Ratio - Avg Motor and Front
            iCase=2; jCase=2;
            curspectra=curspectra_Motort;
            indIn=find((curspectra(:,4)==3)&(ismember(curspectra(:,2),[1]))&(ismember(curspectra(:,3),[2 3 4]))); avgspectracell{iCase,jCase}=mean(curspectra(indIn,5:end),1);
            %------------------
            % Motort- Ratio - Occ
            iCase=1; jCase=2;
            curspectra=curspectra_Motort;
            indIn=find((curspectra(:,4)==3)&(ismember(curspectra(:,2),[2]))&(ismember(curspectra(:,3),[1]))); avgspectracell{iCase,jCase}=mean(curspectra(indIn,5:end),1);
            %------------------
        end
        %==============================================================
        if hasWrkmem
            % Wrkmem- Raw Post - Avg Motor and Front
            iCase=2; jCase=3;
            curspectra=curspectra_Wrkmem;
            indIn=find((curspectra(:,4)==2)&(ismember(curspectra(:,2),[2]))&(ismember(curspectra(:,3),[2 3 4]))); avgspectracell{iCase,jCase}=mean(curspectra(indIn,5:end),1);
            %------------------
            % Wrkmem- Raw Post - Occ
            iCase=1; jCase=3;
            curspectra=curspectra_Wrkmem;
            indIn=find((curspectra(:,4)==2)&(ismember(curspectra(:,2),[1]))&(ismember(curspectra(:,3),[1]))); avgspectracell{iCase,jCase}=mean(curspectra(indIn,5:end),1);
            %------------------
            % Wrkmem- Ratio - Avg Motor and Front
            iCase=2; jCase=4;
            curspectra=curspectra_Wrkmem;
            indIn=find((curspectra(:,4)==3)&(ismember(curspectra(:,2),[2]))&(ismember(curspectra(:,3),[2 3 4]))); avgspectracell{iCase,jCase}=mean(curspectra(indIn,5:end),1);
            %------------------
            % Wrkmem- Ratio - Occ
            iCase=1; jCase=4;
            curspectra=curspectra_Wrkmem;
            indIn=find((curspectra(:,4)==3)&(ismember(curspectra(:,2),[1]))&(ismember(curspectra(:,3),[1]))); avgspectracell{iCase,jCase}=mean(curspectra(indIn,5:end),1);
            %------------------
        end
        %==============================================================
        if hasStoryM,
            % StoryM- Raw Post - Avg Motor and Front
            iCase=2; jCase=5;
            curspectra=curspectra_StoryM;
            indIn=find((curspectra(:,4)==2)&(ismember(curspectra(:,2),[1]))&(ismember(curspectra(:,3),[2 3 4]))); avgspectracell{iCase,jCase}=mean(curspectra(indIn,5:end),1);
            %------------------
            % StoryM- Raw Post - Occ
            iCase=1; jCase=5;
            curspectra=curspectra_StoryM;
            indIn=find((curspectra(:,4)==2)&(ismember(curspectra(:,2),[2]))&(ismember(curspectra(:,3),[1]))); avgspectracell{iCase,jCase}=mean(curspectra(indIn,5:end),1);
            %------------------
            % StoryM- Ratio - Avg Motor and Front
            iCase=2; jCase=6;
            curspectra=curspectra_StoryM;
            indIn=find((curspectra(:,4)==3)&(ismember(curspectra(:,2),[1]))&(ismember(curspectra(:,3),[2 3 4]))); avgspectracell{iCase,jCase}=mean(curspectra(indIn,5:end),1);
            %------------------
            % StoryM- Ratio - Occ
            iCase=1; jCase=6;
            curspectra=curspectra_StoryM;
            indIn=find((curspectra(:,4)==3)&(ismember(curspectra(:,2),[2]))&(ismember(curspectra(:,3),[1]))); avgspectracell{iCase,jCase}=mean(curspectra(indIn,5:end),1);
            %------------------
        end
        %==============================================================
        if hasRestin
            % Restin- Raw Post - Avg Motor and Front
            iCase=2; jCase=7;
            curspectra=curspectra_Restin;
            indIn=find((ismember(curspectra(:,3),[2 3 4]))); avgspectracell{iCase,jCase}=mean(curspectra(indIn,5:end),1);
            %------------------
            % Restin- Raw Post - Occ
            iCase=1; jCase=7;
            curspectra=curspectra_Restin;
            indIn=find((ismember(curspectra(:,3),[1]))); avgspectracell{iCase,jCase}=mean(curspectra(indIn,5:end),1);
            %------------------
        end
        
    end
%===============================================================================
%===============================================================================
%===============================================================================


    function[]=plotavgspectra(plotspectracell)
        for iH=1:Nsubplots_H,
            for iW=1:Nsubplots_W,
                axes(parH.subplotaxes(iH,iW));
                hold on;
                plot(curfreqs,plotspectracell{iH,iW},'k');
            end
        end
        
    end
%===============================================================================
%===============================================================================
%===============================================================================


    function[outavgpeaksraw]=computeallpeaks(inavgspectracell)
        
        outavgpeaksraw=[];
        allFreqs=curfreqs;
        for iW=1:Nsubplots_W,
            for iH=1:2,
                
                outavgpeaksraw{iH,iW}=[];
                
                
                tmpSpec=inavgspectracell{iH,iW};
                if ~isnan(tmpSpec(1))
                    
                    if rem(iW,2)==1,
                        xlog=log(allFreqs([1:3 end-15:end]));
                        ylog=log(tmpSpec([1:3 end-15:end]));
                        linmodel=polyfit(xlog,ylog,1);
                        allx=allFreqs;
                        allxlog=log(allFreqs);
                        allylog=log(tmpSpec);
                        modylog=linmodel(2)+linmodel(1)*allxlog;
                        remy=exp(allylog-modylog);
                        
                    else
                        allx=allFreqs;
                        ally=-(tmpSpec-mean(tmpSpec));%./std(tmpSpec);
                        remy=ally;
                    end
                    
                    meandif=mean(abs(diff(remy)));
                    stddif=std(abs(diff(remy)));
                    thresdif=meandif;
                    [PKS,LOCS]= findpeaks(remy);
                    
                    signPeakLocs=[];
                    countSign=1;
                    for iPeak=1:length(LOCS),
                        if (LOCS(iPeak)>4)&(LOCS(iPeak)<length(allx)-3),
                            curDist=(abs(diff(remy(LOCS(iPeak)+[0 -3])))+abs(diff(remy(LOCS(iPeak)+[0 3]))))./2;
                            if curDist>thresdif,
                                signPeakLocs(countSign)=LOCS(iPeak);
                                countSign=countSign+1;
                            end
                        end;
                    end
                    outavgpeaksraw{iH,iW}=[signPeakLocs; curfreqs(signPeakLocs) ;  tmpSpec(signPeakLocs)  ];
                    
                    
                end
            end
        end
        
        disp('peaks computed...')
    end
%===============================================================================
%===============================================================================
%===============================================================================
    function[]=plotavgrawpeaks(inavgpeaksraw)
        
        for iW=1:Nsubplots_W,
            for iH=1:2,
                axes(parH.subplotaxes(iH,iW));
                
                hold on;
                if ~isempty(inavgpeaksraw{iH,iW})
                    curPeakIndx=inavgpeaksraw{iH,iW}(1,:);
                    curPeakFreq=inavgpeaksraw{iH,iW}(2,:);
                    curPeakVals=inavgpeaksraw{iH,iW}(3,:);
                    
                    NcurPeaks=length(curPeakVals);
                    if NcurPeaks>0
                        
                        for iPeak=1:NcurPeaks,
                            plot(curPeakFreq(iPeak),curPeakVals(iPeak),'.b','markersize',5);
                        end
                        
                    end
                end
            end
        end
        
    end
%===============================================================================
%===============================================================================
%===============================================================================
    function[totalpeakinfo_alpha,totalpeakinfo_beta]=assignalphabetafromrawpeaks(inavgspectracell,inavgpeaksraw)
        
        totalpeakinfo_alpha=[];
        totalpeakinfo_beta=[];
        
        for iW=1:Nsubplots_W,
            
            tmpPeakInfo_A=[nan(9,1) ;0];
            tmpPeakInfo_B=[nan(9,1) ;0];
            tmppeakinfo_alpha=[];
            tmppeakinfo_beta=[];
            curPeakIndx=[];
            curPeakFreq=[];
            curPeakVals=[];
            curPeakValsZ=[];
            NcurPeaks=[];
            
            
            
            
            for iH=1:2
                tmpSpec=inavgspectracell{iH,iW};
                tmpSpecMean=mean(tmpSpec);
                tmpSpecStd=std(tmpSpec);
                if isempty(inavgpeaksraw{iH,iW}),
                    tmppeakinfo_alpha=tmpPeakInfo_A;
                    tmppeakinfo_beta=tmpPeakInfo_B;
                else
                    curPeakIndx=inavgpeaksraw{iH,iW}(1,:);
                    curPeakFreq=inavgpeaksraw{iH,iW}(2,:);
                    curPeakVals=inavgpeaksraw{iH,iW}(3,:);
                    
                    curPeakValsZ=(curPeakVals-tmpSpecMean)./tmpSpecStd;
                    if rem(iW,2)==0
                        curPeakValsZ=-curPeakValsZ;
                    end
                    NcurPeaks=length(curPeakVals);
                    
                    %====================================================
                    
                    %-------------------------------------------------------
                    indRin=find((curPeakFreq>=range_alpha(1))&(curPeakFreq<=range_alpha(2)));
                    Nrangepeaks=length(indRin);
                    if Nrangepeaks==1,
                        tmpPeakInfo_A(iH,1)=curPeakFreq(indRin);
                        tmpPeakInfo_A(iH+3,1)=curPeakVals(indRin);
                        tmpPeakInfo_A(iH+6,1)=curPeakValsZ(indRin);
                    elseif Nrangepeaks>1,
                        multiPfreq_a=curPeakFreq(indRin);
                        multiPvals_a=curPeakVals(indRin);
                        multiPvalsZ_a=curPeakValsZ(indRin);
                        [maxV,maxInd]=max(multiPvalsZ_a);
                        tmpPeakInfo_A(iH,1)=multiPfreq_a(maxInd);
                        tmpPeakInfo_A(iH+3,1)=multiPvals_a(maxInd);
                        tmpPeakInfo_A(iH+6,1)=multiPvalsZ_a(maxInd);
                    end
                    tmppeakinfo_alpha=tmpPeakInfo_A;
                    %-------------------------------------------------------
                    
                    indRin=find((curPeakFreq>=range_beta(1))&(curPeakFreq<=range_beta(2)));
                    Nrangepeaks=length(indRin);
                    if Nrangepeaks==1,
                        tmpPeakInfo_B(iH,1)=curPeakFreq(indRin);
                        tmpPeakInfo_B(iH+3,1)=curPeakVals(indRin);
                        tmpPeakInfo_B(iH+6,1)=curPeakValsZ(indRin);
                    elseif Nrangepeaks>1,
                        multiPfreq_a=curPeakFreq(indRin);
                        multiPvals_a=curPeakVals(indRin);
                        multiPvalsZ_a=curPeakValsZ(indRin);
                        [maxV,maxInd]=max(multiPvalsZ_a);
                        tmpPeakInfo_B(iH,1)=multiPfreq_a(maxInd);
                        tmpPeakInfo_B(iH+3,1)=multiPvals_a(maxInd);
                        tmpPeakInfo_B(iH+6,1)=multiPvalsZ_a(maxInd);
                    end
                    tmppeakinfo_beta=tmpPeakInfo_B;
                end
            end
            %-------------------------
            for iRangeCase=1:2
                if iRangeCase==1,
                    tmpPeakInfo=tmppeakinfo_alpha;
                else
                    tmpPeakInfo=tmppeakinfo_beta;
                end
                % Get the freq for each Case
                allpeakfreqs=tmpPeakInfo(1:2,1);
                allpeakvals=tmpPeakInfo(4:5,1);
                allpeakvalsZ=tmpPeakInfo(7:8,1);
                [maxV,maxInd]=max(allpeakvalsZ);
                tmpPeakInfo(3,1)=allpeakfreqs(maxInd);
                tmpPeakInfo(6,1)=allpeakvals(maxInd);
                tmpPeakInfo(9,1)=allpeakvalsZ(maxInd);
                if iRangeCase==1,
                    tmppeakinfo_alpha=tmpPeakInfo;
                else
                    tmppeakinfo_beta=tmpPeakInfo;
                end
                %-------------------------
                
            end
            
            %----------------------
            totalpeakinfo_alpha(:,iW)=tmppeakinfo_alpha;
            totalpeakinfo_beta(:,iW)=tmppeakinfo_beta;
            
        end
        
        
    end
%===============================================================================
%===============================================================================
%===============================================================================
    function[]=plotalphabetapeaks(inpeakinfo_alpha,inpeakinfo_beta)
        
        for iW=1:Nsubplots_W,
            curSelPeakFreq_A=inpeakinfo_alpha(3,iW);
            curSelPeakFreq_B=inpeakinfo_beta(3,iW);
            
            %curSelPeakVal_A=inpeakinfo_alpha(4+4,iW);
            %curSelPeakVal_B=inpeakinfo_beta(4+4,iW);
            
            parH.subplotselpeakpos{1,iW}=[curSelPeakFreq_A inpeakinfo_alpha(4,iW)];
            parH.subplotselpeakpos{2,iW}=[curSelPeakFreq_A inpeakinfo_alpha(5,iW)];
            
            iH=1;
            curPeakFreq_A=inpeakinfo_alpha(iH,iW);
            curPeakVal_A=inpeakinfo_alpha(iH+3,iW);
            curPeakFreq_B=inpeakinfo_beta(iH,iW);
            curPeakVal_B=inpeakinfo_beta(iH+3,iW);
            
            axes(parH.subplotaxes(iH,iW));hold on;
            plot(curPeakFreq_A,curPeakVal_A,'o','color',[205,92,92]./255,'markersize',10);
            plot(curPeakFreq_B,curPeakVal_B,'o','color',[100,149,237]./255,'markersize',10);
            set(parH.subplotselpeak(1,iW),'XData',curSelPeakFreq_A,'YData',parH.subplotselpeakpos{1,iW}(2));
            text(curPeakFreq_A,curPeakVal_A,num2str(curPeakFreq_A),'color',[0 0 0]./255,'Fontweight','Bold');
            text(curPeakFreq_B,curPeakVal_B,num2str(curPeakFreq_B),'color',[0 0 0]./255,'Fontweight','Bold');
            
            
            iH=2;
            curPeakFreq_A=inpeakinfo_alpha(iH,iW);
            curPeakVal_A=inpeakinfo_alpha(iH+3,iW);
            curPeakFreq_B=inpeakinfo_beta(iH,iW);
            curPeakVal_B=inpeakinfo_beta(iH+3,iW);
            
            axes(parH.subplotaxes(iH,iW));hold on;
            plot(curPeakFreq_A,curPeakVal_A,'o','color',[205,92,92]./255,'markersize',10);
            plot(curPeakFreq_B,curPeakVal_B,'o','color',[100,149,237]./255,'markersize',10);
            set(parH.subplotselpeak(2,iW),'XData',curSelPeakFreq_A,'YData',parH.subplotselpeakpos{2,iW}(2));
            text(curPeakFreq_A,curPeakVal_A,num2str(curPeakFreq_A),'color',[0 0 0]./255,'Fontweight','Bold');
            text(curPeakFreq_B,curPeakVal_B,num2str(curPeakFreq_B),'color',[0 0 0]./255,'Fontweight','Bold');
            
        end
        
        
    end
%===============================================================================
%===============================================================================
%===============================================================================
    function[]=putineditboxalphabetapeaks(inpeakinfo_alpha,inpeakinfo_beta)
        
        for iW=1:Nsubplots_W,
            curSelPeakFreq_A=inpeakinfo_alpha(3,iW);
            curSelPeakFreq_B=inpeakinfo_beta(3,iW);
            set(parH.editfreq_alpha_inlet(1,iW),'String',num2str(curSelPeakFreq_A));
            set(parH.editfreq_beta_inlet(1,iW),'String',num2str(curSelPeakFreq_B));
            
        end
        
        
    end
%===============================================================================
%===============================================================================
%===============================================================================
    function[]=getfromeditboxalphabetapeaks(hObject,~)
        origButColor=get(parH.updatebutton,'Backgroundcolor');
        origButStr=get(parH.updatebutton,'String');
        for iW=1:Nsubplots_W,
            
            peakinfo_alpha(3,iW)=str2num(get(parH.editfreq_alpha_inlet(1,iW),'String'));
            peakinfo_beta(3,iW)=str2num(get(parH.editfreq_beta_inlet(1,iW),'String'));
        end
        set(parH.updatebutton,'String','...getting in vars','BackgroundColor',[1 0 0]);
        pause(1);
        set(parH.updatebutton,'String',origButStr,'BackgroundColor',origButColor);
        
    end
%===============================================================================
%===============================================================================
%===============================================================================
    function getclickposonaxes(hObject,~)
        [indH,indW]=find(ismember(parH.subplotaxes,hObject));
        %axIndx=find(ismember([[parH.panel_plot_L_Axes{:}] [parH.panel_plot_R_Axes{:}] [parH.panel_plot_A_Axes{:}] [parH.panel_plot_S_Axes{:}]],hObject));
        pos=get(hObject,'CurrentPoint');
        disp(['You clicked X:',num2str(pos(1,1)),', Y:',num2str(pos(1,2)), ': axes:', num2str(indH), '-', num2str(indW)]);
        
        parH.subplotXhairpos{indH,indW}=[round(pos(1,1)) pos(1,2)];
        set(parH.subplotXhair(indH,indW),'XData',round(pos(1,1)),'YData',pos(1,2));
        if indH==1
            acindH=2;
        else
            acindH=1;
        end
        set(parH.subplotXhair(acindH,indW),'XData',round(pos(1,1)),'YData',parH.subplotXhairpos{acindH,indW}(1,2));
        parH.subplotXhairpos{acindH,indW}(1)=round(pos(1,1));
    end
%===============================================================================
%===============================================================================
%===============================================================================
    function resetpeakfreqs(hObject,~)
        for iH=1:Nsubplots_H,
            for iW=1:Nsubplots_W,
                axes(parH.subplotaxes(iH,iW));
                hold off;
                parH.subplotselpeakpos{iH,iW}=[0,0];
                parH.subplotXhairpos{iH,iW}=[0,0];
                parH.subplotselpeak(iH,iW)=plot(parH.subplotselpeakpos{iH,iW}(1),parH.subplotselpeakpos{iH,iW}(2),'m+','markersize',30);
                hold on
                parH.subplotXhair(iH,iW)=plot(parH.subplotXhairpos{iH,iW}(1),parH.subplotXhairpos{iH,iW}(2),'k+','markersize',20);
            end
        end
        
        avgspectracell=getavgspectra();
        plotavgspectra(avgspectracell);
        curavgpeaksraw=computeallpeaks(avgspectracell);
        plotavgrawpeaks(curavgpeaksraw);
        [peakinfo_alpha,peakinfo_beta]=assignalphabetafromrawpeaks(avgspectracell,curavgpeaksraw);
        plotalphabetapeaks(peakinfo_alpha,peakinfo_beta);
        putineditboxalphabetapeaks(peakinfo_alpha,peakinfo_beta);
    end
%===============================================================================
%===============================================================================
%===============================================================================
    function getA_fromXhair(hObject,~)
        
        [indH,indW]=find(ismember(parH.getbutton_A,hObject));
        
        if indH==1,
            multiindH=[1 2];
        else
            multiindH=[2 1];
        end
        for iSCase=1:2,
            
            manPeakFreq=parH.subplotXhairpos{multiindH(iSCase),indW}(1);
            indManFreq=nearest(curfreqs,manPeakFreq);
            
            tmpSpec=[];
            tmpSpec(1,:)=avgspectracell{multiindH(iSCase),indW};
            
            manPeakVal=tmpSpec(indManFreq);
            manPeakValZ=(manPeakVal-mean(tmpSpec))./std(tmpSpec);
            if rem(indW,2)
                manPeakValZ=-manPeakValZ;
            end
            
            
            set(parH.subplotselpeak(multiindH(iSCase),indW),'XData',manPeakFreq,'YData',manPeakVal);
            
            
            peakinfo_alpha(3,indW)=manPeakFreq;
            peakinfo_alpha(6,indW)=manPeakVal;
            peakinfo_alpha(9,indW)=manPeakValZ;
            peakinfo_alpha(10,indW)=1;
        end
        %plotalphabetapeaks(peakinfo_alpha,peakinfo_beta);
        putineditboxalphabetapeaks(peakinfo_alpha,peakinfo_beta);
    end
%===============================================================================
%===============================================================================
%===============================================================================
    function getB_fromXhair(hObject,~)
        
        [indH,indW]=find(ismember(parH.getbutton_B,hObject));
        
        if indH==1,
            multiindH=[1 2];
        else
            multiindH=[2 1];
        end
        for iSCase=1:2,
            
            manPeakFreq=parH.subplotXhairpos{multiindH(iSCase),indW}(1);
            indManFreq=nearest(curfreqs,manPeakFreq);
            
            tmpSpec=[];
            tmpSpec(1,:)=avgspectracell{multiindH(iSCase),indW};
            
            manPeakVal=tmpSpec(indManFreq);
            manPeakValZ=(manPeakVal-mean(tmpSpec))./std(tmpSpec);
            if rem(indW,2)
                manPeakValZ=-manPeakValZ;
            end
            
            
            %set(parH.subplotselpeak(multiindH(iSCase),indW),'XData',manPeakFreq,'YData',manPeakVal);
            
            
            peakinfo_beta(3,indW)=manPeakFreq;
            peakinfo_beta(6,indW)=manPeakVal;
            peakinfo_beta(9,indW)=manPeakValZ;
            peakinfo_beta(10,indW)=1;
        end
        %plotalphabetapeaks(peakinfo_alpha,peakinfo_beta);
        putineditboxalphabetapeaks(peakinfo_alpha,peakinfo_beta);
    end

%===============================================================================
%===============================================================================
%===============================================================================
    function setallplotfreqrange(hObject,~)
        tmpXlimLow=str2num(get(parH.editfreqrangeLow,'String'));
        tmpXlimHigh=str2num(get(parH.editfreqrangeHigh,'String'));
        
        for iH=1:Nsubplots_H,
            for iW=1:Nsubplots_W,
                set(parH.subplotaxes(iH,iW),'Xlim',[tmpXlimLow tmpXlimHigh]);
            end
        end
        
    end

%===============================================================================
%===============================================================================
%===============================================================================
    function[]=exitaftersave2file(hObject,~)
        
        alphapeak=str2num(get(parH.editfreqFinalA,'String'));
        betapeak=str2num(get(parH.editfreqFinalB,'String'));

        save(outSaveFile,'alphapeak','betapeak','peakinfo_alpha', 'peakinfo_beta');
        
        origButColor=get(parH.exitbutton,'Backgroundcolor');
        origButStr=get(parH.exitbutton,'String');
        
        set(parH.exitbutton,'String','...saving in file','BackgroundColor',[1 0 0]);
        pause(1);
        set(parH.exitbutton,'String',origButStr,'BackgroundColor',origButColor);
        pause(1);
        close all;
        
        
    end

end

