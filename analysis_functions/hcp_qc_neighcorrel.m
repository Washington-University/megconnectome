function [badchannel_neighcorr, dataclean] = hcp_qc_neighcorrel (dataraw,options_correl)

% This function loads the entire raw MEG data and check for channels that have low
% correlation with their neighbours. The average correlation with neighbours
% within a distance of  0.037m is between 0.6 and 0.8. So the user here can set a
% threshold for identifying such "bad" sensors. The reason for using this approach
% is that "bad" sensor is one that has inherent noise uncorrelated with other
% sensors. Typically the same sensors should show the lowest correlation with
% neighbours for different datasets. In order to avoid cases where the average
% correlation of the entire sensor array is low (when subject's head is away from
% the sensor) an additional threshold is used. This is the mean neighbour
% correlation from all channels minus 3 standard deviations of the neighbour
% correlation of all channels. So in order a bad channel to be selected , its
% neighbour correlation should also lower than this additional threshold.
%
% INPUTS:
%     dataraw : This is the raw data structure. It should contain the
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
%               corrthreshold:  The threshold on correlation with neighbours for each
%                                       sensor, below which channels are
%                                       identified as bad. An
%                                      emprirically sensible value is 0.4.
%
% OUPUTS:
%               badchannel_neighcorr: : These are the labels of bad sensors sorted according
%                                       to the corelation ascending
%               dataclean11  : This is a data structure , identical to the input dataset dataraw but
%                             with the bad channels identified from the neighbour correlation metric discarded

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

resultprefix   = ft_getopt(options_correl, 'resultprefix');
corrthreshold  = ft_getopt(options_correl, 'corrthreshold');        % skipped intervals previously identified

Fsample=dataraw.fsample;


if length(dataraw.trial)==1,
    disp(['input dataset has a single trial - spliting in pseudo trials']);
    
    nSamples=size(dataraw.trial{1},2);
    defTrDur=1; %sec
    defTrSamps=floor(defTrDur.*Fsample);
    Ntrials=floor(nSamples./defTrSamps);
    trlStart=defTrSamps.*((0:Ntrials-1)')+1;
    trlEnd=defTrSamps.*((1:Ntrials)');
    trlDur=repmat(0,Ntrials,1);
    trl=[trlStart trlEnd trlDur];
    cfg          = [];
    cfg.trl=trl;
    cfg.headerformat='4d';
    cfg.dataformat='4d';
    dataclean1=ft_redefinetrial(cfg,dataraw);
    datachannels=ft_channelselection('MEG',dataclean1.label);
    
    
else
    error(['input dataset should have the entire data set in a single trial']);
    return;
end


cfg          = [];
cfg.detrend   = 'yes';
cfg.feedback = 'none';
cfg.headerformat='4d';
cfg.dataformat='4d';
dataclean1=ft_preprocessing(cfg,dataclean1);




Nchannels=length(dataclean1.label);
%==================================================
%--COMPUTE CORRELATION OF ALL CHANNELS WITH ALL CHANNELS
trialsCorr = cellfun(@(x) corrcoef(x'), dataclean1.trial, 'UniformOutput', false);
for iTr=1:Ntrials,
    trialsCorr{iTr}(isnan(trialsCorr{iTr}))=0;
end
%====================================
%-- COMPUTE NEIGHBOURS

cfg = [];
cfg.method='distance';
cfg.grad       = dataraw.grad;
cfg.neighbourdist = 0.037; %second order neighbours
layout='4D248.mat';

neighLabels         = ft_prepare_neighbours(cfg);
%--------------
allCorrelLabels=dataclean1.label;
tmpNeighLabels=struct2cell(neighLabels);
neighIndices=cell(Nchannels,1);
for iChan = 1:Nchannels
    baseChan=dataclean1.label{iChan};
    [indx1,indx2] = match_str(tmpNeighLabels(1,1,:), baseChan);
    tmpNeighChannels=tmpNeighLabels{2,1,indx1};
    [indx1,indx2] = match_str(dataclean1.label,tmpNeighChannels);
    tmpNeighIndx=indx1';
    neighIndices{iChan}=tmpNeighIndx;
end
%=============================================
%-- COMPUTE AVERAGE CORRELATION WITH NEIGHBOURS
[trialsCorrAv]=cellfun(@compute_avercorr, trialsCorr,repmat({neighIndices},1,Ntrials), 'UniformOutput', false);
trialsCorrAv=cell2mat(trialsCorrAv);
avgCorrPerChan=(mean(trialsCorrAv,2));
%=============================================
%----- USE THRESHOLD TO SELECT BAD CHANNELS
tmpThres=mean(avgCorrPerChan(:))-3*std(avgCorrPerChan);
[indBad]=find((avgCorrPerChan<tmpThres)&avgCorrPerChan<corrthreshold);
%====================================================================
%----- ASSIGN OUTPUT VARIABLES
dataclean2= dataclean1; % It seems to be more memory efficient this way;
if ~isempty(indBad)
    badchannel_neighcorr=allCorrelLabels(indBad)';
    badChanCorrelLabels=cellfun(@(x) ['-',x],badchannel_neighcorr,'UniformOutput',false)';
    keepChans=datachannels;
    dataclean1=ft_selectdata(dataraw,'channel',[keepChans; badChanCorrelLabels]);
    
else
    badchannel_neighcorr={};
    dataclean1=ft_selectdata(dataraw,'channel',datachannels);
end
clear dataclean2;
%====================================================================
%===== OUTPUT AND PLOTTING==========================================
disp('PLOTTING');
[indx1,indx2]=match_str(dataclean1.grad.label, allCorrelLabels);
sens3Dpos=dataclean1.grad.chanpos(indx1,:);
cfg=[];
cfg.layout=layout;
sensLay=ft_prepare_layout(cfg);
tmpLay=sensLay;
%--------------------------------
%----- TOPOPLOT of correlation with bad sensors marked
%{
r=0.1;
ang=[0:0.01:2*pi]';
xadd=r*cos(ang);
yadd=r*sin(ang);
x_left= -0.384100425519062+xadd;
y_left=  0.345302231087568+yadd;
x_right=  0.389749219378296+xadd;
y_right=  0.345302231087568+yadd;
tmpLay.mask{2}=[x_left y_left];
tmpLay.mask{3}=[x_right y_right];
%}
%--------------------------------
[indx1topo,indx2topo]=match_str(allCorrelLabels,tmpLay.label);
tmpLay.pos=tmpLay.pos(indx2topo,:);
tmpLay.width=tmpLay.width(indx2topo);
tmpLay.height=tmpLay.height(indx2topo);
tmpLay.label=tmpLay.label(indx2topo);

chanX=tmpLay.pos(:,1);
chanY=tmpLay.pos(:,2);
datavector=avgCorrPerChan;

f1=figure; hold on;
ft_plot_lay(tmpLay, 'box', 'off');
ft_plot_topo(chanX,chanY,datavector,...
    'interpmethod','v4',...
    'interplim','mask',...
    'gridscale',150,...
    'outline',tmpLay.outline,...
    'shading','flat',...
    'isolines',6,...
    'mask',tmpLay.mask,...
    'style','surfiso');
axis equal;
hold on;
plot(tmpLay.pos(indBad,1),tmpLay.pos(indBad,2),'.k','markersize',30);
plot(tmpLay.pos(indBad,1),tmpLay.pos(indBad,2),'.w','markersize',15);
plot(tmpLay.pos(indBad,1),tmpLay.pos(indBad,2),'+r','markersize',20);
%axis([-0.6 0.6 -0.6 0.6]);
%axis off;
caxis([0.3 0.8]);
colorbar
title('Average Cor. with neighbours: cross at bad channels ');
disp('PLOTTED TOPO');

%----- 3D plot
f2=figure;
ft_plot_topo3d(sens3Dpos,avgCorrPerChan,'contourstyle','color');
hold on
plot3(sens3Dpos(indBad,1),sens3Dpos(indBad,2),sens3Dpos(indBad,3),'.k','markersize',30);
plot3(sens3Dpos(indBad,1),sens3Dpos(indBad,2),sens3Dpos(indBad,3),'+r','markersize',20);
caxis([0.3 0.8]);
colorbar
for iB=1:length(indBad),
    text(sens3Dpos(indBad(iB),1),sens3Dpos(indBad(iB),2),sens3Dpos(indBad(iB),3),allCorrelLabels(indBad(iB)))
end
view(270,90)
title('Average Cor. with neighbours 3D: cross at bad channels');
disp('PLOTTED 3D');

%----- SCATTER plot
f3=figure;
plot(avgCorrPerChan,'.')
hold on
plot(indBad,avgCorrPerChan(indBad),'.r','markersize',20);
legend('all','bad');
plot(avgCorrPerChan);
for iB=1:length(indBad),
    text(indBad(iB),avgCorrPerChan(indBad(iB)),allCorrelLabels(indBad(iB)))
end
title('Average Cor. with neighbours - red is bad');
xlabel('# Channel')
%axis([1 Nchannels 0 1])
disp('PLOTTED SIMPLE');

fname1=[resultprefix,'_badchan_cor_topo.png'];
fname2=[resultprefix,'_badchan_cor_topo3D.png'];
fname3=[resultprefix,'_badchan_cor_scatter.png'];

disp('SAVING');
hcp_write_figure(fname3, f3);
hcp_write_figure(fname2, f2);
hcp_write_figure(fname1, f1);
disp('SAVED');

%zip(zipimagefile,{fname1,fname2,fname3});
%delete(fname1)
%delete(fname2)
%delete(fname3)

close(f1);
close(f2);
close(f3);

%{
%--- topoplot
figure; hold on;
ft_plot_lay(sensLay, 'box', 'off');
ft_plot_topo(sensLay.pos(1:248,1),sensLay.pos(1:248,2),stdPerChan,'gridscale',150,'outline',sensLay.outline,'mask',sensLay.mask,'interpmethod','nearest');
title('Variance of Channels');
%}
%--- END PLOTTING
%=====================================================================
% badCorr=zCorrPerChan(indBad);
% [badCorrAscend,indxAscend]=sort(badCorr,1,'descend');

% badCorr=avgCorrPerChan(indBad);
% [badCorrAscend,indxAscend]=sort(badCorr,1,'ascend');
%
% badCh=[];
% badCh.threshold=badCorrThresh;
% badCh.ind=indBad(indxAscend);
% badCh.label=dataclean1.label(indBad);
% badCh.label=badCh.label(indxAscend);
% badCh.corr=badCorrAscend;

dataclean=dataclean1;
clearvars -except badchannel_neighcorr dataclean

end


%====================================================================
%====================================================================
%=====================================================================
%=====================END OF MAIN FUNCTION ====================================
%=====================================================================
%=====================================================================
%====================================================================
function[averCorr]=compute_avercorr(trialCorr,neighChanIndx)

Nchan=size(trialCorr,1);
averCorr=zeros(Nchan,1);
for iCh=1:Nchan,
    tmpCorr=(trialCorr(iCh,neighChanIndx{iCh}));
    tmpCorr(isnan(tmpCorr))=[];
    averCorr(iCh)=mean(tmpCorr);
end

end
