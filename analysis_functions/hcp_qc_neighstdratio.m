function [badchannel_neighstd] = hcp_qc_neighstdratio(datain, options_std)

% This function computes for each channel the ratio [Sensor std- Neighbour
% Std]./Neighbour Std. Channels with high variance relative to the sensors tend to
% cause some Indpendend Components in ICA to tru to describe the single sensor
% higher variance. In order to eliminate such cases a threshold is used
% empirically. A value of 0.5 seems to be adequate for identifying such sensors.
%
% INPUTS:
%     datain : Data structure containing trials
%
%     stdratiothreshold:  The threshold on of the std ratio with naighbours
%
% OUPUTS:
%     badchannel_neighstd : These are the labels of bad sensors according
%                           to std ratio

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

resultprefix   = ft_getopt(options_std, 'resultprefix');
stdratiothreshold     = ft_getopt(options_std, 'stdratiothreshold');        % skipped intervals previously identified


%=======================================
%-- Remove line noise --------------------
tok = tokenize(resultprefix, '_');
subjID=tok{1};
if strcmp(subjID,'CP10128')|strcmp(subjID,'CP10129')
    linefreq=[50 100];
else
    linefreq=[60 120];
end
%{
cfg=[];
cfg.detrend = 'no';
cfg.dftfilter = 'yes';
cfg.dftfreq = linefreq;
datain = ft_preprocessing(cfg,datain);
%}

cfg=[];
cfg.bsfilter = 'yes';
cfg.bsfreq = linefreq(1)+[-1 1];
datain = ft_preprocessing(cfg,datain);
cfg.bsfreq = linefreq(2)+[-1 1];
datain = ft_preprocessing(cfg,datain);

%=======================================


% =======================================
%--- Select MEG data
%channelData=ft_selectdata(datain,'channel',{'MEG','-AF*'});
channelData=datain;
totalData=[channelData.trial{:}];
stdTotal=std(totalData,[],2);
Nchannels=size(totalData,1);
%====================================


% =======================================
%--- COMPUTE NEIGHBOURS ------------
cfg = [];
cfg.method='distance';
cfg.grad       = channelData.grad;

cfg.neighbourdist = 0.037; %second order neighbours
layout='4D248.mat';

neighLabels         = ft_prepare_neighbours(cfg);
%--------------------------------------
tmpNeighLabels=struct2cell(neighLabels);
neighIndices=cell(Nchannels,1);
for iChan = 1:Nchannels
    baseChan=channelData.label{iChan};
    [indx1,indx2] = match_str(tmpNeighLabels(1,1,:), baseChan);
    tmpNeighChannels=tmpNeighLabels{2,1,indx1};
    [indx1,indx2] = match_str(channelData.label,tmpNeighChannels);
    tmpNeighIndx=indx1';
    neighIndices{iChan}=tmpNeighIndx;
end
%====================================

%====================================
%------- COMPUTE STD RATIO
stdRatio=zeros(1,Nchannels);
for iChan=1:Nchannels,
    baseStd=stdTotal(iChan);
    tmpNeighStd=stdTotal(neighIndices{iChan});
    meanNeighStd=mean(tmpNeighStd);
    stdRatio(iChan)=(baseStd-meanNeighStd)./meanNeighStd;
end

meanstdratio=mean(stdRatio);
stdstdratio=std(stdRatio);

%indBad=find(abs(stdRatio)>stdratiothreshold);
indBad=find((abs(stdRatio)>stdratiothreshold)&(abs(stdRatio-meanstdratio)>= 3*stdstdratio));


if ~isempty(indBad)
    badchannel_neighstd=channelData.label(indBad)';
else
    badchannel_neighstd={};
end

% badChan=[];
% badChan.stdratiothreshold=stdratiothreshold;
% if ~isempty(indBad)
%     badChan.label=channelData.label(indBad);
%
% else
%     badChan.label=[];
% end
% badChan.stdRatio=stdRatio;

%==========================================

%==========================================
%----- SAVE AND PLOT
%----- simple plot
f1=figure;
plot(stdRatio,'.')
hold on
plot(indBad,stdRatio(indBad),'.r','markersize',20);
legend('all','bad');
plot(stdRatio);
for iB=1:length(indBad),
    text(indBad(iB),stdRatio(indBad(iB)),channelData.label(indBad(iB)))
end
title('Std Difference Ratio with neighbours - red is bad');
xlabel('# Channel')
%---------------------
%--------------------------------
%----- TOPOPLOT of std ratio with bad sensors marked
[indx1,indx2]=match_str(channelData.grad.label, channelData.label);
sens3Dpos=channelData.grad.chanpos(indx1,:);
cfg=[];
cfg.layout=layout;
sensLay=ft_prepare_layout(cfg);
tmpLay=sensLay;


% r=0.1;
% ang=[0:0.01:2*pi]';
% xadd=r*cos(ang);
% yadd=r*sin(ang);
% x_left= -0.384100425519062+xadd;
% y_left=  0.345302231087568+yadd;
% x_right=  0.389749219378296+xadd;
% y_right=  0.345302231087568+yadd;
% tmpLay.mask{2}=[x_left y_left];
% tmpLay.mask{3}=[x_right y_right];


[indx1topo,indx2topo]=match_str( channelData.label,tmpLay.label);
tmpLay.pos=tmpLay.pos(indx2topo,:);
tmpLay.width=tmpLay.width(indx2topo);
tmpLay.height=tmpLay.height(indx2topo);
tmpLay.label=tmpLay.label(indx2topo);
%---------------------------
chanX=tmpLay.pos(:,1);
chanY=tmpLay.pos(:,2);
datavector=abs(stdRatio);

f2=figure; hold on;
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
%%caxis([0.3 0.8]);
colorbar
title('Std Ratio with neighbours: cross at bad channels ');




fname1=[resultprefix,'_badchan_std_scatter.png'];
fname2=[resultprefix,'_badchan_std_topo.png'];

hcp_write_figure(fname1, f1);
hcp_write_figure(fname2, f2);

% zip(zipimagefile,{fname1,fname2});
% delete(fname1)
% delete(fname2)

close(f1);
close(f2);

