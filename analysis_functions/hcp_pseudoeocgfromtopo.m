function [pseudodataECG,pseudodataEOG] = hcp_pseudoeocgfromtopo(inputCfg, iteration, datain)

% This function is used by hcp_icaclass script. It inputs ICA components
%  and computes vector angle with a template topology for Eye and Heart artifacts.
% Based on
% Yandong Li et al 2006 Physiol. Meas. 27 425 doi:10.1088/0967-3334/27/4/008
%
% INPUTS
%     inputCfg.templateICAFile: This a mat file containing topology
%                               templates for EOG and ECG comopnents. (Currently in sandbox)
%     inputCfg.ouputfile : Filename where to print the EOG and ECG
%                          recognized IC components. At the end of the name is appended ECG and
%                          VEOG respectively
%     iteration:  Structure containing multiple comp structures from
%                 multiple ICA decompositions. Computed in hcp_icaclass
%     datain: The data that was used for the above ICA decompositions.
%             Computed in hcp_icaclass
% OUTPUTS
%     pseudodataECG: data structure containing one virtual channel 'ECG'
%     pseudodataEOG: data structure containing one virtual channel 'VEOG'

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

%=====================================
Nrecurse=length(iteration);
templateICAFile=inputCfg.templateICAFile;
outputDataFile=inputCfg.outputfile;

hcp_read_matlab(templateICAFile,'compeyes','compheart');

%----------------------------------------------------
%=====================================

tmpChanInd=match_str(compeyes.topolabel,datain.label);
eyeUnmix=compeyes.unmixing(tmpChanInd);
eyeTopo=compeyes.topo(tmpChanInd);
tmpChanInd=match_str(compheart.topolabel,datain.label);
heartUnmix=compheart.unmixing(tmpChanInd);
heartTopo=compheart.topo(tmpChanInd);
Ntrials=length(datain.trial);
%=====================================

%%
totMinAngDistEye=[];
totMinAngDistHeart=[];
totMinAngIndxDistEye=[];
totMinAngIndxDistHeart=[];
for iRep=1:Nrecurse,
    
    
    comp=iteration(iRep).comp;
    
    
    NicaComputed=size(comp.topo,2);
    
    angDistEye=[];
    angDistHeart=[];
    for iComp=1:NicaComputed,
        angDistEye(iComp)=acos(abs(comp.topo(:,iComp)'*eyeTopo)./(norm(comp.topo(:,iComp))*norm(eyeTopo)));
        angDistHeart(iComp)=acos(abs(comp.topo(:,iComp)'*heartTopo)./(norm(comp.topo(:,iComp))*norm(heartTopo)));
    end
    
    [tmpMinEyeAng,indMinEye]=min(angDistEye);
    [tmpMinHeartAng,indMinHeart]=min(angDistHeart);
    
    totMinAngDistEye=tmpMinEyeAng;
    totMinAngDistHeart=tmpMinHeartAng;
    totMinAngIndxDistEye=indMinEye;
    totMinAngIndxDistHeart=indMinHeart;
    
    %tmpMinEyeAng
    %tmpMinHeartAng
    figure;
    plot(angDistEye);
    hold on;
    plot(angDistHeart,'r');
    pause(5);
    close;
    
end


[totalMinEye,indRecMinEye]=min(totMinAngDistEye);
[totalMinHeart,indRecMinHeart]=min(totMinAngDistHeart);
recIndxEOG=totMinAngIndxDistEye(indRecMinEye);
recIndxECG=totMinAngIndxDistHeart(indRecMinHeart);





%%
%{
threshAngle=(pi/4);
[indxEOG]=find(angDistEye<=threshAngle);
[indxECG]=find(angDistHeart<=threshAngle);
%}
%-----------
tmpcomp=iteration(indRecMinHeart).comp;
compECG=tmpcomp;
for iTrial=1:Ntrials,
    compECG.trial{iTrial}=tmpcomp.unmixing(recIndxECG,:)*datain.trial{iTrial};
end
compECG.time=datain.time;
compECG.label=compECG.label(recIndxECG);
compECG.topo=tmpcomp.topo(:,recIndxECG);
compECG.unmixing=tmpcomp.unmixing(recIndxECG,:);
pseudodataECG=rmfield(compECG,{'unmixing','grad','sampleinfo','topo','topolabel','order'});
pseudodataECG.sampleinfo=datain.sampleinfo;
pseudodataECG.grad=datain.grad;
pseudodataECG.label={'ECG'};
%-----------
tmpcomp=iteration(indRecMinEye).comp;
compEOG=tmpcomp;
for iTrial=1:Ntrials,
    compEOG.trial{iTrial}=tmpcomp.unmixing(recIndxEOG,:)*datain.trial{iTrial};
end
compEOG.time=datain.time;
compEOG.label=compEOG.label(recIndxEOG);
compEOG.topo=tmpcomp.topo(:,recIndxEOG);
compEOG.unmixing=tmpcomp.unmixing(recIndxEOG,:);
pseudodataEOG=rmfield(compEOG,{'unmixing','grad','sampleinfo','topo','topolabel','order'});
pseudodataEOG.sampleinfo=datain.sampleinfo;
pseudodataEOG.grad=datain.grad;
pseudodataEOG.label={'VEOG'};
%=====================================================
%=============  PLOT AND SAVE =================




%========================================================
%---- PLOT ECG

Naddtrials=Ntrials;
Nsamps=size(compECG.trial{1},2);
Ncomps=length(compECG.label);
tmpConcData=([compECG.trial{1:Naddtrials}]);

cfg = [];
cfg.method = 'mtmfft';
cfg.taper  = 'hanning';
cfg.foi = [0:0.5:40];
freq       = ft_freqanalysis(cfg, compECG);

cfg=[];
cfg.layout='4D248.mat';
cfg.comment=' ';
f1=figure('position',[40   213   870   421]);

ax2=axes('position',[0.1 0.45 0.25 0.5]);
cfg.component=1;
ft_topoplotIC(cfg,compECG);
title(['ECG comp']);
ax3=axes('position',[0.4 0.5 0.18 0.36]);
cfg.component=1;
compheart.time=compECG.time;
compheart.trial=compECG.trial;
ft_topoplotIC(cfg,compheart);
title('template');
ax4=axes('position',[0.65 0.45 0.25 0.5]);
plot(freq.freq,freq.powspctrm,'linewidth',3);
title('Spectrum');
xlabel('freq(Hz)');
tmpTime=(1./compECG.fsample)*[1:length(tmpConcData)];
ax1=axes('position',[0.1 0.1 0.8 0.3]);
plot(tmpTime,tmpConcData);
axis tight;
title('sample interval');
xlabel('time(sec)');

figName=[outputDataFile,'ECG.png'];
hcp_write_figure(figName, f1);


%=========================================================
%=========================================================

%========================================================
%---- PLOT EOG





Naddtrials=Ntrials;
Nsamps=size(compEOG.trial{1},2);
Ncomps=length(compEOG.label);
tmpConcData=([compEOG.trial{1:Naddtrials}]);

cfg = [];
cfg.method = 'mtmfft';
cfg.taper  = 'hanning';
cfg.foi = [0:0.5:40];
freq       = ft_freqanalysis(cfg, compEOG);

cfg=[];
cfg.layout='4D248.mat';
cfg.comment=' ';
f1=figure('position',[40   213   870   421]);

ax2=axes('position',[0.1 0.45 0.25 0.5]);
cfg.component=1;
ft_topoplotIC(cfg,compEOG);
title(['EOG comp.']);
ax3=axes('position',[0.4 0.5 0.18 0.36]);
cfg.component=1;
compeyes.time=compEOG.time;
compeyes.trial=compEOG.trial;
ft_topoplotIC(cfg,compeyes);
title('template');
ax4=axes('position',[0.65 0.45 0.25 0.5]);
plot(freq.freq,freq.powspctrm,'linewidth',3);
title('Spectrum');
xlabel('freq(Hz)');
tmpTime=(1./compEOG.fsample)*[1:length(tmpConcData)];
ax1=axes('position',[0.1 0.1 0.8 0.3]);
plot(tmpTime,tmpConcData);
axis tight;
title('sample interval');
xlabel('time(sec)');

figName=[outputDataFile,'VEOG.png'];
hcp_write_figure(figName, f1);

%=========================================================
%=========================================================

