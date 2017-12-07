function[badTrials,trialsStdMat]=hcp_arti_usestd(inData,stdRatio)

% This function idenitifies but trials by using the ratio of the std of
% single trials over the std of all trials for each channel. The user
% defines a threshold. For each trials all channels are checked to see if
% the std artio is above the threshold. If at least ONE(this is hard coded. should be dynamic?)
% channel is found aboce the threshold then this trial is assumed bad.
% Please notice that this methid tends to remove trials that have an
% unusually high variance. Channels with persistently high variance might
% escape.
%
%
% INPUTS
%  stdRatio: This is the rate of the std of individual trial ./ std of all trials
%  per channel used to find bad trials. Trials in this category are checked and
%  according to the number of channels for which this rate is higher
% (set below Nchannels4BadTrial), they are marked as bad or not
%
% OUPUTS
%  badTrials: The indices of the bad Trials
%  stdMat: Matrix Nchannels x Ntrials with the std for each trial and channel

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


threshStdRate=stdRatio;

Nchannels4BadTrial=1; % This is the least number of channels for which a bad trial must have the std ratio higher than the threshold. TODO: make it an input argument ?
%----
%===========================

trialsStd = cellfun(@(x) std(x,0,2), inData.trial, 'UniformOutput', false);
trialsStdMat=[trialsStd{:}];
% [COEFF, SCORE, LATENT, TSQUARED] = princomp(trialsStdMat);
Nchannels=size(inData.trial{1},1);
NsampPerTrial=size(inData.trial{1},2);
Ntrials=length(inData.trial);
%=========================================================
totalData=[inData.trial{:}];
stdTotalPerChan=std(totalData,0,2);
%=========================================================
stdMatNorm=trialsStdMat./repmat(stdTotalPerChan,1,Ntrials);
[chanIndx,trlIndx]=find(stdMatNorm>threshStdRate); % TODO: This rate should be dynamically set by the user
[uniqueChan,dummyIndx,uniqueChanMap]=unique(chanIndx);
[uniqueTrial,dummyIndx,uniqueTrialMap]=unique(trlIndx);
NchanUniq=length(uniqueChan);
NtrialUniq=length(uniqueTrial);
perChanCounts=zeros(1,NchanUniq);
perTrialCounts=zeros(1,NtrialUniq);
for iCh=1:NchanUniq,
    perChanCounts(iCh)=sum(uniqueChanMap==iCh);
end
for iTrl=1:NtrialUniq,
    perTrialCounts(iTrl)=sum(uniqueTrialMap==iTrl);
end
%=========================================================
indx1=find(perTrialCounts>Nchannels4BadTrial); % TODO: This number of Channels in which a trial appears bad  should be dynamically set by the user
badTrials=uniqueTrial(indx1);

