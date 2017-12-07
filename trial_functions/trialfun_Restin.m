function[trl]=trialfun_Restin(inputCfg)

% This is a generic Trial Definition Function for the experiments in which
% no trigger has been sent to the parallel port.
% In these cases the user defined the duration of each trial and the entire
% scan data is split  in trials with the specified duration.
% 
% 
% Input: 
%       inputCfg: Structure with fields that contain parameters 
%                 for extracting trials 
%                 Required Fields in inputCfg:
%                 ---------------------------------------------------------
%                 .datafile: Char String with the filename of the raw data
%                            from the MEG scan
%                 .trialdef.trialDuration: Desired Trial duration - in seconds.
%                            Must be positive. 
%                 ---------------------------------------------------------
%                
% Output:  
%         trl: Matrix of size  Ntrialsx3.
%              Column 1: Start sample for trials
%              Column 2: End sample for trials
%              Column 3: Offset of beginning of each trial from Trigger
%              onset in Samples. For this specific trial function it is
%              always 0.
% Example:
% 
% cfg=[];
% cfg.dataset='here goes the name of the raw data file from the MEG scan ';
% cfg.trialfun='trialfun_Restin';
% cfg.trialdef.trialDuration=1;
% cfgDefTr = ft_definetrial(cfg);
% cfgDefTr.dataformat='4d';
% cfgDefTr.headerformat='4d';
% dataRaw = ft_preprocessing(cfgDefTr);

% Copyright (C) 2011-2013 by the Human Connectome Project, WU-Minn Consortium (1U54MH091657)
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

%===========================================================
    datafile=inputCfg.dataset;
    trialDuration=inputCfg.trialdef.trialDuration;
    hdr   = ft_read_header(datafile,'headerformat','4d');
    Nsamples=hdr.nSamples;
    Fsample=hdr.Fs;
    trialSamples=floor(trialDuration*Fsample);
    
    tmpTrialRef=1:trialSamples:Nsamples;
    startSamples=tmpTrialRef(1:end-1);
    endSamples=tmpTrialRef(2:end)-1;
    
    trl=[startSamples' endSamples' repmat(0,length(startSamples),1)];
    %===========================================================
    
