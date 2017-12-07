function [data] = hcp_ICA_preprocessing(dataraw, options)

% hcp_ICA_preprocessing reads in and preprocesses data that will be later
% subjected to ICA.
% Output data is a raw data structure where bad channels and segments have
% been removed, filtered and downsampled.
%
% Use as
% [data] = hcp_ICA_preprocessing(filename, options)
%
% where filename is a string that points to a raw data file in the database
% and options is a cell-array specifying the behaviour of the
% algorithm. Cell-arrays need to be organized as sets of key-value pairs,
% i.e. {'key1', 'value1', ...}.
%
% Options needs to contain the following key:
%   resamplefs: downsample frequency
%   or
%   decimation: downsampling factor (only if 'resamplefs' is not used)
%   channels: channel selection
%   knownbadchannels:  bad channels already known to be bad, it should be a row cell % GIORGOS
%   skipped_intervals: Nx2 matrix of skipped intervals (t11 t12 ; t21 t22; ... ; tn1 tn2)
%   bandpass: band pass intervals ([f1 f2])
%   bandstop: band stop intervals ([fs1 fs2])
%
% The following steps are performed:
% -reading in all data from disk
% -band pass and if required band stop filtering of the selected channels
% -downsampling of the data to a specified frequency
%
%
% Example use:
% options = {'resamplefs', 400 , 'channel', 'MEG', 'bandpass', [1 150]};
% data    = hcp_ICA_preprocessing(filename, options);

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

%%%%%%%%% deal with the input options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the options
selchan        = ft_getopt(options, 'channels',  'MEG');
knownbadchan   = ft_getopt(options, 'knownbadchannels',  []); % GIORGOS
resamplefs     = ft_getopt(options, 'resamplefs');
dec            = ft_getopt(options, 'decimation');
bpfreq         = ft_getopt(options, 'bandpass');
bsfreq         = ft_getopt(options, 'bandstop');
bad_segments   = ft_getopt(options, 'skipped_intervals');

%--------------------------
% GIORGOS
if ~isempty(knownbadchan)
    
    selchan=[selchan, knownbadchan];
    
end
%--------------------------

% check the options
doresample = false;
if ~isempty(resamplefs) && ~isempty(dec)
    error('ambiguous definition using both resamplefs and decimation');
elseif ~isempty(resamplefs)
    fsample = resamplefs;
    doresample = true;
elseif ~isempty(dec)
    fsample = dataraw.fsample/dec;
    doresample = true;
else
    fsample = dataraw.fsample; % don't do resampling
end

if isempty(bpfreq)
    if(fsample>375)
        bpfreq = [1 150]; % default passband in Hz
    else
        % ensure the high frequency cut of is sufficiently below the Nyquist freq
        bpfreq = [1 floor(fsample/25)*10];
    end
end

if fsample<2*bpfreq(1,2),
    error('Sampling frequency must be > 2*bpfreq(1,2)');
end

%%%%%%%%% MEG channels preprocessing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reading from file and filtering
% cfg          = [];
% cfg.channel  = selchan;
% cfg.bpfilter = 'yes';
% cfg.bpfreq   = bpfreq;
% if ~isempty(bsfreq)
%   cfg.bsfilter = 'yes';
%   cfg.bsfreq   = bsfreq;
% end
% data     = ft_preprocessing(cfg,dataraw);

cfgin          = [];
cfgin.channel  = selchan;
cfgin.hpfilter = 'yes';
cfgin.hpfreq   = bpfreq(1,1);
if(bpfreq(1,1)<2)
    cfgin.hpfiltord=4;
end
junk     = ft_preprocessing(cfgin,dataraw);

cfgin          = [];
cfgin.lpfilter = 'yes';
cfgin.lpfreq   = bpfreq(1,2);
junk2     = ft_preprocessing(cfgin,junk);

if ~isempty(bsfreq)
    cfgin          = [];
    cfgin.bsfilter = 'yes';
    cfgin.bsfreq   = bsfreq;
    data     = ft_preprocessing(cfgin,junk2);
else
    data=junk2
end

% Remove the bad segments
% if strcmp(dopreproc,'n')
ntmp=size(bad_segments,1);
datalength=size(dataraw.trial{1},2);
cutwin_filt=round(fsample*3);
bad_segments(ntmp+1,:)=[1 cutwin_filt];
bad_segments(ntmp+2,:)=[datalength-cutwin_filt datalength];

% cutwin_filt=round(fsample*50);
% cutwin_filt2=round(fsample*25);
% bad_segments(ntmp+1,:)=[1 cutwin_filt];
% bad_segments(ntmp+2,:)=[datalength-cutwin_filt2 datalength];

cfg = [];
cfg.artfctdef.badsegments.artifact = bad_segments;
cfg.artfctdef.reject               = 'partial';
data = ft_rejectartifact(cfg, data);
% end

% Resampling
if doresample
    cfg            = [];
    cfg.resamplefs = fsample;
    cfg.detrend    = 'no' ;
    cfg.demean     = 'no' ;
    cfg.feedback   = 'text';
    cfg.trials     = 'all' ;
    data           = ft_resampledata(cfg, data);
end
