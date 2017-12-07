function [comp] = hcp_ICA_freq(comp, options, datain)

% HCP_ICA_FREQ allows basic frequency analysis and plots Independent
% Components computed using the FT_COMPONENTANALYSIS function. This
% function estimate the IC Power Spectrum Density and plots the topographic
% distribution (using FT_TOPOPLOTIC), the Power Spectrum Density and the
% Time Course of the Indipendent Components, with two different layouts
% (see cfg.plottype)
%
% Use as
%   hcp_ICA_freq(comp, options, datain)
%
%     where the input comp structure should be obtained from
%     FT_COMPONENTANALYSIS and the datain structure is a fieldtrip data
%     structure contining the channel time courses used as imput of the
%     FT_COMPONENTANALYSIS .
%
% Options needs to contain the following keys:
%   doplot: 'no' or 'yes' to plot the results of the frequency analysis
%           (default='no')
%
% Options used by hcp_plotICARMEG
%   component:  field that contains the list of independent component(s) to
%               be plotted as color (default all)
%   frange:     frequency range of the PSD plot ( default is [1 150] or
%               comp.fsample/2 if comp.fsample<375 Hz)
%   plottype:   'summary' or 'components' (default = 'summary')
%   saveres:    'yes' or 'no' save the resutls in a file (default = 'no')
%   saveformat: image file format (default 'fig')
%   grad:       fieldtrip gradiometer structure
%
% See also HCP_ICA_PLOT HCP_ICA_PLOTCLASSIFICATION HCP_ICA_RMEG_CLASSIFICATION

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

switch nargin
    case 3
        % reconstruct the component time courses from the data and the unmixing
        % matrix
        cfg           = [];
        cfg.unmixing  = comp.unmixing;
        cfg.topolabel = comp.topolabel;
        cfg.order     = comp.order;
        comp          = ft_componentanalysis(cfg, datain);
        comp.order    = cfg.order;
    case 2
        % check whether IC time courses are in the first input
        if ~isfield(comp,'trial') , error('IC time course not defined, data matrix required'); end
    case 1
        error('too few input arguments');
        if ~isfield(comp,'trial') , error('IC time course not defined, data matrix required'); end
    otherwise
        error('too many input arguments');
end

cfgin.doplot   = ft_getopt(options, 'doplot', 'no');
cfgin.grad     = ft_getopt(options, 'grad');

if ~isempty(cfgin.grad),
    comp.grad = cfgin.grad;
end

%----> POWER SPECTRAL DENSITY ESTIMATION <-----
win = 2048/(comp.fsample*2^(round(log2(1024/comp.fsample))));

cfg          = [];
cfg.length   = win;
comp_seg     = ft_redefinetrial(cfg, comp);

cfg          = [];
cfg.method   = 'mtmfft';
cfg.output   = 'pow';
cfg.trials   = 'all';
cfg.taper    = 'hamming';

freq_comp    = ft_freqanalysis(cfg,comp_seg);
comp.freq_comp=freq_comp;

%----> IC POWER TIME COURSE <----
cfg               = [];
cfg.bpfilter      = 'no'; %'no' or 'yes'  bandpass filter (default = 'no')
cfg.bpfreq        = [1 150];
cfg.hilbert       = 'abs';
pow_data          = ft_preprocessing(cfg,comp);

for k = 1:numel(pow_data.trial)
    pIC{k} = pow_data.trial{k}.^2;
end

comp.pIC      = pIC;
comp.pow_data = pow_data;

if strcmp(cfgin.doplot,'yes')
    hcp_ICA_plot(comp,options); % plot results
end
