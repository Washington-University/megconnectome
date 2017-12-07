function hcp_ICA_plot(comp_freq, options)

% hcp_ICA_plot allows summury plots of Independent Components computed using
% the hcp_ICA_freq function.
% This function plots the topographic distribution (using FT_TOPOPLOTIC),
% the Power Spectrum Density and the Time Course of the Indipendent Components, with
% two different layouts (see cfg.plottype)
%
% Use as
%   hcp_ICA_plot(comp_freq, options, datain)
% where the input comp_freq structure should be obtained from FT_COMPONENTANALYSIS and
% the datain structure is a fieldtrip data structure contining the channel time courses
% used as imput of the FT_COMPONENTANALYSIS .
%
% Options needs to contain the following key:
%   component          : field that contains the independent
%                        component(s) to be plotted as color (default all)
%   plottype           : 'summary' or 'components' (default = 'summary')
%   frange             : frequency range of the PSD plot (default is [1 150] or comp_freq.fsample/2 if comp_freq.fsample<375 Hz);
%   saveres            : 'yes' or 'no' save the resutls in a file (default = 'no')
%   saveformat         : image file format (default 'fig')
%   grad               : fieldtrip gradiometer structure
%
% See also HCP_ICA_FREQ HCP_ICA_PLOTCLASSIFICATION HCP_ICA_RMEG_CLASSIFICATION

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

Nc = size(comp_freq.topo,2);

cfgin.plottype   = ft_getopt(options, 'plottype',   'components'); %summary or components
cfgin.saveres    = ft_getopt(options, 'saveres'); 
cfgin.saveformat = ft_getopt(options, 'saveformat', 'fig');
cfgin.fileout    = ft_getopt(options, 'fileout'); 
cfgin.parameter  = ft_getopt(options, 'parameter',  'topo');
cfgin.component  = ft_getopt(options, 'component',  1:Nc);
cfgin.comment    = ft_getopt(options, 'comment',    'no');
cfgin.frange     = ft_getopt(options, 'frange');
cfgin.grad       = ft_getopt(options, 'grad');
cfgin.dataset    = ft_getopt(options, 'dataset');
cfgin.zlim       = ft_getopt(options, 'zlim',       'maxabs');
cfgin.colorbar   = ft_getopt(options, 'colorbar',   'yes');
cfgin.posi       = ft_getopt(options, 'position');
fvis             = ft_getopt(options, 'visible', 'on');
cfgin.modality   = ft_getopt(options, 'modality','MEG');

if strcmp(cfgin.modality,'MEG')
    cfgin.layout='4D248.mat';
elseif strcmp(cfgin.modality,'EEG')
    cfgin.layout='EEG1010.lay';
end

if isempty(cfgin.grad)
    if ~isempty(cfgin.dataset)
        sens           = ft_read_sens(cfgin.dataset) ;
        cfgin.grad     = sens;
        comp_freq.grad = cfgin.grad;
    else
        cfgin.grad = comp_freq.grad;
        warning('myControl:FieldTrip','the gradiometer structure may be undefined (correctly)')
    end
end

% Set the High Frequency limit in the PSD plot
if isempty(cfgin.frange)
    if(comp_freq.fsample>375)
        frange = [1 150]; % default passband if comp_freq.fsample>375 Hz
    else
        frange = [1 (floor(comp_freq.fsample/20)*10)]; % set highest frequency to the Nyquist freq if comp_freq.fsample<375 Hz
    end
else frange = cfgin.frange;
end

% create temporary variable
selcomp = cfgin.component;

freq_comp = comp_freq.freq_comp;
mspec     = sqrt(freq_comp.powspctrm(selcomp,:));
F         = freq_comp.freq;

if strcmp(cfgin.plottype,'summary')
    
    %----> Topological maps of mixing matrix columns <----
    f1 = figure;
    set(f1, 'Position', [50 50 1700 900]);
    ft_topoplotIC(cfgin, comp_freq)
    
    %----> POWER SPECTRAL DENSITY PLOTTING <-----
    
    % allow multiplotting
    nplots = numel(selcomp);
    nyplot = ceil(sqrt(nplots));
    nxplot = ceil(nplots./nyplot);
    
    f2 = figure;
    set(f2, 'Position', [50 50 1700 900]);
    for i = 1:length(selcomp)
        subplot(nxplot, nyplot, i);
        plot(F,mspec(i,:));
        title(['PSD ic ' num2str(selcomp(i))]);
        Prms = mspec(i,:);
        lags = find(F>1 & F<frange(1,2));
        max_val = max(Prms(lags))+0.1*(max(Prms(lags))-min(Prms(lags)));
        axis([frange(1,1) frange(1,2) 0 max_val]);
    end
    
    %----> IC TIME COURSE <----
    
    f3 = figure;
    set(f3, 'Position', [50 50 1700 900]);
    for i = 1:length(selcomp)
        subplot(nxplot, nyplot, i); hold on
        for m =  1:numel(comp_freq.trial)
            plot(comp_freq.time{m}, comp_freq.trial{m}(selcomp(i),:));
        end
        title(['Timecourse ic ' num2str(selcomp(i))]);
    end
    
    if (strcmp(cfgin.saveres,'yes'))
        hcp_write_figure('topo_IC_summary', f1, 'format', cfgin.saveformat);
        hcp_write_figure('PSD_IC_summary',  f2, 'format', cfgin.saveformat);
        hcp_write_figure('TC_IC_summary',   f3, 'format', cfgin.saveformat);
    end
    
elseif strcmp(cfgin.plottype,'components')
    
    leg = char('OneOverF (>0.5)  ' , 'Specflat (<2.0)  ' , 'Kurtosis (>15)   ' , 'elecSig (>0.1)  ', 'elecPow (>0.25) ', 'elecSpe (>0.95) ', 'Brain           ');
    
    for ix = 1:length(selcomp)
        if(~isempty(cfgin.posi))
            f = figure;
            set(f, 'visible', fvis,'on','Position', cfgin.posi);
        else
            f = figure;
            set(f, 'visible', fvis);
        end
        
        subplot(2,2,1)
        cfgin.component = selcomp(ix);
        ft_topoplotIC(cfgin, comp_freq);
        xlabel(['IC ' num2str(selcomp(ix))],'color','k');
        title('IC sensor map','FontSize',12);
        
        subplot(2,2,[3 4])
        plot(F,mspec(ix,:));
        title(['IC ' num2str(selcomp(ix))]);
        Prms = mspec(ix,:);
        lags = find(F>1 & F<frange(1,2));
        max_val = max(Prms(lags))+0.1*(max(Prms(lags))-min(Prms(lags)));
        axis([frange(1,1) frange(1,2) 0 max_val]);
        ylabel(['PSD']);
        xlabel(['Frequency (Hz)'],'color','k');
        if isfield(comp_freq,'class')
            str(1,:)=num2str(abs(comp_freq.class.spectrum_one_over_f(1,selcomp(ix))),'%05.3f');
            str(2,:)=num2str(abs(comp_freq.class.spectrum_flat(1,selcomp(ix))),'%05.1f');
            str(3,:)=num2str(abs(comp_freq.class.time_kurtosis(1,selcomp(ix))),'%05.1f');
            str(4,:)=num2str(abs(comp_freq.class.elc_signal_correlation(1,selcomp(ix))),'%05.3f');
            str(5,:)=num2str(abs(comp_freq.class.elc_power_correlation(1,selcomp(ix))),'%05.3f');
            str(6,:)=num2str(abs(comp_freq.class.elc_spectrum_correlation(1,selcomp(ix))),'%05.3f');
            if(find(comp_freq.class.brain_ic==selcomp(ix)))
                str(7,:)='yes  ';
            else
                str(7,:)='no   ';
            end
            aio=[leg str];
            legend(aio);
            set(legend,'FontSize',8);
        end
        
        subplot(2,2,2);
        timeic=cell2mat(comp_freq.time);
        trialic=cell2mat(comp_freq.trial);
        plot(timeic(1,:),trialic(ix,:));
        sig = trialic(ix,:);
        base = timeic;
        axis([base(1) base(end) min(sig)-0.1*(max(sig)-min(sig)) max(sig)+0.1*(max(sig)-min(sig))]);
        ylabel(['amplitude']);
        xlabel(['time (Sec)']);
        title('IC time course','FontSize',12);
        
        if (strcmp(cfgin.saveres,'yes'))
            if(isempty(cfgin.fileout))
                hcp_write_figure(['ICA_components' num2str(selcomp(ix))], f, 'format', cfgin.saveformat)
            elseif(size(selcomp,2)>1)
                hcp_write_figure([cfgin.fileout '_' num2str(selcomp(ix))], f, 'format', cfgin.saveformat)
            else
                hcp_write_figure(cfgin.fileout, f, 'format', cfgin.saveformat)
            end
            close
            clear f;
        end
    end
    
else error('plot settings must be summary or components');
end
