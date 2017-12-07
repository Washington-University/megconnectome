%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setup the execution environment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opengl software;

% ensure that the time and date of execution are not stored in the provenance information
global ft_default
ft_default.trackcallinfo = 'no';

% allow the user to specify the path where additional data is present, e.g. the channel layout or anatomy files
if exist('path', 'var')
    addpath(path)
end

if ~exist('filename', 'var')
    error('filename should be specified')
end

% the filename is assumed to be something like
% 'rawdatadir/Phase1MEG/Subjects/SUBJECTID/Experiments/SUBJECTID_MEG/Scans/1-Rnoise_MNN_V1/Resources/4D/c,rfDC'
tok = tokenize(filename, '/');

if ~exist('subjectid', 'var')
    subjectid = tok{end-7};
end

if ~exist('experimentid', 'var')
    experimentid = tok{end-5};
end

if ~exist('scanid', 'var')
    scanid = tok{end-3};
end

if ~exist('pipelinedatadir', 'var')
    pipelinedatadir = hcp_pathdef;
end

if ~exist('lfreq', 'var')
    lfreq = 60;
end

if ~exist('coildiameter', 'var')
    % the radius is 9 mm
    coildiameter = 0.018;
end

% print the matlab and megconnectome version to screen for provenance
ver('megconnectome')

% remove temporary local variables
clear tok

% print the value of all local variables to screen for provenance
w = whos;
w = {w.name};
w = setdiff(w, {'w', 'ans'});
for i=1:length(w)
    fprintf(hcp_printstruct(w{i}, eval(w{i})));
end

% change to the location of the processed data (input and output)
cd(pipelinedatadir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% execute the pipeline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

resultprefix = sprintf('%s_%s', experimentid, scanid);
dirname = fileparts(filename);

% Perform a number of sanitychecks on a MEG dataset.
%
% The header and events are read and a textfile is generated providing some basic information about the dataset:
%   - sampling frequency
%   - length in seconds
%   - number of recorded MEG channels
%   - number of recorded REF channels
%   - number of recorded EEG channels
%   - number of Trigger channels
%   - number of Response channels
%   - number of events
%   - additional event information
%   - a quantification of head movement (comparing COH with COH2)
%
% Subsequently the data is read in and some basic quantification is done to assess
% the quality of the data. The following figures will be generated:
%   - coregistration image, showing the hs_file (head surface) in
%     combination with the sensor-array
%   - power spectra of MEG channels (computed after chopping up the data in
%     1024 sample snippets
%   - power spectra of MEGREF channels (computed after chopping up the data
%     in 1024 sample snippets
%   - power spectra of EEG channels (computed after chopping up the data in
%     1024 sample snippets
%   - spatial topography of the power line noise

% add the general pipeline name "datacheck" to the files that are created
output = fullfile(pipelinedatadir, [resultprefix '_datacheck']);

filename = [dirname, '/c,rfDC'];
hsname   = [dirname, '/hs_file'];

hdr   = ft_read_header(filename);
event = ft_read_event(filename);

% provide some textual feedback: this should go into a text file
textfile = [output, '_info.txt'];
fid = hcp_fopen(textfile, 'w');
if fid<0
    error('cannot open output file');
end

hcp_fprintf(fid, '\n%s : %s\n', 'datafile', filename);
hcp_fprintf(fid, '%s : %5.2f %s\n', 'sampling frequency ', hdr.Fs, 'Hz');
hcp_fprintf(fid, '%s : %5.2f %s\n', 'length             ', (hdr.nSamples*hdr.nTrials)/hdr.Fs, 'seconds');
hcp_fprintf(fid, '%s : %d\n', 'number of MEG channels     ', sum(strcmp(hdr.chantype, 'megmag')));
hcp_fprintf(fid, '%s : %d\n', 'number of REF channels     ', sum(strcmp(hdr.chantype, 'megref')));
hcp_fprintf(fid, '%s : %d\n', 'number of EEG channels     ', sum(strcmp(hdr.chantype, 'eeg')));
hcp_fprintf(fid, '%s : %d\n', 'number of ECG channels     ', sum(strcmp(hdr.chantype, 'ecg')));
hcp_fprintf(fid, '%s : %d\n', 'number of EMG channels     ', sum(strcmp(hdr.chantype, 'emg')));
hcp_fprintf(fid, '%s : %d\n', 'number of TRIGGER channels ', sum(strcmp(hdr.chantype, 'trigger')));
hcp_fprintf(fid, '%s : %d\n\n', 'number of events           ', numel(event));

% show event table
if isempty(event)
    hcp_fprintf(fid, 'no events were found in the datafile\n');
else
    % show some details
    eventtype = unique({event.type});
    Neventtype = length(eventtype);
    if Neventtype==0
        hcp_fprintf(fid, 'no events were found in the datafile\n');
    else
        hcp_fprintf(fid, 'the following events were found in the datafile\n');
        for i=1:Neventtype
            sel = find(strcmp(eventtype{i}, {event.type}));
            try
                eventvalue = unique({event(sel).value});            % cell-array with string value
                eventvalue = sprintf('''%s'' ', eventvalue{:});     % translate into a single string
            catch
                eventvalue = unique(cell2mat({event(sel).value}));  % array with numeric values or empty
                eventvalue = num2str(eventvalue);                   % translate into a single string
            end
            hcp_fprintf(fid, 'event type: ''%s'' ', eventtype{i});
            hcp_fprintf(fid, 'with event values: %s', eventvalue);
            hcp_fprintf(fid, '\n');
        end
    end
end % isempty event

% show head localization info
for k = 1:numel(hdr.orig.user_block_data)
    type{k} = hdr.orig.user_block_data{k}.hdr.type;
end
sel  = find(strcmp('B_COH_Points', type));
if numel(sel)==2
    pnt1 = hdr.orig.user_block_data{sel(1)}.pnt; % Coils in dewar space
    pnt2 = hdr.orig.user_block_data{sel(2)}.pnt;
    
    dpnt = sqrt(sum((pnt1-pnt2).^2,2))*1000; % amount of movement in mm
    dpnt = mean(dpnt);
    hcp_fprintf(fid, '\n');
    hcp_fprintf(fid, 'average coil movement: %s mm', num2str(dpnt, '%5.3f'));
    hcp_fprintf(fid, '\n');
else
    hcp_fprintf(fid, '\n');
    hcp_fprintf(fid, 'no head coil information found');
    hcp_fprintf(fid, '\n');
end

hcp_fprintf(fid, '\n');
hcp_fclose(fid);

seleeg = ft_channelselection('E*', hdr.label);     haseeg = ~isempty(seleeg); % FIXME: EEG as heuristic does not work for this set of labels, returning empty
selmeg = ft_channelselection('MEG', hdr.label);    hasmeg = ~isempty(selmeg);
selref = ft_channelselection('MEGREF', hdr.label); hasref = ~isempty(selref);

% get the layout for the MEG topography
cfg = [];
cfg.layout = '4D248.mat';
lay = ft_prepare_layout(cfg);

% plot sensors and headshape
if exist(hsname, 'file')
    shape = ft_read_headshape(hsname);
    f1 = figure;
         a = [-0.15 0.15 -0.10 0.10 -0.10 0.15];
    subplot(2,3,1); ft_plot_headshape(shape); view([-90 90]); axis(a); % top
    subplot(2,3,4); hold on; ft_plot_headshape(shape, 'fidlabel', false, 'fidcolor', 'none'); view([-90 90]); ft_plot_sens(hdr.grad, 'chantype', 'megmag', 'coildiameter', coildiameter);
    subplot(2,3,2); ft_plot_headshape(shape); view([ 00 00]); axis(a); % right
    subplot(2,3,5); hold on; ft_plot_headshape(shape, 'fidlabel', false, 'fidcolor', 'none'); view([0 0]); ft_plot_sens(hdr.grad, 'chantype', 'megmag', 'coildiameter', coildiameter);
    subplot(2,3,3); ft_plot_headshape(shape); view([-90 00]); axis(a); % back
        subplot(2,3,6); hold on; ft_plot_headshape(shape, 'fidlabel', false, 'fidcolor', 'none'); view([90 0]); ft_plot_sens(hdr.grad, 'chantype', 'megmag', 'coildiameter', coildiameter);
    %     subplot(2,2,1); ft_plot_headshape(shape); view([-90 90]); axis(a); % top
    %     subplot(2,2,2); ft_plot_headshape(shape); view([ 00 00]); axis(a); % right
    %     subplot(2,2,3); ft_plot_headshape(shape); view([-90 00]); axis(a); % back
    %     subplot(2,2,4); hold on; ft_plot_headshape(shape, 'fidlabel', false, 'fidcolor', 'none'); view([0 0]); ft_plot_sens(hdr.grad, 'chantype', 'megmag', 'coildiameter', coildiameter); axis(a);
    figurefile = [output, '_headshape'];
    hcp_write_figure(figurefile, f1);
    close(f1);
end

% create a trl matrix
trl      = (1:1024:(hdr.nSamples*hdr.nTrials))';
trl(:,2) = min(trl(:,1)+1023, hdr.nSamples*hdr.nTrials);
trl(:,3) = 0;
% only keep the trials that are precisely 1024 samples, exclude left-overs at the end
trllen = trl(:,2)-trl(:,1)+1;
trl = trl(trllen==1024,:);

% set up the cfgs for the preprocessing and spectral analysis
cfg          = [];
cfg.datafile = filename;
cfg.feedback = 'none';

% some general preprocessing, no demeaning first to get the low freq stuff
cfg1         = [];
cfg1.demean  = 'yes';
cfg1.feedback = 'none';

% for the time course of the line noise
cfg2           = [];
cfg2.dftfilter = 'yes';
cfg2.dftfreq   = lfreq + [-1 0 1].*(hdr.Fs./1024);
cfg2.feedback  = 'none';
cfg2.dftinvert = 'yes';
cfg2.rectify   = 'yes';
cfg2.boxcar    = 0.2;

% for the detection of jumps
cfg3               = [];
cfg3.absdiff       = 'yes';
cfg3.medianfilter  = 'yes';
cfg3.medianfiltord = 9;
cfg3.feedback      = 'none';

% for the analysis of low frequency power
cfg4          = [];
cfg4.lpfilter = 'yes';
cfg4.lpfreq   = 2;
cfg4.feedback = 'none';

% for the spectral analysis
cfgf        = [];
cfgf.method = 'mtmfft';
cfgf.output = 'pow';
cfgf.taper  = 'hanning';
cfgf.feedback = 'none';

% allocate some memory for later
minval = zeros(248,1)+inf;
maxval = zeros(248,1)-inf;
sumval = zeros(248,1);
ssqval = zeros(248,1);
covmat = zeros(248);

if hasmeg,
    for k = 1:200:size(trl,1)
        % chop up the data in 1024 sample snippets; the correlation heuristic
        % is likely to fail when the individual snippets are too long
        % (correlation then mainly reflects the very low frequency fluctuations
        % due to the environment
        cfg.trl     = trl(k:min(size(trl,1),k+199),:);
        cfg.channel = selmeg;
        tmpdata     = ft_preprocessing(cfg);
        tmpfreq     = ft_freqanalysis(cfgf, tmpdata);
        
        % post-process the data for low-freq power <2Hz: before demeaning
        tmp = ft_preprocessing(cfg4, tmpdata);
        if k==1,
            tmptrace0 = zeros(numel(tmp.label),0);
        end
        for m = 1:numel(tmp.trial)
            tmptrace0 = cat(2, tmptrace0, var(tmp.trial{m},[],2));
        end
        
        % demean the data and compute some quantities
        tmpdata = ft_preprocessing(cfg1, tmpdata);
        for m = 1:numel(tmpdata.trial)
            minval  = min(minval, min(tmpdata.trial{m},[],2));
            maxval  = max(maxval, max(tmpdata.trial{m},[],2));
            sumval  = sumval + sum(tmpdata.trial{m},2);
            covmat  = covmat + tmpdata.trial{m}*tmpdata.trial{m}';
            ssqval  = ssqval + sum(tmpdata.trial{m}.^2,2);
        end
        
        % post-process the data line noise
        tmp = ft_preprocessing(cfg2, tmpdata);
        tmptmp = zeros(numel(tmp.label),0);
        for m = 1:numel(tmp.trial)
            if size(tmp.trial{m},2)<512
                continue;
            elseif size(tmp.trial{m},2)<1024
                ix = 256;
            else
                ix = [256 768];
            end
            tmptmp = cat(2, tmptmp, tmp.trial{m}(:,ix)); % sample two points of the smoothed linenoise time courses
        end
        if k==1,
            powline = tmptmp;
        else
            powline = cat(2, powline, tmptmp);
        end
        % FIXME what to do with this?
        
        % post-process the data for squid jumps
        tmp  = ft_preprocessing(cfg3, tmpdata);
        if k==1,
            tmptrace1 = zeros(numel(tmp.label),0);
        end
        for m = 1:numel(tmp.trial)
            tmptrace1 = cat(2, tmptrace1, max(tmp.trial{m},[],2));
        end
        
        % accumulate the power spectra
        if k==1,
            pow  = tmpfreq.powspctrm.*numel(tmpdata.trial);
        else
            pow  = pow  + tmpfreq.powspctrm.*numel(tmpdata.trial);
        end
    end
    pow  = pow./size(trl,1);
    
    % plot for jumps
    f1=figure;
    plot(tmptrace1');
    xlabel('epoch (#)');
    ylabel('amplitude (T)');
    title('MEG squid jumps');
    figurefile = [output, '_jumps'];
    hcp_write_figure(figurefile, f1);
    close(f1);
    
    % plot for low frequency power
    f1=figure;
    plot(tmptrace0');
    xlabel('epoch (#)');
    ylabel('power (T^2/Hz)');
    title('MEG low frequency power');
    figurefile = [output, '_MEG_lowfreq_power'];
    hcp_write_figure(figurefile, f1);
    close(f1);
    
    % make a plot of the powerspectra
    f1=figure;
    loglog(tmpfreq.freq,pow);
    abc = axis;
    axis([1 tmpfreq.freq(end) abc(3:4)]);
    xlabel('log10 frequency (Hz)');
    ylabel('log10 power (T^2/Hz)');
    title('powerspectra MEG');
    figurefile = [output, '_MEG_powspctrm'];
    hcp_write_figure(figurefile, f1);
    close(f1);
    
    % make a topographical plot of the power line noise
    [a,b] = match_str(lay.label, tmpfreq.label);
    c     = nearest(tmpfreq.freq, lfreq);
    
    f1=figure; hold on;
    lay.label{end-1} = 'power (T^2/Hz)';
    ft_plot_lay(lay, 'box', 'off');
    ft_plot_topo(lay.pos(a,1),lay.pos(a,2),pow(b,c),'gridscale',150,'outline',lay.outline,'mask',lay.mask,'interpmethod','nearest');
    axis([-0.6 0.6 -0.6 0.6]);
    axis off;
    abc = caxis;
    caxis([-1 1]*abc(2));
    colorbar
    title('power line noise on MEG sensors');
    figurefile = [output, '_MEG_powerline_noise'];
    hcp_write_figure(figurefile, f1);
    close(f1);
end % if hasmeg

% if haseeg,
%   for k = 1:200:size(trl,1)
%     cfg.trl     = trl(k:min(size(trl,1),k+199),:);
%     cfg.channel = seleeg;
%     tmpeeg      = ft_preprocessing(cfg);
%     tmpfreqe    = ft_freqanalysis(cfgf, tmpeeg);
%
%     if k==1,
%       powe = tmpfreqe.powspctrm.*numel(tmpfreqe.cumtapcnt);
%     else
%       powe = powe + tmpfreqe.powspctrm.*numel(tmpfreqe.cumtapcnt);
%     end
%   end
%   powe = powe./size(trl,1);
%
%   f2=figure;
%   loglog(tmpfreq.freq,powe);
%   abc = axis;
%   axis([1 tmpfreq.freq(end) abc(3:4)]);
%   title('powerspectra EEG');
%   figurefile = [output, '_EEG_powspctrm'];
%   hcp_write_figure(figurefile, f1);
%   close(f1);
% end % if haseeg

if hasref,
    for k = 1:200:size(trl,1)
        cfg.trl     = trl(k:min(size(trl,1),k+199),:);
        cfg.channel = selref;
        tmpref      = ft_preprocessing(cfg);
        tmpfreqr    = ft_freqanalysis(cfgf, tmpref);
        
        if k==1,
            powr = tmpfreqr.powspctrm.*numel(tmpref.trial);
        else
            powr = powr + tmpfreqr.powspctrm.*numel(tmpref.trial);
        end
    end
    powr = powr./size(trl,1);
    
    f1=figure;
    loglog(tmpfreq.freq,powr);
    abc = axis;
    axis([1 tmpfreq.freq(end) abc(3:4)]);
    xlabel('log10 frequency (Hz)');
    ylabel('log10 power (T^2/Hz)');
    title('powerspectra MEGREF');
    figurefile = [output, '_MEGREF_powspctrm'];
    hcp_write_figure(figurefile, f1);
    close(f1);
end % if hasref

% compute mean value (meaningless)
nsmp    = hdr.nSamples*hdr.nTrials;
meanval = sumval./nsmp;

% compute variance and std
varval  = (ssqval - (sumval.^2./nsmp))./(nsmp-1);
stdval  = sqrt(varval);

% compute channel correlation matrix
covmat  = covmat - (sumval*sumval')./nsmp;
corrmat = covmat./sqrt(diag(covmat)*diag(covmat)');

% Set correlation to zeror for zero covariance
corrmat(~isfinite(corrmat)) = 0;

% use Giorgos' strategy to assess noisy channels
cfg = [];
cfg.method     = 'distance';
cfg.neighbourdist = 0.035;
cfg.grad       = hdr.grad;
neighb         = ft_prepare_neighbours(cfg);
Nmat           = zeros(248);
for k = 1:248
    [a,b] = match_str(hdr.grad.label(1:248), neighb(k).neighblabel);
    Nmat(a,k) = 1;
end
corrmat(~Nmat) = nan;
C = nanmean(corrmat);

% Considering different arrangment of labels in "lay.label" and "neighb.label"
% it's probably safe not to assume anything
[a,b] = match_str(lay.label, {neighb.label});

% plot the average neighbour correlation
f1=figure; hold on;
ft_plot_lay(lay, 'box', 'off');
ft_plot_topo(lay.pos(a,1),lay.pos(a,2),C(b),'gridscale',150,'outline',lay.outline,'mask',lay.mask,'interpmethod','nearest');
axis([-0.6 0.6 -0.6 0.6]);
axis off;
caxis([0.4 0.65]);
colorbar
title('average correlation with neighbours');
figurefile = [output, '_neighb_correlation'];
hcp_write_figure(figurefile, f1);
close(f1);

indTrigCh=find(strcmp('TRIGGER', hdr.label));
indRespCh=find(strcmp('RESPONSE', hdr.label));

dat = ft_read_data(filename,'chanindx',indTrigCh);
trig = dat;
dat = ft_read_data(filename,'chanindx',indRespCh);
resp = dat;

f1 = figure;
% Raw trigger, response channels
subplot(2,1,1), plot(trig); ylabel('TRIGGER'); set(gca,'xticklabel',[]);set(gca,'Position',[0.15 0.55 0.8 0.43]); axis tight;
subplot(2,1,2), plot(resp); ylabel('RESPONSE'); xlabel('sample');set(gca,'Position',[0.15 0.1 0.8 0.43]); axis tight;
figurefile = [output, '_triggers'];
hcp_write_figure(figurefile, f1);
close(f1);

%% Process EOG/ECG/EMG and plot

cfg = [];
cfg.dataset = filename;
dataraw = ft_preprocessing(cfg);

montage = hcp_exgmontage(subjectid, experimentid, scanid);

if isfield(montage, 'labelorg') && ~isempty(montage.labelorg) % This check is performed for cases where no Electrodes have been recorded (i.e Glasgow data)
    hasELEC = 1;
    dataELEC = ft_selectdata(dataraw,'channel',montage.labelorg);
    dataELECnew = ft_apply_montage(dataELEC,montage);
    
    % Remove line noise ?
    
    cfg = [];
    cfg.detrend = 'no';
    cfg.bsfilter = 'yes';
    cfg.bsfreq = lfreq+[-1 1];
    dataELECnew2 = ft_preprocessing(cfg, dataELECnew);
    clear dataELECnew
    
    cfg = [];
    cfg.detrend = 'no';
    cfg.hpfilter = 'yes';
    cfg.hpfreq = 0.25;
    cfg.hpfiltord=4;
    dataELECnew2 = ft_preprocessing(cfg, dataELECnew2);
    
    
    cfg = [];
    cfg.detrend = 'no';
    cfg.lpfilter = 'yes';
    cfg.lpfreq = 150;
    dataELECnew2 = ft_preprocessing(cfg, dataELECnew2);
    
    %-freq analysis of elec channels
    
    
    
    %=====================================================
    % The following  code has been modified by francesco
    
    %{ 
      OLD CODE by RObert
    cfgf           = [];
    cfgf.method    = 'mtmfft';
    cfgf.output    = 'pow';
    cfgf.taper     = 'hanning';
    cfgf.feedback  = 'text';
    cfgf.foi       = 1:0.1:100;
    cfgf.tapsmofrq = 0.1;
    freqELEC = ft_freqanalysis(cfgf,dataELECnew2);
    %}
    
    % NEW CODE by Francesco
    win=2048/(dataELECnew2.fsample*2^(round(log2(1024/dataELECnew2.fsample))));
    
    cfg = [];
    cfg.length               = win;
    ref_data_seg             = ft_redefinetrial(cfg, dataELECnew2);
    
    cfg=[];
    cfg.method = 'mtmfft';
    cfg.output = 'pow';
    cfg.trials = 'all';
    cfg.taper = 'hamming';
    
    [freq_ref_data] = ft_freqanalysis(cfg,ref_data_seg);
    
    selec = sqrt(freq_ref_data.powspctrm);
    
    freqELEC =freq_ref_data;
    freqELEC.powspctrm = selec ;
   %=====================================================
    
    %---------------------------------------
    Nelecs=length(dataELECnew2.label);
    for iElec=1:Nelecs,
        hp1=figure;
        set(hp1,'position',[10         100        6*200         200]);
        hax1=subplot(1,2,1);
        hax2=subplot(1,2,2);
        set(hax1,'position',[0.05 0.175 0.7 0.7],'fontsize',8);
        set(hax2,'position',[0.79 0.175 0.2 0.7],'fontsize',8);
        axes(hax1);
        plot(dataELECnew2.time{1}, dataELECnew2.trial{1}(iElec,:));axis tight;
        xlabel('time(sec)');ylabel('signal');title([regexprep(freqELEC.label{iElec},'_','\\_'),'  - Time series']);
        axes(hax2);
%         semilogx(freqELEC.freq, freqELEC.powspctrm(iElec,:));axis tight;
        plot(freqELEC.freq, freqELEC.powspctrm(iElec,:));axis tight;
%         hlabx=xlabel('frequency (Hz)');title(['Spectrum: Power vs log(freq.)']);
         hlabx=xlabel('frequency (Hz)');title(['Spectrum: Power vs freq.']);

        axis([0 70 0 max(freqELEC.powspctrm(iElec,:))])
%         set(hax2,'XTick',[1 100])
        %set(hlabx,'Units','normalized')
        %set(hlabx,'position',[0.15 -0.035 1])
        
        %figpos=get(hp1,'position');
        %pappos=get(hp1,'paperpos');
        %pappos(4)=(figpos(4)./figpos(3))*pappos(3);
        pappos=[0.25 1 6*2 2];
        set(hp1,'paperpos',pappos);
        
        figurefile = [output, '_elecchan_',freqELEC.label{iElec}];
        hcp_write_figure(figurefile, hp1);
        close(hp1);
    end
else
    hasELEC=0;
    display(['NO ELECTRODES WERE PROVIDED BY EXGMONTAGE for ', subjectid, '_', experimentid, '_', scanid]);
end


%% ensure that the expected pipeline output is present
hcp_check_pipelineoutput('datacheck', 'subject', subjectid, 'experiment', experimentid, 'scan', scanid);
