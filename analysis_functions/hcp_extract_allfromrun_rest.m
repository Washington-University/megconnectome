function [data] = hcp_extract_allfromrun_rest(inCfg)
%% This function performs the core processing of the tmegpreproc and rmegpreproc pipelines.
% It extracts  trial data for a give data group and fuses them with the
% results from hcp_baddata.m pipeline so that bad channels and trials that
% coincide with noisy periods are removed. In the case that a trial spans
% a long block of data, in order to avoid removing it entirely, the bad
% segments in the trials are replaced with nan. (This is the case for the Story/Math data groups BSENT and BUN).
% Then the IC components identified as related to heart or eye activity by hcp_icaclass.m pipeline are removed.
% Then the data is resampled to one fourth of the original sampling frequency in order to reduce the size of the
% dataset.
%
%
% INPUT:
%-------------------------------------------------------------------
% inCfg : This is a structure containing required parameters for the
%          analysis
%          Fields:
%                .datafile:  This is the raw data filename for a given scan.
%                .trl:       This is the trials definition matrix. It has 3 columns and number of rows equal to number of trials. The
%                              first column is the start sample, the second is the end sample and time offset of the start of the trial
%                              relative to the 0 reference point. This information is created by the trial definition functions.
%                .badchanfile: This is the file containing information about bad channels. This information is created by
%                                 hcp_baddata.m pipeline.
%                .badsegmfile: This is the file containing information about bad segments. This information is created by
%                                 hcp_baddata.m pipeline.
%                .icainfofile: This is the file containing information about artifactual Independent Components in the data. This
%                                 information is created hcp_icaclass.m pipeline
%                .badsegmode:  This variable defines if trials containing bad segments will be removed in full or the bad segments will be replaced by NANs.
%                              'remfull' for remove full trial or 'repnan'
%                              for replacing with NANs.  This field is set in the alltrialdefparams_*.m scripts.
%
%                .montage:     This is the montage containing the emg channels. This
%                                variable is constructed by the hcp_exgmontage.m function.
%                                In .labelnew subfield , the expected emg channels names are  expected to be
%                                {'EMG_LH','EMG_LF','EMG_RH','EMG_RF'}
%                .lineFreq:    Numerical array that contain the frequencies
%                                of line current to be filtered out i.e. [60 120].
%                .outputfile:   The filename of the data file where the cleaned data will be saved.
%
%
% OUTPUT
%-----------------------------------------------------------
% cleanTrl : %   This is a numerical matrix similar to the input inCfg.trl trials definition matrix described above.
%                It has 3 columns and number of rows equal to number of CLEAN trials that remained after the cleaning performed.
%                The first column is the start sample, the second is the end sample and time offset of the start of the trial
%                relative to the 0 reference point.
%-------------------------------------------------------------

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

datafile = inCfg.datafile;
icainfofile = inCfg.icainfofile;
trl = inCfg.trl;

outputfile = inCfg.outputfile;
badsegmode = inCfg.badsegmode; % tmpCfg.badsegmode = 'remfull'; % This is used by tmeg_preproc to know if trials should be fully removed, or bad segments should be replaced by nans(i.e. when each trial contins an entire block)
montage = inCfg.montage;       % EEG and or EMG (for both one can use E*)
lineFreq = inCfg.lineFreq;

badchanfile = inCfg.badchanfile;
badsegmfile = inCfg.badsegmfile;

bandpass = [1.3 150]; % band pass frequency
bandstop = [59 61 ; 119 121]; % band stop frequency


useTrlStd = 0;

cfg = [];
cfg.dataset = datafile;
dataRaw=ft_preprocessing(cfg);


% -------------------------------------------------------------------
% remove the Supine balancing coefficients, because it's likely incorrect.
% the call to ft_datatype_sens is needed because ft_apply_montage strips
% the chanunit and chantype for some reason
dataRaw.grad = ft_datatype_sens(ft_apply_montage(dataRaw.grad, dataRaw.grad.balance.Supine, 'inverse', 'yes', 'keepunused', 'yes'));



hcp_read_matlab(icainfofile, 'comp_class');
hcp_read_ascii(badchanfile);% badchannel
hcp_read_ascii(badsegmfile); % badsegment

badChannels = badchannel.all;

disp('bad segments concatenation')
% collect the results in a structure
maxsmp2 = max(badsegment.ica(:));
maxsmp3 = max(badsegment.manual(:));
maxsmp  = max([maxsmp2, maxsmp3]);
badsample = zeros(1,maxsmp);
if ~isempty(badsegment.ica)
    for j=1:size(badsegment.ica,1)
        badsample(badsegment.ica(j,1):badsegment.ica(j,2)) = 1;
    end
end
if ~isempty(badsegment.manual)
    for j=1:size(badsegment.manual,1)
        badsample(badsegment.manual(j,1):badsegment.manual(j,2)) = 1;
    end
end

if ~isempty(badsample)
    flank  = diff([0 badsample 0]);
    badsegment_ica=[];
    badsegment_ica(:,1) = find(flank== 1);
    badsegment_ica(:,2) = find(flank==-1) - 1;
    
    badsegment.ica       = badsegment_ica;
end

badSegments = badsegment.ica;


sel_channels={'MEG' 'MEGREF'};
if ~(isempty([badChannels{:}])),
    for ich=1:size(badChannels,2)
        sel_channels(1,ich+2)={['-' badChannels{1,ich}]};
    end
end

options = {'channels', sel_channels, 'skipped_intervals', badSegments, 'bandpass', bandpass, 'bandstop', bandstop}; % GIORGOS
dataNEURO = hcp_ICA_preprocessing(dataRaw, options);
dataNEURO_old = dataNEURO;
refdata = ft_selectdata(dataNEURO, 'channel', 'MEGREF'); 
dataNEURO = ft_selectdata(dataNEURO, 'channel', 'MEG'); 
badicacomps= [1:comp_class.class.total_ic_number];
badicacomps= setdiff(badicacomps,comp_class.class.brain_ic_vs);

for itrl=1:numel(comp_class.trial)
dataNEURO.trial{itrl}(:,:)=dataNEURO.trial{itrl}(:,:)-comp_class.topo(:,badicacomps)*comp_class.trial{itrl}(badicacomps,:);
end

% % % % % junk=dataNEURO; clear dataNEURO
% % % % % % -- Denoise with Reference Sensors --
% % % % % cfg = [];
% % % % % dataNEURO = ft_denoise_pca(cfg, junk, refdata);
% % % % % clear junk ;


data=dataNEURO;


clear dataNEURO origDataNEURO; % clear data;

% =======================================
% % if strcmp(badsegmode, 'repnan')
% %     addTrlColmn = trlNanFlag;
% % elseif strcmp(badsegmode, 'remfull')|strcmp(badsegmode, 'partial')
% %     addTrlColmn = zeros(length(data.trial), 1);
% % end
% % data.trialinfo = [data.trialinfo addTrlColmn];
% % % =======================================
% % data.grad = dataNEUROres.grad;
% % clear dataNEUROres;
% % 
% % % =========================
% % cleanTrl = data.trialinfo;
% % varinfo = whos('data');
% % if (varinfo.bytes/1000000000)>2,
% %     hcp_write_matlab(outputfile, 'data', '-v7.3');
% % else
% %     hcp_write_matlab(outputfile, 'data');
% % end
% % clear data;

end % function