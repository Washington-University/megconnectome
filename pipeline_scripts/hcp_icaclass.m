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

% print the matlab and megconnectome version to screen for provenance
ver('megconnectome')

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

% ensure that the expected input data exists
hcp_check_pipelineoutput('baddata', 'subject', subjectid, 'experiment', experimentid, 'scan', scanid);

cfg = [];
cfg.dataset = filename;
dataraw=ft_preprocessing(cfg);

%%%%%%%%%%% do ICA %%%%%%%%%%%%%%%%%
badsegments = hcp_read_ascii([resultprefix '_baddata_badsegments.txt']);
%     badsegments = badsegments.badsegment.all;

disp('bad segments concatenation')
% collect the results in a structure
maxsmp2 = max(badsegments.badsegment.ica(:));
maxsmp3 = max(badsegments.badsegment.manual(:));
maxsmp  = max([maxsmp2, maxsmp3]);
badsample = zeros(1,maxsmp);
if ~isempty(badsegments.badsegment.ica)
    for j=1:size(badsegments.badsegment.ica,1)
        badsample(badsegments.badsegment.ica(j,1):badsegments.badsegment.ica(j,2)) = 1;
    end
end
if ~isempty(badsegments.badsegment.manual)
    for j=1:size(badsegments.badsegment.manual,1)
        badsample(badsegments.badsegment.manual(j,1):badsegments.badsegment.manual(j,2)) = 1;
    end
end

if ~isempty(badsample)
    flank  = diff([0 badsample 0]);
    badsegment_ica=[];
    badsegment_ica(:,1) = find(flank== 1);
    badsegment_ica(:,2) = find(flank==-1) - 1;
    
    badsegments.badsegment.ica       = badsegment_ica;
end

badsegments = badsegments.badsegment.ica;

badchannels = hcp_read_ascii([resultprefix '_baddata_badchannels.txt']);
badchannels = badchannels.badchannel.all;

sel_channels={'MEG'};
grad=dataraw.grad;

if ~(isempty([badchannels{:}])),
    for ich=1:size(badchannels,2)
        sel_channels(1,ich+1)={['-' badchannels{1,ich}]};
    end
end
bandpass = [1.3 150]; % band pass frequency
if strcmp(subjectid,'CP10128') || strcmp(subjectid,'CP10129')
    bandstop = [49 51 ; 99 101]; % band stop frequency
else
    bandstop = [59 61 ; 119 121]; % band stop frequency
end

if ~exist('ica_iter', 'var')
    if ~isempty(regexp(scanid,'Rest'))
        ica_iter=20;
    else
        ica_iter=20;
    end
end
ica_iter

options = {'resultprefix', resultprefix, 'channels', sel_channels, 'ica_iterations', ica_iter, 'skipped_intervals', badsegments, 'bandpass', bandpass, 'bandstop', bandstop}; % GIORGOS
data_meg = hcp_ICA_preprocessing(dataraw, options);

% do the computation
[iteration datain]= hcp_ICA_unmix(data_meg, options);
%provvisory
hcp_write_matlab([resultprefix '_icaiteration'],'iteration','options')

%try
%----------------------------------------------
% The hcp_exgmontage.m function should return the EOG ECG reference
% channels
montage       = hcp_exgmontage(subjectid, experimentid, scanid);

if ~isempty(montage.labelorg) % This check is performed for cases where no Electrodes have been recorded (i.e Glasgow data)
    hasELEC=1;
    dataELEC=ft_selectdata(dataraw,'channel',montage.labelorg);
    dataELECnew=ft_apply_montage(dataELEC,montage);
else
    hasELEC=0;
    dataELECnew=[];
    display(['NO ELECTRODES WERE PROVIDED BY EXGMONTAGE for ', subjectid, '_', experimentid, '_', scanid]);
end



%ref     = montage.labelorg';
%====================================================
%====================================================
%====================================================
%====================================================
%======================================================
% ---- Check if the data is from MOTOR TASK for which no ECG EOG channels
% are available and construct pseudo ECG EOG channels based on topology
% templates. As we do not use EEG anymore in Phase 2 , MOTOR task data
% should also have EO/CG channels and this part should not be necessary
oldMotorExperimentIds={'CP10128_MEG_v1',... % Glasgow scan
    'CP10129_MEG_v1',... % Glasgow scan
    'CP10141_MEG_v1',... % There are 2 scans with 2 difference EMG electrode locations. For each scan there are 3 blocks for each hand and foot with 12 trials in each block. So each condition has only 36 trials.
    'CP10138_MEG',...    % There are 6 different scan with variable ISI and time jittering in some. For each scan there are 3 blocks for each hand and foot with beween 8 to 12 trials in each block. So each condition has only 24 to 36 trials.
    'CP10167_MEG',...    % 1 scan with the latest protocol. 8 blocks per hand and foot with 10 trials each. 1200 msec ISI. However there is no MRI
    'CP10113_MEG_v2'...  % There are 2 different scans. NOt clear what the difference is between the two. There are 8 blocks for each hand and foot, with 10 trials per block. Unfortunately EMG signal is bad for both Hands and trials cannot be extracted based on the EMG signal
    };

isMotorTask= ~isempty(regexp(scanid,'Motor'));
isScanWithOldMotor= ismember(experimentid,oldMotorExperimentIds);
if (isMotorTask && isScanWithOldMotor) || (hasELEC==0)
    % TODO: In the case that there are no electrodes at all(i.e. Glasgow) or
    % there are no ECOG(i.e. Motor Tasks) pseudo reference channels are
    % extracted from the ICA decomposition in the variable iteration.
    % A decision need to be made how these cases (mainly the Motor tasks)
    % will be handled. Probably in these cases the classification in Brain
    % and artifact components should be performed without using any reference channels
    % and the ICs should be examined manually in order to assign heart and
    % eye artifact related components.
    templateICAFile='templatecomps.mat';% This fle contains the template topologies for EOCG ICs and is locates in the sandbox folder
    
    cfg=[];
    cfg.templateICAFile=templateICAFile;
    cfg.outputfile=[experimentid,'_',scanid,'_icaclass_virtchan'];
    %------------------------------------------------------------
    
    
    [pseudodataECG,pseudodataEOG] = hcp_pseudoeocgfromtopo(cfg,iteration,datain);
    %------------------------
    icaref=[pseudodataECG.label,pseudodataEOG.label];
    ref_data=ft_appenddata([],pseudodataECG,pseudodataEOG);
    options     = ft_setopt(options, 'ref_channels', icaref);
    options     = ft_setopt(options, 'subject', resultprefix);
    options_ref = ft_setopt(options, 'channels',icaref);
    options     = ft_setopt(options, 'grad', grad); % GIORGOS
    clear templateICAFile pseudodataECG pseudodataEOG
else
    
    [tmpIndx1,tmpIndx2]=match_str(dataELECnew.label,{'ECG' 'VEOG' 'HEOG'});
    if ~isempty(tmpIndx1)
        icaref=dataELECnew.label(tmpIndx1);
    else
        icaref=[];
        error(['NO ECG VEOG ot HEOG found in the dataset . Check data and/or montage for  ',experimentid,'_',scanid]);
        return;
    end
    clear tmpIndx1 tmpIndx2;
    
    options     = ft_setopt(options, 'ref_channels', icaref);
    options     = ft_setopt(options, 'subject', resultprefix);
    options_ref = ft_setopt(options, 'channels',icaref);
    options     = ft_setopt(options, 'grad', grad); % GIORGOS
    ref_data    = hcp_ICA_preprocessing(dataELECnew, options_ref);
end
%========================================================
clear icaref options_ref;

comp_class = hcp_ICA_RMEG_classification(ref_data,options,iteration,datain);

hcp_ICA_plotclassification(comp_class,options,datain);

hcp_write_matlab([resultprefix '_icaclass'],'comp_class','options')

vs=[];
vs.good=[];
vs.bad=[];
vs.total_ic_number=comp_class.class.total_ic_number;
vs.brain_ic_number=comp_class.class.brain_ic_number;
vs.brain_ic=comp_class.class.brain_ic;
vs.ecg_eog_ic= comp_class.class.ecg_eog_ic;
vs.flag=0;
vs.physio=[];
hcp_write_ascii(sprintf('%s_icaclass.txt', resultprefix), 'vs');

clear iteration datain comp_class options;
clear resultprefix dataraw badsegments badchannels sel_channels grad bandpass bandstop ica_iter;
clear data_meg montage dataELEC dataELECnew ref_data;
clear oldMotorExperimentIds isMotorTask isScanWithOldMotor;
clear cfg;

% ensure that the expected output files were created
hcp_check_pipelineoutput('icaclass', 'subject', subjectid, 'experiment', experimentid, 'scan', scanid);
%catch me
%    save([resultprefix '_montageerror'],'resultprefix')
%end
