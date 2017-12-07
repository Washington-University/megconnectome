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

if ~exist('aband', 'var')
    aband=[1:8];
end

resultprefix = sprintf('%s_%s', experimentid, scanid);

% change to the location of the processed data (input and output)
cd(pipelinedatadir)

% ensure that the output from the previous pipelines is present
hcp_check_pipelineoutput('baddata', 'subject', subjectid, 'experiment', experimentid, 'scan', scanid);
hcp_check_pipelineoutput('icaclass', 'subject', subjectid, 'experiment', experimentid, 'scan', scanid);
hcp_check_pipelineoutput('icaclass_qc', 'subject', subjectid, 'experiment', experimentid, 'scan', scanid);
hcp_check_pipelineoutput('anatomy', 'subject', subjectid);

% print the matlab and megconnectome version to screen for provenance
ver('megconnectome')

% print the value of all local variables to screen for provenance
w = whos;
w = {w.name};
w = setdiff(w, {'w', 'ans'});
for i=1:length(w)
    fprintf(hcp_printstruct(w{i}, eval(w{i})));
end

scanmnem='Restin';
fname=filename;
trialFunName=['trialfun_',scanmnem];

hdr = ft_read_header(fname);
Fsample=hdr.Fs;

montage         = hcp_exgmontage(subjectid, experimentid, scanid);

tmpCfg=[];
tmpCfg.lockMnem= 'BREST';
tmpCfg.cutMode= 'blocks';
tmpCfg.prestimTime= 0;
tmpCfg.poststimTime= 0;
tmpCfg.montage= montage;
tmpCfg.trialDuration=(hdr.nSamples-0.2).*(1/Fsample);

trlCfg=[];

trlCfg.trialdef        = tmpCfg;
trlCfg.dataset        = fname;
trlCfg.trialfun        = trialFunName;

eval(['[trl]=',trialFunName,'(trlCfg);']);


%%
badchansuffix='baddata_badchannels'; % mnemonic for file that contains bad channel and segment info
badsegmsuffix='baddata_badsegments'; % mnemonic for file that contains bad channel and segment info
icainfosuffix='icaclass_vs'; % mnemonic to save results

iCaseGroup=1;
trlInfo=[];
trlInfo.lockNames{iCaseGroup}=trlCfg.trialdef.lockMnem;
trlInfo.lockTrl{iCaseGroup}=[trl 1+0*trl(:,1)];

Ncasegroups=length(trlInfo.lockTrl);
cleanTrlInfo=trlInfo;
datagroupid=trlInfo.lockNames{iCaseGroup};
savesuffix_general='testpartial';

cfg  =  [];
cfg.badchanfile     = [resultprefix,'_',badchansuffix,'.txt'];
cfg.badsegmfile     = [resultprefix,'_',badsegmsuffix,'.txt'];
cfg.icainfofile     = [resultprefix,'_',icainfosuffix];
cfg.datafile        = fname;
cfg.badsegmode      = 'partial'; %'repnan'; %tmpCfg.badsegmode='remfull'; % This is used by tmeg_preproc to know if trials should be fully removed, or bad segments should be replaced by nans(i.e. 'BU' case)
cfg.montage         = montage; %hcp_exgmontage(subjectid, experimentid, scanid);
cfg.outputfile      = [experimentid,'_',scanid,'_',savesuffix_general,'_',datagroupid];
cfg.lineFreq        = [60 120];

cfg.trl=trlInfo.lockTrl{iCaseGroup};

[data]=hcp_extract_allfromrun_rest(cfg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% execute the pipeline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
multiscanid='Restin';
datagroupid='BREST';
savesuffix_general='testpartial';
datafile      = [experimentid,'_',scanid,'_',savesuffix_general,'_',datagroupid];
gridtype='2D';
%--------------------------------------

% load(datafile,'data');

gridAllLF=[];

Fsample=data.fsample;
if isempty(gridAllLF),
    grad=ft_convert_units(data.grad,'cm');
    fileGrid2D=[experimentid,'_anatomy_sourcemodel_2d.mat'];     % Subject's 2D cortical sheet source model
    fileVol=[experimentid,'_anatomy_headmodel.mat'];             % Subject's brain volume
    
    
    hcp_read_matlab(fileGrid2D,'sourcemodel2d');
    grid=ft_convert_units(sourcemodel2d,'cm');
    
    hcp_read_matlab(fileVol,'headmodel');
    vol=ft_convert_units(headmodel,'cm');
    
    %----- Create Leadfields ---------------------------------------
    allChansMEG=ft_channelselection({'MEG'},grad.label);
    cfg=[];
    cfg.grid=grid; % Grid for Individual's Brain in MEG sensor space
    cfg.vol=vol;
    cfg.grad=grad;
    cfg.reducerank      = 2; %(default = 3 for EEG, 2 for MEG)
    cfg.normalize       = 'yes' ; %Normalise Leadfield: 'yes' for beamformer
    cfg.normalizeparam  = 1;      %depth normalization parameter (default = 0.5).
    cfg.feedback='no';
    cfg.channel=allChansMEG;
    if strcmp(gridtype,'2D') %JUST TO MAKE SURE ENTIRE CORTICAL SHEET IN SOURCE SPACE
        cfg.inwardshift=-1;
    end
    gridAllLF= ft_prepare_leadfield(cfg);
    gridAllLF.label=allChansMEG;
end

%------------------------
% declare the output file
outputfile = [resultprefix,'_bfblpenv'];

blp_bands = [ 1.3 4.5 ; 3 9.5 ; 6.3 16.5 ; 12.5 29 ; 22.5 39 ; 30 55 ;  45 82 ; 70 125 ; 1.3 150];
band_prefix={
    'delta'
    'theta'
    'alpha'
    'betalow'
    'betahigh'
    'gammalow'
    'gammamid'
    'gammahigh'
    'whole'
    };

blp_step=20; % in ms
blp_window=400; % in ms

lambda=75;
lambdastr=[num2str(lambda) '%'];

for ib=aband
    
    options_blp  = {'dataprefix', resultprefix, 'band_prefix', band_prefix{ib}, ...
        'blp_band', blp_bands(ib,:), 'blp_step', blp_step, 'blp_window', blp_window,'bf_lambda',lambdastr};
    
    
    source_blp = hcp_bf_blp(gridAllLF,data,vol,grad,options_blp);
    disp('saving icapowenv results')
    
    % save it as a cifti file
    hcp_write_cifti([outputfile,'_' band_prefix{ib} '.power'],source_blp, 'parameter', 'power', 'type', 'dtseries');
    
    % save it as a matlab file as well
    hcp_write_matlab([outputfile '_' band_prefix{ib}],'source_blp');
    
    hcp_check_pipelineoutput('bfblpenv', 'subject', subjectid, 'experiment', experimentid, 'scan', scanid, 'band', band_prefix{ib});

    clear source_blp;
end
