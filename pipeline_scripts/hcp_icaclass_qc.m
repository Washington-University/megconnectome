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
hcp_check_pipelineoutput('icaclass', 'subject', subjectid, 'experiment', experimentid, 'scan', scanid);

%%%%%%%%%%% read vs components %%%%%%%%%%%%%%%%%
if exist([resultprefix '_icaclass_vs.txt'],'file');
    hcp_read_ascii([resultprefix '_icaclass_vs.txt']);
elseif exist([resultprefix '_icaclass.txt'],'file');
    hcp_read_ascii([resultprefix '_icaclass.txt']);
else
    vs=[];
    vs.good=[];
    vs.bad=[];
    vs.physio=[];
    hcp_write_ascii(sprintf('%s_icaclass_vs.txt', resultprefix), 'vs');
end

good = vs.good;
bad = vs.bad;
if(isfield(vs,'physio'))
    physio=vs.physio;
else
    physio=[];
    vs.physio=[];
end
vs.flag=1;

disp(['loading run ' resultprefix])
load([resultprefix '_icaclass'])

disp(['brain ics -> ' num2str(comp_class.class.brain_ic)])
comp_class.class.brain_ic_vs=comp_class.class.brain_ic;
for i=1:size(good,2)
    comp_class.class.brain_ic_vs(1,end+1)=good(1,i);
    [junk indxb]=find(comp_class.class.ecg_eog_ic==good(1,i));
    comp_class.class.ecg_eog_ic(indxb)=[];
end
comp_class.class.brain_ic_vs=unique(comp_class.class.brain_ic_vs);

for i=1:size(bad,2)
    [junk indxb]=find(comp_class.class.brain_ic_vs==bad(1,i))
    comp_class.class.brain_ic_vs(indxb)=[];
end
comp_class.class.brain_ic_vs_number=size(comp_class.class.brain_ic_vs,2);

for i=1:size(physio,2)
    comp_class.class.ecg_eog_ic(1,end+1)=physio(1,i);
end
comp_class.class.ecg_eog_ic=unique(comp_class.class.ecg_eog_ic);
comp_class.class.physio=physio;

disp(['total ics number -> ' num2str(comp_class.class.total_ic_number)])
vs.total_ic_number=comp_class.class.total_ic_number;
disp(['brain ics vs number -> ' num2str(comp_class.class.brain_ic_vs_number)])
vs.brain_ic_vs_number=comp_class.class.brain_ic_vs_number;
disp(['brain ics vs -> ' num2str(comp_class.class.brain_ic_vs)])
vs.brain_ic_vs=comp_class.class.brain_ic_vs;

vs.ecg_eog_ic= comp_class.class.ecg_eog_ic;

hcp_write_ascii(sprintf('%s_icaclass_vs.txt', resultprefix), 'vs');

disp(['modified vs-brain = ' num2str(good) ' vs-art = ' num2str(bad)])

disp(['saving ' resultprefix '_icaclass_vs'])
hcp_write_matlab([resultprefix '_icaclass_vs'],'comp_class','options')


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

options = {'subject',subjectid,'resultprefix', resultprefix, 'channels', sel_channels, 'skipped_intervals', badsegments, ...
    'grad', grad, 'bandpass', bandpass, 'bandstop', bandstop,'textsize',24};
data_meg = hcp_ICA_preprocessing(dataraw, options);

hcp_ica_plotreport(comp_class,options,data_meg);

