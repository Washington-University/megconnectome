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

clear tok

if ~exist('pipelinedatadir', 'var')
    pipelinedatadir = hcp_pathdef;
end

% print the value of all local variables to screen for provenance
w = whos;
w = {w.name};
w = setdiff(w, {'w', 'ans'});
for i=1:length(w)
    % fprintf(hcp_printstruct(w{i}, eval(w{i})));
end

% change to the location of the processed data (input and output)
cd(pipelinedatadir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% execute the pipeline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hcp_check_pipelineoutput('rmegpreproc', 'subject', subjectid, 'experiment', experimentid, 'scan', scanid);

resultprefix  = sprintf('%s_%s', experimentid, scanid);
inputfile     = [resultprefix,'_rmegpreproc'];
outputfile    = [resultprefix,'_powavg'];

hcp_read_matlab(inputfile);

cfg = [];
cfg.foilim = [0 150];
cfg.method = 'mtmfft';
cfg.taper = 'hanning';
cfg.keeptrials = 'no';
freq = ft_freqanalysis(cfg, data);

hcp_write_matlab(outputfile, 'freq');

cfg = [];
cfg.parameter = 'powspctrm';
cfg.operation = 'log10';
log10freq = ft_math(cfg, freq);

cfg = [];
cfg.layout = '4D248.mat';
cfg.xlim = [5 55];
cfg.axes = 'no';
cfg.box = 'yes';

figure
ft_multiplotER(cfg, log10freq);
h = title(resultprefix); set(h, 'interpreter', 'none');
hcp_write_figure([outputfile '_multiplot']);
close

% cfg = [];
% cfg.xlim = [5 50];
% figure
% ft_singleplotER(cfg, log10freq);

figure
plot(log10freq.freq, log10freq.powspctrm);
axis([5 55 -30 -26]);
xlabel('frequency (Hz)');
ylabel('log10(power)');
grid on
h = title(resultprefix); set(h, 'interpreter', 'none');
hcp_write_figure([outputfile '_singleplot']);
close

hcp_check_pipelineoutput('powavg', 'subject', subjectid, 'experiment', experimentid, 'scan', scanid);
