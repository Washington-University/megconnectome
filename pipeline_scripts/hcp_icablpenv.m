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
% 'rawdatadir/Phase1MEG/Subjects/CP10018/Experiments/CP10018_MEG/Scans/1-Rnoise_MNN_V1/Resources/4D/c,rfDC'
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
    aband=[];
end


resultprefix = sprintf('%s_%s', experimentid, scanid);

% change to the location of the processed data (input and output)
cd(pipelinedatadir)

% ensure that the output from the previous pipelines is present
% hcp_check_pipelineoutput('baddata', 'subject', subjectid, 'experiment', experimentid, 'scan', scanid);
% hcp_check_pipelineoutput('icaclass', 'subject', subjectid, 'experiment', experimentid, 'scan', scanid);
%     hcp_check_pipelineoutput('icamne', 'subject', subjectid, 'experiment', experimentid, 'scan', scanid);

% print the matlab and megconnectome version to screen for provenance
ver('megconnectome')

% print the value of all local variables to screen for provenance
w = whos;
w = {w.name};
w = setdiff(w, {'w', 'ans'});
for i=1:length(w)
    fprintf(hcp_printstruct(w{i}, eval(w{i})));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% execute the pipeline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%------------------------
% declare the input files
inputfile1 = fullfile([resultprefix,'_baddata_badsegments.txt']);
inputfile2 = fullfile([resultprefix,'_baddata_badchannels.txt']);
inputfile4 = fullfile([resultprefix,'_icaclass_vs.mat']);
inputfile5 = fullfile([resultprefix,'_icamne.mat']);

%------------------------
% declare the output file
outputfile = [resultprefix,'_icablpenv'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the channel level representation of the data,
% using the information from the qualitycheck pipeline


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the component level data (and information with respect
% to which components are the actual brain components) from
% the classification pipeline
% FIXME it seems that the output fo the classification pipeline contains
% the options: I suggest to create a text file just like the one specifying
% the badchannels/segments

hcp_read_matlab(inputfile4);
hcp_read_matlab(inputfile5, 'source');

comp=comp_class;
mixing = comp_class.topo;
if(max(size(source.val))>2)
for i = 1:size(mixing, 2)
    mixing(:, i) = mixing(:, i)/source.val(i);
    for jic=1:size(comp.trial,2)
        comp.trial{jic}(i,:)=comp.trial{jic}(i,:)*source.val(i);
    end
end
end
comp.topo=mixing;

% adjust the time axis to avoid memory problems during resampling:
% exact time information is discarded anyway

% old_time=comp.time;
% for k = 1:numel(comp.time)
%     comp.time{k} = comp.time{k} - comp.time{k}(1);
% end

if isempty(aband)
    aband=[1 2 3 4 5 6 7];
end

blp_bands = [ 1.3 4 ; 3 8 ; 6 15 ; 12.5 28.5 ; 30 75 ; 70 150 ; 1 150];
band_prefix={
    'delta'
    'theta'
    'alpha'
    'beta'
    'lowgamma'
    'highgamma'
    'whole'
    };

blp_step=20; % in ms
blp_window=400; % in ms

for ib=aband
    
    options_blp  = {'dataprefix', resultprefix, 'band_prefix', band_prefix{ib}, ...
        'blp_band', blp_bands(ib,:), 'blp_step', blp_step, 'blp_window', blp_window};
    
    source_blp = hcp_ica_blp(source,comp,options_blp);
    disp('saving icapowenv results')
    
      hcp_write_matlab([outputfile '_' band_prefix{ib}],'source_blp');
%       save([outputfile '_' band_prefix{ib}],'source_blp','-v7.3');

    % hcp_check_pipelineoutput('icapowenv', 'subject', subjectid, 'experiment', experimentid, 'scan', scanid, 'sourcemodel', sourcemodel_type, 'band', band_prefix{ib});
    
    clear source_blp;
end