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

smodel_type = {'2D';'3D'}; dimindx=1; % 1 for 2D cortical sheet and 2 for 3D gird
griddim = {'4mm';'6mm';'8mm'}; gridindx=1; % $ 1,2,3 for 3D 4mm,6mm and 8mm grid

if(strcmp(smodel_type{dimindx},'2D'))
  sourcemodel_type=smodel_type{dimindx};
elseif(strcmp(smodel_type{dimindx},'3D'))
  gridname = 'sourcemodel3d';
  sourcemodel_type=[smodel_type{dimindx} griddim{gridindx}];
end

%------------------------
% declare the output file
resultprefix = sprintf('%s_%s', experimentid, scanid);
outputfile = [resultprefix,'_icaimagcoh_connect_' sourcemodel_type];

freque = [20];

for f=freque
  
  hcp_read_matlab([resultprefix,'_icaimagcoh_' sourcemodel_type,'_freq',num2str(f)]);
  imagesc(imagcoh.mimspctrm)
  caxis([0 0.1])
  h1 = gcf;
  % set(h1, 'paperposition', [1 1 10 7]);
  % hcp_write_figure(imgname, h0, 'format', 'png')
  % hcp_write_figure(imgname, h1)
  hcp_write_figure([outputfile '_freq',num2str(f)], h1, 'format', 'png')
  %  hcp_write_figure(outputfile, h1, 'format', 'fig')
  close(h1)
  
  clear imagcoh
end % nfreq

