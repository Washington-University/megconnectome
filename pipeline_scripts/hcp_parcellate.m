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

if ~exist('subjectid', 'var')
    subjectid = tok{end-7};
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

% ensure that the anatomy data exists
hcp_check_pipelineoutput('anatomy', 'subject', subjectid);

% first make a selection of the files that contain dense connectomes: this
% will be more straightforward once the names of the files have become
% standardized.
fnames = {'177746_MEG_3-Restin_icaimagcoh_2D_freq1.mat'};

% load the cortical sheet description for this subject
anatomydatadir = strrep(pipelinedatadir, 'analysis', 'anatomy');
fname_corticalsheet = fullfile(anatomydatadir, [subjectid, '_MEG_anatomy_sourcemodel_2d']);
hcp_read_matlab(fname_corticalsheet);

% specify the parcellations and load them
parcellations = {'expFM200AtlasFile_8k.mat', ...
                 'expFM400AtlasFile_8k.mat', ...
                 'expFM800AtlasFile_8k.mat', ...
                 'parcellations_VGD11b_8k.mat', ...
                 'RSN-networks_8k.mat'}';

parceldir = '/home/language/jansch/matlab/megconnectome/template/parcellations';
atlas     = cell(numel(parcellations),1);
for k = 1:numel(parcellations)
  parcellations{k} = fullfile(parceldir, parcellations{k});
  tmp = load(parcellations{k});
  atlas{k} = tmp.atlas;
end
  
% do the parcellation and save the results
cfg              = [];
cfg.method       = 'mean';
cfg.parcellation = 'parcellation';
for k = 1:numel(fnames)
  
  % read in the dense connectome and get it represented with a uniform
  % variable name
  tmp = hcp_read_matlab(fnames{k});
  fn  = fieldnames(tmp);
  if numel(fn)>1
    error('more than one variable per file is not supported')
  else
    data_dense = tmp.(fn{1});
    clear tmp;
  end
  
  
  for m = 1:numel(atlas)
    
    % select the potential parcellations (there can be more than one
    % parcellation per atlas
    fn  = fieldnames(atlas{m});
    sel = ~cellfun('isempty', strfind(fn, 'label'));
    
    labfn = fn(sel);
    fn    = fn(~sel);
    sel   = false(numel(fn),1);
    for p = 1:numel(labfn)
      sel = sel | strcmp(fn, labfn{p}(1:end-5));
    end
    parcelparam = fn(sel);   
    
    for p = 1:numel(parcelparam)
      sourcemodel2d.parcellation      = atlas{m}.(parcelparam{p});
      sourcemodel2d.parcellationlabel = atlas{m}.([parcelparam{p},'label']);
      data = ft_sourceparcellate(cfg, data_dense, sourcemodel2d);
      
      [pn,f,e]  = fileparts(fnames{k});
      [pn,f2,e] = fileparts(parcellations{m}); 
      outputfilename = fullfile(pipelinedatadir, [f,'_',f2,'_',parcelparam{p}]);
      hcp_write_matlab(outputfilename, 'data');
      clear data;
    end % for number of parcellations in a given atlas
  end % for atlas
end % for file

%hcp_check_pipelineoutput('tmegpreproc', 'subject', subjectid, 'experiment', experimentid, 'scan', scanid,'datagroup',trlInfo.lockNames);
