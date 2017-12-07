function atlas = hcp_read_atlas(inputfile)

% HCP_READ_ATLAS converts a set of gii label-files containing a parcellation into a
% FieldTrip compatible parcellation structure that can be saved as a mat file. The
% parcellations corresponding to the left and right hemisphere are concatenated, LEFT
% HEMISPHERE FIRST, and the indices for the labels for the right hemisphere parcels
% are adjusted accordingly. The parcel labels are moreover prefixed with L_ and R_.
%
% The input should be a cell-array with the file names for the left and right
% hemisphere. For example
%
% labels = hcp_read_atlas({'parcellations_VGD11b.L.8k_fs_LR.label.gii', 'parcellations_VGD11b.R.8k_fs_LR.label.gii'})
%
%          parcellation1: [15684x1 double]
%     parcellation1label: {109x1 cell}
%          parcellation2: [15684x1 double]
%     parcellation2label: {86x1 cell}
%          parcellation3: [15684x1 double]
%     parcellation3label: {4x1 cell}
%             hemisphere: [15684x1 double]
%        hemispherelabel: {2x1 cell}
%
% This can be combined with
%
% geometry = ft_read_headshape({'Conte69.L.midthickness.8k_fs_LR.surf.gii', 'Conte69.R.midthickness.8k_fs_LR.surf.gii'})
%
%                 pnt: [15684x3 double]
%                 tri: [31360x3 double]
%          hemisphere: [15684x1 double]
%     hemispherelabel: {2x1 cell}
%
% and subsequently combined into a parcellation as
%
% atlas = hcp_mergestruct(labels, geometry)
% 
%          parcellation1: [15684x1 double]
%     parcellation1label: {109x1 cell}
%          parcellation2: [15684x1 double]
%     parcellation2label: {86x1 cell}
%          parcellation3: [15684x1 double]
%     parcellation3label: {4x1 cell}
%               coordsys: 'unknown'
%        hemispherelabel: {2x1 cell}
%             hemisphere: [15684x1 double]
%                    pnt: [15684x3 double]
%                    tri: [31360x3 double]
%                    
% See also FT_READ_ATLAS, FT_READ_HEADSHAPE, FT_DATATYPE_PARCELLATION

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

atlasleft  = ft_read_atlas(inputfile{1});
atlasright = ft_read_atlas(inputfile{2});

fnames = fieldnames(atlasleft);
sel = false(size(fnames));
for k = 1:numel(fnames)
  sel(k) = isfield(atlasleft, [fnames{k} 'label']);
end
fnames = fnames(sel);

% update the right hemisphere indices
for k = 1:numel(fnames)
  atlasright.(fnames{k}) = atlasright.(fnames{k})+max(atlasleft.(fnames{k})(:));
end

% prefix the labels
for k = 1:numel(fnames)
  pname = [fnames{k},'label'];
  if all(strncmp(atlasleft.(pname), 'L_', 2))
    % don't prefix
  else
    for m = 1:numel(atlasleft.(pname))
      atlasleft.(pname){m} = ['L_',atlasleft.(pname){m}];
    end
  end
  if all(strncmp(atlasright.(pname), 'R_', 2))
    % don't prefix
  else
    for m = 1:numel(atlasright.(pname))
      atlasright.(pname){m} = ['R_',atlasright.(pname){m}];
    end
  end
end

% concatenate the parcellations
atlas = atlasleft;
for k = 1:numel(fnames)
  atlas.(fnames{k}) = cat(1,atlasleft.(fnames{k}),atlasright.(fnames{k}));
  atlas.([fnames{k},'label']) = cat(1,atlasleft.([fnames{k},'label']),atlasright.([fnames{k},'label']));
end

% add the indices and labels for the hemispheres, just like ft_read_headshape
atlas.hemispherelabel = inputfile(:);
atlas.hemisphere      = cat(1, 1*ones(length(atlasleft.(fnames{1})),1), 2*ones(length(atlasleft.(fnames{1})),1));
