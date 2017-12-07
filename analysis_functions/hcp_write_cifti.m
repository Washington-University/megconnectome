function hcp_write_cifti(filename, source, varargin)

% HCP_WRITE_CIFTI writes a MATLAB source estimate of activity or connectivity to a
% CIFTI format 2 file. It works just like ft_write_cifti, but also stores the
% corresponding providence information in an XML file.
%
% Use as
%   hcp_write_cifti(filename, source, ...)
% where the input source should be a valid FieldTrip source-level data structure.
%
% Additional input arguments should be specified as key-value pairs
% and may include
%   'parameter'      = string, fieldname that contains the functional data
%   'brainstructure' = string, fieldname that describes the brain structures (default = 'brainstructure')
%   'parcellation'   = string, fieldname that describes the parcellation (default = 'parcellation')
%   'precision'      = string, can be 'single', 'double', 'int32', etc. (default ='single')
%   'writesurface'   = boolean, can be false or true (default = true)
%
% See also HCP_WRITE_ASCII, HCP_WRITE_MATLAB, HCP_WRITE_PROVENANCE

% Copyright (C) 2014 by the Human Connectome Project, WU-Minn Consortium (1U54MH091657)
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

parameter       = ft_getopt(varargin, 'parameter');
brainstructure  = ft_getopt(varargin, 'brainstructure');
parcellation    = ft_getopt(varargin, 'parcellation');
precision       = ft_getopt(varargin, 'precision');
writesurface    = ft_getopt(varargin, 'writesurface');

%% some data handling, this part is copied from ft_sourcewrite

% keep the transformation matrix
if isfield(source, 'transform')
  transform = source.transform;
elseif isfield(source, 'brainordinate') && isfield(source.brainordinate, 'transform')
  transform = source.brainordinate.transform;
else
  transform = [];
end

if isfield(source, 'brainordinate')
  % it is a parcellated source representation, i.e. the main structure one channel for each parcel
  brainordinate = source.brainordinate;
  source = rmfield(source, 'brainordinate');
  
  % split them and check individually
  source        = ft_checkdata(source, 'datatype', {'timelock', 'freq', 'chan'}, 'feedback', 'yes');
  brainordinate = ft_checkdata(brainordinate, 'datatype', 'parcellation', 'parcellationstyle', 'indexed', 'hasunit', 'yes');
  
  % merge them again
  source = copyfields(brainordinate, source, setdiff(fieldnames(brainordinate), {'cfg'}));
else
  source = ft_checkdata(source, 'datatype', 'source', 'hasunit', true, 'feedback', 'yes');
end

% keep the transformation matrix
if ~isempty(transform)
  source.transform = transform;
end

%% fix the filename, this section is consistent with ft_write_cifti

switch getdimord(source, parameter)
  case {'pos' 'pos_scalar'}
    % NIFTI_INTENT_CONNECTIVITY_DENSE_SCALARS
    extension   = '.dscalar.nii';
  case 'pos_pos'
    % NIFTI_INTENT_CONNECTIVITY_DENSE
    extension = '.dconn.nii';
  case 'pos_time'
    % NIFTI_INTENT_CONNECTIVITY_DENSE_SERIES
    extension = '.dtseries.nii';
  case 'pos_freq'
    % NIFTI_INTENT_CONNECTIVITY_DENSE_SERIES
    extension = '.dtseries.nii';
  case {'chan' 'chan_scalar'}
    % NIFTI_INTENT_CONNECTIVITY_PARCELLATED_SCALARS
    extension   = '.pscalar.nii';
  case 'chan_chan'
    % NIFTI_INTENT_CONNECTIVITY_PARCELLATED
    extension = '.pconn.nii';
  case 'chan_time'
    % NIFTI_INTENT_CONNECTIVITY_PARCELLATED_SERIES
    extension = '.ptseries.nii';
  case 'chan_freq'
    % NIFTI_INTENT_CONNECTIVITY_PARCELLATED_SERIES
    extension = '.ptseries.nii';
  case 'pos_pos_time'
    % NIFTI_INTENT_CONNECTIVITY_PARCELLATED_SERIES
    extension = '.dconnseries.nii';    
  otherwise
    error('unsupported dimord "%s"', dimord);
end % switch

[p, f, x] = fileparts(filename);
if isequal(x, '.nii')
  filename = fullfile(p, f); % strip the extension
end

[p, f, x] = fileparts(filename);
if any(isequal(x, {'.dtseries', '.ptseries', '.dconn', '.pconn', '.dscalar', '.pscalar', '.dconnseries'}))
  filename = fullfile(p, f); % strip the extension
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the extension does not need to be added for ft_write_cifti

%% add the full cifti extension to the filename
%[p, f, x] = fileparts(filename);
%filename = fullfile(p, [f x extension]);

%% write the nii and xml file

% brainstructure should represent the global anatomical structure, such as CortexLeft, Thalamus, etc.
% parcellation should represent the detailled parcellation, such as BA1, BA2, BA3, etc.
ft_write_cifti(filename, source, 'parameter', parameter, 'brainstructure', brainstructure, 'parcellation', parcellation, 'precision', precision, 'writesurface', writesurface);

% by now the filename has some additional
% dscalar,dtseries,dconn,pscalar,ptseries,or pconn, just before the .nii

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the extension needs to be added for provenance writing

% add the full cifti extension to the filename
[p, f, x] = fileparts(filename);
filename = fullfile(p, [f x extension]);

% write the corresponding provenance information
hcp_write_provenance(filename);

