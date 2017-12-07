function value = hcp_read_matlab(filename, varargin)

% HCP_READ_MATLAB reads one or multiple MATLAB variable from a *.mat file.
% It works just like the standard MATLAB load command. Besides reading the
% variables, it also prints prints the MD5 checksum to screen.
%
% Use as
%   hcp_read_matlab(filename)
% to load all variables from the text file into the local workspace, or as
%    struct = hcp_read_matlab(filename)
% to load multiple variables from the text file into a single MATLAB structure.
%
% See also LOAD, HCP_WRITE_MATLAB, HCP_READ_ASCII, HCP_WRITE_ASCII

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

filewithpath = hcp_which(filename);

if isempty(filewithpath)
  error('file "%s" does not exist', filename);
else
  filename = filewithpath;
  clear filewithpath
end

[p, f, x] = fileparts(filename);
if ~strcmp(x, '.mat')
  warning('appending the extension .mat');
  x = '.mat';
  filename = fullfile(p, [f x]);
end
fprintf('reading file "%s" with md5sum %s\n', [f x], hcp_md5sum(filename));
clear p f x

value = load(filename, varargin{:});

if ~nargout
  % copy the variables into the caller workspace
  fn = fieldnames(value);
  for i=1:length(fn)
    assignin('caller', fn{i}, value.(fn{i}));
  end
  clear value
else
  % return the variables as a structure
end
