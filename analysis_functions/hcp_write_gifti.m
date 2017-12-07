function hcp_write_gifti(filename, data)

% HCP_WRITE_GIFTI writes a MATLAB variable to a *.gii file. 
% It only works when the input data variable is a structure that contains the fields
% pnt and tri, describing a mesh. It also stores the correspoding provenance 
% information in an XML file.
%
% Use as
%   hcp_write_gifti(filename, data)
% where the input consists of the filename and the data that should be saved.
%
% See also FT_WRITE_HEADSHAPE, HCP_WRITE_ASCII, HCP_WRITE_FIGURE, HCP_WRITE_PROVENANCE

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

% this function uses some local variables with a cryptical name to prevent
% variable name colisions with the local copy of the input variables.
filename_c9e61b166b = filename; 
clear filename

[p, f, x] = fileparts(filename_c9e61b166b);
if ~strcmp(x, '.gii')
  warning('appending the extension .gii');
  x = '.mat';
  filename_c9e61b166b = fullfile(p, [f x]);
end
clear p f x

if ~isfield(data, 'pnt') || ~isfield(data, 'tri')
    error('the input variable should contain both a ''pnt'' and a ''tri'' field');
end
ft_write_headshape(filename_c9e61b166b, data, 'format', 'gifti');

% write the corresponding provenance information
hcp_write_provenance(filename_c9e61b166b);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to prevent an eval(sprintf(...))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setvariable(name, val)
assignin('caller', name, val)
