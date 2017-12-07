function hcp_write_matlab(filename, varargin)

% HCP_WRITE_MATLAB writes one or multiple MATLAB variables to a *.mat file. 
% It works just like the standard MATLAB save command, but also stores the
% correspoding providence information in an XML file.
%
% Use as
%   hcp_write_matlab filename variable
%   hcp_write_matlab filename variable1 variable2 ...
% or in the functional form as
%   hcp_write_matlab(filename, name1, name2, ...)
% where the input consists of strings pointing to the variables that should be saved.
%
% See also SAVE, HCP_READ_MATLAB, HCP_WRITE_ASCII, HCP_WRITE_FIGURE, HCP_WRITE_PROVENANCE

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
if ~strcmp(x, '.mat')
  warning('appending the extension .mat');
  x = '.mat';
  filename_c9e61b166b = fullfile(p, [f x]);
end
clear p f x

if isempty(varargin)
  varargin = evalin('caller', 'whos');
  varargin = {varargin.name};
end

% get the anonymous input variables into the local workspace
for index_c9e61b166b=1:length(varargin)
    if ~strcmp(varargin{index_c9e61b166b}(1),'-')
        setvariable(varargin{index_c9e61b166b}, hcp_scrubcallinfo(evalin('caller', varargin{index_c9e61b166b})));
    end
end
clear index_c9e61b166b

save(filename_c9e61b166b, '-v7', varargin{:});

% write the corresponding provenance information
hcp_write_provenance(filename_c9e61b166b);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to prevent an eval(sprintf(...))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setvariable(name, val)
assignin('caller', name, val)
