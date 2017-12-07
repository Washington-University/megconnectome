function dir = hcp_select_datadir(dirlist)

% HCP_SELECT_DATADIR checks for the presence of each directory from a list
% and returns the first from the list that is present. It helps in managing
% the HCP data on multiple computers.
%
% Use as
%   dir = hcp_check_datadir(dirlist)
% where the input dirlist is a cell array with strings and the output is a
% single string.
%
% Example:
%
% dataprefix = hcp_check_datadir({ ...
%   '/data/jansch/HCP/database' ...
%   '/data/roboos/webdav' ...
%   '/scratch/r.oostenveld/webdav' ...
%   });

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

dir = [];
for i=1:length(dirlist)
  if isdir(dirlist{i})
    dir = dirlist{i};
    return
  end
end
