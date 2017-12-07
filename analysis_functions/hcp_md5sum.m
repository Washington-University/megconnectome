function hash = hcp_md5sum(filename)

% HCP_MD5SUM finds the file on the path and computes the MD5 sum
%
% Use as
%   hash = hcp_md5sum(filename)
%
% See also HCP_READ_ASCII, HCP_WRITE_ASCII, HCP_READ_MATLAB,
% HCP_WRITE_MATLAB

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

% the MATLAB exists function will look for the file anywhere on the path
% which returns the actual full file with path

filewithpath = hcp_which(filename);

if isempty(filewithpath)
  error('file "%s" does not exist', filename);
else
  filename = filewithpath;
  clear filewithpath
end

% status will be 0 once the md5sum is known
status = -1;

if status~=0
  [status, str] = system(sprintf('md5sum %s', filename));
  % this looks like
  % ee338cf30b7e10f249b61dade77339d1 filename
  hash = strtok(str);
end

if status~=0
  [status, str] = system(sprintf('md5 %s', filename));
  % this looks like
  % MD5 (filename) = ee338cf30b7e10f249b61dade77339d1
  [tok, rem] = strtok(str, '=');
  hash = strtrim(rem(2:end));
end

if status~=0
  warning('failed to compute MD5 hash for file "%s', filename);
  hash = 'unknown';
end
