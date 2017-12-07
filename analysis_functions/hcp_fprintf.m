function fid = hcp_fprintf(fid, varargin)

% HCP_FPRINTF works just like FOPEN/FPRINTF/FCLOSE but also writes an XML
% provenance file as soon as the text file is closed.
%
% Example
%  fid = hcp_fopen(filename, 'w')
%  hcp_fprintf(fid, 'the number is %d\n', 42);
%  hcp_fclose(fid)
%
% See also HCP_FOPEN, HCP_FPRINTF, HCP_FCLOSE, HCP_WRITE_PROVENANCE

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

persistent remember
if isempty(remember)
  remember = {};
end

if ischar(fid) && strcmp(fid, 'open')
  fid = fopen(varargin{:});
  remember{fid} = varargin; % remember the file
  
elseif ischar(fid) && strcmp(fid, 'close')
  fid = varargin{1};
  hcp_write_provenance(remember{fid}{1}); % write the provenance
  remember(fid) = []; % forget about the file
  fclose(fid);
  
else
  fprintf(fid, varargin{:});
  
end
