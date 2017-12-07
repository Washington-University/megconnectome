function hcp_write_provenance(filename)

% HCP_WRITE_PROVENANCE writes an XML file with provenance information to
% complement a MATLAB or an ascii file. Thit function is called automatically
% immediately following the writing of a MATLAB or ascii file.
%
% Use as
%   hcp_write_provenance(filename)
%
% The filename should refer to the *.mat or *.txt file. The provenance
% information will be saved in a file with corresponding name but with the
% extension *.prov.xml
%
% See also HCP_WRITE_MATLAB, HCP_WRITE_ASCII

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

global hcp_default

% the expected input filename extension is *.txt, *.mat or *.png
[p, f, x] = fileparts(filename);

% get the general provenance information (e.g. user and hostname)
prov = hcp_provenance;

% when called from the megconnectome application, it passes some information about the script that is executed
if isfield(hcp_default, 'prov') && isfield(hcp_default.prov, 'script') && isfield(hcp_default.prov.script, 'filename')
  prov.script.filename = hcp_default.prov.script.filename;
end
if isfield(hcp_default, 'prov') && isfield(hcp_default.prov, 'script') && isfield(hcp_default.prov.script, 'md5sum')
  prov.script.md5sum = hcp_default.prov.script.md5sum;
end

% get the specific provenance information for the file that is written
prov.filename = [f x]; % filename with extension, exclude path
prov.md5sum   = hcp_md5sum(filename);

fprintf('writing file "%s" with md5sum %s\n', prov.filename, prov.md5sum);

% huisdier.soort = 'kanarie';
% huisdier.naam  = 'geeltje';
% huisdier.geboortedatum = '2012-10-11';
% xml = struct2xml(struct('huisdier', huisdier));

fileexchange_struct2xml(struct('megconnectome', prov), filename);

