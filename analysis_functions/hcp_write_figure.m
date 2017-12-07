function hcp_write_figure(filename, varargin)

% HCP_WRITE_FIGURE prints a MATLAB figure to a *.png file. It works just
% like the standard "print -dpng" command, but also stores the corresponding
% providence information in an XML file..
%
% Use as
%   hcp_write_figure(filename, handle, ...)
% where the input is a figure handle, or
%   hcp_write_figure(filename, ...)
% to export the current figure to a bitmap file.
%
% Optional input arguments should be specified as key-value pairs
% and may include
%  'resolution'  = number (default = 600)
%  'format'      = string, can be png, fig, svg (default is all formats)
%
% It is also possible to specify the format as a cell-array to export it
% to multiple formats.
%
% See also HCP_WRITE_ASCII, HCP_WRITE_MATLAB, HCP_WRITE_PROVENANCE

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

if isempty(varargin)
  handle = gcf;
elseif ~ishandle(varargin{1})
  handle = gcf;
elseif ishandle(varargin{1})
  handle = varargin{1};
  varargin = varargin(2:end);
end

[p, f, x] = fileparts(filename);

resolution = ft_getopt(varargin, 'resolution', 600);
if ~isempty(x)
  format   = ft_getopt(varargin, 'format', x(2:end)); % without the '.'
else
  % set the default if nothing otherwise is specified
  % format = ft_getopt(varargin, 'format', {'png', 'svg', 'fig'});
  format   = ft_getopt(varargin, 'format', 'png');
end

if ~iscell(format)
  format = {format};
end

if any(strcmp(format, 'png'))
  filename = fullfile(p, [f '.png']); % ensure that it has the right extension
  % first write it to a temporary file
  tempfile = [tempname '.png'];
  print(handle, '-dpng', sprintf('-r%d', resolution), tempfile);
  % the MATLAB written file has the date of creating in it as metadata
  % see https://github.com/Washington-University/megconnectome/issues/126
  % the solution is to read the RGB values and write them again without metadata
  fprintf('scrubbing metadata from "%s"\n', tempfile);
  rgb = imread(tempfile);
  delete(tempfile);
  savepng(rgb, filename, 4095); % use maximal compression
  % write the corresponding provenance information
  hcp_write_provenance(filename);
end

if any(strcmp(format, 'fig'))
  filename = fullfile(p, [f '.fig']); % ensure that it has the right extension
  saveas(handle, filename, 'fig')
  % write the corresponding provenance information
  hcp_write_provenance(filename);
end

if any(strcmp(format, 'svg'))
  ft_hastoolbox('plot2svg', 1);
  filename = fullfile(p, [f '.svg']); % ensure that it has the right extension
  plot2svg(filename, handle);
  % write the corresponding provenance information
  hcp_write_provenance(filename);
end
