function varargout = hcp_printstruct(varargin)

% HCP_PRINTSTRUCT is a robust ewrapper around the corresponding FieldTrip
% function.

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

try
  varargout{:} = printstruct(varargin{:});
catch
  warning(lasterr)
  varargout{1} = sprintf('%s could not be printed to the screen due to the error above\n', varargin{1});
end
