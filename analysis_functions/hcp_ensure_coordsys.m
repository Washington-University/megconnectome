function object = hcp_ensure_coordsys(object, transform, desired)

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

if ~isfield(object, 'coordsys')
  % this one returns an updated structure
  object = ft_determine_coordsys(object);
end

if strcmp(object.coordsys, desired)
  % nothing to do
else
  T = transform.(sprintf('%s2%s', object.coordsys, desired));
  object = ft_transform_geometry(T, object);
  object.coordsys = desired;
end
