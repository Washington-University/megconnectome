function data = hcp_scrubcallinfo(data)

% HCP_SCRUBCALLINFO removes the data.s.callinfo field from data.s all previous
% fields. The callinfo contains the date and time at which the analysis was
% performed. Within the HCP we don't want to share details on when a specific
% analysis was performed, especially not if the analysis is done shortly after
% acquisition of a subject.
%
% Use as
%   data = hcp_scrubcallinfo(data)
%
% This function is called by HCP_WRITE_MATLAB and HCP_WRITE_ASCII

% Copyright (C) 2013-2014 by the Human Connectome Project, WU-Minn Consortium (1U54MH091657)
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

switch class(data)
  case 'cell'
    for i=1:numel(data)
      if isstruct(data{i})
        data{i} = scrubcallinfo(data{i});
      end
    end
  case 'struct'
    for i=1:numel(data)
      data(i) = scrubcallinfo(data(i));
    end
  case 'double'
    if ishandle(data)
      scrubfigure(data);
    end
  otherwise
    % do nothing
end

end % hcp_scrubcallinfo

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function s = scrubcallinfo(s)
% this function expects a structure or structure-array as input
if isfield(s, 'callinfo')
  warning('removing callinfo');
  s = rmfield(s, 'callinfo');
end
% s.previous is expected to be either a structure or a cell-array, use recursion to deal with it
% it might also be that a fieldtrip structure is nested in another non-fieldtrip structure, hence also test all other fields
for k=1:numel(s)
  fn = fieldnames(s);
  for j=1:numel(fn)
    f = s(k).(fn{j});
    if isstruct(f)
      f = scrubcallinfo(f);
    elseif iscell(f)
      for i=1:numel(f)
        if isstruct(f{i})
          f{i} = scrubcallinfo(f{i});
        end % isstruct
      end % for i, loop over cell-array
    end % if struct or cell
    s(k).(fn{j}) = f;
  end % for j, loop over fieldnames
end % for k, loop over structure-array

end % scrubcallinfo
