function value = hcp_read_ascii(filename)

% HCP_READ_ASCII reads one or multiple MATLAB variable from an ASCII file.
% It works just like the standard MATLAB load command. Besides reading the
% variables, it also prints the MD5 checksum to screen.
%
% Use as
%   hcp_read_ascii(filename)
% to load all variables from the text file into the local workspace, or as
%    struct = hcp_read_ascii(filename)
% to load multiple variables from the text file into a single MATLAB structure.
%
% See also LOAD, HCP_WRITE_ASCII, HCP_READ_MATLAB, HCP_WRITE_MATLAB

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

filewithpath = hcp_which(filename);

if isempty(filewithpath)
  error('file "%s" does not exist', filename);
else
  filename = filewithpath;
  clear filewithpath
end

[p, f, x] = fileparts(filename);
fprintf('reading file "%s" with md5sum %s\n', [f x], hcp_md5sum(filename));
clear p f x

% this function uses some local variables with a cryptical name to prevent
% variable name colisions with the local copy of the input variables.
filename_c9e61b166b = filename;
clear filename

fid_c9e61b166b = fopen(filename_c9e61b166b, 'rt');

if fid_c9e61b166b<0
  error('cannot open file "%s"', filename_c9e61b166b);
  
else
  str_c9e61b166b = fread(fid_c9e61b166b, [1 inf], 'char=>char');
  fclose(fid_c9e61b166b);
  
  w1_c9e61b166b = []; % this variable has to exist just prior to "whos" to ensure that it gets included
  w1_c9e61b166b = whos;
  
  evalc(str_c9e61b166b);
  
  % determine the new variables that were created by the evalc command
  w2_c9e61b166b   = whos;
  name_c9e61b166b = setdiff({w2_c9e61b166b.name}, {w1_c9e61b166b.name});
  
  if nargout==0
    for indx_c9e61b166b=1:length(name_c9e61b166b)
      assignin('caller', name_c9e61b166b{indx_c9e61b166b}, eval(name_c9e61b166b{indx_c9e61b166b}));
    end
    
  elseif nargout==1
    for indx_c9e61b166b=1:length(name_c9e61b166b)
      value.(name_c9e61b166b{indx_c9e61b166b}) = eval(name_c9e61b166b{indx_c9e61b166b});
    end
    
  else
    error('multiple output arguments not supported')
  end
  
end

