function[contrastList]=contrast_Wrkmem(trialinfo)

% This function produces a default list of contrasts for the Working Memory
% paradigm. It uses the 'trialinfo' input matrix which contains all the
% trial information for a SINGLE run (FIXME:append runs). 
% Each item in the list is a structure with fields 
%   .mnemonic  % encoded mnemonic
%   .description= % cell array with detailed description of the contrast.
%   .selection - indices of contrast trials from input trialinfo.
%   .operation : type of comparison . can be 'average difference'  or 'ratio'

% Copyright (C) 2011-2013 by the Human Connectome Project, WU-Minn Consortium (1U54MH091657)
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

contrastList=contrast_Wrkmem_Basic(trialinfo);
