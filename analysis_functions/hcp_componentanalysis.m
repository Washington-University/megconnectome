function comp = hcp_componentanalysis(cfg, data)

% HCP_COMPONENTANALYSIS performs a fastica decomposition of the
% input data using the FT_COMPONENTANALYSIS. In case fastica fails
% to converge for the first component, it will restart.
%
% Use as
%    comp = hcp_componentanalysis(cfg, data)
% where
%   cfg.method  = string describing the method (default = 'fastica')
%   cfg.restart = how many times to restart the decomposition in case of convergence problems (default = 10)
%
% See FT_COMPONENTANALYSIS for the other options
%
% See also http://bugzilla.fcdonders.nl/show_bug.cgi?id=1519

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

% get the defaults
cfg.method  = ft_getopt(cfg, 'method', 'fastica');
cfg.restart = ft_getopt(cfg, 'restart', 10);

% separate the options for this function and the ft function
restart = cfg.restart;
tmpcfg  = rmfield(cfg, 'restart');

while restart>0
  try
    comp    = ft_componentanalysis(tmpcfg, data);
    restart = 0;
  catch
    disp(lasterr);
    restart = restart - 1;
    
    randn('state',restart*100);
    tmpcfg.fastica.initGuess = randn(numel(data.label));
    
    if restart>0
      warning('problem with convergence, restarting the decomposiiton');
    else
      error('problem with convergence, not restarting any more');
    end
  end % try
end % while
