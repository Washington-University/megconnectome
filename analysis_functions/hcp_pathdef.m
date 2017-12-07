function [pipelinedatadir, rawdatadir] = hcp_pathdef

% HCP_PATHDEF returns the directories where the raw and processed data
% can be found on this computer, together with a list of all raw datasets.
%
% Use as
%   [pipelinedatadir, rawdatadir] = hcp_pathdef
% which returns
%   pipelinedatadir = string
%   rawdatadir      = string
%   dataset         = cell-array with strings
%
% See also HCP_DB_QUERY, HCP_DB_GETFILENAME

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

% this is where the pipeline output is written and where a pipeline can find the output from the preceding steps
pipelinedatadir = hcp_select_datadir({
  '/home/mrphys/roboos/data/hcp/processedRO' % this should go before jansch
  '/home/language/jansch/public/hcp/processedJM'
  '/data/jansch/HCP/database/processedJM'
  '/data/roboos/pipeline'
  '/scratch/r.oostenveld/pipeline'
  '/mnt/hps/slurm/michalareasg/pipeline'
  % 'N:\pipeline_testing'
  % 'N:\PhaseII_testing'
  % 'N:\PhaseII_testing_x-nat'
  % 'N:\PhaseII_testing_newfit'
  % 'N:\PhaseII_testing_newfit\bump_problem_analysis'
  'N:\PhaseII_testing_newfit\new_win_filter'
  'W:\PhaseII_testing_newfit' %%%point to W:\PhaseII_testing\laura dir for elaborated data
  % 'N:\PhaseII_testing_MEG_EEG'
  
  pwd
  });

% specify where the raw data can be found
rawdatadir = hcp_select_datadir({
  '/data/HCP/webdav'
  '/data/jansch/HCP/database'
  '/data/roboos/webdav'
  '/home/mrphys/roboos/data/hcp/webdav' % this should go before jansch
  '/home/language/jansch/public/hcp/webdav'
  '/mnt/hps/home/michalareasg/data/hcp/database'
  '/scratch/r.oostenveld/webdav'
  'D:\webdav'
  % 'E:\HCP\webdav'
  'I:\webdav'
  % 'N:\Phase1MEG_fra\webdav'
  'N:\webdav'
  'W:\webdav\'
  });

