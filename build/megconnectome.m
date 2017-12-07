function megconnectome(varargin)

% MEGCONNECTOME is the entry function of the compiled "megconnectome.exe"
% application that includes fieldtrip and the megconnectome/analysis_functions.
% The compiled application can be used to execute the scripts that are found
% in megconnectome/analysis_scripts.
%
% This function can be started on the MATLAB command line as
%   megconnectome  scriptname.m
%   megconnectome  script1.m script2.m ...
%   megconnectome  jobfile.mat
% or after compilation on the Linux command line as
%   megconnectome.sh <MATLABROOT>
%   megconnectome.sh <MATLABROOT>  scriptname.m
%   megconnectome.sh <MATLABROOT>  script1.m script2.m ...
%   megconnectome.sh <MATLABROOT>  jobfile.mat
%
% It is possible to pass additional options on the MATLAB command line like
% this on the MATLAB command line
%   megconnectome --option value scriptname.m
% or on the Linux command line
%   megconnectome.sh <MATLABROOT> --option value scriptname.m
% The options and their values are automaticallly made available as local
% variables in the script execution environment.
%
% More information is available on
%   https://wiki.humanconnectome.org/display/FieldTrip/FieldTrip
%
% See also COMPILE_MEGCONNECTOME

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

% this function uses assignin/evalin. The alternative is to use assign/eval
% but then the script execution might collide with local variables inside this
% megconnectome function
workspace = 'caller';

% separate the --options from the filenames
[options, varargin] = hcp_getopt(varargin{:});

if any(strcmp(options(1:2:end), 'version'))
  % only show the provenance information, do not run anything
  varargin = {};
end

% get the general provenance information (e.g. user and hostname)
prov = hcp_provenance;

for i=1:length(varargin)
  fname = varargin{i};
  
  % get the specific provenance information for the script
  prov.script{i}.filename = fname;
  prov.script{i}.md5sum   = hcp_md5sum(fname);
  
  % The following section with "which" does not work in the compiled application
  %
  %   if isempty(which(fname))
  %     error('The file %s cannot be found\n',fname);
  %   end
  %   % ensure that the script name includes the full path and extension
  %   fname = which(fname);
  
  [p, f, ext] = fileparts(fname);
  
  switch ext
    case '.m'
      if ~exist(fname, 'file')
        error('The script %s cannot be found\n',fname);
      end
      
      fid = fopen(fname);
      if fid == -1, error('Cannot open %s',fname); end
      S = fscanf(fid,'%c');
      fclose(fid);
      
      prov.script{i}.content  = sprintf('\n%s\n', S);
      
      % capture all screen output in a diary file
      diary off
      diaryfile = tempname;
      diary(diaryfile);
      
      try
        % ensure that subsequent scripts do not interfere with each other
        % keep the hcp_default global variable
        global hcp_default
        localcopy = hcp_default;
        evalin(workspace,'clear all');
        global hcp_default
        hcp_default = localcopy;
        
        % make the options available as local variables
        for j=1:2:length(options)
          key = options{j};
          val = options{j+1};
          
          if ismember(val(1), {'[', '{'})
            % the value is a string that represents an array
            evalin(workspace, sprintf('%s = %s;', key, val));
          elseif ismember(val(1), {'1', '2', '3', '4', '5', '6', '7', '8', '9', '0', '-', '+'}) && all(ismember(cellstr(val(:)), {'1', '2', '3', '4', '5', '6', '7', '8', '9', '0', '-', '+', '.', 'e'}))
            % the value is a string that represents a number
            evalin(workspace, sprintf('%s = %s;', key, val));
          else
            assignin(workspace, key, val);
          end
        end
        
        % add some information about the script to the global hcp_default variable
        % this will be picked up by hcp_write_provenance
        hcp_default.prov.script = prov.script{i};
        
        % evaluate this script
        evalin(workspace,S);
        
        % remove the provenance information about the script
        hcp_default.prov.script = [];
        
      catch err
        % remove the provenance information about the script
        hcp_default.prov.script = [];
        
        fprintf('Execution failed: %s\n', fname);
        rethrow(err);
      end
      
    case '.mat'
      % The job input consists of a file like this
      %   /Users/robert/tmp/robert_mbp_p36978_b1_j001_input.mat
      % which after execution results in an output file
      %   /Users/robert/tmp/robert_mbp_p36978_b1_j001_output.mat
      
      if length(fname)<10
        error('The jobfile %s is incorrect\n',fname);
      elseif ~strcmp(fname(end-9:end), '_input.mat')
        error('The jobfile %s is incorrect\n',fname);
      end
      
      jobid = fname(1:end-10);
      [p, jobid] = fileparts(jobid);
      cd(p);
      
      % read the job input arguments from jobid_input.mat, delete the input file and execute the job
      qsubexec(jobid);
      
      if isempty(getenv('PBS_ENVIRONMENT'))
        % this suggests that the job is running outside of the torque environment, e.g. for local testing
        % execution of qsubexec is expected to result in a stdout and stderr log file from torque
        logout       = fullfile(p, sprintf('%s.o', jobid));
        logerr       = fullfile(p, sprintf('%s.e', jobid));
        fclose(fopen(logout, 'w'));
        fclose(fopen(logerr, 'w'));
      end
      
      % read the job results from jobid_output.mat, delete the output file
      % FIXME what to do with the output variable?
      argout = qsubget(jobid);
      
    otherwise
      error('Unknown input format for %s, should be *.m or *.mat\n',fname);
  end % switch type of input argument
  
  % force close all figures that were created
  close all
  
  diary off
  fid = fopen(diaryfile);
  if fid == -1, error('Cannot open %s',diaryfile); end
  S = fscanf(fid,'%c');
  fclose(fid);
  delete(diaryfile);
  
  % add the screen output to the provenance structure
  prov.script{i}.stdout = sprintf('\n%s\n', S);
  
end % for each of the input arguments

% display the provenance structure on screen as an XML stream
disp(fileexchange_struct2xml(struct('megconnectome', prov)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function assign(key, val)
assignin('caller', key, val);
