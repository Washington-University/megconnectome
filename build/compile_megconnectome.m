function compile_megconnectome(fieldtriproot, hcproot)

% COMPILE_MEGCONNECTOME compiles the hcp and fieldtrip functions along with
% the "megconnectome" entry function into a stand-alone compiled executable.
%
% The compiled executable includes
%  - all main fieldtrip m-files
%  - all main fieldtrip's m-files dependencies for as long as these
%    dependencies are in the fieldtrip modules on the path, matlab
%    built-in, or toolbox/(stats/images/signal) functions
%  - the megconnectome/analysis_functions and its subdirectories
%
% See also MEGCONNECTOME

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

% clear all variables, globals, functions and MEX links
clear all

fname = 'megconnectome';

v = ver('MATLAB');
if ~strcmp(v.Version, '8.0')
    error('the megconnectome application should be compiled with MATLAB 2012b (8.0)');
end

if nargin<1 || isempty(fieldtriproot)
    fieldtriproot = fileparts(which('ft_defaults'));
end

if nargin<2 || isempty(hcproot)
    % this script is in hcproot/build
    hcproot = fileparts(which(mfilename));
    hcproot = fileparts(hcproot);
end

origdir  = pwd;
builddir = fullfile(hcproot, 'build');
bindir   = fullfile(hcproot, 'bin');

cd(builddir);

% create a file with the timestamp of the compilation
fid = fopen('buildtimestamp.m', 'wt');
fprintf(fid, 'function s = buildtimestamp\n');
fprintf(fid, 's = ''%s'';\n', datestr(now));
fclose(fid);

cd(bindir);

fprintf('Using fieldtrip     from "%s"\n', fieldtriproot);
fprintf('Using megconnectome from "%s"\n', hcproot);

% clean the path
restoredefaultpath;

%--------------------------------
% FIELDTRIP RELATED PATH SETTINGS

% add the path to fieldtrip
addpath(fieldtriproot);

% ensure that the path to the default modules is specified
clear ft_defaults;
ft_defaults;

% do not use my personal defaults, but rather FieldTrip standard defaults
global ft_default
ft_default = [];

% do not use my personal defaults, but rather HCP standard defaults
global hcp_default
hcp_default = [];

% ensure that these special modules are also added
ft_hastoolbox('qsub', 1);
ft_hastoolbox('engine', 1);

% ensure that all external toolboxes are added to the path
% excluding spm2 (leading to more than one spm on the path -> confusion in FT)
exclude = {
    '.'
    '..'
    '.svn'
    'dipoli'
    'dmlt'
    'iso2mesh'
    'simbio'
    'spm2'
    'sqdproject'
    'yokogawa'
    };
extd = dir([fieldtriproot,'/external']);
extd = setdiff({extd([extd.isdir]).name}, exclude);
for k = 1:numel(extd)
    addpath(fullfile(fieldtriproot,'external',extd{k}));
end

%--------------------------
% HCP RELATED PATH SETTINGS

addpath(hcproot);
addpath(fullfile(hcproot, 'external'));
addpath(fullfile(hcproot, 'analysis_functions'));
addpath(fullfile(hcproot, 'trial_functions'));

%-------------------
% DO THE COMPILATION

cmd = ['mcc -R -singleCompThread -N -o ' fname ' -m ' fname ...
    ' -a ' fieldtriproot '/*.m' ...
    ' -a ' fieldtriproot '/utilities/*.m' ...
    ' -a ' fieldtriproot '/fileio/*.m' ...
    ' -a ' fieldtriproot '/forward/*.m' ...
    ' -a ' fieldtriproot '/inverse/*.m' ...
    ' -a ' fieldtriproot '/plotting/*.m' ...
    ' -a ' fieldtriproot '/statfun/*.m' ...
    ' -a ' fieldtriproot '/trialfun/*.m' ...
    ' -a ' fieldtriproot '/preproc/*.m' ...
    ' -a ' fieldtriproot '/qsub/*.m' ...
    ' -a ' fieldtriproot '/engine/*.m' ...
    ' -a ' fieldtriproot '/external/plot2svg/*.m' ...
    ' -a ' hcproot       '/Contents.m' ...
    ' -a ' hcproot       '/build/buildtimestamp.m' ...
    ' -a ' hcproot       '/build/megconnectome.m' ...
    ' -a ' hcproot       '/external/*.m' ...
    ' -a ' hcproot       '/external/*.mex*' ...
    ' -a ' hcproot       '/analysis_functions/*.m' ...
    ' -a ' hcproot       '/trial_functions/*.m' ...
    ' -p ' matlabroot    '/toolbox/signal' ...
    ' -p ' matlabroot    '/toolbox/images' ...
    ' -p ' matlabroot    '/toolbox/stats' ...
    ' -p ' matlabroot    '/toolbox/optim' ...
    ' -p ' matlabroot    '/toolbox/curvefit' ...
    ];
eval(cmd);

% somehow I don't manage to get this going with more than one directory to be added when calling mcc directly:
% mcc('-N', '-a', '/home/common/matlab/fieldtrip/*.m', '-o', fname, '-m', fname, '-p', [matlabroot,'/toolbox/signal:', matlabroot,'/toolbox/images:', matlabroot,'/toolbox/stats']);

% remove the additional files that were created during compilation
delete mccExcludedFiles.log
delete readme.txt
delete run_megconnectome.sh

fprintf('Finished compilation\n');
cd(origdir);
