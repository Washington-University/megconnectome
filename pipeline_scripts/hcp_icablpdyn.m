%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setup the execution environment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opengl software;

% ensure that the time and date of execution are not stored in the provenance information
global ft_default
ft_default.trackcallinfo = 'no';

% allow the user to specify the path where additional data is present, e.g. the channel layout or anatomy files
if exist('path', 'var')
    addpath(path)
end

if ~exist('filename', 'var')
  %error('filename should be specified')
  tok = {};
else
  % this is old-style and unnecessary. it is a bit strange to having to define a raw data file name, while 
  % the data itself is not used, for the sole purpose of tokenizing the string.
  % the filename is assumed to be something like
  % 'rawdatadir/Phase1MEG/Subjects/CP10018/Experiments/CP10018_MEG/Scans/1-Rnoise_MNN_V1/Resources/4D/c,rfDC'
  tok = tokenize(filename, '/');
end

if ~exist('subjectid', 'var') && ~isempty(tok)
  subjectid = tok{end-7}; % hard-coded assumption
elseif ~exist('subjectid', 'var') && isempty(tok)
  error('the subjectid needs to be specified');
elseif isnumeric(subjectid)
  % convert to string
  subjectid = num2str(subjectid);
end

if ~exist('experimentid', 'var') && ~isempty(tok)
  experimentid = tok{end-5};
elseif ~exist('experimentid', 'var')
  experimentid = [subjectid,'_MEG'];
end

% scanid should be something like 3-Restin
if ~exist('scanid', 'var') && ~isempty(tok)
  scanid = tok{end-3}; % hard coded assumption
elseif ~exist('scanid', 'var') && isempty(tok) 
  error('the scanid needs to be specified');
end

% this is the directory where the results will be saved
if ~exist('pipelinedatadir', 'var')
  pipelinedatadir = hcp_pathdef;
end

% this is where the anatomical results can be found
if ~exist('datadir_anatomy', 'var')
  % make empty to keep Francesco's style of working operational
  datadir_anatomy = pipelinedatadir;
end

% this is where the source-level band-limited power can be found
if ~exist('datadir_analysis', 'var')
  % make empty to keep Francesco's style of working operational
  datadir_analysis = pipelinedatadir;
end

if ~exist('aband', 'var')
  aband =4;
end 

% the following defines the time step and the time window, default is 0.5 seconds and 10, respectively
if ~exist('timestep', 'var')
  timestep = 0.5; % in seconds
end
if ~exist('timewindow', 'var')
  timewindow = 10; % in seconds
end

% look whether a parcellation is requested
if ~exist('parcellationfile', 'var')
  parcellationfile = '';
end

% flag for keeping the individual time slice dconnfiles
if ~exist('keepdense', 'var')
  keepdense = true;
end

if ~exist('dofig', 'var')        
  dofig = 'no';
end

% print the matlab and megconnectome version to screen for provenance
ver('megconnectome')

% print the value of all local variables to screen for provenance
w = whos;
w = {w.name};
w = setdiff(w, {'w', 'ans'});
for i=1:length(w)
  fprintf(hcp_printstruct(w{i}, eval(w{i})));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% execute the pipeline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the variable aband indexes into the following cell-array
band_prefix = {
    'delta'
    'theta'
    'alpha'
    'betalow'
    'betahigh'
    'gammalow'
    'gammamid'
    'gammahigh'
    'whole'
    };

% convert to samples
Fs      = 50;    % in Hz, this should be information that can be recovered from the band-limited timecourses, right?
step    = round(Fs*timestep);
window  = round(Fs*timewindow); % in points;
for ib = aband
  nwin_tot = 0;
    
  % check whether the band-limited power envelope time courses exist
  hcp_check_pipelineoutput('icablpenv', 'subject', subjectid, 'experiment', fullfile(datadir_analysis,experimentid), 'scan', scanid, 'band', band_prefix{ib});
   
  % load the data
  outstr = sprintf('%s_%s_icablpenv_%s', experimentid, scanid, band_prefix{ib});
  disp(['loading blp file ' outstr])
  hcp_read_matlab(fullfile(datadir_analysis,outstr))
  ntp  = size(source_blp.power,2);
  nwin = fix((ntp-window)/step);
    
  % compute a running sum
  for k=1:nwin
    vect1        = [(k-1)*step+1:(k-1)*step+window];
    time_corr(k) = (1/Fs)*mean(vect1); % in seconds
     
    tmp          = find(vect1 >= 1 & vect1 <= ntp);
    connect_stat = corr(source_blp.power(:,vect1(tmp))');
     
    connect          = [];
    connect.pos      = source_blp.pos;
    %connect.time     = time_corr(k);
    connect.blpcorr  = connect_stat;
    connect.blpcorrdimord   = 'pos_pos';%_time';
    if isfield(source_blp, 'tri'),                 connect.tri                 = source_blp.tri;                 end
    if isfield(source_blp, 'brainstructure'),      connect.brainstructure      = source_blp.brainstructure;      end
    if isfield(source_blp, 'brainstructurelabel'), connect.brainstructurelabel = source_blp.brainstructurelabel; end
       
    % save it as a cifti
    outputfile = fullfile(pipelinedatadir,[experimentid '_' scanid '_icablpdyn_' band_prefix{ib} '_windowlength' num2str(timewindow,'%3.1f') 's_timepoint' num2str(time_corr(k),'%05.1f') 's']);
       
    %hcp_write_cifti(outputfile, connect, 'parameter', 'blpcorr', 'type', 'dconnseries');
    hcp_write_cifti(outputfile, connect, 'parameter', 'blpcorr', 'type', 'dconn', 'writesurface', false);
     
    if ~isempty(parcellationfile)
      % parcellate it with wb_command, using system call, using the specified parcellation
      inputfile = [outputfile, '.dconn.nii'];
      tmpfile   = [inputfile,'.temp'];
      pconnfile = strrep(inputfile, 'dconn', 'pconn');
      systemcall1 = ['wb_command -cifti-parcellate ',inputfile,' ',parcellationfile,' ROW ',tmpfile];
      systemcall2 = ['wb_command -cifti-parcellate ',tmpfile,' ',parcellationfile,' COLUMN ',pconnfile];
      system(systemcall1);
      system(systemcall2);
      delete(tmpfile);

      if ~keepdense
        delete(inputfile);
        delete([inputfile,'.xml']);
      end
    end
 
    if strcmp(dofig,'yes')
      % check for existence and load in the cortical sheet
      hcp_check_pipelineoutput('anatomy', 'subject', fullfile(datadir_anatomy,subjectid));
      hcp_read_matlab(fullfile(datadir_anatomy,[subjectid '_MEG_anatomy_sourcemodel_2d']));
      
      imgname    = outputfile;
      options_pl = {'outputfile',imgname,'sorting','yes','parcel_type','Yeo',...
                    'mask','yes','mask_edist','yes','edist_radius',r_dist,'color_extr',[0 1]};
      hcp_icaplotconnectome_yeo(connect_stat, sourcemodel2d, options_pl)
            
      % extract the index of a vertex that belongs to a node in the DMN (PCC)
      net_seeds = hcp_mcw_netdef('DMN', []);
      z_s       = connect_stat(net_seeds.cortex_index(1),:);
            
      imgname    = [outputfile '_view_' net_seeds.label{1} '.png'];
      color_extr = [0 1];
      options_pl = {'outputfile',imgname, 'mask','yes','color_extr',color_extr,'color_map','jet'};
      hcp_icaplotcortex(z_s, subjectid, options_pl)
    end

    if ~keepdense
      % remove the dconn file
    end    

  end
  clear source_blp connect_stat;
    
end % for each frequency band
