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
    error('filename should be specified')
end

% the filename is assumed to be something like
% 'rawdatadir/Phase1MEG/Subjects/CP10018/Experiments/CP10018_MEG/Scans/1-Rnoise_MNN_V1/Resources/4D/c,rfDC'
tok = tokenize(filename, '/');

if ~exist('subjectid', 'var')
    subjectid = tok{end-7};
end

if ~exist('experimentid', 'var')
    experimentid = tok{end-5};
end

scanid={'3-Restin' ; '4-Restin' ; '5-Restin'};

if ~exist('pipelinedatadir', 'var')
    pipelinedatadir = hcp_pathdef;
end

if ~exist('aband', 'var')
    aband=[1 2 3 4 5 6 7 8 9];
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

% change to the location of the processed data (input and output)
cd(pipelinedatadir)

hcp_check_pipelineoutput('anatomy', 'subject', subjectid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% execute the pipeline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hcp_read_matlab([subjectid '_MEG_anatomy_sourcemodel_2d']);
sourcemodelsubj = sourcemodel2d;
head=ft_read_headshape([subjectid '.L.midthickness.4k_fs_LR.surf.gii']);
head2=ft_read_headshape([subjectid '.R.midthickness.4k_fs_LR.surf.gii']);
sourcemodel=head;
sourcemodel.pnt=[head.pnt ; head2.pnt];
sourcemodel.tri=[head.tri ; head2.tri];


band_prefix={
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


window_corr=25; % window for stationary corr in sec
Fs=50; % BLP sampling frequency in Hz
window_corr2=window_corr*Fs;
source_blp=[];
r_dist=3.5; % radius of the mask related to the vertices euclidean distance
for ib=aband
    nwin_tot=0;
    for i_scan=1:3
        clear source_blp
        
        % check whether the band-limited power envelope time courses exist
        hcp_check_pipelineoutput('icablpenv', 'subject', subjectid, 'experiment', experimentid, 'scan', scanid{i_scan}, 'band', band_prefix{ib});
        
        % load the data
        outstr = sprintf('%s_%s_icablpenv_%s', experimentid, scanid{i_scan}, band_prefix{ib});
        disp(['loading blp file ' outstr])
        hcp_read_matlab(outstr)
        
        if (i_scan==1)
            % allocate memory for the connectivity matrix in the first iteration for the first scan
            connect_stat=zeros(size(source_blp.power,1));
        end
        ntp=size(source_blp.power,2);
        nwin=floor(ntp/(Fs*window_corr));
        
        % compute a running sum
        for i=1:nwin
            vect=[((i-1)*window_corr2)+1:window_corr2*i];
            connect_stat=corr(source_blp.power(:,vect)')+connect_stat;
        end % for each window
        nwin_tot=nwin_tot+nwin;
    end % for each resting-state scan
    
    % normalise with the number of windows and cast to single precision
    connect_stat = connect_stat/nwin_tot;
    
    % create the output structure
    connect.pos      = source_blp.pos;
    connect.dimord   = 'pos_pos';
    connect.blpcorr  = connect_stat;
    connect.blp_band = source_blp.blp_band;
    if isfield(source_blp, 'tri'), connect.tri = source_blp.tri; end
    if isfield(source_blp, 'brainstructure'), connect.brainstructure = source_blp.brainstructure; end
    if isfield(source_blp, 'brainstructurelabel'), connect.brainstructurelabel = source_blp.brainstructurelabel; end
    
    % save it as a cifti
    outputfile   = [experimentid '_Restin_icablpcorr_' band_prefix{ib} '.blpcorr'];
    hcp_write_cifti(outputfile, connect, 'parameter', 'blpcorr', 'type', 'dconn');
    
    % also save it as a matlab file
    hcp_write_matlab(outputfile,'connect');
    
    dofig='yes';
    if strcmp(dofig,'yes')
        
        imgname = outputfile;
        options_pl={'outputfile',imgname, 'sorting', 'yes','parcel_type','RSN','mask','yes','mask_edist','yes','edist_radius',r_dist,'color_extr',[0 1]};
        hcp_icaplotconnectome(connect_stat, sourcemodel2d, options_pl)
        
        % extract the index of a vertex that belongs to a node in the DMN (PCC)
        net_seeds=hcp_mcw_netdef('DMN', []);
        z_s=connect.blpcorr(net_seeds.cortex_index(1),:);
        
        imgname = [outputfile '_view_' net_seeds.label{1} '.png'];
        color_extr= [0 1];
        options_pl={'outputfile',imgname, 'mask','yes','color_extr',color_extr,'color_map','jet'};
        
        hcp_icaplotcortex(z_s, subjectid, options_pl)
        
    end
    clear source_blp connect connect_stat;
    
end % for each frequency band
