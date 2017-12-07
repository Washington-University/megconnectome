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
ft_default.checksize     = inf; % we need to retain the shapemri, will be otherwise cleared
ft_default.trackconfig   = 'silent';

% allow the user to specify the path where additional data is present, e.g. the channel layout or anatomy files
if exist('path', 'var')
    addpath(path)
end

if ~exist('subjectid', 'var')
    error('subjectid should be specified')
elseif isnumeric(subjectid)
    subjectid = num2str(subjectid);
end

fprintf('executing the anatomy pipeline for subject %s\n', subjectid);

if ~exist('outputdir', 'var')
    outputdir = pwd;
    fprintf('using %s as the directory to save the results\n', outputdir);
end

if ~exist('structuralpreprocdir', 'var')
    % we cannot use the strucural preprocessing results, revert to old-style pipeline
    fprintf('not using the high quality structural preprocessing results\n');
else
    fprintf('using the structural preprocessing results from %s\n', structuralpreprocdir);
    hrmrifile      = fullfile(structuralpreprocdir, 'T1w', 'T1w_acpc_dc_restore.nii.gz');
    inputsurffile  = fullfile(structuralpreprocdir, 'T1w', 'fsaverage_LR32k', [subjectid,'.L.midthickness.32k_fs_LR.surf.gii']);
    inputsphere    = fullfile(structuralpreprocdir, 'MNINonLinear', 'fsaverage_LR32k', [subjectid,'.L.sphere.32k_fs_LR.surf.gii']);
    outputsphere   = fullfile(outputdir, 'Sphere.4k.L.surf.gii');
    outputsurffile = fullfile(outputdir, [subjectid,'.L.midthickness.4k_fs_LR.surf.gii']);
end

% the following flags pertain to the three main parts of the pipeline
if ~exist('dopipeinteractive', 'var'), dopipeinteractive = 0; end
if ~exist('dopipeautomatic', 'var'),   dopipeautomatic   = 0; end
if ~exist('doqualitycheck', 'var'),    doqualitycheck    = 0; end

if dopipeinteractive,
    % set flags to facilitate debugging or running parts of the pipeline
    if ~exist('dofiducials',      'var'),                dofiducials      = 1; end
    if ~exist('dolandmarks',      'var'),                dolandmarks      = 1; end
    if ~exist('docoregistration', 'var'),                docoregistration = 1; end
    if docoregistration && ~exist('docoreg_spm', 'var'), docoreg_spm      = 1; end
    if docoregistration && ~exist('docoreg_bti', 'var'), docoreg_bti      = 1; end
    
    % specify some default parameters
    if docoreg_bti && ~exist('docoreg_bti_icp', 'var'),         docoreg_bti_icp         = 1; end
    if docoreg_bti && ~exist('docoreg_bti_interactive', 'var'), docoreg_bti_interactive = 1; end
    
    % perform some checks on conditional presence of some variables
    if (dofiducials || dolandmarks) && ~exist('dicomfile', 'var')
        % we need a pointer to a file in the dicom-series that contains the T1-weighted anatomical image
        error('for the interactive part of the anatomy pipeline, a pointer to a file from the dicom series is needed');
    end
end

if dopipeautomatic,
    % set flags to facilitate debugging or running parts of the pipeline
    if ~exist('doheadmodel', 'var'),        doheadmodel     = 1; end
    if ~exist('dosourcemodel3d', 'var'),    dosourcemodel3d = 1; end
    if ~exist('dosourcemodel2d', 'var'),    dosourcemodel2d = 1; end
    if ~exist('dofreesurfer', 'var'),       dofreesurfer    = 0; end % this functionality can probably be removed
    if ~exist('domnesuite', 'var'),         domnesuite      = 0; end % this functionality can probably be removed
    
    % specify some default parameters
    if dosourcemodel3d && ~exist('gridresolution', 'var'),  gridresolution  = [4 6 8]; end
    if doheadmodel     && ~exist('headmodelthr', 'var'),    headmodelthr    = 0.5;     end
    if doheadmodel     && ~exist('headmodelsmooth', 'var'), headmodelsmooth = 5;       end
    
    % perform some checks on conditional presence of some variables
    if dosourcemodel2d && ~exist('structuralpreprocdir', 'var') && ~exist('mnepath', 'var')
        error('when computing the cortical sheet based source model the path to MNE-suite software needs to be specified as ''mnepath''');
    end
    if dosourcemodel2d && ~exist('structuralpreprocdir', 'var') && ~exist('surffile', 'var')
        error('when computing the cortical sheet based source model a file pointer to the high resolution cortical sheets needs to be provided as ''surffile''');
    end
    if dosourcemodel2d && ~exist('structuralpreprocdir', 'var') && ~exist('surflabeldir', 'var')
        error('when computing the cortical sheet based source model a pointer to the directory where the labeling of the surfaces is stored needs to be provided as ''surflabeldir''');
    end
    
    if docoregistration && ~exist('hsfile', 'var')
        % we need a pointer to an hs-file
        error('for the automatic coregistration a pointer to an hsfile is needed');
    end
end

%--------------------
% declare outputfiles
outputprefix       = fullfile(outputdir, [subjectid,'_MEG_anatomy']);
nifti_anatomical   = [outputprefix,'_anatomical.nii'];
textfile_fiducials = [outputprefix,'_fiducials.txt'];
textfile_landmarks = [outputprefix,'_landmarks.txt'];
textfile_transform = [outputprefix,'_transform.txt'];
%--------------------------


% print the value of all local variables to screen for provenance
w = whos;
w = {w.name};
w = setdiff(w, {'w', 'ans'});
for i=1:length(w)
    fprintf(hcp_printstruct(w{i}, eval(w{i})));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HERE THE COMPUTATIONS START

%% The following is the interactive part of the pipeline
if dopipeinteractive,
    fprintf('\n');
    fprintf('-------------------------------------------------------------------------\n');
    fprintf('Running the interactive part of the anatomy pipeline\n');
    fprintf('\n');
    
    %---------------------------------------------
    % execute the interactive coregistration steps
    if dolandmarks,
        mriorig = ft_read_mri(dicomfile);
        mri     = mriorig;
        
        % coregister to anterior commissure based RAS space
        fprintf('Please identify the Anterior Commissure, Posterior Commissure, a point on the positive Z and X axes, and a point on the right part of the head\n');
        cfg             = [];
        cfg.interactive = 'yes';
        mri             = ft_volumerealign(cfg, mriorig);
        landmarks       = mri.cfg.landmark;
        landmarks.coordsys = 'vox';
        hcp_write_ascii(textfile_landmarks, 'landmarks');
    end
    
    if dofiducials,
        mriorig = ft_read_mri(dicomfile);
        mri     = mriorig;
        
        % coregister to subject headspace
        fprintf('Please identify the LPA, RPA, nasion, and a point on the positive Z-axis\n');
        cfg             = [];
        cfg.interactive = 'yes';
        mri             = ft_volumerealign(cfg, mriorig);
        fiducials       = mri.cfg.fiducial;
        fiducials.coordsys = 'vox';
        hcp_write_ascii(textfile_fiducials, 'fiducials');
        
        % save the pre-registered anatomical to a nifti file
        cfg           = [];
        cfg.parameter = 'anatomy';
        cfg.filename  = nifti_anatomical;
        cfg.filetype  = 'nifti';
        ft_volumewrite(cfg, mri);
    end
    
    %--------------------------------
    % coregistration
    if docoregistration
        
        % check that at least the nifti with the anatomy exists
        if exist(nifti_anatomical, 'file')
            mri = ft_read_mri(nifti_anatomical);
        else
            error('for the automatic part of the anatomy pipeline a nifti file containing the anatomy is needed: run the interactive part of the pipeline first and/or ensure that the nifti file created in this step is in the output directory');
        end
        mriorig = mri;
        % mriorig.transform = eye(4); % FIXME does this need to be done?
        % what when voxels are not isotropic
        
        if exist(textfile_fiducials, 'file')
            hcp_read_ascii(textfile_fiducials);
        else
            error('for the automatic coregistration a textfile with fiducial locations needs to exist');
        end
        if exist(textfile_landmarks, 'file')
            hcp_read_ascii(textfile_landmarks);
        else
            error('for the automatic coregistration a textfile with landmark locations needs to exist');
        end
        
        if exist(textfile_transform, 'file')
            hcp_read_ascii(textfile_transform);
        end
        
        if docoreg_spm
          if exist('transform', 'var')
            % if the transform already exists
            % scrub the fields that have 'spm' in them
            fn = fieldnames(transform);
            removefields = ~cellfun('isempty', strfind(fn, 'spm'));
            transform    = rmfield(transform, fn(removefields));
          end
          
          % do the coregistration to MNI space
          fprintf('\n');
          fprintf('Coregistering the anatomy to the axes of the MNI coordinate system\n');
          fprintf('\n');
          
          cfg          = [];
          cfg.landmark = landmarks;
          mri          = ft_volumerealign(cfg, mriorig);
          
          if exist('hrmrifile', 'var') && exist(hrmrifile, 'file')
            % coregister the low-resolution MRI to the high resolution in
            % ACPC-space; landmarks are not needed.
            fprintf('\n');
            fprintf('Coregistering the low-resolution anatomy to the high-resolution anatomy\n');
            fprintf('\n');
            
            % do the coregistration to MNI space
            targetmri = ft_read_mri(hrmrifile);
            targetmri.coordsys = 'spm';
            
            cfg             = [];
            cfg.method      = 'spm';
            cfg.spm.regtype = 'rigid';
            mri             = ft_volumerealign(cfg, mri, targetmri);
            
            % here make a control figure to check the coregistration
            % reconstruct the scalp surface from the mris, and overlay them
            tmpcfg = [];
            tmpcfg.output = 'scalp';
            seg           = ft_volumesegment(tmpcfg, mri);
            seg_target    = ft_volumesegment(tmpcfg, targetmri);
            
            tmpcfg = [];
            tmpcfg             = [];
            tmpcfg.tissue      = 'scalp';
            tmpcfg.method      = 'projectmesh';
            tmpcfg.numvertices = 20000;
            scalp              = ft_prepare_mesh(tmpcfg, seg);
            scalp_target       = ft_prepare_mesh(tmpcfg, seg_target);
            
            figure;
            subplot('position',[0.01 0.51 0.48 0.48]);hold on;
            ft_plot_mesh(scalp,       'edgecolor','none','facecolor','r','facealpha',0.3);
            ft_plot_mesh(scalp_target,'edgecolor','none','facecolor','b','facealpha',0.3); view(180,-90);
            plot3([-130 130],[0 0],[0 0],'k');plot3([0 0],[-120 120],[0 0],'k');plot3([0 0],[0 0],[-100 150],'k');
            subplot('position',[0.51 0.51 0.48 0.48]);hold on;
            ft_plot_mesh(scalp,       'edgecolor','none','facecolor','r','facealpha',0.3);
            ft_plot_mesh(scalp_target,'edgecolor','none','facecolor','b','facealpha',0.3); view(0,90);
            plot3([-130 130],[0 0],[0 0],'k');plot3([0 0],[-120 120],[0 0],'k');plot3([0 0],[0 0],[-100 150],'k');
            subplot('position',[0.01 0.01 0.48 0.48]);hold on;
            ft_plot_mesh(scalp,       'edgecolor','none','facecolor','r','facealpha',0.3);
            ft_plot_mesh(scalp_target,'edgecolor','none','facecolor','b','facealpha',0.3); view(90,0);
            plot3([-130 130],[0 0],[0 0],'k');plot3([0 0],[-120 120],[0 0],'k');plot3([0 0],[0 0],[-100 150],'k');
            subplot('position',[0.51 0.01 0.48 0.48]);hold on;
            ft_plot_mesh(scalp,       'edgecolor','none','facecolor','r','facealpha',0.3);
            ft_plot_mesh(scalp_target,'edgecolor','none','facecolor','b','facealpha',0.3); view(0,0);
            plot3([-130 130],[0 0],[0 0],'k');plot3([0 0],[-120 120],[0 0],'k');plot3([0 0],[0 0],[-100 150],'k');
            axis on;
            grid on;
            set(gcf,'color','w')
            hcp_write_figure([outputprefix,'_coregistration_lowres2hires.png'], gcf, 'resolution', 300); close all;
            
            transform.vox2spm_interactive = mri.transformorig;
            transform.vox2spm_registered  = mri.transform;
          end
          
          % set the transformation matrix
          transform.vox2spm = mri.transform;
          transform.spm2vox = inv(transform.vox2spm);
        end
        
        if docoreg_bti
          if exist('transform', 'var'),          
            % scrub the fields that have 'bti' in them
            fn = fieldnames(transform);
            removefields = ~cellfun('isempty', strfind(fn, 'bti'));
            transform    = rmfield(transform, fn(removefields));
          end

          fprintf('\n');
          fprintf('-------------------------------------------------------------------------\n');
          fprintf('\n');
          fprintf('Coregistering the anatomy to the axes of the MEG coordinate system\n');
          fprintf('\n');
          
          % do an initial coregistration to BTI space
          cfg          = [];
          cfg.fiducial = fiducials;
          mri          = ft_volumerealign(cfg, mriorig);
          
          transform.vox2bti_interactive = mri.transform;
          
          % refine the coregistration by doing a icp-based coregistration using
          % the hs_file and the scalp surface reconstructed from the 1mm anatomy
          fprintf('\n');
          fprintf('Refining the coregistration using the headshape file\n');
          fprintf('\n');
          
          cfg           = [];
          cfg.method    = 'headshape';
          cfg.headshape.headshape = ft_read_headshape(hsfile);
          cfg.headshape.icp       = docoreg_bti_icp;
          cfg.headshape.interactive = docoreg_bti_interactive;
          % weight the points below the xy-plane and on the forehead more
          cfg.weights   = ones(size(cfg.headshape.headshape.pnt,1),1);
          %cfg.weights(cfg.headshape.headshape.pnt(:,3)<0)    = 1.5;
          %cfg.weights(cfg.headshape.headshape.pnt(:,3)>0.08 & cfg.headshape.headshape.pnt(:,3)<0.1) = 1.5;
          %cfg.weights(cfg.headshape.headshape.pnt(:,1)>0.05)  = 2;
          %cfg.weights(cfg.headshape.headshape.pnt(:,1)<-0.05) = 2;
          cfg.headshape.scalpsmooth = 1;%'no';
          cfg.headshape.scalpthreshold = 0.08;
          mri           = ft_volumerealign(cfg, mri);
          
          if isfield(mri.cfg.headshape, 'headshape')
            headshape     = struct(mri.cfg.headshape.headshape); % convert back from config object
            headshapemri  = struct(mri.cfg.headshape.headshapemri);
          else
            % backward compatibility
            headshape     = struct(mri.cfg.headshape); % convert back from config object
            headshapemri  = struct(mri.cfg.headshapemri);
          end
          headshape.coordsys    = 'bti';
          headshapemri.coordsys = 'bti';
            
          mrifid.pnt   = ft_warp_apply(mri.transform, [fiducials.nas;fiducials.lpa;fiducials.rpa]);
          mrifid.label = {'NZinteractive';'Linteractive';'Rinteractive'};
          headshapemri.fid = mrifid;
          
          % write the headshapes
          hcp_write_matlab([outputprefix,'_headshape'],     'headshape');
          hcp_write_matlab([outputprefix,'_headshapemri'],  'headshapemri');
          
          % quality check for the coregistration between headshape and mri
          headshape    = hcp_ensure_units(headshape,    'mm');
          headshapemri = hcp_ensure_units(headshapemri, 'mm');
          
          figure;
          subplot('position',[0.01 0.51 0.48 0.48]);hold on;
          ft_plot_mesh(headshapemri,'edgecolor','none','facecolor',[0.7 0.7 0.7],'fidcolor','y','facealpha',0.3);
          ft_plot_mesh(headshapemri,'edgecolor','none','vertexcolor',headshapemri.distance); view(180,-90);
          plot3([-130 130],[0 0],[0 0],'k');plot3([0 0],[-120 120],[0 0],'k');plot3([0 0],[0 0],[-100 150],'k');
          subplot('position',[0.51 0.51 0.48 0.48]);hold on;
          ft_plot_mesh(headshapemri,'edgecolor','none','facecolor',[0.7 0.7 0.7],'fidcolor','y','facealpha',0.3);
          ft_plot_mesh(headshapemri,'edgecolor','none','vertexcolor',headshapemri.distance); view(0,90);
          plot3([-130 130],[0 0],[0 0],'k');plot3([0 0],[-120 120],[0 0],'k');plot3([0 0],[0 0],[-100 150],'k');
          subplot('position',[0.01 0.01 0.48 0.48]);hold on;
          ft_plot_mesh(headshapemri,'edgecolor','none','facecolor',[0.7 0.7 0.7],'fidcolor','y','facealpha',0.3);
          ft_plot_mesh(headshapemri,'edgecolor','none','vertexcolor',headshapemri.distance); view(90,0);
          plot3([-130 130],[0 0],[0 0],'k');plot3([0 0],[-120 120],[0 0],'k');plot3([0 0],[0 0],[-100 150],'k');
          subplot('position',[0.51 0.01 0.48 0.48]);hold on;
          ft_plot_mesh(headshapemri,'edgecolor','none','facecolor',[0.7 0.7 0.7],'fidcolor','y','facealpha',0.3);
          ft_plot_mesh(headshapemri,'edgecolor','none','vertexcolor',headshapemri.distance);  view(0,0); colorbar('east');
          plot3([-130 130],[0 0],[0 0],'k');plot3([0 0],[-120 120],[0 0],'k');plot3([0 0],[0 0],[-100 150],'k');
          axis on;
          grid on;
          set(gcf,'color','w')
          hcp_write_figure([outputprefix,'_headshape_distance.png'], gcf, 'resolution', 300); close all;
          
          v = headshapemri.pnt;
          f = headshapemri.tri;
          [f,v]=reducepatch(f,v, 0.2);
          headshapemri.pnt = v;
          headshapemri.tri = f;
          
          figure;
          subplot('position',[0.01 0.51 0.48 0.48]);hold on;
          ft_plot_mesh(headshapemri,'edgecolor','none','facecolor',[0.5 0.6 0.8],'fidcolor','y','facealpha',0.3);
          ft_plot_headshape(headshape,'vertexsize',3); view(180,-90);
          plot3([-130 130],[0 0],[0 0],'k');plot3([0 0],[-120 120],[0 0],'k');plot3([0 0],[0 0],[-100 150],'k');
          subplot('position',[0.51 0.51 0.48 0.48]);hold on;
          ft_plot_mesh(headshapemri,'edgecolor','none','facecolor',[0.5 0.6 0.8],'fidcolor','y','facealpha',0.3);
          ft_plot_headshape(headshape,'vertexsize',3); view(0,90);
          plot3([-130 130],[0 0],[0 0],'k');plot3([0 0],[-120 120],[0 0],'k');plot3([0 0],[0 0],[-100 150],'k');
          subplot('position',[0.01 0.01 0.48 0.48]);hold on;
          ft_plot_mesh(headshapemri,'edgecolor','none','facecolor',[0.5 0.6 0.8],'fidcolor','y','facealpha',0.3);
          ft_plot_headshape(headshape,'vertexsize',3); view(90,0);
          plot3([-130 130],[0 0],[0 0],'k');plot3([0 0],[-120 120],[0 0],'k');plot3([0 0],[0 0],[-100 150],'k');
          subplot('position',[0.51 0.01 0.48 0.48]);hold on;
          ft_plot_mesh(headshapemri,'edgecolor','none','facecolor',[0.5 0.6 0.8],'fidcolor','y','facealpha',0.3);
          ft_plot_headshape(headshape,'vertexsize',3); view(0,0);
          plot3([-130 130],[0 0],[0 0],'k');plot3([0 0],[-120 120],[0 0],'k');plot3([0 0],[0 0],[-100 150],'k');
          axis on;
          grid on;
          set(gcf,'color','w')
          hcp_write_figure([outputprefix,'_headshape.png'], gcf, 'resolution', 300); close all;
          
          % create figures at the landmarks' locations, for QC
          crosshair=@(pos)plot([-100 100],pos(2),'y',pos(1)*[1 1],[-100 100],'y');
          cfg = [];
          cfg.locationcoordinates = 'voxel';
          cfg.location    = landmarks.ac;
          cfg.interactive = 'no';
          figure;ft_sourceplot(cfg, mri);
          hcp_write_figure([outputprefix,'_landmarks_ac.png'], gcf, 'resolution', 500); close;
          cfg.location = landmarks.pc;
          figure;ft_sourceplot(cfg, mri);
          hcp_write_figure([outputprefix,'_landmarks_pc.png'], gcf, 'resolution', 500); close;
          cfg.location = landmarks.xzpoint;
          figure;ft_sourceplot(cfg, mri);
          hcp_write_figure([outputprefix,'_landmarks_xzpoint.png'], gcf, 'resolution', 500); close;
          cfg.location = landmarks.rpoint;
          figure;ft_sourceplot(cfg, mri);
          hcp_write_figure([outputprefix,'_landmarks_rpoint.png'], gcf, 'resolution', 500); close;
          
          % create figures at the fiducials' location, for QC
          cfg = [];
          cfg.locationcoordinates = 'voxel';
          cfg.location    = fiducials.lpa;
          cfg.interactive = 'no';
          figure;ft_sourceplot(cfg, mri);
          hcp_write_figure([outputprefix,'_fiducials_lpa.png'], gcf, 'resolution', 500); close;
          cfg.location = fiducials.rpa;
          figure;ft_sourceplot(cfg, mri);
          hcp_write_figure([outputprefix,'_fiducials_rpa.png'], gcf, 'resolution', 500); close;
          cfg.location = fiducials.nas;
          figure;ft_sourceplot(cfg, mri);
          hcp_write_figure([outputprefix,'_fiducials_nas.png'], gcf, 'resolution', 500); close;
          cfg.location = fiducials.zpoint;
          figure;ft_sourceplot(cfg, mri);
          hcp_write_figure([outputprefix,'_fiducials_zpoint.png'], gcf, 'resolution', 500); close;
          
          transform.vox2bti = mri.transform;
          transform.bti2vox = inv(transform.vox2bti);
        end
        
        % add additional transformations
        transform.spm2bti = transform.vox2bti/transform.vox2spm;
        transform.bti2spm = transform.vox2spm/transform.vox2bti;
        [p,f,e] = fileparts(nifti_anatomical);
        transform.vox_filename = [f,e];
        
        if exist('hrmrifile', 'var') && exist(hrmrifile, 'file')
          targetmri = ft_read_mri(hrmrifile);
          targetmri.coordsys = 'spm';  
                        
          % create some additional transformation matrices
          transform.vox07mm2spm = targetmri.transform;
          transform.vox07mm2bti = transform.spm2bti*transform.vox07mm2spm;
          
          transform.bti2vox07mm = inv(transform.vox07mm2bti);
          transform.spm2vox07mm = inv(transform.vox07mm2spm);
          
          [p,f,e] = fileparts(hrmrifile);
          transform.vox07mm_filename = [f,e];
        end
        
        try, transform = rmfield(transform, 'vox2spm_interactive'); end
        try, transform = rmfield(transform, 'vox2bti_interactive'); end
        try, transform = rmfield(transform, 'vox2spm_registered'); end
        
        % write the transform structure
        hcp_write_ascii(textfile_transform, 'transform', transform, '%% FIXME put some comment here');
        
        fprintf('\n');
        fprintf('-------------------------------------------------------------------------\n');
        fprintf('\n');
    end
    
    fprintf('Here ends the interactive part of the anatomy pipeline\n');
    fprintf('-------------------------------------------------------------------------\n');
    fprintf('\n');
end

%% The following is the non-interactive/automatic part of the pipeline
if dopipeautomatic,
    fprintf('\n');
    fprintf('-------------------------------------------------------------------------\n');
    fprintf('Running the non-interactive part of the anatomy pipeline\n');
    fprintf('\n');
    
    % check that at least the nifti with the anatomy exists
    if exist(nifti_anatomical, 'file')
        mri = ft_read_mri(nifti_anatomical);
    else
        error('for the automatic part of the anatomy pipeline a nifti file containing the anatomy is needed: run the interactive part of the pipeline first and/or ensure that the nifti file created in this step is in the output directory');
    end
    
    if exist(textfile_transform, 'file')
        hcp_read_ascii(textfile_transform);
    else
        error('when the automatic part of the anatomy pipeline is run without coregistration, a textfile with coregistration information needs to exist');
    end
    
    if isfield(transform, 'vox2bti'),
      % old style transform structure
      mri4D     = mri; mri4D.transform  = transform.vox2bti; mri4D.coordsys  = 'bti';
      mriRAS    = mri; mriRAS.transform = transform.vox2spm; mriRAS.coordsys = 'spm';
    else 
      mri4D     = mri; mri4D.transform  = transform.vox2bti; mri4D.coordsys  = 'bti';
      mriRAS    = mri; mriRAS.transform = transform.vox2spm; mriRAS.coordsys = 'spm';
    end
    
    %---------------------------------------------
    % execute the non-interactive part of the pipeline
    
    % specify some variables and options
    mnitemplate = fullfile(outputdir, 'T1.nii'); % external file assumed to be present in the outputdirectory
    tpm         = {fullfile(outputdir, 'grey.nii') fullfile(outputdir, 'white.nii') fullfile(outputdir, 'csf.nii')};
    
    % ensure the environment variables to be defined
    if dofreesurfer
        % obsolete
    end
    
    if dosourcemodel2d && domnesuite
        p = getenv('MNE_ROOT');
        if isempty(p)
            setenv('MNE_ROOT', mnepath);
        end
    end
    
    
    %------------------------------
    % create volume conductor model
    if doheadmodel
        % single shell for MEG
        fprintf('\n');
        fprintf('-------------------------------------------------------------------------\n');
        fprintf('Creating a single shell volume conduction model for MEG forward modelling\n');
        fprintf('\n');
        
        % perform the segmentation
        cfg          = [];
        cfg.tpm      = tpm;
        cfg.template = mnitemplate;
        cfg.output   = 'brain';
        cfg.brainthreshold = headmodelthr;
        cfg.brainsmooth    = headmodelsmooth;
        seg          = ft_volumesegment(cfg, mri4D);
        
        % create the boundary
        cfg = [];
        cfg.method = 'singleshell';
        cfg.numvertices = 5000;
        headmodel  = ft_prepare_headmodel(cfg, seg);
        headmodel.coordsys = 'bti';
        headmodel  = ft_convert_units(headmodel);
        
        % write the headmodel
        hcp_write_matlab([outputprefix,'_headmodel'], 'headmodel');
        
        % create figures for qualitycheck
        hcp_read_ascii(textfile_fiducials);
        hcp_read_ascii(textfile_transform);
        fiducials = hcp_ensure_units(fiducials, 'mm');
        fiducials = hcp_ensure_coordsys(fiducials, transform, 'bti');
        xyz       = [fiducials.nas;fiducials.lpa;fiducials.rpa];
        headmodel = hcp_ensure_units(headmodel, 'mm');
        
        figure;
        subplot(2,2,1); hold on; view(180,-90); ft_plot_vol(headmodel); plot3(xyz(:,1),xyz(:,2),xyz(:,3),'b*');
        subplot(2,2,2); hold on; view(0,90);    ft_plot_vol(headmodel); plot3(xyz(:,1),xyz(:,2),xyz(:,3),'b*');
        subplot(2,2,3); hold on; view(90,  0);  ft_plot_vol(headmodel); plot3(xyz(:,1),xyz(:,2),xyz(:,3),'b*');
        subplot(2,2,4); hold on; view(0, 0);    ft_plot_vol(headmodel); plot3(xyz(:,1),xyz(:,2),xyz(:,3),'b*');
        axis on; grid on;
        hcp_write_figure([outputprefix,'_headmodel.png'], 'format', 'png');
        hcp_write_figure([outputprefix,'_headmodel.fig'], 'format', 'fig');
        
        fprintf('\n');
        fprintf('-------------------------------------------------------------------------\n');
        fprintf('\n');
        
        % qualitycheck for the headmodel
        
    end
    
    %-----------------------
    % create source model(s)
    if dosourcemodel3d
        % 3D grid
        fprintf('\n');
        fprintf('-------------------------------------------------------------------------\n');
        fprintf('Creating 3-dimensional dipole grid for source modelling with beamformers\n');
        fprintf('\n');
        
        for k = 1:numel(gridresolution)
            fprintf('Grid resolution: %s\n\n',num2str(gridresolution(k)));
            
            templategridname = fullfile(outputdir,['standard_sourcemodel3d',num2str(gridresolution(k)),'mm']);
            
            % create the sourcemodel
            hcp_read_matlab(templategridname);
            cfg                 = [];
            cfg.grid.warpmni    = 'yes';
            cfg.grid.nonlinear  = 'yes';
            cfg.grid.templatemri = mnitemplate;
            cfg.grid.template   = sourcemodel;
            
            cfg.mri        = mri4D;
            sourcemodel3d    = ft_prepare_sourcemodel(cfg);
            
            % remove the mri-structure from grid.cfg
            sourcemodel3d.cfg = rmfield(sourcemodel3d.cfg, 'mri');
            
            % ensure correct units and coordinate system
            sourcemodel3d.coordsys = 'bti';
            sourcemodel3d          = ft_convert_units(sourcemodel3d);
            
            % write the sourcemodels
            hcp_write_matlab([outputprefix,'_sourcemodel_3d',num2str(gridresolution(k)),'mm'], 'sourcemodel3d');
            
            % create figure for qualitycheck
            hcp_read_ascii(textfile_fiducials);
            hcp_read_ascii(textfile_transform);
            fiducials     = hcp_ensure_units(fiducials, 'mm');
            fiducials     = hcp_ensure_coordsys(fiducials, transform, 'bti');
            xyz           = [fiducials.nas;fiducials.lpa;fiducials.rpa];
            sourcemodel3d = hcp_ensure_units(sourcemodel3d, 'mm');
            
            figure;
            subplot(2,2,1); hold on; view(180,-90); ft_plot_mesh(sourcemodel3d.pos(sourcemodel3d.inside,:)); plot3(xyz(:,1),xyz(:,2),xyz(:,3),'b*');
            subplot(2,2,2); hold on; view(0,90);    ft_plot_mesh(sourcemodel3d.pos(sourcemodel3d.inside,:)); plot3(xyz(:,1),xyz(:,2),xyz(:,3),'b*');
            subplot(2,2,3); hold on; view(90,0);    ft_plot_mesh(sourcemodel3d.pos(sourcemodel3d.inside,:)); plot3(xyz(:,1),xyz(:,2),xyz(:,3),'b*');
            subplot(2,2,4); hold on; view(0,0);     ft_plot_mesh(sourcemodel3d.pos(sourcemodel3d.inside,:)); plot3(xyz(:,1),xyz(:,2),xyz(:,3),'b*');
            axis on; grid on;
            hcp_write_figure([outputprefix,'_sourcemodel_3d',num2str(gridresolution(k)),'mm.png'], 'format', 'png');
            hcp_write_figure([outputprefix,'_sourcemodel_3d',num2str(gridresolution(k)),'mm.fig'], 'format', 'fig');
            
        end
        fprintf('\n');
        fprintf('-------------------------------------------------------------------------\n');
        fprintf('\n');
    end
    
    if dosourcemodel2d
        % 2D cortical sheet
        fprintf('\n');
        fprintf('-------------------------------------------------------------------------\n');
        fprintf('Creating dipole grid based on the cortical sheet\n');
        fprintf('\n');
        
        % if structural_preproc output exists, use this
        if exist('structuralpreprocdir', 'var')
            % use wb_command to resample the 32k_LR registered meshes to a lower
            % resolution
            fprintf('\n');
            fprintf('-------------------------------------------------------------------------\n');
            fprintf('Using workbench to downsample subject specific registered surfaces\n');
            fprintf('\n');
            
            str = ['wb_command -surface-resample ',inputsurffile,' ',inputsphere,' ',outputsphere,' BARYCENTRIC ',outputsurffile];
            fprintf('executing %s\n', str);
            system(str);
            
            str = strrep(str, '.L', '.R');
            fprintf('executing %s\n', str);
            system(str);
            
            % also do the inflated surfaces
            str = strrep(str, 'midthickness', 'inflated');
            fprintf('executing %s\n', str);
            system(str);
            
            str = strrep(str, '.R', '.L');
            fprintf('executing %s\n', str);
            system(str);
            
            % also do the sulc and curvature (just in case we need it)
            %       str = strrep(str, '.inflated', '.curv');
            %       str = strrep(str, '.surf',     '.shape');
            %       str = strrep(str, 'resample-surface', 'resample-metric');
            %       fprintf('executing %s\n', str);
            %       system(str);
            %
            %       str = strrep(str, '.L', '.R');
            %       fprintf('executing %s\n', str);
            %       system(str);
            %
            %       str = strrep(str, 'curv', 'sulc');
            %       fprintf('executing %s\n', str);
            %       system(str);
            %
            %       str = strrep(str, '.R', '.L');
            %       fprintf('executing %s\n', str);
            %       system(str);
            
            % also do the aparc files
            inputlabelfile = strrep(inputsurffile, 'T1w', 'MNINonLinear');
            inputlabelfile = strrep(inputlabelfile, '.surf', '.label');
            inputlabelfile = strrep(inputlabelfile, '.midthickness', '.aparc');
            
            outputlabelfile = strrep(outputsurffile, '.surf', '.label');
            outputlabelfile = strrep(outputlabelfile, '.midthickness', '.aparc');
            
            str = ['wb_command -label-resample ',inputlabelfile,' ',inputsphere,' ',outputsphere,' BARYCENTRIC ',outputlabelfile];
            fprintf('executing %s\n', str);
            system(str);
            
            str = strrep(str, '.L', '.R');
            fprintf('executing %s\n', str);
            system(str);
            
            str = strrep(str, '.aparc', '.aparc.a2009s');
            fprintf('executing %s\n', str);
            system(str);
            
            str = strrep(str, '.R', '.L');
            fprintf('executing %s\n', str);
            system(str);
            
            str = strrep(str, '.aparc.a2009s', '.BA');
            fprintf('executing %s\n', str);
            system(str);
            
            str = strrep(str, '.L', '.R');
            fprintf('executing %s\n', str);
            system(str);
            
            sourcemodel2d = ft_read_headshape({outputsurffile strrep(outputsurffile,'.L','.R')});
            if isfield(sourcemodel2d, 'hemisphere')
              sourcemodel2d.brainstructure = sourcemodel2d.hemisphere;
              sourcemodel2d.brainstructurelabel = {'CORTEX_LEFT' 'CORTEX_RIGHT'};
              sourcemodel2d = rmfield(sourcemodel2d, {'hemisphere', 'hemispherelabel'});
            end

            % convert to mm units and change the coordinate system
            sourcemodel2d = ft_convert_units(sourcemodel2d, 'mm');
            sourcemodel2d = ft_transform_geometry(transform.spm2bti, sourcemodel2d);
            
            % convert back to cm units
            sourcemodel2d = ft_convert_units(sourcemodel2d, 'cm');
            
            sourcemodel2d.pos      = double(sourcemodel2d.pnt);
            sourcemodel2d.tri      = double(sourcemodel2d.tri);
            sourcemodel2d          = rmfield(sourcemodel2d, 'pnt');
            sourcemodel2d.coordsys = 'bti';
            
            rootstr      = 'midthickness'; % hack for backward compatibility
            surflabeldir = outputdir;      % these two are needed to get the parcellation info
        else
            %       % we need to create the source space with MNE, using these mesh in
            %       % freesurfer format. These are in a different coordinate system, so we
            %       % need to adjust this.
            %
            %       % get the transformation matrix
            %       bnd1 = ft_read_headshape(fullfile(pathname_hrfsdir,'surf','lh.white.surf.gii'));
            %       bnd2 = ft_read_headshape(fullfile(pathname_hrfsdir,'surf','lh.white'));
            %
            %       pnt1 = double(bnd1.pnt(1:150:end,:));
            %       pnt2 = double(bnd2.pnt(1:150:end,:));
            %       T    = [pnt2 ones(size(pnt2,1),1)]\[pnt1 ones(size(pnt1,1),1)];
            %       T    = T'; %transforms from 2 to 1 (i.e. FS to native)
            %       % The previous snippet of code is only needed when not working with the gifties
            %       % The gifties are already appropriately registered
            
            % only execute the lengthy creation of the high resolution 2D-grid if it has not yet been done.
            if dofreesurfer
                % segment the mri
                cfg        = [];
                cfg.output = 'skullstrip';
                cfg.smooth = 5;
                cfg.tpm    = tpm;
                cfg.template = mnitemplate;
                seg        = ft_volumesegment(cfg, mriRAS);
                
                % the mri and seg variables need to be written to disk in order for
                % freesurfer to work out
                cfg = [];
                cfg.filetype  = 'nifti';
                cfg.parameter = 'anatomy';
                cfg.filename  = fullfile(pwd,subjectid);
                ft_volumewrite(cfg, mriRAS);
                cfg.filename  = fullfile(pwd,[subjectid 'skullstrip']);
                ft_volumewrite(cfg, seg);
                
                str     = which('freesurferscript.sh');
                [p,f,e] = fileparts(str);
                
                % run the first part of the freesurfer pipeline
                system([p,'/freesurferscript.sh ', pwd,' ', subjectid]);
            end
            
            if domnesuite
                % convert gifties into freesurfer format
                if ft_filetype(surffile, 'caret_surf')
                    tmpdir = fullfile(outputdir, subjectid);
                    mkdir([tmpdir filesep 'surf']);
                    
                    [p,f,e] = fileparts(surffile);
                    tok     = tokenize(f, '.');
                    rootstr = tok{3};
                    
                    tmpdir2 = fullfile(tmpdir, 'surf');
                    tmpbnd  = ft_read_headshape(surffile);
                    ft_write_headshape(fullfile(tmpdir2,['lh.',rootstr]), tmpbnd, 'format', 'freesurfer');
                    tmpbnd  = ft_read_headshape(strrep(surffile, '.L.', '.R.'));
                    ft_write_headshape(fullfile(tmpdir2,['rh.',rootstr]), tmpbnd, 'format', 'freesurfer');
                    tmpbnd  = ft_read_headshape(strrep(surffile, rootstr, 'sphere'));
                    ft_write_headshape(fullfile(tmpdir2,'lh.sphere'), tmpbnd, 'format', 'freesurfer');
                    tmpbnd  = ft_read_headshape(strrep(surffile, ['.L.',rootstr], ['.R.sphere']));
                    ft_write_headshape(fullfile(tmpdir2,'rh.sphere'), tmpbnd, 'format', 'freesurfer');
                    
                else
                    % get the files that are needed temporarily at a location that
                    % mnesuite can deal with it.
                    tmpdir = fullfile(outputdir, subjectid);
                    mkdir(tmpdir);
                    [p,f,e] = fileparts(surffile);
                    system(['cp ' p filesep '* ' tmpdir filesep 'surf']);
                    rootstr = 'orig';
                end
                
                str     = which('mnescript.sh');
                [p,f,e] = fileparts(str);
                system([p,'/mnescript.sh ', outputdir,' ', subjectid,' ', rootstr]);
            end
            
            % read in the source model
            fname         = fullfile(tmpdir,'bem',[subjectid,'-oct-6-src.fif']);
            sourcemodel2d = ft_read_headshape(fname, 'format', 'mne_source');
            
            % convert to mm units and change the coordinate system
            sourcemodel2d = ft_convert_units(sourcemodel2d, 'mm');
            sourcemodel2d = ft_transform_geometry(transform.spm2bti, sourcemodel2d);
            
            % convert back to cm units
            sourcemodel2d = ft_convert_units(sourcemodel2d, 'cm');
            
            sourcemodel2d.pos      = sourcemodel2d.pnt;
            sourcemodel2d          = rmfield(sourcemodel2d, 'pnt');
            sourcemodel2d.coordsys = 'bti';
            
            % hack to let the following work
            outputsurffile = surffile;
        end
        
        % write the sourcemodel
        hcp_write_matlab([outputprefix,'_sourcemodel_2d'], 'sourcemodel2d');
        
        % qualitycheck figures, for this we need the mri
        if exist('hrmrifile', 'var')
          mri           = ft_read_mri(hrmrifile);
          mri.coordsys  = 'bti';
          mri.transform = transform.vox07mm2bti;
        else
          mri           = ft_read_mri(nifti_anatomical);
          mri.coordsys  = 'bti';
          mri.transform = transform.vox2bti;
        end
        
        % also read in the headmodel: this will provide a crash if it does
        % not exist.
        hcp_read_matlab([outputprefix,'_headmodel']);
        
        % and the sourcemodel2d should be in mm
        sourcemodel2d = ft_convert_units(sourcemodel2d, 'mm');
        headmodel     = ft_convert_units(headmodel,     'mm');
        
        figure;
        options = {'transform',mri.transform,'intersectmesh',{sourcemodel2d headmodel.bnd}};
        subplot(2,2,1); hold on; hcp_plot_slice(mri.anatomy, 'location', [0  0 60], 'orientation', [0 0 1], options{:}); view(0,90);
        subplot(2,2,2); hold on; hcp_plot_slice(mri.anatomy, 'location', [0  0 20], 'orientation', [0 0 1], options{:}); view(0,90);
        subplot(2,2,3); hold on; hcp_plot_slice(mri.anatomy, 'location', [0 20  0], 'orientation', [1 0 0], options{:}); view(90,0);
        subplot(2,2,4); hold on; hcp_plot_slice(mri.anatomy, 'location', [0 20  0], 'orientation', [0 1 0], options{:}); view(0,0);
        set(gcf, 'Renderer', 'zbuffer')
        hcp_write_figure([outputprefix,'_sourcemodel_2d.png']);

        options = {'transform', mri.transform,'nslice',16,'intersectmesh',{sourcemodel2d headmodel.bnd},...
           'intersectlinewidth',1,'slicesize',[300 300]};

        figure;
        hcp_plot_montage(mri.anatomy, 'location', [0 0 0], 'orientation', [0 0 1], 'slicerange', [-20 120], options{:});
        set(gcf, 'Renderer', 'zbuffer');   
        hcp_write_figure([outputprefix,'_slice1.png'], 'resolution', 300);

        figure;
        hcp_plot_montage(mri.anatomy, 'location', [0 0 0], 'orientation', [0 1 0], 'slicerange', [-60 60], options{:});
        set(gcf, 'Renderer', 'zbuffer');
        hcp_write_figure([outputprefix,'_slice2.png'], 'resolution', 300);

        figure;
        hcp_plot_montage(mri.anatomy, 'location', [0 0 0], 'orientation', [1 0 0], 'slicerange', [-70 110], options{:});
        set(gcf, 'Renderer', 'zbuffer');
        hcp_write_figure([outputprefix,'_slice3.png'], 'resolution', 300);
        
        fprintf('\n');
        fprintf('-------------------------------------------------------------------------\n');
        fprintf('\n');
    end
    
    fprintf('Here ends the non-interactive part of the anatomy pipeline\n');
    fprintf('-------------------------------------------------------------------------\n');
    fprintf('\n');
end % dopipeautomatic

%% The following is checking the quality of the output
if doqualitycheck,
    %     hcp_read_ascii(textfile_fiducials); % fiducials
    %     hcp_read_ascii(textfile_landmarks); % landmarks
    %     hcp_read_ascii(textfile_transform); % transform
    %     hcp_read_matlab([outputprefix,'_headshape']);
    %     hcp_read_matlab([outputprefix,'_headshapemri']);
    %     hcp_read_matlab([outputprefix,'_headmodel']);
    %     hcp_read_matlab([outputprefix,'_sourcemodel2d']);
    %     hcp_read_matlab([outputprefix,'_sourcemodel3d8mm']);
    %     mri           = ft_read_mri(nifti_anatomical);
    %     mri.coordsys  = 'bti';
    %     mri.transform = transform.vox2bti;
    %     hcp_anatomy_qualitycheck(outputprefix, fiducials, landmarks, transform, mri, headshape, headshapemri, headmodel, sourcemodel2d, sourcemodel3d);
    
    % check whether the files exist FIXME
    
    %% The following is checking whether the output is present as expected
    % ensure that the expected pipeline output is present
    % see https://wiki.humanconnectome.org/display/EEG/MEG+anatomy
    cd(outputdir);
    hcp_check_pipelineoutput('anatomy', 'subject', subjectid);
    
end
