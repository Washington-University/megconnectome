function hcp_check_pipelineoutput(pipeline, varargin)

% HCP_CHECK_PIPELINEOUTPUT checks that the expected output files from a
% certain pipeline exist. The expected files are determined based on the
% documentation on the HCP wiki. The purpose of this script is to call it
% at the end of the respective pipeline, and at the beginning of any other
% pipeline that expects the output of the respective pipeline.
%
% Use as
%   hcp_check_pipelineoutput(pipeline, ...)
% where the input is a string with the pipeline name (see below).
%
% The additional input arguments should be specified as key-value pairs and
% may include
%   subject     = string with subject ID,    e.g. CP10128
%   experiment  = string with experiment ID, e.g. CP10128_MEG_v2
%   scan        = string with scan ID,       e.g. 2-Wrkmem_MBP_VG1
%
% In general for the pipeline described on
%   https://wiki.humanconnectome.org/display/EEG/xxxx
% you would call this function as
%   hcp_check_pipelineoutput('xxxx', ...)
%
% The pipelines included in this function are documented ar
%    https://wiki.humanconnectome.org/display/EEG/
%
% Depending on the pipeline name, subject, experiment, and scan identifier,
% the expected files according to the documentation are determined and their
% presence is checked. Please refer to the code of this function and the
% online documentation of the wiki for the actual output file names of each
% function.

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

% get the optional input arguments
subject     = ft_getopt(varargin, 'subject', '');
experiment  = ft_getopt(varargin, 'experiment', '');
scan        = ft_getopt(varargin, 'scan', '');
band        = ft_getopt(varargin, 'band', '');
net_seed    = ft_getopt(varargin, 'net_seed', '');
freq        = ft_getopt(varargin, 'freq', '');
datagroup   = ft_getopt(varargin, 'datagroup', '');
contrasts   = ft_getopt(varargin, 'contrasts', '');
avgmode     = ft_getopt(varargin, 'avgmode', '');
gridtype    = ft_getopt(varargin, 'gridtype', '');
savedir     = ft_getopt(varargin, 'savedir', '');
sourcemodel_type = ft_getopt(varargin, 'sourcemodel', '');

switch pipeline
    
    case 'baddata'
        assert_file('%s_%s_baddata_badchannels.txt', experiment, scan);
        assert_file('%s_%s_baddata_badsegments.txt', experiment, scan);
        assert_file('%s_%s_baddata_manual_badsegments.txt', experiment, scan);
        assert_file('%s_%s_baddata_manual_badchannels.txt', experiment, scan);
        
        assert_file('%s_%s_baddata_badchan_cor_scatter.png', experiment, scan);
        assert_file('%s_%s_baddata_badchan_cor_topo3D.png', experiment, scan);
        assert_file('%s_%s_baddata_badchan_cor_topo.png', experiment, scan);
        assert_file('%s_%s_baddata_badchan_std_scatter.png', experiment, scan);
        assert_file('%s_%s_baddata_badchan_std_topo.png', experiment, scan);
        
        tok=tokenize(scan,'-');
        scanmnem=tok{2};
        isTask=strcmp(scanmnem,'Wrkmem')|strcmp(scanmnem,'Motort')|strcmp(scanmnem,'StoryM');
        if isTask,
            assert_file('%s_%s_baddata_rawtrialinfo_QC.txt', experiment, scan);
        end
        
    case 'icaclass'
        assert_file('%s_%s_icaclass.mat', experiment, scan);
        assert_file('%s_%s_icaclass_1.png', experiment, scan); % there can be more figures, but at least one is expected
        
        %-------------------------------------------------------
        % ---- Check if the data is from MOTOR TASK for which no ECG EOG channels
        % are available and construct pseudo ECG EOG channels based on topology
        % templates. As we do not use EEG anymore in Phase 2 , MOTOR task data
        % should also have EO/CG channels and this part should not be necessary
        oldMotorExperimentIds = {
            'CP10128_MEG_v1'  % Glasgow scan
            'CP10129_MEG_v1'  % Glasgow scan
            'CP10141_MEG_v1'  % There are 2 scans with 2 difference EMG electrode locations. For each scan there are 3 blocks for each hand and foot with 12 trials in each block. So each condition has only 36 trials.
            'CP10138_MEG'     % There are 6 different scan with variable ISI and time jittering in some. For each scan there are 3 blocks for each hand and foot with beween 8 to 12 trials in each block. So each condition has only 24 to 36 trials.
            'CP10167_MEG'     % 1 scan with the latest protocol. 8 blocks per hand and foot with 10 trials each. 1200 msec ISI. However there is no MRI
            'CP10113_MEG_v2'  % There are 2 different scans. NOt clear what the difference is between the two. There are 8 blocks for each hand and foot, with 10 trials per block. Unfortunately EMG signal is bad for both Hands and trials cannot be extracted based on the EMG signal
            };
        
        isMotorTask = ~isempty(regexp(scan,'Motor', 'once'));
        isScanWithOldMotor = ismember(experiment, oldMotorExperimentIds);
        if isMotorTask && isScanWithOldMotor
            assert_file('%s_%s_icaclass_virtchanVEOG.png', experiment, scan); % there can be more figures, but at least one is expected
            assert_file('%s_%s_icaclass_virtchanECG.png', experiment, scan); % there can be more figures, but at least one is expected
        end
        %-------------------------------------------------------
        
    case 'icaclass_qc'
        assert_file('%s_%s_icaclass_vs.mat', experiment, scan);
        assert_file('%s_%s_icaclass_vs.txt', experiment, scan);
        assert_file('%s_%s_icaclass_vs_1.png', experiment, scan); % there can be more figures, but at least one is expected
        
    case 'icaqc'
        error('the "icaqc" pipeline has been deprecated, please use the "baddata" pipeline instead');
        
    case 'rmegpreproc'
        assert_file('%s_%s_rmegpreproc.mat', experiment, scan);
        
    case 'tmegpreproc'
        %----------------
        tmpscan = scan;
        isInCell = iscell(tmpscan);
        if isInCell,
            Nfiles = length(tmpscan);
        else
            Nfiles = 1;
            tmpscan = {tmpscan};
        end
        for iScan = 1:Nfiles,
            assert_file('%s_%s_tmegpreproc_trialinfo.mat', experiment, tmpscan{iScan});
        end
        %----------------
        if ~isempty(datagroup)
            tmpgroup = datagroup;
            isInCell = iscell(tmpgroup);
            if isInCell,
                Ngroups = length(tmpgroup);
            else
                Ngroups = 1;
                tmpgroup = {tmpgroup};
            end
            %----------------
            for iFile = 1:Nfiles
                for iGroup = 1:Ngroups
                    assert_file('%s_%s_tmegpreproc_%s.mat', experiment, tmpscan{iFile}, tmpgroup{iGroup});
                end
            end
        else
            disp('No datagroup names provided - Checking only for trialinfo information');
        end
        %----------------
        
    case 'eravg'
        %----------------
        tmpscan = scan;
        isInCell = iscell(tmpscan);
        if isInCell,
            Nfiles = length(tmpscan);
        else
            Nfiles = 1;
            tmpscan = {tmpscan};
        end
        
        %----------------
        tmpcontr = contrasts;
        isInCell = iscell(tmpcontr);
        if isInCell,
            Ncontr = length(tmpcontr);
        else
            Ncontr = 1;
            tmpcontr = {tmpcontr};
        end
        %----------------
        tmpavgmode = avgmode;
        isInCell = iscell(tmpavgmode);
        if isInCell,
            Navgmodes = length(tmpavgmode);
        else
            Navgmodes = 1;
            tmpavgmode = {tmpavgmode};
        end
        %----------------
        for iFile = 1:Nfiles
            for iCntr = 1:Ncontr
                for iAvgMode = 1:Navgmodes
                    if strcmp(tmpavgmode{iAvgMode},'mag')
                        assert_file('%s_%s_%s_[MODE-mag].mat', experiment, tmpscan{iFile}, tmpcontr{iCntr});
                    elseif strcmp(tmpavgmode{iAvgMode},'planar')
                        assert_file('%s_%s_%s_[MODE-planar].mat', experiment, tmpscan{iFile}, tmpcontr{iCntr});
                    end
                end
            end
        end
        
    case 'tfavg'
        %----------------
        tmpscan = scan;
        isInCell = iscell(tmpscan);
        if isInCell,
            Nfiles = length(tmpscan);
        else
            Nfiles = 1;
            tmpscan = {tmpscan};
        end
        %----------------
        tmpcontr = contrasts;
        isInCell = iscell(tmpcontr);
        if isInCell,
            Ncontr = length(tmpcontr);
        else
            Ncontr = 1;
            tmpcontr = {tmpcontr};
        end
        %----------------
        tmpavgmode = avgmode;
        isInCell = iscell(tmpavgmode);
        if isInCell,
            Navgmodes = length(tmpavgmode);
        else
            Navgmodes = 1;
            tmpavgmode = {tmpavgmode};
        end
        %----------------
        for iFile = 1:Nfiles
            for iCntr = 1:Ncontr
                for iAvgMode = 1:Navgmodes
                    if strcmp(tmpavgmode{iAvgMode},'mag')
                        assert_file('%s_%s_%s_[MODE-mag].mat', experiment, tmpscan{iFile}, tmpcontr{iCntr});
                    elseif strcmp(tmpavgmode{iAvgMode},'planar')
                        assert_file('%s_%s_%s_[MODE-planar].mat', experiment, tmpscan{iFile}, tmpcontr{iCntr});
                    end
                end
            end
        end
        
    case 'tfsens'
        %----------------
        tmpscan = scan;
        isInCell = iscell(tmpscan);
        if isInCell,
            Nfiles = length(tmpscan);
        else
            Nfiles = 1;
            tmpscan = {tmpscan};
        end
        %----------------
        tmpgroup = datagroup;
        isInCell = iscell(tmpgroup);
        if isInCell,
            Ngroups = length(tmpgroup);
        else
            Ngroups = 1;
            tmpgroup = {tmpgroup};
        end
        %----------------
        tmpband = band;
        isInCell = iscell(tmpband);
        if isInCell,
            Nbands = length(tmpband);
        else
            Nbands = 1;
            tmpband = {tmpband};
        end
        %----------------
        for iFile = 1:Nfiles
            for iGroup = 1:Ngroups
                for iBand = 1:Nbands
                    assert_file('%s_%s_tfsens_%s_%s.mat', experiment, tmpscan{iFile}, tmpgroup{iGroup}, tmpband{iBand});
                end
            end
        end
        
    case 'srcavglcmv'
        %----------------
        tmpscan = scan;
        isInCell = iscell(tmpscan);
        if isInCell,
            Nfiles = length(tmpscan);
        else
            Nfiles = 1;
            tmpscan = {tmpscan};
        end
        %----------------
        tmpcontr = contrasts;
        isInCell = iscell(tmpcontr);
        if isInCell,
            Ncontr = length(tmpcontr);
        else
            Ncontr = 1;
            tmpcontr = {tmpcontr};
        end
        %----------------
        if ~isempty(savedir)
            experiment=[savedir,experiment];
        end
        %----------------
        for iFile = 1:Nfiles
            for iCntr = 1:Ncontr
                assert_file('%s_%s_%s.mat', experiment, tmpscan{iFile}, tmpcontr{iCntr});
            end
        end
        
    case 'anatomy'
        assert_file('%s_MEG_anatomy_anatomical.nii',       subject);
        assert_file('%s_MEG_anatomy_fiducials.txt',        subject);
        assert_file('%s_MEG_anatomy_landmarks.txt',        subject);
        assert_file('%s_MEG_anatomy_transform.txt',        subject);
        assert_file('%s_MEG_anatomy_sourcemodel_2d.mat',    subject);
        assert_file('%s_MEG_anatomy_sourcemodel_3d4mm.mat', subject);
        assert_file('%s_MEG_anatomy_sourcemodel_3d6mm.mat', subject);
        assert_file('%s_MEG_anatomy_sourcemodel_3d8mm.mat', subject);
        assert_file('%s_MEG_anatomy_headshape.mat',        subject);
        assert_file('%s_MEG_anatomy_headshapemri.mat',     subject);
        assert_file('%s_MEG_anatomy_headmodel.mat',        subject);
        assert_file('%s_MEG_anatomy_headmodel.png',        subject);
        assert_file('%s_MEG_anatomy_sourcemodel_2d.png',   subject);
        assert_file('%s_MEG_anatomy_sourcemodel_3d4mm.png', subject);
        assert_file('%s_MEG_anatomy_slice1.png',           subject);
        assert_file('%s_MEG_anatomy_slice2.png',           subject);
        assert_file('%s_MEG_anatomy_slice3.png',           subject);
        assert_file('%s_MEG_anatomy_headshape.png',        subject);
        
    case 'datacheck'
        assert_file('%s_%s_datacheck_info.txt', experiment, scan);
        % assert_file('%s_%s_datacheck_headshape.png', experiment, scan); % this one does not exist for empty room recordings
        assert_file('%s_%s_datacheck_jumps.png', experiment, scan);
        assert_file('%s_%s_datacheck_MEG_lowfreq_power.png', experiment, scan);
        assert_file('%s_%s_datacheck_MEG_powerline_noise.png',experiment, scan);
        assert_file('%s_%s_datacheck_MEG_powspctrm.png', experiment, scan);
        % assert_file('%s_%s_datacheck_EEG_powspctrm.png', experiment, scan); % not relevant any more in the Phase 2 protocol
        assert_file('%s_%s_datacheck_MEGREF_powspctrm.png', experiment, scan);
        assert_file('%s_%s_datacheck_neighb_correlation.png', experiment, scan);
        
    case 'icamne'
        assert_file('%s_%s_icamne.mat', experiment, scan);
        assert_file('%s_%s_icamne_1.png', experiment, scan);% there can be more figures, but at least one is expected
        
    case 'icaimagcoh'
        assert_file('%s_%s_icaimagcoh_%s.mat', experiment, scan, band);
        
    case 'icablpenv'
        assert_file('%s_%s_icablpenv_%s.power.dtseries.nii', experiment, scan, band);
    
    case 'bfblpenv'
        assert_file('%s_%s_bfblpenv_%s.power.dtseries.nii', experiment, scan, band);
        
    case 'icamcw'
        assert_file('%s_%s_mcw_%s_%s_%s.mat', experiment, scan, sourcemodel_type, net_seed, band);
        
    case 'powavg'
        assert_file('%s_%s_powavg.mat', experiment, scan);
        assert_file('%s_%s_powavg_multiplot.png', experiment, scan);
        assert_file('%s_%s_powavg_singleplot.png', experiment, scan);
        
    otherwise
        error('unknown pipeline "%s"', pipeline);
        
end % switch pipeline

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTION to keep the code above more clean
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function assert_file(varargin)
% the following will check for the file being present anywhere on the path
filename = sprintf(varargin{:});
if ~exist(filename, 'file')
    error('the file "%s" does not exist', filename);
else
    fprintf('found "%s"\n', filename);
end
