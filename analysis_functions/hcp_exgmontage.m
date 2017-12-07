function [montage] = hcp_exgmontage(subjectid, experimentid, scanid)

% HCP_EXGMONTAGE produces a montage for the EEG channel data.
%
% Use as
%   montage = hcp_exgmontage(subjectid)
% or
%   montage = hcp_exgmontage(subjectid, experimentid, scanid)

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

if nargin < 3
    % let the sessions with cap prevail
    scanid = '_B';
end

if ~isempty(strfind(scanid, 'Rnoise')) || ~isempty(strfind(scanid, 'Pnoise'))
    montage = [];
else
    
    switch subjectid
        case {'CP10018', 'CP10019', 'CP10033', 'CP10053', 'CP10066', 'CP10141', 'CP10138', 'CP10167'}
            % scans without EEG cap ('CP10066' has unuseful EEG data)
            montage = [];
            montage.labelorg = { 'E29' 'E30' 'E31' 'E32'}';
            montage.labelnew = {'HEOG' 'ECG' 'VEOG'}';
            montage.tra = [1 -1 0 0; 0 0 1 0; 0 0 0 1];
            
        case 'CP10040'
            % scans with EEG cap
            montage = [];
            montage.labelorg = {'E31' 'E32' 'E64'}';
            montage.labelnew = {'ECG' 'VEOG' 'HEOG'}';
            montage.tra = eye(3);
            
        case {'CP10106'}
            % scans with EEG cap
            montage = [];
            montage.labelorg = { 'E29' 'E30' 'E31' 'E32'}';
            montage.labelnew = {'HEOG' 'ECG' 'VEOG'}';
            montage.tra = [1 -1 0 0; 0 0 1 0; 0 0 0 1];
            
        case {'CP10113'}
            % scans with EEG cap
            montage = [];
            montage.labelorg = { 'E29' 'E30' 'E31' 'E32'}';
            montage.labelnew = {'HEOG' 'ECG' 'VEOG'}';
            montage.tra = [1 -1 0 0; 0 0 1 0; 0 0 0 1];
            
        case {'CP10168'}
            % scans with EEG cap
            montage = [];
            montage.labelorg = {'HEOG' 'ECG' 'VEOG'}';
            montage.labelnew = {'HEOG' 'ECG' 'VEOG'}';
            montage.tra = eye(numel(montage.labelnew));
            
        case {'CP10128', 'CP10129'}
            % these are pilot subjects recorded in Glasgow
            error('Subject %s is not supported', subjectid);
            
            
        case {'CP10172', 'CP10173', 'CP10174', 'CP10188', 'CP10192', 'CP10195', 'CP10197', 'CP10198'};
            % These are the latest scans with finalized protocols.
            % According to Abbas email the channels should be
            % ECG: monopolar electrodes 'E1' and 'E2'
            % VEOG: monopolar electrodes 'E3' and 'E4'
            % HEOG: monopolar electrodes 'E5' and 'E6'
            % EMG for hands and feet: bipolar electrodes 'E31', 'E32', 'E63', 'E64'
            montage = [];
            montage.labelorg = { 'ECG+' 'ECG-' 'VEOG+' 'VEOG-' 'HEOG+' 'HEOG-'}';
            montage.labelnew = {'ECG' 'VEOG' 'HEOG'}';
            montage.tra = [1 -1 0 0 0 0; 0 0 1 -1 0 0;0 0 0 0 1 -1];
            
            
        otherwise
            % ASSUMING THAT ALL SUBJECTS THAT END UP HERE ARE FROM PHASE 2
            % WHERE THE ELECTRODE CONFIGURATION IS STANDARDIZED.
            montage = [];
            montage.labelorg = { 'ECG+' 'ECG-' 'VEOG+' 'VEOG-' 'HEOG+' 'HEOG-'}';
            montage.labelnew = {'ECG' 'VEOG' 'HEOG'}';
            montage.tra = [1 -1 0 0 0 0; 0 0 1 -1 0 0;0 0 0 0 1 -1];
            
    end % switch
    
    
    if ~isempty(regexp(scanid, 'Motort$', 'once'))
        % VERY QUICK FIX for the latest 2 scans so that EOC ECG and EMG channels all exist in the data
        montagenew = [];
        montagenew.labelorg = {'ECG+' 'ECG-' 'VEOG+' 'VEOG-' 'HEOG+' 'HEOG-' 'EMG_LH' 'EMG_RH' 'EMG_LF' 'EMG_RF'}';
        montagenew.labelnew = {'ECG' 'VEOG' 'HEOG' 'EMG_LH' 'EMG_RH' 'EMG_LF' 'EMG_RF' }';
        montagenew.tra = zeros(7, 10);
        montagenew.tra(1:3, 1:6) = montage.tra;
        montagenew.tra(4:7, 7:10) = eye(4);
        
        montage = montagenew;
        
    elseif ~isempty(strfind(scanid, 'Motort')) && isempty(strfind(scanid, '_VG'))
        % Motor task scans collected in SLU
        Motorlabelorg = { 'E31' 'E32' 'E63' 'E64'}';
        Motorlabelnew = { 'EMG_LH' 'EMG_RH' 'EMG_LF' 'EMG_RF' }';
        
        idxcol = [];
        for i=1:length(Motorlabelorg)
            sel = find(strcmp(Motorlabelorg(i), montage.labelorg));
            if ~isempty(sel)
                idxcol = [idxcol sel];
            end
        end
        tra = montage.tra;
        labelorg = montage.labelorg;
        labelnew = montage.labelnew;
        idxrow = find(any(tra(:, idxcol)')');
        tra(:, idxcol) = [];
        tra(idxrow, :) = [];
        labelorg(idxcol) = '';
        labelnew(idxrow) = '';
        
        montage.labelorg = [labelorg; Motorlabelorg];
        montage.labelnew = [labelnew; Motorlabelnew];
        montage.tra = [tra, zeros(size(tra, 1), 4); zeros(4, size(tra, 2)), eye(4)];
        
        % sort "montage.labelorg"
        lab = char(montage.labelorg);
        [sx, si] = sort(str2double(cellstr(lab(:, 2:end))));
        montage.labelorg = montage.labelorg(si);
        montage.tra = montage.tra(:, si);
        
    elseif ~isempty(strfind(scanid, 'Motort')) && ~isempty(strfind(scanid, '_VG'))
        error('Subject %s is not supported for MOTOR task', subjectid);
    end
    
end
