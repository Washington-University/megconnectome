function[trl,trlInfoColDescr,trialSummary,scanStartSamp,scanEndSamp,warninfo] = trialfun_Motort_BaseExtractAll( cfg )
%% This is the Trial Definition Function for Motor experiment.
% It extracts the trial definition for ALL trials within each of different
% datagroup.
% There are 2 different data groups for Motor task.
%
% DATA GROUPS
%----------------------------------
% 1. mnemonic: TFLA.     description: Onset of flashing cross that instructs 
%                                       subject to perform movement by hand or foot.
% 2. mnemonic: TEMG. 	 description: Onset of the emg signal from hand or 
%                                       foot recorded muscles.
%
% 
%
% INPUT VARIABLE
%----------------------------------
% cfg : This is a structure containing information required for extracting
%       the trials for either of the 2 data groups described above.
%       Fields:
%               .datafile: This is the filename of the raw data file.
%               .trialdef: This is a structure containing the parameters
%                          required to split the data into the trials of either data
%                          group.
%                          Fields:
%                          .trialdef.TrigBasedOnEMG = 'yes' or 'no';   % Defines the 0 reference time of eahc trial
%                                                                          according to the desired data group.
%                                                                          'yes' for TEMG,'no' for TFLA data group
%                                                                      
%                          .trialdef.cutMode = 'trials';               % This representes that data will be cut in trials
%                                                                      %  and not in blocks. The same for both data groups
%                          .trialdef.preStimTime = 1.2;                % Time interval prior to 0 reference point for each trial 
%                          .trialdef.postStimTime = 1.2;               % Time interval after the 0 reference point for each trial 
%                          .trialdef.montage                % This is the montage containing the emg channels. This
%                                                               variable is constructed by the hcp_exgmontage.m function. 
%                                                               In .labelnew subfield , the expected 
%                                                               emg channels names are  expected to be 
%                                                               {'EMG_LH','EMG_LF','EMG_RH','EMG_RF'}         
%
%
%
% OUTPUT VARIABLES
%-----------------------------------
% trl: This is a numerical matrix.  Each column corresponds to a trial and
%      each column to a specific condition of piece of information, quantified numerically , regarding
%      each trial. The first 3 Columns describe the trial start sample, end sample and time offset 
%      of the start of the trial relative to the 0 reference point. These 3
%      columns are used by fieldtrip as the necessary information required
%      to extract the data for each trial.
%      The rest of the columns encode various types of information about
%      each trial. This part of the trl Matrix , from column 4 to the last
%      column, should be refered to as 'trialinfo' part of the trl matrix.
%
% trlInfoColDescr: This is a cell array. Each element is a string
%                  describing the type of information encoded by the corresponding column of
%                  the 'trialinfo' part of trl matrix described above. 
%
% trialSummary:    This is a structure containing an overview of the breakdown of the data
%                  in trials of the main conditions for BOTH data groups.
%                  This is used for Quality Control purposes in order to
%                  ensure that the correct number of trials for each of
%                  the main conditions is present in the data and identifying any
%                  discrepances due to problems in the stimulus presentation protocol.
%
% scanStartSamp:   This is the sample number of the onset of the first
%                  stimulus block.
% scanEndSamp:     This is the sample number of the offset of the last
%                  stimulus block.
% warninfo:        This is cell variable used for Quality Control , in
%                  order to catch problems with unexpected trigger values.
%
%
% Below is presented, for convenience of the reader, the description of
% each column of the 'trialinfo' part of trl matrix described above and 
% contained in the trlInfoColDescr output variable. 
%
%
% Trialinfo Column Description as presented in trlInfoColDescr output variable
%-------------------------------------------------------------------------------
%========================================================================
% data group: TFLA
%------------------
%  ' 1. Block Index  '
%  ' 2. Block Stim Code:    1-Left Hand,  2 - Left Foot, 4 - Right Hand. 5 - Right Foot, 6 - Fixation'
%  ' 3. Trial Index in Block 
%  ' 4. Trial Onset Sample'
%  ' 5. prev. Block Stim Code'};
%========================================================================
% data group: TEMG
%------------------
%  ' 1. Block Index  '
%  ' 2. Block Stim Code:    1-Left Hand,  2 - Left Foot, 4 - Right Hand. 5 - Right Foot, 6 - Fixation'
%  ' 3. Trial Index in Block (This is derived by finding the flash cross onset just before the EMG onset)'
%  ' 4. Trial EMG Onset Sample'
%  ' 5. prev. Block Stim Code'
%  ' 6. Time from EMG onset to the previous Flashing Cross'
%========================================================================

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

%__________________________________________________________________________

%%
if ~isfield(cfg.trialdef, 'TrigBasedOnEMG'),    cfg.trialdef.TrigBasedOnEMG = 'yes'; end
if ~isfield(cfg.trialdef, 'cutMode'),           cfg.trialdef.cutMode = 'trials'; end
if ~isfield(cfg.trialdef, 'plotresults'),       cfg.trialdef.plotresults = 'no'; end

datafile = cfg.datafile;
%stimulCode = [1 2 4 5]; %These define motor block codes; for Fixation the code is 6 which is hard coded in the code processing fixation

TrigBasedOnEMG = cfg.trialdef.TrigBasedOnEMG;
cutMode = cfg.trialdef.cutMode;% 'trials' or 'blocks'
prestimTime = cfg.trialdef.prestimTime;
poststimTime = cfg.trialdef.poststimTime;
montage=cfg.trialdef.montage;

% summaryfile=cfg.trialdef.summaryfile;
if strcmpi(cfg.trialdef.plotresults , 'yes' ),
    plot_flag = 1;
    %1->Left Hand, 2->Left Foot, 3->Tongue, 4->Right Hand, 5->Right Foot, 6->Fixation
    %titletxt = {'Left Hand', 'Left Foot', 'Tongue', 'Right Hand', 'Right Foot', 'Fixation' };
else
    plot_flag = 0;
end


%{
% Trigger information
%-----------------------------------------------------------------
According to E-prime files the Trigger codes are the following:
%----------
LeftHandCueProcedure :  18
CrossLeftPROC: 22
%----------
LeftFootCuePROC:  34
CrossLeftPROC : 38
%----------
RightHandCuePROC: 66
CrossRightPROC: 70
%----------
RightFoottCuePROC: 130
CrossRightPROC: 134
%----------
These can be modelled as:

%-----------------------
Photodiode:
The photodiode is activated on the onset of the Cues and the onset of the
flashing Cross.
%}

%==========================================================================
%%  READ TRIGGER FROM DATA FILES
hdr = ft_read_header(datafile);
Fsample = hdr.Fs;
detailedTrigger = ft_read_data(datafile,'chanindx',1,'header',hdr,'eventformat','4d','dataformat','4d');
Nsamples=length(detailedTrigger);
prestimSamples = floor(prestimTime*Fsample);
poststimSamples = floor(poststimTime*Fsample);

% SPLIT TRIGGERS FROM PARALLEL PORT AND TRIGGERS FROM PHOTODIODE
trigComb=detailedTrigger;
trigParal=trigComb-256*(floor(trigComb./256));
trigPhoto=256*(floor(trigComb./256));

%----------------------------------------------------------
% the following contains informatin about e-prime triggers for cues and
% crosses. COL1: trigger calue. COL2: 1:cue 2:cross. COL3:1:Left Hand
% 2:Left Foot 4: Right Hand  5: Right Foot 0:Fixation
% The trigger value of 2 at then end corresponds to fixation block

epTrigInfo=[18  1 1
    22  2 1
    34  1 2
    38  2 2
    66  1 4
    70  2 4
    130 1 5
    134 2 5
    2   0 6];

%-----------------------------------------------------------
% Find the onsets of the photodiode and the onsets of cues and crosses ni
% eprime. These 2 sets should have the same number of events and their time
% difference should be in the order of milliseconds jsut due to the
% projector delay.

epTrigsCueAndCross=epTrigInfo(ismember(epTrigInfo(:,2),[1 2]),1)';

%---
difTrigParal=diff(trigParal);
fwdTrigParal=trigParal(2:end);
bwdTrigParal=trigParal(1:end-1);
%---
difTrigPhoto=diff(trigPhoto);
fwdTrigPhoto=trigPhoto(2:end);
bwdTrigPhoto=trigPhoto(1:end-1);
%---
indParalUP=find(ismember(fwdTrigParal,epTrigsCueAndCross) & (fwdTrigParal~=bwdTrigParal) )+1; % Find onset of parallel trigger for cues and crosses
indPhotoUP=find(difTrigPhoto>0)+1;
%%
%=========================================
%-- Check correspondence between photodiode and parallel port triggered
%--  cue and cross events.
% Create stimMap variable with Columns:
% COL1: 1:cue 2:cross 0: fixation.
% COL2: 1:Left Hand, 2:Left Foot 4: Right Hand  5: Right Foot 0:Fixation
% COL3: Event onset based on Parallel port trigger
% COL4: Event onset based on Parallel port trigger
% COL5: Block Number in Run
% COL6: Trial Sequence number within a block of Flashing Cross events. For cues and fixation this is 0.
% COL7: Parallel Port Trigger Code
warninfo=[];

badThresh=4/60; %This is the maximum time difference threshold between a trigger in photodiode and the corresponding trigger in parallel port;
badThreshSamps=floor(badThresh.*hdr.Fs);
stimMap=[];
if length(indParalUP)~=length(indPhotoUP)
    %error('The number of onsets from photodiode is not the same with the number of cue and cross onset derived from e-prime trigger values');
    warning('The number of onsets from photodiode is not the same with the number of cue and cross onset derived from e-prime trigger values');
    warninfo.dif_ph2pp.difNevents=length(indPhotoUP)-length(indParalUP);
    
    countStims=1;
    for iStim=1:length(indParalUP)
        indIn=find(abs(indPhotoUP-indParalUP(iStim))<=badThreshSamps,1,'first');
        if isempty(indIn)
            error(['No close-by photodiode trigger found for stimulus at sample:',num2str(indParalUP(iStim))]);
        else
            addinfo=epTrigInfo(find(ismember(epTrigInfo(:,1),trigParal(indParalUP(iStim)))),2:end);
            stimMap(countStims,1:4)=[addinfo indParalUP(iStim) indPhotoUP(indIn)];
            stimMap(countStims,7)=trigParal(indParalUP(iStim));
        end
        countStims=countStims+1;
    end
else
    tmpTimeDif=(indPhotoUP-indParalUP)./hdr.Fs;
    indbad=find(abs(tmpTimeDif)>badThresh);
    if ~isempty(indbad),
        warning('there seem to be time difference > 4 refresh frames between the parallel port and photodiode triggers for some stimulus')
        warninfo.dif_ph2pp.samp_longdelay=indPhotoUP(indbad);
        
        stimMap=[];
        countStims=1;
        for iStim=1:length(indParalUP)
            indIn=find(abs(indPhotoUP-indParalUP(iStim))<=badThreshSamps,1,'first');
            if ~isempty(indIn)
                error(['No close-by photodiode trigger found for stimulus at sample:',num2str(indParalUP(iStim))]);
            else
                addinfo=epTrigInfo(find(ismember(epTrigInfo(:,1),trigParal(indParalUP(iStim)))),2:end);
                stimMap(countStims,1:4)=[addinfo indParalUP(iStim) indPhotoUP(indIn)];
                stimMap(countStims,7)=trigParal(indParalUP(iStim));
                
            end
            countStims=countStims+1;
        end
    else
        stimMap=[];
        for iStim=1:length(indParalUP)
            addinfo=epTrigInfo(find(ismember(epTrigInfo(:,1),trigParal(indParalUP(iStim)))),2:end);
            stimMap(iStim,1:4)=[addinfo indParalUP(iStim) indPhotoUP(iStim)];
            stimMap(iStim,7)=trigParal(indParalUP(iStim));
        end
    end
end

tmpDifPH2PP=stimMap(:,4)-stimMap(:,3);
if  nanmean(abs(tmpDifPH2PP)./Fsample) > (2/60)
    warning('The avergae difference between Photodiode and Parallel Port Triggers is more than 2 monitor refresh rates.'); 
    warninfo.dif_ph2pp.avgdelay=nanmean(abs(tmpDifPH2PP)./Fsample);    
end
%-------------------------------------------------------------------
%% Add fixation bock information in stimMap variable.
%  Fixation blocks can only be defined by the parallel port

indParalFix=find((fwdTrigParal==2)&(bwdTrigParal~=2))+1;

Nfix=length(indParalFix);
addinfo=epTrigInfo(find(ismember(epTrigInfo(:,1),2)),2:end);
stimMap=[stimMap; [repmat(addinfo,Nfix,1)  indParalFix' indParalFix' nan(Nfix,2,1) repmat(2,Nfix,1)]]; % add fixation block information to stimMap.
stimMap=sortrows(stimMap,3);    % sort all events so that fixation blocks end up at the right sequence
%% Derive the sequence of cross onsets within each block of trials and the sequence of Blocks
blockNum=0;
seqTrlNum=1;
for iEv=1:size(stimMap,1)
    
    if stimMap(iEv,1)==1, % if a cue
        blockNum=blockNum+1;
        stimMap(iEv,5)=blockNum;
        stimMap(iEv,6)=0;  % This is a cue so trial sequence number set to 0
        seqTrlNum=1;
    elseif stimMap(iEv,1)==2, % cross case
        stimMap(iEv,5)=blockNum;
        stimMap(iEv,6)=seqTrlNum;
        seqTrlNum=seqTrlNum+1;
    elseif stimMap(iEv,1)==0,   % fixatin case
        blockNum=blockNum+1;
        stimMap(iEv,5)=blockNum;
        stimMap(iEv,6)=0;  % This is fixation so trial sequence number set to 0
        seqTrlNum=1;
    end
end
%%
%-------------------------------
%  Fixation blocks  cannot be defined by the photodiode directly as there is
%  not a photrodiode trigger for this. But it can be inferred by adding the
%  inter stimulus interval of flashing cross , 1200 msec, to the time of
%  the last photodiode trigger from flashing cross within a block.

ISI=1.2; %inter stimulus interval in sec
ISIsamps=ceil(ISI.*hdr.Fs);
tmpIndFix=find(stimMap(:,1)==0);
indPhotoFix=stimMap(tmpIndFix-1,4)+ISIsamps;
stimMap(tmpIndFix,4)=indPhotoFix;

%=================================================================
%% FORM MOTOR TRIAL INFO based on Flashing Cross
% This  creates a trial information matrix only for the flashing cross
% trials with the following columns:
%Columns:
%-------------------
% 1. Block Index
% 2. Block Stim Code: 1:Left Hand, 2:Left Foot 4: Right Hand  5: Right Foot
% 3. Trial Index in Block
% 4. Trial Onset Trigger Sample
% 5. prev. Block Stim Code

motorTrialInfo=[];
prevBlockStimcode=nan;
countCrosses=1;
for iEv=1:length(stimMap)
    if (stimMap(iEv,1)==1)|(stimMap(iEv,1)==0), % cue or fixatoin
        curBlockNum=stimMap(iEv,5);
        if iEv>1, prevBlockStimcode=stimMap(iEv-1,2); end;
    else
        motorTrialInfo(countCrosses,:)=[stimMap(iEv,[5 2 6 4]) prevBlockStimcode];
        countCrosses=countCrosses+1;
    end
end

%===============================================
%% FROM BLOCK INFO
% This creates a matrix with information about all blocks.
% allBlockInfo
%===Columns===
% 1. Block Index
% 2. Stim Code
% 3. Cue Onset sample . For fixation this is the start sample of the fixation block
% 4. Cue Offset sample . For fixation nan
% 5. Cross 1 Onset Sample. For fixation nan
% 6. Cross 2 Onset Sample. For fixation nan
% 7. Cross 3 Onset Sample. For fixation nan
% 8. Cross 4 Onset Sample. For fixation nan
% 9. Cross 5 Onset Sample. For fixation nan
% 10. Cross 6 Onset Sample. For fixation nan
% 11. Cross 7 Onset Sample. For fixation nan
% 12. Cross 8 Onset Sample. For fixation nan
% 13. Cross 9 Onset Sample. For fixation nan
% 14. Cross 10 Onset Sample. For fixation nan
% 15. Block End Sample
% 16. Previous Block Code
fixBlockDur=15; % This is the expected duration of fixation blocks in sec
fixBlockSamps=floor(fixBlockDur*hdr.Fs);

indNonCross=find(stimMap(:,1)~=2);
noncrossMap=stimMap(indNonCross,:);
allBlockInfo=[];
countBlocks=1;
for iEv=1:length(noncrossMap)
    if (noncrossMap(iEv,1)==1), % cue
        allBlockInfo(countBlocks,1:3)=stimMap(indNonCross(iEv),[5 2 4]);
        % Find the sample of cue off from the photodiode trigger
        indCueOff=find(([1:Nsamples-1]>stimMap(indNonCross(iEv),4))&(fwdTrigPhoto~=trigPhoto(stimMap(indNonCross(iEv),4)))&(bwdTrigPhoto==trigPhoto(stimMap(indNonCross(iEv),4))),1,'first')+1;
        % % Find the sample of cue off from the parallel port trigger
        % indCueOff=find(([1:Nsamples-1]>stimMap(indNonCross(iEv),3))&(fwdTrigParal~=stimMap(indNonCross(iEv),7))&(bwdTrigParal==stimMap(indNonCross(iEv),7)),1,'first');
        allBlockInfo(countBlocks,4)=indCueOff;
        allBlockInfo(countBlocks,5:14)=stimMap(indNonCross(iEv)+1:indNonCross(iEv)+10,4)';
        allBlockInfo(countBlocks,15)=stimMap(indNonCross(iEv+1),4)-1;%Set the end of the block to the sample just before the onset of the next block
    elseif (noncrossMap(iEv,1)==0), % fixation
        allBlockInfo(countBlocks,1:3)=stimMap(indNonCross(iEv),[5 2 4]);
        allBlockInfo(countBlocks,4:14)=nan;
        if iEv<length(noncrossMap)
            allBlockInfo(countBlocks,15)=stimMap(indNonCross(iEv+1),4)-1;%Set the end of the block to the sample just before the onset of the next block
        else
            allBlockInfo(countBlocks,15)=stimMap(indNonCross(iEv),4)+fixBlockSamps; % This is the last fixation block which is also the last block in the scan
        end
    end
    if iEv>1
        allBlockInfo(countBlocks,16)=noncrossMap(iEv-1,2); % prev Block Stim Type
    else
        allBlockInfo(countBlocks,16)=nan; % prev Block Stim Type
    end
    countBlocks=countBlocks+1;
end

%-------------------------------------------
%======================================================
%% Identify the start and end of scan as the onset of the first cue and the end of the last block
% A default temporal padding of 3 second is used in either end
defEndPad=floor(3*Fsample);
scanStartSamp=allBlockInfo(1,3)-defEndPad;
if scanStartSamp<1 
    scanStartSamp=1;
end
scanEndSamp=allBlockInfo(end,15)+defEndPad;
if scanEndSamp>Nsamples,
    scanEndSamp=Nsamples;
end
%======================================================
%%
%==================================================
%% FORM FIXATION PSEUDO TRIALS
%  This code takes the fixation blocks and divides them in pseudo trials of
%  the same length as the Inter Stimulus Interval.
%-------------------
% 1. Block Index
% 2. Block Stim Code: 6: Fixation
% 3. pseudo Trial Index in Block
% 4. pseudo Trial Onset Trigger Sample
% 5. prev. Block Stim Code
fixBlockInfo=allBlockInfo((allBlockInfo(:,2)==6),:);

halfISIsamps=floor(ISIsamps./2);
fixTrialInfo=[];
for iBlock=1:size(fixBlockInfo,1),
    indxStart   = fixBlockInfo(iBlock,3)+halfISIsamps; %Exclude the first 600msec from fixation
    indxEnd     = fixBlockInfo(iBlock,15)-ISIsamps;
    pseudoTrlStartSamp=[indxStart:ISIsamps:indxEnd];
    
    tmpNtrials=length(pseudoTrlStartSamp);
    fixTrialInfo=[fixTrialInfo;
        repmat(fixBlockInfo(iBlock,1),tmpNtrials,1),...%1
        repmat(6,tmpNtrials,1),...%2
        [1:tmpNtrials]',...%3
        pseudoTrlStartSamp',...%4
        repmat(fixBlockInfo(iBlock,16),tmpNtrials,1)];%5
end

%==================================================
%==================================================
%==================================================
%==================================================
%==================================================
%%  LOAD EMG DATA

emgExpChans={'EMG_LH'
    'EMG_RH'
    'EMG_LF'
    'EMG_RF'};

hasEMG=1;
[ind1,ind2]=match_str(emgExpChans,montage.labelnew);
if isempty(ind1)
    hasEMG=0;
    disp('WARNING ! No EMG channels with the standard naming convention were found');
    disp('Assuming no EMG is available. Trials will be extracted based only on the flashing cross.');
end

if hasEMG
    if (length(ind1)~=4)
        error('Not all EMG channels are present or there are duplicates.');
        return;
    end
    if length(unique(emgExpChans(ind1)))~=4
        error('There are EMG channel duplicates.');
        return;
    end
    tmpTra=montage.tra(ind2,:);
    [k1,l1]=find(tmpTra~=0);
    indOrig=unique(l1);
    chanOrig=montage.labelorg(indOrig);
    chanNew=emgExpChans(ind1);
    newTra=tmpTra(:,indOrig);
    
    newmontage=montage;
    newmontage.labelorg=chanOrig;
    newmontage.labelnew=chanNew;
    newmontage.tra=newTra;
    
    [ind1,ind2]=match_str(newmontage.labelorg,hdr.label);
    if (length(ind2)~=4)
        error('Not all EMG channels were found in the dataset.');
        return;
    end
    
    emgdataorig= ft_read_data(datafile,'chanindx',ind2,'header',hdr,'eventformat','4d','dataformat','4d');
    emgdatanew=newmontage.tra*emgdataorig;
    
    [ind1,ind2]=match_str(chanNew,'EMG_LH');
    emgdata_LH=emgdatanew(ind1,:);
    [ind1,ind2]=match_str(chanNew,'EMG_LF');
    emgdata_LF=emgdatanew(ind1,:);
    [ind1,ind2]=match_str(chanNew,'EMG_RH');
    emgdata_RH=emgdatanew(ind1,:);
    [ind1,ind2]=match_str(chanNew,'EMG_RF');
    emgdata_RF=emgdatanew(ind1,:);
    clear emgdataorig emgdatanew;
else
    emgdata_LH=[];
    emgdata_LF=[];
    emgdata_RH=[];
    emgdata_RF=[];
end

%==================================================================
%%  MOTOR TRIAL INFO BASED ON EMG
%Columns:
%-------------------
% 1. Block Index
% 2. Block Stim Code
% 3. Trial Index in Block (This is derived by finding the flash cross onset just before the EMG onset)
% 4. Trial Onset EMG Sample
% 5. prev. Block Stim Code
% 6. Time from EMG onset to the previous Flashing Cross

motorBlockColumnDescription={''}; % IF BLOCKS ARE ANALYZED AS A WHOLE ADD DESCRIPTION
motorTrlColumnDescriptionEMG={' 1. Block Index  '
    ' 2. Block Stim Code:    1-Left Hand,  2 - Left Foot, 4 - Right Hand. 5 - Right Foot, 6 - Fixation'
    ' 3. Trial Index in Block (This is derived by finding the flash cross onset just before the EMG onset)'
    ' 4. Trial EMG Onset Sample'
    ' 5. prev. Block Stim Code'
    ' 6. Time from EMG onset to the previous Flashing Cross'};
motorTrlColumnDescription={' 1. Block Index  '
    ' 2. Block Stim Code:    1-Left Hand,  2 - Left Foot, 4 - Right Hand. 5 - Right Foot, 6 - Fixation'
    ' 3. Trial Index in Block '
    ' 4. Trial Onset Sample'
    ' 5. prev. Block Stim Code'};

motorBlockInfo=allBlockInfo((allBlockInfo(:,2)~=6),:);

tmpInfoEMG=[];
if hasEMG
    for iBlock = 1 : size(motorBlockInfo,1)
        indBlockStart   = motorBlockInfo(iBlock,4); % Cue Off point
        indBlockEnd     = motorBlockInfo(iBlock,15); % End of Block
        tmpBlockStimcode   = motorBlockInfo(iBlock,2);
        
        
        %  trg1 = bitand( trg(ix1) , trgmsk(tmpBlockCode) ) .* (bitand( trg(ix1) , trgcode(tmpBlockCode) ) == trgcode(tmpBlockCode) );
        indArrayBlock=indBlockStart:indBlockEnd;
        trigBlock=trigPhoto(indArrayBlock);
        
        if tmpBlockStimcode==1,
            emgdata=emgdata_LH;
        elseif tmpBlockStimcode==2,
            emgdata=emgdata_LF;
        elseif tmpBlockStimcode==4,
            emgdata=emgdata_RH;
        elseif tmpBlockStimcode==5,
            emgdata=emgdata_RF;
        end
        if isempty(emgdata)
            continue;
        end
        
        trigEMG = Extract_EMG_Trigger( emgdata(indBlockStart:indBlockEnd) , hdr.Fs , trigBlock );
        
        indTrlEMG = indArrayBlock( find( diff(trigEMG) > 0 ) + 1 );
        tmpNtrials=length(indTrlEMG);
        tmpInfoEMG=[tmpInfoEMG;
            repmat(motorBlockInfo(iBlock,1),tmpNtrials,1),...%1
            repmat(tmpBlockStimcode,tmpNtrials,1),...%2
            nan(tmpNtrials,1),...%3
            indTrlEMG',...%4
            repmat(motorBlockInfo(iBlock,16),tmpNtrials,1),...%5
            nan(tmpNtrials,1)];%6
        
        
    end
    
    %     if isempty(tmpInfoEMG)
    %         disp( [ 'There is no trial in "' cfg.datafile '"'])
    %         disp( 'based on parameters of following cfg:')
    %         disp( cfg.trialdef )
    %         trl = [];
    %         return;
    %     end
    %==================================================
    %-- Fuse with trial info from flashing cross
    %-------------------------------------------------------------------------
    NemgTrials=size(tmpInfoEMG,1);
    indPreEmpty=[];
    for iTrial=1:NemgTrials,
        iTrial
        tmpIndx=find(motorTrialInfo(:,4)<tmpInfoEMG(iTrial,4),1,'last');
        if isempty(tmpIndx),
            indPreEmpty=[indPreEmpty ; iTrial];
        end
    end
    if ~isempty(indPreEmpty)
       tmpInfoEMG(indPreEmpty,:)=[]; 
       NemgTrials=size(tmpInfoEMG,1);
    end
    %-------------------------------------------------------------------------
    for iTrial=1:NemgTrials,
        tmpIndx=find(motorTrialInfo(:,4)<tmpInfoEMG(iTrial,4),1,'last');
        tmpInfoEMG(iTrial,3)=motorTrialInfo(tmpIndx,3);
        tmpInfoEMG(iTrial,6)=(1/Fsample)*(tmpInfoEMG(iTrial,4) -motorTrialInfo(tmpIndx,4));
    end
    %-------------------------------------------------------------------------
    
end
motorTrialInfoEMG=tmpInfoEMG;

%end
%==================================================
%====================================================
%==========================================================
%================================================================
%% FORM OUTPUT TRL
if strcmp(cutMode,'blocks')
    trl=[allBlockInfo(:,3)-prestimSamples allBlockInfo(:,15)+poststimSamples -repmat(prestimSamples,Nblocks,1)  allBlockInfo];
    [i1,j1]=find(trl(:,2)>Nsamples);
    trl(i1,:)=Nsamples;
    psfixtrl=[];
    titletxt = ['cutMode:blocks'];
    
    trlInfoColDescr=motorBlockColumnDescription;
    
elseif strcmp(cutMode,'trials')
    
    titletxt = ['cutMode:trials'];
    if strcmpi(TrigBasedOnEMG, 'yes')
        if ~isempty(motorTrialInfoEMG)
            trl=[motorTrialInfoEMG(:,4)-prestimSamples  motorTrialInfoEMG(:,4)+poststimSamples  -repmat(prestimSamples,size(motorTrialInfoEMG,1),1)  motorTrialInfoEMG];
            psfixtrl=[fixTrialInfo(:,4)  fixTrialInfo(:,4)+ISIsamps  repmat(0,size(fixTrialInfo,1),1) fixTrialInfo nan(size(fixTrialInfo,1),1)];
            titletxt=[titletxt,'  :  ','lockedon:EMG'];
        else
            trl=[];
            psfixtrl=[];
            titletxt=[titletxt,'  :  ','lockedon:EMG no trials found'];
        end
    elseif strcmpi(TrigBasedOnEMG, 'no')
        trl=[motorTrialInfo(:,4)-prestimSamples  motorTrialInfo(:,4)+poststimSamples  -repmat(prestimSamples,size(motorTrialInfo,1),1) motorTrialInfo];
        psfixtrl=[fixTrialInfo(:,4)  fixTrialInfo(:,4)+ISIsamps  repmat(0,size(fixTrialInfo,1),1) fixTrialInfo];
        titletxt=[titletxt,'  :  ','lockedon:Flash'];
    end
    if ~isempty(trl)
        [i1,j1]=find(trl(:,2)>Nsamples);
        trl(i1,:)=[];
    end
    if ~isempty(psfixtrl)
        [i1,j1]=find(psfixtrl(:,2)>Nsamples);
        psfixtrl(i1,:)=[];
    end
    
    
    if strcmpi(TrigBasedOnEMG, 'yes'),
        trlInfoColDescr=motorTrlColumnDescriptionEMG;
    elseif strcmpi(TrigBasedOnEMG, 'no'),
        trlInfoColDescr=motorTrlColumnDescription;
    end
end




%=========================================================
%% -- MERGE TRL with FIXATION PSEUDO TRL (makes a difference only when data is cut in trials)
trl=[trl; psfixtrl];


trialSummary=[];
%===========================================
%% --- Compute summary of trial information to print is summary file
trialSummary.Nblocks=size(allBlockInfo,1);
trialSummary.Nblocks_Fixation=length(find(allBlockInfo(:,2)==6));
trialSummary.Nblocks_LH=length(find(allBlockInfo(:,2)==1));
trialSummary.Nblocks_LF=length(find(allBlockInfo(:,2)==2));
trialSummary.Nblocks_RH=length(find(allBlockInfo(:,2)==4));
trialSummary.Nblocks_RF=length(find(allBlockInfo(:,2)==5));
trialSummary.Ntrials=size(motorTrialInfo,1);
trialSummary.Ntrials_LH=length(find(motorTrialInfo(:,2)==1));
trialSummary.Ntrials_LF=length(find(motorTrialInfo(:,2)==2));
trialSummary.Ntrials_RH=length(find(motorTrialInfo(:,2)==4));
trialSummary.Ntrials_RF=length(find(motorTrialInfo(:,2)==5));
if ~isempty(motorTrialInfoEMG)
    trialSummary.NtrialsEMG_LH=length(find(motorTrialInfoEMG(:,2)==1));
    trialSummary.NtrialsEMG_LF=length(find(motorTrialInfoEMG(:,2)==2));
    trialSummary.NtrialsEMG_RH=length(find(motorTrialInfoEMG(:,2)==4));
    trialSummary.NtrialsEMG_RF=length(find(motorTrialInfoEMG(:,2)==5));
    trialSummary.RatioTrlEMG_LH=trialSummary.NtrialsEMG_LH./trialSummary.Ntrials_LH;
    trialSummary.RatioTrlEMG_LF=trialSummary.NtrialsEMG_LF./trialSummary.Ntrials_LF;
    trialSummary.RatioTrlEMG_RH=trialSummary.NtrialsEMG_RH./trialSummary.Ntrials_RH;
    trialSummary.RatioTrlEMG_RF=trialSummary.NtrialsEMG_RF./trialSummary.Ntrials_RF;
else
    trialSummary.NtrialsEMG_LH=0;
    trialSummary.NtrialsEMG_LF=0;
    trialSummary.NtrialsEMG_RH=0;
    trialSummary.NtrialsEMG_RF=0;
    trialSummary.RatioTrlEMG_LH=0;
    trialSummary.RatioTrlEMG_LF=0;
    trialSummary.RatioTrlEMG_RH=0;
    trialSummary.RatioTrlEMG_RF=0;
end
%===========================================

end
%%  END OF MAIN FUNCTION

%%=========================================================================
%%=========================================================================
%%=========================================================================
%%
function [ emg_trg, Trg_EMG_delay, emgflt] = Extract_EMG_Trigger( emg , Fs , trg , PLT_FLG1)
EMG_STD_THRESHOLD = 0; % Thereshold for z-transform
if nargin < 4
    PLT_FLG1 = 0;% falg for end results plot
end
PLT_FLG2 = 0; % falg for very detailed plots

emg = emg(:)';
mx = max( trg );
mn = min( trg );
ay = abs(emg);
trg = mean( ay( ay > 5*mean(ay) ) ) * (trg - mn )/(mx - mn);

emgflt = ft_preproc_highpassfilter(emg, Fs, 10);       % highpassfilter
emghlb = abs(hilbert(emgflt')');                       % hilbert transform
emgcnv = conv(emghlb , ones(1,round(Fs/5)), 'same');   % smooth using convolution
emgstd = (emgcnv - mean(emgcnv))/std(emgcnv);          % z-transform, i.e. mean=0 and stdev=1
emgtrl = emgstd > EMG_STD_THRESHOLD;                   % detect the muscle activity

ixtrg  = find( diff([trg(1);trg(:)]) > 0 );
ixemgp = find( diff([emgtrl(1);emgtrl(:)]) > 0 );
ixemgn = find( diff([emgtrl(1);emgtrl(:)]) < 0 );

dt = mean( diff(ixtrg) );
emg_trg = zeros( size(trg) );
emgstd2 = emgstd;
Trg_EMG_delay = [];
for i = 1 : length(ixtrg)
    ixtrg1 = ixtrg(i) + round( -0.38*dt:0.58*dt );
    [mvtrg mitrg1] = max( emgstd2(ixtrg1) );
    if mvtrg > EMG_STD_THRESHOLD
        ix1megU = ixemgp(max( find( ixemgp < ( mitrg1 + ixtrg1(1) ) ) ));% Trigger based on EMG starting time
        if ~isempty(ix1megU)
            ix1megD = ixemgn(min( find( ixemgn > ix1megU ) ));
            if ~isempty(ix1megD)
                %========= Trigger based on EMG pick time =================================
                %               tmp = emgstd(ix1megU:ix1megD);
                %               tmpix = find( tmp > .75*max(tmp) );
                %               tmpix = round( sum(tmpix .* tmp(tmpix)) / sum( tmp(tmpix) ) );
                %               ix1megU = ix1megU + tmpix - 1;
                %==========================================================================
                %========= Trigger based on thresold of median of EMG trial ===============
                tmp = abs(emgflt(ix1megU:ix1megD));
                mtmp = mean( tmp( tmp >= median(tmp) ) );
                tmpix = find( tmp >= mtmp );
                ix1megU = ix1megU + min(tmpix) - 1;
                %==========================================================================
                emg_trg(ix1megU:ix1megD)=1;
                emgstd2(ix1megU:ix1megD)=0;
                emgstd2(1:ix1megD)=0;
                Trg_EMG_delay(i) = ix1megU - ixtrg(i);
            else
                Trg_EMG_delay(i) = NaN;
            end
        else
            Trg_EMG_delay(i) = NaN;
        end
    else
        Trg_EMG_delay(i) = NaN;
    end
    % plot results for each trial
    if PLT_FLG2
        figure;clf;%if i== 1, clf, end
        plot( emgstd , ':' )
        hold on
        plot( emgstd2 )
        plot( emgtrl ,'g')
        plot( emg_trg ,'r')
        plot(2*trg/max(trg),'k')
        plot(ixtrg1,2*trg(ixtrg1)/max(trg),'c')
        plot( abs(emgflt)/max(abs(emgflt)),'b' )
        grid on
    end
end

% plot results for whole block
if PLT_FLG1
    figure;clf
    plot( 1.0*abs(emgflt)/max(emgflt)*max(emgstd), 'g' );
    hold on;grid on
    plot( emgstd);
    plot( 1.5*trg/max(trg) , 'r' )
    plot( emg_trg , 'k','linewidth' , 2)
    legend('abs of EMG' , 'EMG after HighPass, Hilbert, Smooth & Z-Trans' , 'Flash on', 'EMG trigger')
    xlim( [ 0 length(trg)] )
    
    %     figure;clf
    %     hist(emgstd,1000)
    %     hold on
    %     plot( [ mean(emgstd) mean(emgstd)], [0 150] , '-r' , 'linewidth' , 2 )
    %     title( 'Histogram of z-transform' )
    
    %     figure;clf
    %     subplot(311)
    %     plot( emg ); hold on;
    %     plot( .75*trg , 'r')
    %     legend('EMG' , 'Flash on')
    %     subplot(312)
    %     plot( emghlb);hold on;
    %     plot( .75*trg , 'r')
    %     legend('Hilbert EMG' , 'Flash on')
    %     subplot(313)
    %     plot( emgcnv); hold on;
    %     plot( .25*max(emgcnv)* trg/max(trg) , 'r')
    %     legend('Hilbert EMG' , 'Flash on')
    
end
end


