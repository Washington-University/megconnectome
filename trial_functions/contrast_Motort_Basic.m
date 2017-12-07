function[contrastList]=contrast_Motort_Basic(trlInfo)

% This function produces a default list of contrasts for the Motor
% paradigm. It uses the 'trialinfo' input matrix which contains all the
% trial information for a SINGLE run (FIXME:append runs).
% Each item in the list is a structure with fields
%   .mnemonic  % encoded mnemonic
%   .description= % cell array with detailed description of the contrast.
%   .selection - indices of contrast trials from input trialinfo.
%   .operation : type of comparison . can be 'average difference'  or 'ratio'

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

lockNames=trlInfo.lockNames;
lockTrl=trlInfo.lockTrl;

if length(lockNames)~=length(lockTrl)
    error(['The number of groups does not match the number of trialinfo matrices'])
end



% INDIVIDUAL CONDITIONS - extracting condition trial indices and mnemonic labels -
%==============================================================================
%--- LOCKED ON IMAGE ONSET -------------------------------

iTFLA=find(strcmp(lockNames,'TFLA'));
trialinfo=lockTrl{iTFLA};

% Case - Fixation
[seltfla_FIX,labtfla_FIX]=trialsplit_Motort(trialinfo, 6,[]);
% Case - Left Hand
[seltfla_LH,labtfla_LH]=trialsplit_Motort(trialinfo, 1,[]);
% Case - Right Foot
[seltfla_LF,labtfla_LF]=trialsplit_Motort(trialinfo, 2,[]);
% Case - Right Hand
[seltfla_RH,labtfla_RH]=trialsplit_Motort(trialinfo, 4,[]);
% Case - Right Foot
[seltfla_RF,labtfla_RF]=trialsplit_Motort(trialinfo, 5,[]);


iTEMG=find(strcmp(lockNames,'TEMG'));
if ~isempty(iTEMG)
    hasEMG=1;
    trialinfo=lockTrl{iTEMG};
    if ~isempty(trialinfo)
        % Case - Fixation
        [seltemg_FIX,labtemg_FIX]=trialsplit_Motort(trialinfo, 6,[]);
        % Case - Left Hand
        [seltemg_LH,labtemg_LH]=trialsplit_Motort(trialinfo, 1,[]);
        % Case - Right Foot
        [seltemg_LF,labtemg_LF]=trialsplit_Motort(trialinfo, 2,[]);
        % Case - Right Hand
        [seltemg_RH,labtemg_RH]=trialsplit_Motort(trialinfo, 4,[]);
        % Case - Right Foot
        [seltemg_RF,labtemg_RF]=trialsplit_Motort(trialinfo, 5,[]);
    else
        hasEMG=0;
        disp(['WARNING!!! - trialinfo for TEMG lock mode is empty - Creating contrasts only for TFLA']);
    end
else
    hasEMG=0;
end

%=============================================================================
%=============================================================================
%------- CREATING CONTRASTS ===========================================
% -- Here all the contrasts to be investigated for this specific paradigm
% -- are hard coded . Any new contrasts to be added, should added in the
% -- sequence below after the individual condition trial and label selection
% -- has been peformed above
%--------------------------------------------------
%-----------------------------


cntrst = [];
%============================================================
%============================================================
%{
 %Main structure of a cntrst

        contr =[];
        contr.pipeline    = [];
        contr.lockmode   = [];
        contr.mnemtrl   = [];
        contr.freqband   =  [];
        contr.freq       =  [];
        contr.operation   = [];
        contr.timeperiods = [];
        contr.timedef     = [];
        contr.baseline    = [];
        contr.baselinetype= [];
        contr.connemetric = [];
        contr.invfiltertype =[];
        contr.selection   = [];
        contr.description = [];
        contr.mnemprint    = [];
%}
%============================================================
if hasEMG
    LockModes={'TFLA'    % Lock modes represent different pools of trial data.
        'TEMG'};
else
    LockModes={'TFLA'};
end
Nlm=length(LockModes);

%--- LockMode  'TFLASH'
tfavgTimes{1}=[-1.2:0.025:1.2]';
srcTimesFine{1}=[fliplr(0:-0.02:-0.8) 0.02:0.02:0.8]';
srcTimesCoarseSing{1}=[0 0.3
    0.3 0.6
    0.6 0.9];
srcTimesCoarseComp{1}=[-0.3 0
    0 0.3
    0.3 0.6
    0.6 0.9];
srcTimesCoarseCompFIX{1}=[0 1.2];
basePeriod{1}=[-0.3 0];

if hasEMG
    %--- LockMode  'TEMG'
    tfavgTimes{2}=[-1.2:0.025:1.2]';
    srcTimesFine{2}=[fliplr(0:-0.02:-0.8) 0.02:0.02:0.8]';
    srcTimesCoarseSing{2}=[0 0.3
        0.3 0.6
        0.6 0.9];
    srcTimesCoarseComp{2}=[-0.3 0
        0 0.3
        0.3 0.6
        0.6 0.9];
    srcTimesCoarseCompFIX{2}=[0 1.2];
    basePeriod{2}=[-0.3 0];
end
%=======================================================================
%=======================================================================
%=======================================================================
%=======================================================================
tfavgFreqs=[1:100];

freqBandsCoarse={'D','TH','A','Blow','Bhigh','Glow','Gmid','Ghigh'};
freqBandsFine={'D1','D2','TH1','TH2','A1','A2','A3','B1','B2','B3','B4','G1','G2','G3','G4','G5','G6','G7'};
tmpSingCases={'LH','RH','LF','RF','FIX'};
tmpCompCases={ 'LH' 'FIX'
    'RH' 'FIX'
    'LF' 'FIX'
    'RF' 'FIX'};

for iLM=1:Nlm,
    
    lmstr=lower(LockModes{iLM});
    
    for iCase=1:length(tmpSingCases)
        %-- All Stimuli in Trials cut in n
        %---------------------------------------------ind----------------
        cntrst{end+1}=newcntrst(...
            'pipeline'     , 'eravg',...
            'lockmode'     ,  {LockModes{iLM}},...
            'mnemtrl'      ,  tmpSingCases(iCase),...
            'baseline'     , {basePeriod{iLM}},...
            'baselinetype' , 'diff',...
            'selection'    , {eval(['sel',lmstr,'_',tmpSingCases{iCase}])},...
            'description'  , {eval(['lab',lmstr,'_',tmpSingCases{iCase}])}...
            );
        cntrst{end}.mnemprint=createcontrmnem(cntrst{end});
        if strcmp(tmpSingCases{iCase},'FIX'),
            cntrst(end)=[];
        end
        %-------------------------------------------------------------
        cntrst{end+1}=newcntrst(...
            'pipeline'     , 'tfavg',...
            'lockmode'     ,  {LockModes{iLM}},...
            'mnemtrl'      ,  tmpSingCases(iCase),...
            'freq'         , tfavgFreqs,...
            'timeperiods'  , {tfavgTimes{iLM}},...
            'baseline'     , {basePeriod{iLM}},... % In tfavg the Baseline is used ONLY FOR PLOTTING
            'baselinetype' , 'diff',...            % In tfavg the Baseline is used ONLY FOR PLOTTING
            'selection'    , {eval(['sel',lmstr,'_',tmpSingCases{iCase}])},...
            'description'  , {eval(['lab',lmstr,'_',tmpSingCases{iCase}])}...
            );
        cntrst{end}.mnemprint=createcontrmnem(cntrst{end});
        if strcmp(tmpSingCases{iCase},'FIX'),
            cntrst(end)=[];
        end        
        %-------------------------------------------------------------
        if hasEMG
            % Quick way of defining the computation of corticomuscular
            % coherence at sensor level. Can be usefull to examine the
            % quality of EMG and to identify in which frequencies there is
            % the highest coherence.
            cntrst{end+1}=newcntrst(...
                'pipeline'     , 'tfavg',...
                'connemetric' , 'emgcoh',... % SPECIAL CASE - in tfavg only coherence with emg is computed to get an idea of the band of coherency
                'lockmode'     ,  {LockModes{iLM}},...
                'mnemtrl'      ,  tmpSingCases(iCase),...
                'freq'         , tfavgFreqs,...
                'timeperiods' ,  {tfavgTimes{iLM}},...
                'selection'    , {eval(['sel',lmstr,'_',tmpSingCases{iCase}])},...
                'description'  , {eval(['lab',lmstr,'_',tmpSingCases{iCase}])}...
                );
            cntrst{end}.mnemprint=createcontrmnem(cntrst{end});
        if strcmp(tmpSingCases{iCase},'FIX'),
            cntrst(end)=[];
        end               
        end
        %-------------------------------------------------------------
        cntrst{end+1}=newcntrst(...
            'pipeline'     , 'srcavglcmv',...
            'lockmode'     , {LockModes{iLM}},...
            'mnemtrl'      , tmpSingCases(iCase),...
            'timeperiods'  , {srcTimesFine{iLM}},...
            'baseline'     , {basePeriod{iLM}},...
            'baselinetype' , 'diff',...
            'invfiltertype', 'avg',...
            'selection'    , {eval(['sel',lmstr,'_',tmpSingCases{iCase}])},...
            'description'  , {eval(['lab',lmstr,'_',tmpSingCases{iCase}])}...
            );
        if strcmp(tmpSingCases(iCase),'FIX'),
            cntrst{end}.baseline={[]};
            cntrst{end}.baselinetype=[];
            cntrst{end}.invfiltertype='all';
            %cntrst{end}.timeperiods={srcTimesCoarseCompFIX{1}(1,:)};
            cntrst{end}.timeperiods={[]};
        end
        
        cntrst{end}.mnemprint=createcontrmnem(cntrst{end});
        % ---For fixation compute also power from pseudo average ERF
        if strcmp(tmpSingCases(iCase),'FIX'),
         cntrst{end+1}=cntrst{end};
         cntrst{end}.invfiltertype='avg';
        end
        cntrst{end}.mnemprint=createcontrmnem(cntrst{end});
        %-------------------------------------------------------------
        %{
        cntrst{end+1}=cntrst{end};
        cntrst{end}.baselinetype= 'relch';
        cntrst{end}.mnemprint=createcontrmnem(cntrst{end});
        %}
        %-------------------------------------------------------------
        
        
        %-------------------------------------------------------------
        cntrst{end+1}=newcntrst(...
            'pipeline'     , 'srcavgdics',...
            'lockmode'     , {LockModes{iLM}},...
            'mnemtrl'      , tmpSingCases(iCase),...
            'timeperiods'  , {srcTimesFine{iLM}},...
            'baseline'     , {basePeriod{iLM}},...
            'baselinetype' , 'diff',...
            'selection'    , {eval(['sel',lmstr,'_',tmpSingCases{iCase}])},...
            'description'  , {eval(['lab',lmstr,'_',tmpSingCases{iCase}])}...
            );
        %'freqband'     , freqBandsCoarse{iBand},...
        if strcmp(tmpSingCases{iCase},'FIX'),        
            cntrst{end}.baseline={[]};
            cntrst{end}.baselinetype=[];
            %cntrst{end}.timeperiods={srcTimesCoarseCompFIX{1}(1,:)};
            cntrst{end}.timeperiods={[]};
            
        end
        cntrst{end}.mnemprint=createcontrmnem(cntrst{end});

        %-------------------------------------------------------------
        if hasEMG
            cntrst{end+1}=newcntrst(...
                'pipeline'     , 'srcavgdics',...
                'connemetric' , 'emgcoh',...
                'lockmode'     , {LockModes{iLM}},...
                'mnemtrl'      , tmpSingCases(iCase),...
                'timeperiods'  , {srcTimesFine{iLM}},...
                'baseline'     , {[]},...
                'baselinetype' , [],...
                'selection'    , {eval(['sel',lmstr,'_',tmpSingCases{iCase}])},...
                'description'  , {eval(['lab',lmstr,'_',tmpSingCases{iCase}])}...
                );
            %'freqband'     , freqBandsCoarse{iBand},...
            cntrst{end}.mnemprint=createcontrmnem(cntrst{end});
        if strcmp(tmpSingCases{iCase},'FIX'),
           cntrst(end)=[];            
        end              
        end
     
        
        
        %-------------------------------------------------------------
        %{
            cntrst{end+1}=cntrst{end};
            cntrst{end}.baselinetype= 'relch';
            cntrst{end}.mnemprint=createcontrmnem(cntrst{end});
        %}
        
%             %-------------------------------------------------------------
%             if hasEMG
%                 cntrst{end+1}=newcntrst(...
%                     'pipeline'     , 'tmegconne',...
%                     'connemetric' , 'emgcoh',...
%                     'lockmode'     , {LockModes{iLM}},...
%                     'mnemtrl'      , tmpSingCases(iCase),...
%                     'freqband'     , freqBandsCoarse{iBand},...
%                     'timeperiods'  , {tfavgTimes{iLM}},...
%                     'timedef'      , 'concat',...
%                     'selection'    , {eval(['sel',lmstr,'_',tmpSingCases{iCase}])},...
%                     'description'  , {eval(['lab',lmstr,'_',tmpSingCases{iCase}])}...
%                     );
%                 cntrst{end}.mnemprint=createcontrmnem(cntrst{end});
%             end
%             %-------------------------------------------------------------
            cntrst{end+1}=newcntrst(...
                'pipeline'     , 'tmegconne',...
                'connemetric' , 'coh',...
                'lockmode'     , {LockModes{iLM}},...
                'mnemtrl'      , tmpSingCases(iCase),...
                'timeperiods'  , {srcTimesCoarseComp{iLM}},...
                'selection'    , {eval(['sel',lmstr,'_',tmpSingCases{iCase}])},...
                'description'  , {eval(['lab',lmstr,'_',tmpSingCases{iCase}])}...
                );
                            %'timedef'      , 'concat',...
                %'baseline'     , {basePeriod{iLM}},...
                %'baselinetype' , 'diff',...
            if strcmp(tmpSingCases{iCase},'FIX'),        
                cntrst{end}.timeperiods={srcTimesCoarseCompFIX{iLM}};
            end
             cntrst{end}.mnemprint=createcontrmnem(cntrst{end});
            %-------------------------------------------------------------
            cntrst{end+1}=cntrst{end};
            cntrst{end}.connemetric= 'imcoh';
            cntrst{end}.mnemprint=createcontrmnem(cntrst{end});
            %-------------------------------------------------------------
            cntrst{end+1}=cntrst{end};
            cntrst{end}.connemetric= 'plv';
            cntrst{end}.mnemprint=createcontrmnem(cntrst{end});
            %-------------------------------------------------------------
            %cntrst{end+1}=cntrst{end};
            %cntrst{end}.connemetric= 'psi';
            %cntrst{end}.mnemprint=createcontrmnem(cntrst{end});
            %-------------------------------------------------------------
            cntrst{end+1}=cntrst{end};
            cntrst{end}.connemetric= 'powc';
            cntrst{end}.mnemprint=createcontrmnem(cntrst{end});
            %-------------------------------------------------------------
            cntrst{end+1}=cntrst{end};
            cntrst{end}.connemetric= 'orthopowc';
            cntrst{end}.mnemprint=createcontrmnem(cntrst{end});
            %-------------------------------------------------------------
        
        
    end
    
    
    
    %============================================================
    %============================================================
    %========CONSTRASTS FOR COMPARISONS BETWEEN CONDITIONS ======
    %============================================================
    
    %for iGroup=1:size(tmpCompCases,1)
        %-------------------------------------------------------------
        % Not sure if it makes sens to compute eravg and favg for comparison
        % between Feet or Hands and Fixation. It makes more sense to
        % reference the power in the souce level to the power during
        % fixation. Now srcavglcmv from these comparison also does not make
        % sense with the comparison being applied at each time point.
        % What makes more sense is to average the power during fixation and
        % then take the average of each condition in time windows and
        % reference it to the fixation average power. Same for srcavgdics
        
        %{
        cntrst{end+1}=newcntrst(...
            'pipeline'     , 'eravg',...
            'operation'    , 'diff',...
            'lockmode'     , {LockModes{iLM} LockModes{iLM}},...
            'mnemtrl'      , tmpCompCases(iGroup,:),...
            'baseline'     , {basePeriod{iLM} basePeriod{iLM}},...
            'baselinetype' , 'diff',...
            'selection'    , {eval(['sel',lmstr,'_',tmpCompCases{iGroup,1}]) eval(['sel',lmstr,'_',tmpCompCases{iGroup,2}])},...
            'description'  , {eval(['lab',lmstr,'_',tmpCompCases{iGroup,1}]) eval(['lab',lmstr,'_',tmpCompCases{iGroup,2}])}...
            );
        cntrst{end}.mnemprint=createcontrmnem(cntrst{end});
        %-------------------------------------------------------------
        cntrst{end+1}=newcntrst(...
            'pipeline'     , 'tfavg',...
            'operation'    , 'diff',...
            'lockmode'     , {LockModes{iLM} LockModes{iLM}},...
            'mnemtrl'      , tmpCompCases(iGroup,:),...
            'freq'         , tfavgFreqs,...
            'timeperiods'  , {tfavgTimes{iLM} tfavgTimes{iLM}},...
            'selection'    , {eval(['sel',lmstr,'_',tmpCompCases{iGroup,1}]) eval(['sel',lmstr,'_',tmpCompCases{iGroup,2}])},...
            'description'  , {eval(['lab',lmstr,'_',tmpCompCases{iGroup,1}]) eval(['lab',lmstr,'_',tmpCompCases{iGroup,2}])}...
            );
        cntrst{end}.mnemprint=createcontrmnem(cntrst{end});
        %}
        %-------------------------------------------------------------
        %{
        cntrst{end+1}=newcntrst(...
            'pipeline'     , 'srcavglcmv',...
            'operation'    , 'diff',...
            'lockmode'     , {LockModes{iLM} LockModes{iLM}},...
            'mnemtrl'      , tmpCompCases(iGroup,:),...
            'timeperiods'  , { srcTimesCoarseComp{iLM}    srcTimesCoarseCompFIX{iLM}},...
            'invfiltertype', 'ind',...
            'selection'    , {eval(['sel',lmstr,'_',tmpCompCases{iGroup,1}]) eval(['sel',lmstr,'_',tmpCompCases{iGroup,2}])},...
            'description'  , {eval(['lab',lmstr,'_',tmpCompCases{iGroup,1}]) eval(['lab',lmstr,'_',tmpCompCases{iGroup,2}])}...
            );
        cntrst{end}.mnemprint=createcontrmnem(cntrst{end});
        %-------------------------------------------------------------
        cntrst{end+1}=cntrst{end};
        cntrst{end}.operation= 'relch';
        cntrst{end}.mnemprint=createcontrmnem(cntrst{end});
        %-------------------------------------------------------------
        %}
        
        %for iBand=1:length(freqBandsCoarse)
            %{
            %-------------------------------------------------------------
            cntrst{end+1}=newcntrst(...
                'pipeline'     , 'srcavgdics',...
                'operation'    , 'diff',...
                'lockmode'     , {LockModes{iLM} LockModes{iLM}},...
                'mnemtrl'      , tmpCompCases(iGroup,:),...
                'freqband'     , freqBandsCoarse{iBand},...
                'timeperiods'  , {srcTimesCoarseComp{iLM}    srcTimesCoarseCompFIX{iLM}},...
                'invfiltertype', 'ind',...
                'selection'    , {eval(['sel',lmstr,'_',tmpCompCases{iGroup,1}]) eval(['sel',lmstr,'_',tmpCompCases{iGroup,2}])},...
                'description'  , {eval(['lab',lmstr,'_',tmpCompCases{iGroup,1}]) eval(['lab',lmstr,'_',tmpCompCases{iGroup,2}])}...
                );
            cntrst{end}.mnemprint=createcontrmnem(cntrst{end});
            %-------------------------------------------------------------
            cntrst{end+1}=cntrst{end};
            cntrst{end}.operation= 'relch';
            cntrst{end}.mnemprint=createcontrmnem(cntrst{end});c
            %}
            %-------------------------------------------------------------
            %{
            cntrst{end+1}=newcntrst(...
                'pipeline'     , 'tmegconne',...
                'connemetric'  ,  'coh',...
                'operation'    , 'diff',...
                'lockmode'     , {LockModes{iLM} LockModes{iLM}},...
                'mnemtrl'      , tmpCompCases(iGroup,:),...
                'freqband'     , freqBandsCoarse{iBand},...
                'timeperiods'  , {srcTimesCoarseComp{iLM}    srcTimesCoarseCompFIX{iLM}},...
                'timedef'      , 'concat',...
                'invfiltertype', 'ind',...
                'selection'    , {eval(['sel',lmstr,'_',tmpCompCases{iGroup,1}]) eval(['sel',lmstr,'_',tmpCompCases{iGroup,2}])},...
                'description'  , {eval(['lab',lmstr,'_',tmpCompCases{iGroup,1}]) eval(['lab',lmstr,'_',tmpCompCases{iGroup,2}])}...
                );
            cntrst{end}.mnemprint=createcontrmnem(cntrst{end});
            %-------------------------------------------------------------
            cntrst{end+1}=cntrst{end};
            cntrst{end}.connemetric= 'imcoh';
            cntrst{end}.mnemprint=createcontrmnem(cntrst{end});
            %-------------------------------com------------------------------
            cntrst{end+1}=cntrst{end};
            cntrst{end}.connemetric= 'plv';
            cntrst{end}.mnemprint=createcontrmnem(cntrst{end});
            %-------------------------------------------------------------
            cntrst{end+1}=cntrst{end};
            cntrst{end}.connemetric= 'psi';
            cntrst{end}.mnemprint=createcontrmnem(cntrst{end});
            %-------------------------------------------------------------
            cntrst{end+1}=cntrst{end};
            cntrst{end}.connemetric= 'powc';
            cntrst{end}.mnemprint=createcontrmnem(cntrst{end});
            %-------------------------------------------------------------
            cntrst{end+1}=cntrst{end};
            cntrst{end}.connemetric= 'orthopowc';
            cntrst{end}.mnemprint=createcontrmnem(cntrst{end});
            %-------------------------------------------------------------
            
            %====================================================================
        end
            
        
    end
            %}
    
    
    
    
end

%============================================================================
%============================================================================
%============================================================================
contrastList=cntrst;
end
%**************************************************************************
% ******** END OF MAIN FUNCTION
%**************************************************************************
%**************************************************************************
%**************************************************************************
%**************************************************************************
%**************************************************************************
%**************************************************************************
%**************************************************************************
%**************************************************************************

function[selection,description]=trialsplit_Motort(trialinfo, memberTypes,prevBlockMemberTypes)
% INPUTS
% trialinfo : matrix with information about each trial
% memberTypes=[1 2 3 4 5 6];            % Param A: 1.Left Hand  2.Left Foot  3.Tongue 4.Right Hand 5.Right Foot 6.Fixation
% prevBlockMemberTypes=[1 2 3 4 5 6];   % Param B: 1.Left Hand  2.Left Foot  3.Tongue 4.Right Hand 5.Right Foot 6.Fixation



%==============================================================
%========  BELOW ARE THE COLUMNS DESCRIPTION for trialinfo matrix produced
%           by the trial function.
%  MOTOR TRIAL INFO
%Columns:
%-------------------
% 1. Block Index
% 2. Block Stim Code
% 3. Trial Index in Block (This is derived by finding the flash cross onset just before the EMG onset)
% 4. Trial Onset EMG Sample
% 5. prev. Block Stim Code
% If trials are cut based on EMG then there is one more extra column
% 6. Time from EMG onset to the previous Flashing Cross




condcfg=[];
condcfg.memberTypes=[1 2 3 4 5 6];            % Param A: 1.Left Hand  2.Left Foot  3.Tongue 4.Right Hand 5.Right Foot 6.Fixation
condcfg.prevBlockMemberTypes=[1 2 3 4 5 6];   % Param B: 1.Left Hand  2.Left Foot  3.Tongue 4.Right Hand 5.Right Foot 6.Fixation


descript=[];
descript.memberTypes          ={ 'Left_Hand' 'Left_Foot' 'Tongue' 'Right Hand' 'Right Foot' 'Fixation'};
descript.prevBlockMemberTypes ={ 'Left_Hand' 'Left_Foot' 'Tongue' 'Right Hand' 'Right Foot' 'Fixation'};

if ~isempty(memberTypes),    condcfg.memberTypes=memberTypes;  end
if ~isempty(prevBlockMemberTypes),    condcfg.prevBlockMemberTypes=prevBlockMemberTypes;  end


isin_memberTypes=ismember(trialinfo(:,2),condcfg.memberTypes);
isin_prevBlockMemberTypes=ismember(trialinfo(:,5),condcfg.prevBlockMemberTypes);

%{
tmpBlockSeq=trialinfo(:,1);
tmpBlockSeqIndx=unique(tmpBlockSeq(ismember(trialinfo(:,2),condcfg.prevBlockMemberTypes),1));
tmpBlockSeqIndx=tmpBlockSeqIndx+1;
isin_prevBlockMemberTypes=ismember(tmpBlockSeq,tmpBlockSeqIndx);

%}
isInCase=(isin_memberTypes&isin_prevBlockMemberTypes);
isInCaseCurrent=isin_memberTypes;
selection=find(isInCase);

if  isequal(condcfg.prevBlockMemberTypes,[1 2 3 4 5 6])
    includeFirstOfAll=1;
else
    includeFirstOfAll=0;
end

if includeFirstOfAll
    if (sum(selection==1)==0)&(isInCaseCurrent(1)==1),
        selection=[1; selection];
    end
end

description=[];

if ~isequal(condcfg.memberTypes,[1 2 3 4 5])
    description=[description, 'memberTypes: '];
    for iTmp=1:length(condcfg.memberTypes)
        if iTmp==1, tmpSymbol=''; else tmpSymbol='-and-'; end;
        description=[description,tmpSymbol,descript.memberTypes{condcfg.memberTypes(iTmp)}];
    end
    description=[description,'\n'];
end


if ~isequal(condcfg.prevBlockMemberTypes,[1 2 3 4 5 6])
    description=[description, 'prevBlockMemberTypes: '];
    for iTmp=1:length(condcfg.prevBlockMemberTypes)
        if iTmp==1, tmpSymbol=''; else tmpSymbol='-and-'; end;
        description=[description,tmpSymbol,descript.prevBlockMemberTypes{condcfg.prevBlockMemberTypes(iTmp)+1}];
    end
    description=[description,'\n'];
end


if isempty(description),
    description='allstimuli';
end

%==========================================
end



function [contr]=newcntrst(varargin)
% This subfunction creates a default cntrst function

contr =[];
contr.pipeline    = ft_getopt(varargin,'pipeline',[]);
%pipelines={'eravg','tfavg','srcavglcmv','srcavgdics','conne'};
%----------------------------------------------
contr.lockmode    = ft_getopt(varargin,'lockmode',{[]});
%for WM datagroups={'TRLIM','TRLRESP'};
%----------------------------------------------
contr.mnemtrl     = ft_getopt(varargin,'mnemtrl', {[]});
% Here go mnemonics for each trial selection i.e. {'0B'}
%----------------------------------------------
contr.freqband    = ft_getopt(varargin,'freqband',[]);
% freqBandsCoarse={'D','TH','A','Blow','Bhigh','Glow','Gmid','Ghigh'};
%----------------------------------------------
contr.freq       =  ft_getopt(varargin,'freq',[]);
% frequencies to be analysed . Currently only used for tfavg pipeline;
%----------------------------------------------
contr.operation   = ft_getopt(varargin,'operation',[]);
% 'diff','rel', or 'relch' This is used when 2 conditions are compared
%----------------------------------------------
contr.timeperiods = ft_getopt(varargin,'timeperiods',{[]});
% This defined the times to be used. If 'all' then the result is
% computed for each time point. For each condition the rows indicate
% The time points for which the analysis is
% performed. If a second column is also
% present then the first column is the
% lower and the second the upper time
% limits within which the analysis should
% be performed. i.e. {[-1 -0.5
%                       -0.5 0
%                       0 0.5]}
% Means that the analysis should be
% performed for these 3 time windows
% How data is integrated in each window is
% defined by the .timedef field of the
% cntrst
%----------------------------------------------
contr.timedef     =  ft_getopt(varargin,'timedef',[]);
% This how data with a time window should be integrated
% It can be 'avg' or 'concat'. If 'avg' then
% the data within a window is averaged
% before computing the result. If 'concat'
% then all the points within the window are
% used for the computation
%----------------------------------------------
contr.baseline    =  ft_getopt(varargin,'baseline',{[]});
% This defines the time period to be used as baseline
%----------------------------------------------
contr.baselinetype= ft_getopt(varargin,'baselinetype',[]);
% It can be 'diff','rel' or 'relch' and defines how
% the baseline will be used on the rest of the data
% When 2 conditions are compared this
% defines how baseline has been used in
% each condition. In the case of 2
% conditions in source space this defines
% how baseline has been used in sensor
% space to derive the inverse solution.
%----------------------------------------------
contr.connemetric =  ft_getopt(varargin,'connemetric',[]);
% conneCases={'coh','plv','imcoh','psi','powc','orthopowc'}; % Next to be implemented ,'xpowc','xpowphase','bigranger'
%----------------------------------------------
contr.invfiltertype = ft_getopt(varargin,'invfiltertype',[]);
% It can be 'com' or 'ind'. This is used when 2 conditions are compared in source space. If 'com' then a common filter
% is derived from both conditions. .If 'ind' the a filter is derived for each condition
%----------------------------------------------
contr.selection   =  ft_getopt(varargin,'selection',[]);
% This is a cell with the indices of the trials to be used
%----------------------------------------------
contr.description =  ft_getopt(varargin,'description',[]);
% This is a cell with basic description of the cntrst (needs to be updated)
%----------------------------------------------
contr.mnemprint    = ft_getopt(varargin,'mnemprint',[]);
% This is the mnemonic that will be used to distinguish the cntrst. Used for saving
%======================== ========================================================================
end


function[printMnem]= createcontrmnem(incontr)

pflags=[];
%pflags.pipeline      = 'PI-';
pflags.lockmode      = 'LM-';
pflags.freqband      = 'FB-';
pflags.operation     = 'OP-';
pflags.timedef       = 'TD-';
pflags.baselinetype  = 'BT-';
pflags.connemetric   = 'CM-';
pflags.invfiltertype  ='IT-';


printMnem=[];
if isempty(incontr.pipeline)
    error('Pipeline must be defined in a constrast')
    return;
end

printMnem=[printMnem,incontr.pipeline];

Nlms=length(incontr.lockmode);

if Nlms==1
    printMnem=[printMnem,'_[',pflags.lockmode,incontr.lockmode{1}];
    printMnem=[printMnem,'-',incontr.mnemtrl{1},']']; %Here a hyphen is used instead of underscore. The next underscore shoudl be after the trial mnemonic
    if ~isempty(incontr.freqband)
        printMnem=[printMnem,'_[',pflags.freqband,incontr.freqband,']'];
    end
    if ~isempty(incontr.operation)
        printMnem=[printMnem,'_[',pflags.operation,incontr.operation,']'];
    end
    if ~isempty(incontr.timedef)
        printMnem=[printMnem,'_[',pflags.timedef,incontr.timedef,']'];
    end
    if ~isempty(incontr.baselinetype)
        if (~strcmp(incontr.pipeline,'tfavg'))&(~strcmp(incontr.pipeline,'srcavglcmv'))&(~strcmp(incontr.pipeline,'srcavgdics'))
            printMnem=[printMnem,'_[',pflags.baselinetype,incontr.baselinetype,']'];
        end
    end
    if ~isempty(incontr.connemetric)
        printMnem=[printMnem,'_[',pflags.connemetric,incontr.connemetric,']'];
    end
    if ~isempty(incontr.invfiltertype)
        printMnem=[printMnem,'_[',pflags.invfiltertype,incontr.invfiltertype,']'];
    end
elseif Nlms==2
    printMnem=[printMnem,'_[',pflags.lockmode,incontr.lockmode{1}];
    printMnem=[printMnem,'-',incontr.mnemtrl{1}]; %Here a hyphen is used instead of underscore. The next underscore shoudl be after the trial mnemonic
    printMnem=[printMnem,'-versus'];
    %printMnem=[printMnem,pflags.lockmode,incontr.lockmode{2}]; This was removed as
    %only contrasts from one
    %datagroup are
    %supported
    printMnem=[printMnem,'-',incontr.mnemtrl{2},']']; %Here a hyphen is used instead of underscore. The next underscore shoudl be after the trial mnemonic
    if ~isempty(incontr.freqband)
        printMnem=[printMnem,'_[',pflags.freqband,incontr.freqband,']'];
    end
    if ~isempty(incontr.operation)
        printMnem=[printMnem,'_[',pflags.operation,incontr.operation,']'];
    end
    if ~isempty(incontr.timedef)
        printMnem=[printMnem,'_[',pflags.timedef,incontr.timedef,']'];
    end
    if ~isempty(incontr.baselinetype)
        if ~strcmp(incontr.pipeline,'tfavg')
            printMnem=[printMnem,'_[',pflags.baselinetype,incontr.baselinetype,']'];
        end
    end
    if ~isempty(incontr.invfiltertype)
        printMnem=[printMnem,'_[',pflags.invfiltertype,incontr.invfiltertype,']'];
    end
    if ~isempty(incontr.connemetric)
        printMnem=[printMnem,'_[',pflags.connemetric,incontr.connemetric,']'];
    end
    
else
    error('1 or 2 different data Lock Modes are supported');
    return;
end
%=================================================================

end



