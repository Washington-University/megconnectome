function[contrastList]=contrast_StoryM_Basic(trlInfo)

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
%--- LOCKED ON EVENTS - ONSET OF SENSTENCE, ONSET OF NUMBER in question or options, ONSET OF calculation WORD , ONSET OF OPTION WORD -------------------------------
% Fixed size trials

iTEV=find(strcmp(lockNames,'TEV'));
trialinfo=lockTrl{iTEV};
datagroupflag=1;

% INDIVIDUAL CONDITIONS - extracting condition trial indices and mnemonic labels -
%--Case - Onset of Sentence in Story
%- c21== 10
[seltev_storsentnon,labtev_storsentnon]  = trialsplit_StoryM(trialinfo,datagroupflag, [10] , [] , []);
%--Case - Onset of Sentence in Math
%- c21 ==20 and choose the first one in the unit
[seltev_mathsentnon,labtev_mathsentnon]  = trialsplit_StoryM(trialinfo,datagroupflag, [20] , 1 , []);
%--Case - Onset of Numbers in Questions
%- c21 ==20
[seltev_mathnumque,labtev_mathnumque]  = trialsplit_StoryM(trialinfo,datagroupflag, [20] , [] , []);
%--Case - Onset of Late Numbers in Questions (It seems that there always 4 numbers in quetions. If this is the case pick the last 2 as late and first two as early)
%- c21 ==20 and choose the ones for which the subunit number in unit is
%higher than half
[seltev_mathnumquelate,labtev_mathnumquelate]   = trialsplit_StoryM(trialinfo,datagroupflag, [20] , 3 , []);
%--Case - Onset of Early Numbers in Questions (It seems that there always 4 numbers in questions. If this is the case pick the last 2 as late and first two as early)
%- c21 ==20 and choose the ones for which the subunit number in unit is
%lower than half
[seltev_mathnumqueearly,labtev_mathnumqueearly]   = trialsplit_StoryM(trialinfo,datagroupflag, [20] , 2 , []);
%--Case - Onset of Numbers in Options
% -c21==23 or 25
[seltev_mathnumopt,labtev_mathnumopt]   = trialsplit_StoryM(trialinfo,datagroupflag, [23 25] , [] , []);
%--Case - Onset of Correct Numbers in Options
% -c21==23 or 25 and fuse with c14 which shows which option is the correct
% one
[seltev_mathnumoptcor,labtev_mathnumoptcor]   = trialsplit_StoryM(trialinfo,datagroupflag, [23 25] , [] , 1);
%--Case - Onset of Wrong Numbers in Options
% -c21==23 or 25 and fuse with c14 which shows which option is the correct
% one
[seltev_mathnumoptwro,labtev_mathnumoptwro]   = trialsplit_StoryM(trialinfo,datagroupflag, [23 25] , [] , 0);
%--Case - Onset of Correct Story Answer in Options
% -c21==13 or 15 and fuse with c14 which shows which option is the correct
% one
[seltev_storoptcor,labtev_storoptcor]   = trialsplit_StoryM(trialinfo,datagroupflag, [13 15] , [] , 1);
%--Case - Onset of Wrong Story Answer in Options
% -c21==13 or 15 and fuse with c14 which shows which option is the correct
% one
[seltev_storoptwro,labtev_storoptwro]   = trialsplit_StoryM(trialinfo,datagroupflag, [13 15] , [] , 0);
%--Case - Onset of Math operations words in Questions (in order to compare with numbers in question the first number should not be included so that the number of trials match - hmmmm)
%-c21 = 21
[seltev_mathoper,labtev_mathoper]   = trialsplit_StoryM(trialinfo,datagroupflag, 21 , [] , []);




%--- LOCKED ON RESPONSE
% fixed size trials
iTRESP=find(strcmp(lockNames,'TRESP'));
trialinfo=lockTrl{iTRESP};
datagroupflag=2;
%--Case - All responses
%- check c21 to get only trials that there was a response
[seltresp_all,labtresp_all]   = trialsplit_StoryM(trialinfo,datagroupflag, [] , [] , []);
%--Case - responses in stories
%- check c21 to get only trials that there was a response and then get
%- c2==1
[seltresp_stor,labtresp_stor]   = trialsplit_StoryM(trialinfo,datagroupflag, 1 , [] , []);
%--Case - responses in maths
%- check c21 to get only trials that there was a response and then get
%- c2==2
[seltresp_math,labtresp_math]   = trialsplit_StoryM(trialinfo,datagroupflag, 2 , [] , []);


%--- LOCKED ON SENTENCE (STORY SENTENCE or MATH SENTENCE excluding option )
% Variable size trials
iBSENT=find(strcmp(lockNames,'BSENT'));
trialinfo=lockTrl{iBSENT};
%--Case - all story sentences
%-c2==1
[selbsent_stor,labbsent_stor]   = trialsplit_StoryM(trialinfo,datagroupflag, 1 , [] , []);
%--Case - all story sentences late (the late half of sentences in a story)
%-c2==1 and use c8 and c11 in order to select the later half sentences in
%each story.
[selbsent_storlate,labbsent_storlate]    = trialsplit_StoryM(trialinfo,datagroupflag, 1 , 3 , []);
%--Case - all story sentences early (the early half of sentences in a story)
%-c2==1 and use c8 and c11 in order to select the earlier half sentences in
%each story.
[selbsent_storearly,labbsent_storearly]    = trialsplit_StoryM(trialinfo,datagroupflag, 1 , 2 , []);
%--Case - all math sentences
%-c2==2
[selbsent_math,labbsent_math]    = trialsplit_StoryM(trialinfo,datagroupflag, 2 , [] , []);


%--- LOCKED ON UNIT (STORY BLOCK or MATH BLOCK including option )
% Variable size trials
iBU=find(strcmp(lockNames,'BU'));
trialinfo=lockTrl{iBU};
%--Case - Story Units
%-c2==1
[selbu_stor,labbu_stor]    = trialsplit_StoryM(trialinfo,datagroupflag, 1 , [] , []);
%--Case - Math Units
%-c2==2
[selbu_math,labbu_math]    = trialsplit_StoryM(trialinfo,datagroupflag, 2 , [] , []);



%============================================================
%============================================================
%============================================================
%============================================================
%============================================================

LockModes={'TEV'    % Lock modes represent different pools of trial data.
    'TRESP'
    'BSENT'
    'BU'};

Nlm=length(LockModes);

%--- LockMode  'TEV'
tfavgTimes{1}=[-1:0.025:3.5]';
srcTimesFine{1}=[-1:0.025:3.5]';
srcTimesCoarseSing{1}=[0 0.5
    0.5 1
    1 1.5
    1.5 2];
srcTimesCoarseComp{1}=[-0.5 0
    0 0.5
    0.5 1
    1 1.5
    1.5 2];
basePeriod{1}=[-0.5 0];

%--- LockMode  'TRESP'
tfavgTimes{2}=[-1.25:0.025:1.25]';
srcTimesFine{2}=[-1.25:0.025:1.25]';
srcTimesCoarseSing{2}=[-0.25 0.25
    0.25 0.75
    0.75 1.25];
srcTimesCoarseComp{2}=[-0.75 -0.25
    -0.25 0.25
    0.25 0.75
    0.75 1.25];
basePeriod{2}=[-0.75 -0.25];

%--- LockMode  'BSENT'
tfavgTimes{3}=[]; % No contrast will be computed at sensor level
srcTimesFine{3}=[]; % As there is variable trial length no fine time resolution will be computed in source level
srcTimesCoarseSing{3}=[]; % Same as above
srcTimesCoarseComp{3}=[]; % In the comparisons the metrics will be computed from the entire trial length
basePeriod{3}=[]; % No baseline required here

%--- LockMode  'BU'
tfavgTimes{4}=[]; % No contrast will be computed at sensor level
srcTimesFine{4}=[]; % As there is variable trial length no fine time resolution will be computed in source level
srcTimesCoarseSing{4}=[]; % Same as above
srcTimesCoarseComp{4}=[]; % In the comparisons the metrics will be computed from the entire trial length
basePeriod{4}=[]; % No baseline required here
%=======================================================================
%=======================================================================
tfavgFreqs=[1:100];
freqBands={'D','TH','A','Blow','Bhigh','Glow','Gmid','Ghigh'};



%{
SENTENCE IN STORY vs ONSET OF SENTENCE IN MATH
ONSET OF LATE NUMBERS IN QUESTIONS vs ONSET OF EARLY NUMBER IS QUESTIONS
ONSET OF CORRECT NUMBERS IN OPTIONS vs ONSET OF WRONG NUMBERS IN OPTIONS
ONSET OF CORRECT ANSWERS IN STORY OPTIONS vs ONSET OF WRONG STORY ANSWERS
ONSET OF NUMBER WORDS vs ONSET OF OPERAND WORDS IN QUESTIONS

TRESP
ONSET OF ALL REPONSES

BSENT:
LATE STORY SENTENCES VS EARLY STORY SENTENCES

BU:
STORIES vs MATHS



%}
cntrst=[];
%===========CONTRASTS FOR TEV
% SINGLE CASES ====================
%=========================================================================
tmpSingCases={'mathnumque','storsentnon','mathsentnon','mathnumopt','mathoper'};
iLM=1;
lmstr=lower(LockModes{iLM});

for iCase=1:length(tmpSingCases)
    
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
    %-------------------------------------------------------------
    cntrst{end+1}=newcntrst(...
        'pipeline'     , 'tfavg',...
        'lockmode'     ,  {LockModes{iLM}},...
        'mnemtrl'      ,  tmpSingCases(iCase),...
        'freq'         , tfavgFreqs,...
        'timeperiods'  , {tfavgTimes{iLM}},...
        'baseline'     , {basePeriod{iLM}},...  % In tfavg the Baseline is used ONLY FOR PLOTTING
        'baselinetype' , 'diff',...             % In tfavg the Baseline is used ONLY FOR PLOTTING
        'selection'    , {eval(['sel',lmstr,'_',tmpSingCases{iCase}])},...
        'description'  , {eval(['lab',lmstr,'_',tmpSingCases{iCase}])}...
        );
    cntrst{end}.mnemprint=createcontrmnem(cntrst{end});
    %-------------------------------------------------------------
    %-------------------------------------------------------------
    cntrst{end+1}=newcntrst(...
        'pipeline'     , 'srcavglcmv',...
        'lockmode'     , {LockModes{iLM}},...
        'mnemtrl'      , tmpSingCases(iCase),...
        'timeperiods'  , {srcTimesFine{iLM}},...
        'baseline'     , {basePeriod{iLM}},...
        'baselinetype' , 'diff',...
        'selection'    , {eval(['sel',lmstr,'_',tmpSingCases{iCase}])},...
        'description'  , {eval(['lab',lmstr,'_',tmpSingCases{iCase}])}...
        );
    cntrst{end}.mnemprint=createcontrmnem(cntrst{end});
    %-------------------------------------------------------------
    cntrst{end+1}=cntrst{end};
    cntrst{end}.baselinetype= 'relch';
    cntrst{end}.mnemprint=createcontrmnem(cntrst{end});
    %-------------------------------------------------------------
    for iBand=1:length(freqBands)
        
        %-------------------------------------------------------------
        cntrst{end+1}=newcntrst(...
            'pipeline'     , 'srcavgdics',...
            'lockmode'     , {LockModes{iLM}},...
            'mnemtrl'      , tmpSingCases(iCase),...
            'freqband'     , freqBands{iBand},...
            'timeperiods'  , {srcTimesFine{iLM}},...
            'baseline'     , {basePeriod{iLM}},...
            'baselinetype' , 'diff',...
            'selection'    , {eval(['sel',lmstr,'_',tmpSingCases{iCase}])},...
            'description'  , {eval(['lab',lmstr,'_',tmpSingCases{iCase}])}...
            );
        cntrst{end}.mnemprint=createcontrmnem(cntrst{end});
        %-------------------------------------------------------------
        cntrst{end+1}=cntrst{end};
        cntrst{end}.baselinetype= 'relch';
        cntrst{end}.mnemprint=createcontrmnem(cntrst{end});
        %-------------------------------------------------------------
        cntrst{end+1}=newcntrst(...
            'pipeline'     , 'tmegconne',...
            'connemetric ' , 'coh',...
            'lockmode'     , {LockModes{iLM}},...
            'mnemtrl'      , tmpSingCases(iCase),...
            'freqband'     , freqBands{iBand},...
            'timeperiods'  , {srcTimesCoarseSing{iLM}},...
            'timedef'      , 'concat',...
            'baseline'     , {basePeriod{iLM}},...
            'baselinetype' , 'diff',...
            'selection'    , {eval(['sel',lmstr,'_',tmpSingCases{iCase}])},...
            'description'  , {eval(['lab',lmstr,'_',tmpSingCases{iCase}])}...
            );
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
    end
end

%====================================================================
% COMPARISONS
tmpCompCases={ 'storsentnon' 'mathsentnon'
    'mathnumquelate' 'mathnumqueearly'
    'mathnumoptcor' 'mathnumoptwro'
    'storoptcor' 'storoptwro'
    'mathnumque' 'mathoper'};


%============================================================
%============================================================
%========CONSTRASTS FOR COMPARISONS BETWEEN CONDITIONS ======
%============================================================

for iGroup=1:size(tmpCompCases,1)
    %-------------------------------------------------------------
    
    
    
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
    %-------------------------------------------------------------
    cntrst{end+1}=newcntrst(...
        'pipeline'     , 'srcavglcmv',...
        'operation'    , 'diff',...
        'lockmode'     , {LockModes{iLM} LockModes{iLM}},...
        'mnemtrl'      , tmpCompCases(iGroup,:),...
        'timeperiods'  , { srcTimesCoarseComp{iLM}    srcTimesCoarseComp{iLM}},...
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
    
    
    for iBand=1:length(freqBands)
        
        %-------------------------------------------------------------
        cntrst{end+1}=newcntrst(...
            'pipeline'     , 'srcavgdics',...
            'operation'    , 'diff',...
            'lockmode'     , {LockModes{iLM} LockModes{iLM}},...
            'mnemtrl'      , tmpCompCases(iGroup,:),...
            'freqband'     , freqBands{iBand},...
            'timeperiods'  , {srcTimesCoarseComp{iLM}    srcTimesCoarseComp{iLM}},...
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
        cntrst{end+1}=newcntrst(...
            'pipeline'     , 'tmegconne',...
            'connemetric'  ,  'coh',...
            'operation'    , 'diff',...
            'lockmode'     , {LockModes{iLM} LockModes{iLM}},...
            'mnemtrl'      , tmpCompCases(iGroup,:),...
            'freqband'     , freqBands{iBand},...
            'timeperiods'  , {srcTimesCoarseComp{iLM}    srcTimesCoarseComp{iLM}},...
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

%====================================================================
%==================================================================
%====================================================================
%==================================================================
%====================================================================
%==================================================================
%===========CONTRASTS FOR TRESP
% SINGLE CASES ====================
%=========================================================================

iLM=2;
lmstr=lower(LockModes{iLM});
tmpSingCases={'all'};

for iCase=1:length(tmpSingCases)
    
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
    %-------------------------------------------------------------
    cntrst{end+1}=newcntrst(...
        'pipeline'     , 'tfavg',...
        'lockmode'     ,  {LockModes{iLM}},...
        'mnemtrl'      ,  tmpSingCases(iCase),...
        'freq'         , tfavgFreqs,...
        'timeperiods'  , {tfavgTimes{iLM}},...
        'baseline'     , {basePeriod{iLM}},...
        'baselinetype' , 'diff',...
        'selection'    , {eval(['sel',lmstr,'_',tmpSingCases{iCase}])},...
        'description'  , {eval(['lab',lmstr,'_',tmpSingCases{iCase}])}...
        );
    cntrst{end}.mnemprint=createcontrmnem(cntrst{end});
    %-------------------------------------------------------------
    %-------------------------------------------------------------
    cntrst{end+1}=newcntrst(...
        'pipeline'     , 'srcavglcmv',...
        'lockmode'     , {LockModes{iLM}},...
        'mnemtrl'      , tmpSingCases(iCase),...
        'timeperiods'  , {srcTimesFine{iLM}},...
        'baseline'     , {basePeriod{iLM}},...
        'baselinetype' , 'diff',...
        'selection'    , {eval(['sel',lmstr,'_',tmpSingCases{iCase}])},...
        'description'  , {eval(['lab',lmstr,'_',tmpSingCases{iCase}])}...
        );
    cntrst{end}.mnemprint=createcontrmnem(cntrst{end});
    %-------------------------------------------------------------
    cntrst{end+1}=cntrst{end};
    cntrst{end}.baselinetype= 'relch';
    cntrst{end}.mnemprint=createcontrmnem(cntrst{end});
    %-------------------------------------------------------------
    for iBand=1:length(freqBands)
        
        %-------------------------------------------------------------
        cntrst{end+1}=newcntrst(...
            'pipeline'     , 'srcavgdics',...
            'lockmode'     , {LockModes{iLM}},...
            'mnemtrl'      , tmpSingCases(iCase),...
            'freqband'     , freqBands{iBand},...
            'timeperiods'  , {srcTimesFine{iLM}},...
            'baseline'     , {basePeriod{iLM}},...
            'baselinetype' , 'diff',...
            'selection'    , {eval(['sel',lmstr,'_',tmpSingCases{iCase}])},...
            'description'  , {eval(['lab',lmstr,'_',tmpSingCases{iCase}])}...
            );
        cntrst{end}.mnemprint=createcontrmnem(cntrst{end});
        %-------------------------------------------------------------
        cntrst{end+1}=cntrst{end};
        cntrst{end}.baselinetype= 'relch';
        cntrst{end}.mnemprint=createcontrmnem(cntrst{end});
        %-------------------------------------------------------------
        cntrst{end+1}=newcntrst(...
            'pipeline'     , 'tmegconne',...
            'connemetric ' , 'coh',...
            'lockmode'     , {LockModes{iLM}},...
            'mnemtrl'      , tmpSingCases(iCase),...
            'freqband'     , freqBands{iBand},...
            'timeperiods'  , {srcTimesCoarseSing{iLM}},...
            'timedef'      , 'concat',...
            'baseline'     , {basePeriod{iLM}},...
            'baselinetype' , 'diff',...
            'selection'    , {eval(['sel',lmstr,'_',tmpSingCases{iCase}])},...
            'description'  , {eval(['lab',lmstr,'_',tmpSingCases{iCase}])}...
            );
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
    end
end

%====================================================================
%====================================================================
%==================================================================
%====================================================================
%==================================================================
%====================================================================
%==================================================================

%====================================================================
%- datagroup BSENT
%----  No single conditions
%----  only comparisons
iLM=3;
lmstr=lower(LockModes{iLM});

% COMPARISONS
tmpCompCases={ 'storlate' 'storearly'};


%============================================================
%============================================================
%========CONSTRASTS FOR COMPARISONS BETWEEN CONDITIONS ======
%============================================================

for iGroup=1:size(tmpCompCases,1)
    %-------------------------------------------------------------
    %-------------------------------------------------------------
    cntrst{end+1}=newcntrst(...
        'pipeline'     , 'srcavglcmv',...
        'operation'    , 'diff',...
        'lockmode'     , {LockModes{iLM} LockModes{iLM}},...
        'mnemtrl'      , tmpCompCases(iGroup,:),...
        'timeperiods'  , { srcTimesCoarseComp{iLM}    srcTimesCoarseComp{iLM}},...
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
    
    
    for iBand=1:length(freqBands)
        
        %-------------------------------------------------------------
        cntrst{end+1}=newcntrst(...
            'pipeline'     , 'srcavgdics',...
            'operation'    , 'diff',...
            'lockmode'     , {LockModes{iLM} LockModes{iLM}},...
            'mnemtrl'      , tmpCompCases(iGroup,:),...
            'freqband'     , freqBands{iBand},...
            'timeperiods'  , {srcTimesCoarseComp{iLM}    srcTimesCoarseComp{iLM}},...
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
        cntrst{end+1}=newcntrst(...
            'pipeline'     , 'tmegconne',...
            'connemetric'  ,  'coh',...
            'operation'    , 'diff',...
            'lockmode'     , {LockModes{iLM} LockModes{iLM}},...
            'mnemtrl'      , tmpCompCases(iGroup,:),...
            'freqband'     , freqBands{iBand},...
            'timeperiods'  , {srcTimesCoarseComp{iLM}    srcTimesCoarseComp{iLM}},...
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

%====================================================================
%==================================================================
%====================================================================
%==================================================================
%====================================================================
%==================================================================

%====================================================================
%- datagroup BU
%----  No single conditions
%----  only comparisons
iLM=4;
lmstr=lower(LockModes{iLM});

% COMPARISONS
tmpCompCases={'stor' 'math'};


%============================================================
%============================================================
%========CONSTRASTS FOR COMPARISONS BETWEEN CONDITIONS ======
%============================================================

for iGroup=1:size(tmpCompCases,1)
    %-------------------------------------------------------------
    cntrst{end+1}=newcntrst(...
        'pipeline'     , 'srcavglcmv',...
        'operation'    , 'diff',...
        'lockmode'     , {LockModes{iLM} LockModes{iLM}},...
        'mnemtrl'      , tmpCompCases(iGroup,:),...
        'timeperiods'  , { srcTimesCoarseComp{iLM}    srcTimesCoarseComp{iLM}},...
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
    
    
    for iBand=1:length(freqBands)
        
        %-------------------------------------------------------------
        cntrst{end+1}=newcntrst(...
            'pipeline'     , 'srcavgdics',...
            'operation'    , 'diff',...
            'lockmode'     , {LockModes{iLM} LockModes{iLM}},...
            'mnemtrl'      , tmpCompCases(iGroup,:),...
            'freqband'     , freqBands{iBand},...
            'timeperiods'  , {srcTimesCoarseComp{iLM}    srcTimesCoarseComp{iLM}},...
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
        cntrst{end+1}=newcntrst(...
            'pipeline'     , 'tmegconne',...
            'connemetric'  ,  'coh',...
            'operation'    , 'diff',...
            'lockmode'     , {LockModes{iLM} LockModes{iLM}},...
            'mnemtrl'      , tmpCompCases(iGroup,:),...
            'freqband'     , freqBands{iBand},...
            'timeperiods'  , {srcTimesCoarseComp{iLM}    srcTimesCoarseComp{iLM}},...
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

%==================================================================







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

function[selection,description]=trialsplit_StoryM(trialinfo,datagroupflag, evType,orderType,corrOptType)
% INPUTS
% trialinfo : matrix with information about each trial
% datagroupflag : 1: TEV locked on the most detailed events . This is
%                     sentences for the Story , NUmbers and operator words
%                     for maths , option 1 -  FIXED length trials
%                    , OR word, opton 2.
%                 2:  TRESP: locked on the response -  FIXED length trials
%                 3:  BSENT: each trial is a sentence. Either a sentence of a story of a Math question sentence - Variable trial length
%                 4:  BU:   each trials is a unit , meaning one story
%                     (including option) or one math problem(including option)
%
%evType = For TEV Column 21 of trialinfo. See below for description. For
%the rest of datagroups is column 2.
%orderType = 1: The first one,  2: The early ones  3 : The late ones
%corrType =  0: Wrong Options  1: Correct Options
%==============================================================
%========  BELOW ARE THE COLUMNS DESCRIPTIONS
% TEV
%---------------------------------------------
%detailed trlAllEventInfo
% Columns:
%----------------------------------------------
% 1. Block Number within Run
% 2. Unit Type : 1.Story 2.Math
% 3. Unit Number within Run
% 4. Total Number of units (N of Stories or N of Math problems ) in same Run
% 5. Unit Number within Block
% 6. Total Number of units (N of Stories or N of Math problems ) in same Block
% 7. Attribute1: For story this is the story number. For Math is the
%    difficulty level
% 8. Unit Narration interval Start Sample - Start with the onset of the
%         first word trigger or the beginning of the first sentence
% 9. Unit Narration interval End Sample
% 10. N subunits within Narration interval
% 11. Unit Option interval Start Sample
% 12. Unit Option interval End Sample
% 13. N subunits within Option interval
% 14 Correct Option- 1 or 2
% 15 Unit Response interval start sample
% 16 Unit Response interval end sample
% 17 Unit Response sample
% 18 is Response Correct
% 19 is Response Early
% 20 Event Sample
% 21 Event Type - 20.Math Narration number word, 21.Math Narration operand
%       word, 22.Math Option Intro Word 23.Math Option 1 Word  24. Math Option
%       OR Word 25.Math Option 2 Word
%       10.Story Sentence, 12.Story Option Intro 13.Story Option 1 Word  14. Story Option
%       OR Word 15.Story Option 2 Word
% 22. Narration Event Number in Narration interval  (Applies only to Sentence or number or operation in Narration interval)
%==============================================================
% TRESP:
%-----------------------------
% trlUnitInfo
% Columns:
%----------------------------------------------
% 1. Block Number within Run
% 2. Unit Type : 1.Story 2.Math
% 3. Unit Number within Run
% 4. Total Number of units (N of Stories or N of Math problems ) in same Run
% 5. Unit Number within Block
% 6. Total Number of units (N of Stories or N of Math problems ) in same Block
% 7. Attribute1: For story this is the story number. For Math is the
%    difficulty level
% 8. Unit Narration interval Start Sample - Start with the onset of the
%         first word trigger or the beginning of the first sentence
% 9.  Unit Narration interval End Sample
% 10. N subunits within Narration interval
% 11. Unit Option interval Start Sample
% 12. Unit Option interval End Sample
% 13. N subunits within Option interval
% 14. Option Intro Sample - "equals" or "That was about"
% 15. Option1 onset sample
% 16. OR onset sample
% 17. Option2 onset sample
% 18. Correct Option- 1 or 2
% 19. Unit Response interval start sample
% 20. Unit Response interval end sample
% 21. Unit Response sample
% 22. is Response Correct
% 23. is Response Early
%=====================================================================
% BSENT
% trlSentInfo
% Each trials corresponds to a Narration Sentence , either in story or in
% math (option blocks not included)
% Columns:
%----------------------------------------------
% 1. Block Number within Run
% 2. Unit Type : 1.Story 2.Math
% 3. Unit Number within Run
% 4. Total Number of units (N of Stories or N of Math problems ) in same Run
% 5. Unit Number within Block
% 6. Total Number of units (N of Stories or N of Math problems ) in same Block
% 7. Attribute1: For story this is the story number. For Math is the
%    difficulty level
% 8. N sentences within Narration interval (always one for math as in math one narration interval corresponds to one math sentence)
% 9. is Response Correct (For Story this refers to the response at the very end of the sentence)
% 10. is Response Early
% 11. Narration Sentence Number in Narration interval  (For math always equal to one)
% 12. Narration Sentence Start Sample
% 13. Narration Sentence End Sample
%-----------------------------------------
if datagroupflag==1 %TEV
    condcfg=[];
    condcfg.evTypes=[10 20 13 15 23 25 21];
    condcfg.orderTypes=[1 2 3];
    condcfg.corrTypes=[0 1];
    
    descript=[];
    descript.evTypes={'Story_Sentence_On' , 'Math_Number_InQuestion' ,'Story_Option_1' , 'Story_Option_2' , 'Math_Option_1' , 'Math_Option_2', 'Math_Operand'};
    descript.orderTypes={'First' , 'Early' , 'Late'};
    descript.corrTypes={'Wrong','Early','Late'};
    
    
    tmptrialinfo=trialinfo;
    tmpDescr=[];
    
    
    
    isin_evType=ismember(tmptrialinfo(:,21),evType);
    totIndx=find(isin_evType);
    tmptrialinfo=tmptrialinfo(isin_evType,:);
    if length(evType)==1,
        tmpIndx=ismember(condcfg.evTypes,evType);
        tmpDescr= descript.evTypes{tmpIndx};
    else
        if (evType==[13 15])
            tmpDescr=['Story_Options'];
            if corrOptType==0
                tmpIndx1=find((tmptrialinfo(:,14)==1)&(tmptrialinfo(:,21)==15));
                tmpIndx2=find((tmptrialinfo(:,14)==2)&(tmptrialinfo(:,21)==13));
                tmpIndx=[tmpIndx1 ; tmpIndx2];
                tmpIndx=unique(tmpIndx); % Check if there are duplicate entries. There shouldnt
                totIndx=totIndx(tmpIndx);
                tmptrialinfo=tmptrialinfo(tmpIndx,:);
                tmpDescr=[tmpDescr,'_Wrong'];
            elseif corrOptType==1
                tmpIndx1=find((tmptrialinfo(:,14)==1)&(tmptrialinfo(:,21)==13));
                tmpIndx2=find((tmptrialinfo(:,14)==2)&(tmptrialinfo(:,21)==15));
                tmpIndx=[tmpIndx1 ; tmpIndx2];
                tmpIndx=unique(tmpIndx); % Check if there are duplicate entries. There shouldnt
                totIndx=totIndx(tmpIndx);
                tmptrialinfo=tmptrialinfo(tmpIndx,:);
                tmpDescr=[tmpDescr,'_Correct'];
            end
            
            
        elseif (evType==[23 25])
            tmpDescr=['Math_Options'];
            if corrOptType==0
                tmpIndx1=find((tmptrialinfo(:,14)==1)&(tmptrialinfo(:,21)==25));
                tmpIndx2=find((tmptrialinfo(:,14)==2)&(tmptrialinfo(:,21)==23));
                tmpIndx=[tmpIndx1 ; tmpIndx2];
                tmpIndx=unique(tmpIndx); % Check if there are duplicate entries. There shouldnt
                totIndx=totIndx(tmpIndx);
                tmptrialinfo=tmptrialinfo(tmpIndx,:);
                tmpDescr=[tmpDescr,'_Wrong'];
            elseif corrOptType==1
                tmpIndx1=find((tmptrialinfo(:,14)==1)&(tmptrialinfo(:,21)==23));
                tmpIndx2=find((tmptrialinfo(:,14)==2)&(tmptrialinfo(:,21)==25));
                tmpIndx=[tmpIndx1 ; tmpIndx2];
                tmpIndx=unique(tmpIndx); % Check if there are duplicate entries. There shouldnt
                totIndx=totIndx(tmpIndx);
                tmptrialinfo=tmptrialinfo(tmpIndx,:);
                tmpDescr=[tmpDescr,'_Correct'];
            end
            
            
        elseif (evType==[13 15 23 25]);
            tmpDescr=['Options'];
            if corrOptType==0
                tmpIndx1=find((tmptrialinfo(:,14)==1)&ismember(tmptrialinfo(:,21),[15 25]));
                tmpIndx2=find((tmptrialinfo(:,14)==2)&ismember(tmptrialinfo(:,21),[13 23]));
                tmpIndx=[tmpIndx1 ; tmpIndx2];
                tmpIndx=unique(tmpIndx); % Check if there are duplicate entries. There shouldnt
                totIndx=totIndx(tmpIndx);
                tmptrialinfo=tmptrialinfo(tmpIndx,:);
                tmpDescr=[tmpDescr,'_Wrong'];
            elseif corrOptType==1
                tmpIndx1=find((tmptrialinfo(:,14)==1)&ismember(tmptrialinfo(:,21),[13 23]));
                tmpIndx2=find((tmptrialinfo(:,14)==2)&ismember(tmptrialinfo(:,21),[15 25]));
                tmpIndx=[tmpIndx1 ; tmpIndx2];
                tmpIndx=unique(tmpIndx); % Check if there are duplicate entries. There shouldnt
                totIndx=totIndx(tmpIndx);
                tmptrialinfo=tmptrialinfo(tmpIndx,:);
                tmpDescr=[tmpDescr,'_Correct'];
            end
            
        end
        
    end
    
    if ~isempty(orderType) % It only
        if orderType==1,
            isin_orderType=(tmptrialinfo(:,22)==1);
            tmpDescr=[tmpDescr,'_First'];
        elseif orderType==2,
            isin_orderType=(tmptrialinfo(:,22)<=floor(tmptrialinfo(:,10)./2));   % Early
            tmpDescr=[tmpDescr,'_Early'];
        elseif orderType==3,
            isin_orderType=(tmptrialinfo(:,22)>floor(tmptrialinfo(:,10)./2));   % Late
            tmpDescr=[tmpDescr,'_Late'];
        end
        tmptrialinfo=tmptrialinfo(isin_orderType,:);
        totIndx=totIndx(isin_orderType);
    end
    
    selection=totIndx;
    description=tmpDescr;
    
    
    
elseif datagroupflag==2 %TRESP
    
    tmptrialinfo=trialinfo;
    tmpDescr=[];
    totIndx=[];
    
    if ~isempty(evType)
        isin_evType=ismember(tmptrialinfo(:,2),evType);
        totIndx=find(isin_evType);
        tmptrialinfo=tmptrialinfo(isin_evType,:);
        
        if evType==1,
            tmpDescr='Story';
        elseif evType==2,
            tmpDescr='Math';
        end
        
    else
        totIndx=[1:size(tmptrialinfo,1)];
    end
    %-----
    if isempty(tmpDescr)
        tmpDescr='all';
    end
    
    
    selection=totIndx;
    description=tmpDescr;
    
    
    
elseif datagroupflag==3 %BSENT
    
    tmptrialinfo=trialinfo;
    tmpDescr=[];
    totIndx=[];
    
    if ~isempty(evType)
        isin_evType=ismember(tmptrialinfo(:,2),evType);
        totIndx=find(isin_evType);
        tmptrialinfo=tmptrialinfo(isin_evType,:);
        
        if evType==1,
            tmpDescr='Story';
        elseif evType==2,
            tmpDescr='Math';
        end
        %------------------------------------
        if ~isempty(orderType)
            if orderType==2,  %early
                isin_orderType=(tmptrialinfo(:,11)<=floor(tmptrialinfo(:,8)./2));   % Early
                tmpDescr=[tmpDescr,'_Early'];
            elseif orderType==3, % late
                isin_orderType=(tmptrialinfo(:,11)>floor(tmptrialinfo(:,8)./2));   % Late
                tmpDescr=[tmpDescr,'_Late'];
            end
            tmptrialinfo=tmptrialinfo(isin_orderType,:);
            totIndx=totIndx(isin_orderType);
        end
        %------------------------------------
    else
        totIndx=[1:size(tmptrialinfo,1)];
    end
    %-----
    if isempty(tmpDescr)
        tmpDescr='all';
    end
    
    
    selection=totIndx;
    description=tmpDescr;
    
    
    
elseif datagroupflag==4 %BU
    
    tmptrialinfo=trialinfo;
    tmpDescr=[];
    totIndx=[];
    
    if ~isempty(evType)
        isin_evType=ismember(tmptrialinfo(:,2),evType);
        totIndx=find(isin_evType);
        tmptrialinfo=tmptrialinfo(isin_evType,:);
        
        if evType==1,
            tmpDescr='Story';
        elseif evType==2,
            tmpDescr='Math';
        end
        
    else
        totIndx=[1:size(tmptrialinfo,1)];
    end
    %-----
    if isempty(tmpDescr)
        tmpDescr='all';
    end
    
    
    selection=totIndx;
    description=tmpDescr;
    
end



end
%==========================================
%==========================================
%==========================================
%==========================================
%==========================================
%==========================================



function [contr]=newcntrst(varargin)
% This subfunction creates a default cntrst function

contr =[];
contr.pipeline    = ft_getopt(varargin,'pipeline',[]);
%pipelines={'eravg','tfavg','srcavglcmv','srcavgdics','conne'};
%----------------------------------------------
contr.lockmode    = ft_getopt(varargin,'lockmode',{[]});
%for WM datagroups={'TRLIM','TRLRESP'};
%----------------------------------------------
contr.mnemtrl     = ft_getopt(varargin,'mnemtrl',{[]});
% Here go mnemonics for each trial selection i.e. {'0B'}
%----------------------------------------------
contr.freqband    = ft_getopt(varargin,'freqband',[]);
% freqBands={'D','TH','A','Blow','Bhigh','Glow','Gmid','Ghigh'};
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
%================================================================================================
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
         if ~strcmp(incontr.pipeline,'tfavg')
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
