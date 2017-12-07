function[cntrstList]=contrast_Wrkmem_Basic(trlInfo)

% This function produces a default list of cntrsts for the Working Memory
% paradigm. It uses the 'trialinfo' input matrix which contains all the
% trial information for a SINGLE run (FIXME:append runs).
% Each item in the list is a structure with fields
%   .mnemonic  % encoded mnemonic
%   .description= % cell array with detailed description of the cntrst.
%   .selection - indices of cntrst trials from input trialinfo.
%   .operation : type of comparison . can be 'average difference'  or 'ratio'
%
% List of cntrsts and comparisons for working memory
%-------------------------------------------------------------
%--Case 1 - fixation
%--Case 2 - all images , memory, target types (no settings , default)
%--Case 3 - OBack
%--Case 4 - 2Back
%--Case 5 - faces
%--Case 6 - tools
%--Case 7 - targets
%--Case 8 - non targets and lures
%--Case 9 - 0-Back Targets
%--Case 10 - 0-Back Non-Targets and Lures
%--Case 11 - 2-Back Targets
%--Case 12 - 2-Back Non-Targets and Lures
%--Case 13 - 2-Back faces
%--Case 14 - 0-Back faces
%-----   Comparisons
% 15. 0-Back VS 2-Back     'average difference'
% 16. 0-Back VS 2-Back       'average ratio'
% 17. faces VS tools       'average difference'
% 18. faces VS tools            'average ratio'
% 19. targets VS nontargets and lures       'average difference'
% 20. targets VS nontargets and lures       'average ratio'
% 21. 0-Back Targets VS 0-Back Non-Targets and Lures       'average difference'
% 22. 0-Back Targets VS 0-Back Non-Targets and Lures    'average ratio'
% 23. 2-Back Targets VS 2-Back Non-Targets and Lures        'average difference'
% 24.  2-Back Targets VS 2-Back Non-Targets and Lures     'average ratio'
% 25.  2-Back Targets VS 0-Back Targets        'average difference'
% 26.  2-Back Targets VS 0-Back Targets    'average ratio'
% 27.  2-Back faces VS 0-Back faces          'average difference'
% 28.  2-Back faces VS 0-Back faces     'average ratio'

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


%===========================================================================
%{ The main selection of cntrsts is based on the input paramters (which should correspond
%  to columns of the trialinfo matrix created by the corresponding
%  trialun). Within this function at the bottom there is a subfunction
%  called trialsplit_<caseprefix>, which is the trial index selection
%  engine. This function takes the folloing inputs which should be
%  extracted directly from the trialinfo matrix.
%{
function[selection,mnemonic]=trialsplit_Wrkmem(...
    trialinfo,...  % arg. 1
    imageTypes,...  % arg. 2
    memoryTypes,...  % arg. 3
    targetTypes,...  % arg. 4
    hasResponded,...  % arg. 5
    responseTypes,...  % arg. 6
    prevBlockImageTypes,...  % arg. 7
    prevBlockMemoryTypes,...  % arg. 8
    prevTrialTargetTypes,...  % arg. 9
    prevTrialHasResponded,...  % arg. 10
    prevTrialResponseTypes)  % arg. 11
%}



% === The following are the main general division of analysis
%freqBands={'D','TH','A','Blow','Bhigh','Glow','Gmid','Ghigh'};
%pipelines={'eravg','tfavg','srcavglcmv','srcavgdics','conne'};
%conneCases={'coh','plv','imcoh','psi','powc','orthopowc','xpowc','xpowphase','bigranger'};

%== The following are the main subsets of data to apply analysis
%datagroups={'TIM','TRESP'};
%tmpStimCases={'allstim','0B','2B','face','tool','targ','nontarg','lure'};


lockNames=trlInfo.lockNames;
lockTrl=trlInfo.lockTrl;

if length(lockNames)~=length(lockTrl)
    error(['The number of groups does not match the number of trialinfo matrices']);
end

%==============================================================================
%--- LOCKED ON IMAGE ONSET -------------------------------

iTIM=find(strcmp(lockNames,'TIM'));
trialinfo=lockTrl{iTIM};
%=====================================================================================
% INDIVIDUAL CONDITIONS - extracting condition trial indices and mnemonic labels -
%--Case - fixation
[seltim_fix,labtim_fix]         = trialsplit_Wrkmem(trialinfo, 0,[],[],[],[],[],[],[],[],[]);
%--Case - all images , memory, target types (no settings , default)
[seltim_allstim,labtim_allstim] = trialsplit_Wrkmem(trialinfo,[],[],[],[],[],[],[],[],[],[]);
%--Case - OBack
[seltim_0B,labtim_0B]           = trialsplit_Wrkmem(trialinfo,[], 1,[],[],[],[],[],[],[],[]);
%--Case - 2Back
[seltim_2B,labtim_2B]           = trialsplit_Wrkmem(trialinfo,[], 2,[],[],[],[],[],[],[],[]);
%--Case - faces
[seltim_face,labtim_face]       = trialsplit_Wrkmem(trialinfo, 1,[],[],[],[],[],[],[],[],[]);
%--Case - tools
[seltim_tool,labtim_tool]       = trialsplit_Wrkmem(trialinfo, 2,[],[],[],[],[],[],[],[],[]);
%--Case - targets
[seltim_targ,labtim_targ]       = trialsplit_Wrkmem(trialinfo,[],[], 1,[],[],[],[],[],[],[]);
%--Case - non targets
[seltim_nontarg,labtim_nontarg]       = trialsplit_Wrkmem(trialinfo,[],[],[2],[],[],[],[],[],[],[]);
%--Case - lure
[seltim_lure,labtim_lure]       = trialsplit_Wrkmem(trialinfo,[],[],[3],[],[],[],[],[],[],[]);
%--Case - non targets and lures
[seltim_nontargandlure,labtim_nontargandlure]       = trialsplit_Wrkmem(trialinfo,[],[],[2 3],[],[],[],[],[],[],[]);
%--Case - 0-Back Targets
[seltim_0Btarg,labtim_0Btarg]=trialsplit_Wrkmem(trialinfo,[],1,1,[],[],[],[],[],[],[]);
%--Case - 0-Back Non-Targets and Lures
[seltim_0Bnontargandlure,labtim_0Bnontargandlure]=trialsplit_Wrkmem(trialinfo,[],[1],[2 3],[],[],[],[],[],[],[]);
%--Case 11 - 2-Back Targets
[seltim_2Btarg,labtim_2Btarg]=trialsplit_Wrkmem(trialinfo,[],[2],[1],[],[],[],[],[],[],[]);
%--Case 12 - 2-Back Non-Targets and Lures
[seltim_2Bnontargandlure,labtim_2Bnontargandlure]=trialsplit_Wrkmem(trialinfo,[],[2],[2 3],[],[],[],[],[],[],[]);
%--Case 13 - 2-Back faces
[seltim_2Bface,labtim_2Bface]=trialsplit_Wrkmem(trialinfo,[1],[2],[],[],[],[],[],[],[],[]);
%--Case 14 - 0-Back faces
[seltim_0Bface,labtim_0Bface]=trialsplit_Wrkmem(trialinfo,[1],[1],[],[],[],[],[],[],[],[]);
%==============================================================================

% LOCKED ON RESPONSE

iTRESP=find(strcmp(lockNames,'TRESP'));
trialinfo=lockTrl{iTRESP};
%=====================================================================================
% INDIVIDUAL CONDITIONS - extracting condition trial indices and mnemonic labels -
%--Case - fixation
[seltresp_fix,labtresp_fix]         = trialsplit_Wrkmem(trialinfo, 0,[],[],[],[],[],[],[],[],[]);
%--Case - all images , memory, target types (no settings , default)
[seltresp_allstim,labtresp_allstim] = trialsplit_Wrkmem(trialinfo,[],[],[],[],[],[],[],[],[],[]);
%--Case - OBack
[seltresp_0B,labtresp_0B]           = trialsplit_Wrkmem(trialinfo,[], 1,[],[],[],[],[],[],[],[]);
%--Case - 2Back
[seltresp_2B,labtresp_2B]           = trialsplit_Wrkmem(trialinfo,[], 2,[],[],[],[],[],[],[],[]);
%--Case - faces
[seltresp_face,labtresp_face]       = trialsplit_Wrkmem(trialinfo, 1,[],[],[],[],[],[],[],[],[]);
%--Case - tools
[seltresp_tool,labtresp_tool]       = trialsplit_Wrkmem(trialinfo, 2,[],[],[],[],[],[],[],[],[]);
%--Case - targets
[seltresp_targ,labtresp_targ]       = trialsplit_Wrkmem(trialinfo,[],[], 1,[],[],[],[],[],[],[]);
%--Case - non targets
[seltresp_nontarg,labtresp_nontarg]       = trialsplit_Wrkmem(trialinfo,[],[],[2],[],[],[],[],[],[],[]);
%--Case - lure
[seltresp_lure,labtresp_lure]       = trialsplit_Wrkmem(trialinfo,[],[],[3],[],[],[],[],[],[],[]);
%--Case - non targets and lures
[seltresp_nontargandlure,labtresp_nontargandlure]       = trialsplit_Wrkmem(trialinfo,[],[],[2 3],[],[],[],[],[],[],[]);
%--Case - 0-Back Targets
[seltresp_0Btarg,labtresp_0Btarg]=trialsplit_Wrkmem(trialinfo,[],1,1,[],[],[],[],[],[],[]);
%--Case - 0-Back Non-Targets and Lures
[seltresp_0Bnontargandlure,labtresp_0Bnontargandlure]=trialsplit_Wrkmem(trialinfo,[],[1],[2 3],[],[],[],[],[],[],[]);
%--Case 11 - 2-Back Targets
[seltresp_2Btarg,labtresp_2Btarg]=trialsplit_Wrkmem(trialinfo,[],[2],[1],[],[],[],[],[],[],[]);
%--Case 12 - 2-Back Non-Targets and Lures
[seltresp_2Bnontargandlure,labtresp_2Bnontargandlure]=trialsplit_Wrkmem(trialinfo,[],[2],[2 3],[],[],[],[],[],[],[]);
%--Case 13 - 2-Back faces
[seltresp_2Bface,labtresp_2Bface]=trialsplit_Wrkmem(trialinfo,[1],[2],[],[],[],[],[],[],[],[]);
%--Case 14 - 0-Back faces
[seltresp_0Bface,labtresp_0Bface]=trialsplit_Wrkmem(trialinfo,[1],[1],[],[],[],[],[],[],[],[]);
%==============================================================================




%=============================================================================
%=============================================================================
%------- CREATING CONTRASTS ===========================================
% -- Here all the cntrsts to be investigated for this specific paradigm
% -- are hard coded . Any new cntrsts to be added, should added in the
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
LockModes={'TIM'    % Lock modes represent different pools of trial data.
    'TRESP'};
Nlm=length(LockModes);

%--- LockMode  'TIM'
tfavgTimes{1}=[-1:0.025:2]';
srcTimesFine{1}=[-1:0.025:2]';
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
tfavgTimes{2}=[-1.25:0.025:1.75]';
srcTimesFine{2}=[-1.25:0.025:1.75]';
srcTimesCoarseSing{2}=[-0.25 0.25
    0.25 0.75
    0.75 1.25
    1.25 1.75];
srcTimesCoarseComp{2}=[-0.75 -0.25
    -0.25 0.25
    0.25 0.75
    0.75 1.25
    1.25 1.75];
basePeriod{2}=[-0.75 -0.25];
%=======================================================================
%=======================================================================
%=======================================================================
%=======================================================================
tfavgFreqs=[1:100];
freqBands={'D','TH','A','Blow','Bhigh','Glow','Gmid','Ghigh'};
%tmpSingCases={'allstim','0B','2B','face','tool','targ','nontarg','lure'};
%tmpCompCases={ '0B' '2B'
%    'face' 'tool'
%    'targ' 'nontarg'
%    'targ' 'lure'
%    'targ' 'nontarg'};
tmpSingCases={'0B','2B','face','tool'};
tmpCompCases={ '0B' '2B'
    'face' 'tool'};


for iLM=1:Nlm,

       lmstr=lower(LockModes{iLM}); 

    for iCase=1:length(tmpSingCases)
        %-- All Stimuli in Trials cut in n
        %-------------------------------------------------------------
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
            'baseline'     , {basePeriod{iLM}},... % In tfavg the Baseline is used ONLY FOR PLOTTING
            'baselinetype' , 'diff',...            % In tfavg the Baseline is used ONLY FOR PLOTTING 
            'timeperiods' , {tfavgTimes{iLM}},...
            'selection'    , {eval(['sel',lmstr,'_',tmpSingCases{iCase}])},...
            'description'  , {eval(['lab',lmstr,'_',tmpSingCases{iCase}])}...
            );
        cntrst{end}.mnemprint=createcontrmnem(cntrst{end});
        %-------------------------------------------------------------
        cntrst{end+1}=newcntrst(...
            'pipeline'     , 'srcavglcmv',...
            'lockmode'     , {LockModes{iLM}},...
            'mnemtrl'      , tmpSingCases(iCase),...
            'timeperiods' , {srcTimesFine{iLM}},...
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
                'connemetric'  , 'coh',...
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
            'timeperiods'  , { srcTimesFine{iLM}    srcTimesFine{iLM}},...
            'invfiltertype', 'com',...
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
                'timeperiods'  , {srcTimesFine{iLM}, srcTimesFine{iLM}},...
                'invfiltertype', 'com',...
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
                'timeperiods'  , {srcTimesCoarseComp{iLM}, srcTimesCoarseComp{iLM}},...
                'timedef'      , 'concat',...
                'invfiltertype', 'com',...
                'selection'    , {eval(['sel',lmstr,'_',tmpCompCases{iGroup,1}]) eval(['sel',lmstr,'_',tmpCompCases{iGroup,2}])},...
                'description'  , {eval(['lab',lmstr,'_',tmpCompCases{iGroup,1}]) eval(['lab',lmstr,'_',tmpCompCases{iGroup,2}])}...
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
            
            %====================================================================
        end
        
    end
    
    
end
%============================================================
%============================================================
%============================================================
%============================================================================

cntrstList=cntrst;


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

function[selection,description]=trialsplit_Wrkmem(...
    trialinfo,...
    imageTypes,...
    memoryTypes,...
    targetTypes,...
    hasResponded,...
    responseTypes,...
    prevBlockImageTypes,...
    prevBlockMemoryTypes,...
    prevTrialTargetTypes,...
    prevTrialHasResponded,...
    prevTrialResponseTypes)


%{
    
    
% INPUT PARAMETERS
imageTypes=[1 2];            %param A % 1. 'face', 2. 'place', 3.'body', 4.'tools' 0.fixation
memoryTypes=[1 2];                %param B % 1.'0Back' 2.'2Back'
targetTypes=[1 2 3];              %param C % 1: target', 2: nontarget', 3: lure ???'
hasResponded=[0 1];               %param D % 0. No  1. Yes
responseTypes=[0 1 2];            %param E % 0.Mistaken 1. Correct  2.Nan (Non responded or pressed both)
prevBlockImageTypes=[1 2 0];    %param F % 1. 'face', 2. 'place', 3.'body', 4.'tools' (As there is always fixation block prior to image block this referes to the previous image block)
prevBlockMemoryTypes=[1 2];       %param G % 1.'0Back' 2.'2Back'
prevTrialTargetTypes=[1 2 3];     %param H % 1: target', 2: nontarget', 3: lure ???'
prevTrialHasResponded=[0 1];      %param I % 0. No  1. Yes
prevTrialResponseTypes=[0 1 2];   %param J % 0.Mistaken 1. Correct  2.Nan (Non responded or pressed both)
    %}
    
    
    
    %==============================================================
    %========  BELOW ARE THE COLUMNS DESCRIPTION for trialinfo matrix produced
    %           by the trial function.
    %trialInfo matrix contains for each presented image the following
    %======== info (Each row corresponds to an image):
    %===== COLUMNS:
    % 1. Run Number
    % 2. Block Number
    % 3. Image Index in the list of all images. This can be found in variable
    %    'encodedRunImageInfo.allImagesList' in the file
    %    'WM_runStimInfoDataFile.mat'.The actual image can be directly
    %    displayed in matlab from the field 'encodedRunImageInfo.allImagesList.allImageData'
    %
    % 4. imgType :'
    %    '                    1: Face'
    %    '                    2: Tools'
    %
    % 5. memoryType :'
    %    '                    1: 0-Back'
    %    '                    2: 2-Back'
    %
    % 6. targetType           1: target'
    %    '                    2: nontarget'
    %    '                    3: lure ???'
    %
    % 7. Trial start Sample
    % 8. Trial end Sample
    % 9. Parallel Port Code used for this image : In the stimuli presentation code this was computed as
    %        curPPCode=20+25*(imgType-1)+5*(blockType-1); see descr. above for
    %        imgType and blockType
    %
    % 10. isPressed :
    %                  0: user did not press any response button
    %                  1: user pressed a response button
    %
    % 11. isPressedLate:
    %                  1: If subject responded after the 2 seconds that the
    %                  image is at the longest displayed and before the next
    %                  trial
    %                  0: If pressed within the presentation time of the image
    %                  NaN: Otherwise
    %
    % 12. isDoubleResponse
    %                  1: user pressed two response buttons in the same trial
    %                  0: user DID NOT press two response buttons in the same trial
    %
    % 13. pressedCode: Code of the pressed button (If not pressed NaN)
    %
    % 14. isCorrect:
    %                  1: If subject has responded that saw a target when a actual
    %                     target was on or that saw a nontarget when a actual
    %                     nontarget was on
    %                  0:  The opposite of the above.
    %                  NaN:  When subject has not responded or has pressed two
    %                       buttons
    %
    % 15.isLureAsCorrect:
    %                  1: If subject has responded that saw a target when a lure image of actual
    %                     target was on
    %                  0: In all other cases that a subject has responded
    %                  NaN:  When subject has not responded or has pressed two
    %                       buttons
    %
    % 16. respTime: The number of samples from onset of Image to response
    %
    % 17. respDuration: Duration of button press in seconds
    %========================================================================
    %++++++++++++++++++++++++++
    % 18. isFirstInBlock
    % 19. isLastInBlockk
    % -----------
    % 20. prev. Trial: Run Number
    % 21. prev. Trial: Block Number
    % 22. prev. Trial: Image Index in the list of all images. This can be found in variable
    %    'encodedRunImageInfo.allImagesList' in the file
    %    'WM_runStimInfoDataFile.mat'.The actual image can be directly
    %    displayed in matlab from the field 'encodedRunImageInfo.allImagesList.allImageData'
    %
    % 23. prev. Trial: imgType :'
    %    '                    1: Face'
    %    '                    2: Tools'
    %
    % 24. prev. Trial: memoryType :'
    %    '                    1: 0-Back'
    %    '                    2: 2-Back'
    %
    % 25. prev. Trial: targetType           1: target'
    %    '                    2: nontarget'
    %    '                    3: lure ???'
    %
    % 26. prev. Trial: Trial start Sample
    % 27. prev. Trial: Trial end Sample
    % 28. prev. Trial: Parallel Port Code used for this image : In the stimuli presentation code this was computed as
    %        curPPCode=20+25*(imgType-1)+5*(blockType-1); see descr. above for
    %        imgType and blockType
    %
    % 29. prev. Trial: isPressed :
    %                  0: user did not press any response button
    %                  1: user pressed a response button
    %
    % 30. prev. Trial: isPressedLate:
    %                  1: If subject responded after the 2 seconds that the
    %                  image is at the longest displayed and before the next
    %                  trial
    %                  0: If pressed within the presentation time of the image
    %                  NaN: Otherwise
    %
    % 31. prev. Trial: isDoubleResponse
    %                  1: user pressed two response buttons in the same trial
    %                  0: user DID NOT press two response buttons in the same trial
    %
    % 32. prev. Trial: pressedCode: Code of the pressed button (If not pressed NaN)
    %
    % 33. prev. Trial: isCorrect:
    %                  1: If subject has responded that saw a target when a actual
    %                     target was on or that saw a nontarget when a actual
    %                     nontarget was on
    %                  0:  The opposite of the above.
    %                  NaN:  When subject has not responded or has pressed two
    %                       buttons
    %
    % 34. prev. Trial: isLureAsCorrect:
    %                  1: If subject has responded that saw a target when a lure image of actual
    %                     target was on
    %                  0: In all other cases that a subject has responded
    %                  NaN:  When subject has not responded or has pressed two
    %                       buttons
    %
    % 35. prev. Trial: respTime: The number of samples from onset of Image to response
    %
    % 36. prev. Trial: respDuration: Duration of button press in seconds
    % 37. prev. Trial: isFirstInBlock
    % 38. prev. Trial: isLastInBlockk
    % 39. Is button pressed during onset of the stimulus (New field - not in glasgow scans)
    %==========================================================================
    
    
    
    condcfg=[];
    condcfg.imageTypes=[1 2];               %param A % 1. 'face', 2.'tools' 0.fixation
    condcfg.memoryTypes=[1 2];                %param B % 1.'0Back' 2.'2Back'
    condcfg.targetTypes=[1 2 3];              %param C % 1: target', 2: nontarget', 3: lure ???'
    condcfg.hasResponded=[0 1];               %param D % 0. No  1. Yes
    condcfg.responseTypes=[0 1 2];            %param E % 0.Mistaken 1. Correct  2.Nan (Non responded or pressed both)
    condcfg.prevBlockImageTypes=[1 2 0];    %param F % 1. 'face', 2. 'place', 3.'body', 4.'tools' (As there is always fixation block prior to image block this referes to the previous image block)
    condcfg.prevBlockMemoryTypes=[1 2 0];       %param G % 1.'0Back' 2.'2Back'
    condcfg.prevTrialTargetTypes=[1 2 3];     %param H % 1: target', 2: nontarget', 3: lure ???'
    condcfg.prevTrialHasResponded=[0 1];      %param I % 0. No  1. Yes
    condcfg.prevTrialResponseTypes=[0 1 2];   %param J % 0.Mistaken 1. Correct  2.Nan (Non responded or pressed both)
    
    descript=[];
    descript.imageTypes={'fixation','face', 'tool'};             %param A %  1. 'face', 2.'tools' 0.fixation
    descript.memoryTypes={'0Back' '2Back'};                            %param B % 1.'0Back' 2.'2Back'
    descript.targetTypes={'target', 'non-target', 'lure'};              %param C % 1: target', 2: nontarget', 3: lure ???'
    descript.hasResponded={'responded' 'noTRESPonded'};                                %param D % 0. No  1. Ys
    descript.responseTypes={'wrong' 'correct'  'other'};            %param E % 0.Mistaken 1. Correct  2.Nan (Non responded or pressed both)
    descript.prevBlockImageTypes={'face', 'tool','fix'};    %param F % 1. 'face', 2. 'place', 3.'body', 4.'tools'
    descript.prevBlockMemoryTypes={'0Back' '2Back','fix'};                         %param G % 1.'0Back' 2.'2Back'
    descript.prevTrialTargetTypes={'targ', 'nontarg', 'lure'};      %param H % 1: target', 2: nontarget', 3: lure ???'
    descript.prevTrialHasResponded={'responded' 'noTRESPonded'};      %param I % 0. No  1. Yes
    descript.prevTrialResponseTypes={'wrong' 'correct'  'other'};    %param J % 0.Mistaken 1. Correct  2.Nan (Non responded or pressed both)
    
    
    if ismember(condcfg.imageTypes,5)&(numel(condcfg.imageTypes)>1)
        error('fixation trials can only be separately extracted from image types - check the image types ');
    end
    
    
    if ~isempty(imageTypes),    condcfg.imageTypes=imageTypes;  end
    if ~isempty(memoryTypes),    condcfg.memoryTypes=memoryTypes;  end
    if ~isempty(targetTypes),    condcfg.targetTypes=targetTypes;  end
    if ~isempty(hasResponded),    condcfg.hasResponded=hasResponded;  end
    if ~isempty(responseTypes),    condcfg.responseTypes=responseTypes;  end
    if ~isempty(prevBlockImageTypes), condcfg.prevBlockImageTypes=prevBlockImageTypes;  end
    if ~isempty(prevBlockMemoryTypes),    condcfg.prevBlockMemoryTypes=prevBlockMemoryTypes;  end
    if ~isempty(prevTrialTargetTypes),    condcfg.prevTrialTargetTypes=prevTrialTargetTypes;  end
    if ~isempty(prevTrialHasResponded),    condcfg.prevTrialHasResponded=prevTrialHasResponded;  end
    if ~isempty(prevTrialResponseTypes),    condcfg.prevTrialResponseType=prevTrialResponseTypes;  end
    
    isin_imageTypes=ismember(trialinfo(:,4),condcfg.imageTypes);
    isin_memoryTypes=ismember(trialinfo(:,5),condcfg.memoryTypes);
    isin_targetTypes=ismember(trialinfo(:,6),condcfg.targetTypes);
    isin_hasResponded=ismember(trialinfo(:,10),condcfg.hasResponded);
    isin_responseTypes=ismember(trialinfo(:,14),condcfg.responseTypes);
    if sum(ismember(condcfg.responseTypes,2))>=1
        isin_responseTypes=(isin_responseTypes|isnan(trialinfo(:,14)));
    end
    
    
    isin_FirstBlock=(trialinfo(:,2)==1);
    if  isequal(condcfg.prevBlockImageTypes,[1 2 0])&isequal(condcfg.prevBlockMemoryTypes,[1 2 0])
        includeFirstBlock=1;
    else
        includeFirstBlock=0;
    end
    %---
    isin_prevBlockImageTypes=ismember(trialinfo(:,23),condcfg.prevBlockImageTypes);
    if includeFirstBlock,
        isin_prevBlockImageTypes=(isin_prevBlockImageTypes|isin_FirstBlock);
    end
    %---
    isin_prevBlockMemoryTypes=ismember(trialinfo(:,24),condcfg.prevBlockMemoryTypes);
    if includeFirstBlock,
        isin_prevBlockMemoryTypes=(isin_prevBlockMemoryTypes|isin_FirstBlock);
    end
    isin_prevTrialTargetTypes=ismember(trialinfo(:,25),condcfg.prevTrialTargetTypes);
    if  isequal(condcfg.prevTrialTargetTypes,[1 2 3])
        isin_prevTrialTargetTypes(1)=1;
    end
    isin_prevTrialHasResponded=ismember(trialinfo(:,29),condcfg.prevTrialHasResponded);
    if  isequal(condcfg.prevTrialHasResponded,[0 1])
        isin_prevTrialHasResponded(1)=1;
    end
    isin_prevTrialResponseTypes=ismember(trialinfo(:,33),condcfg.prevTrialResponseTypes);
    if  isequal(condcfg.prevTrialResponseTypes,[0 1])
        isin_prevTrialResponseTypes(1)=1;
    end
    if sum(ismember(condcfg.prevTrialResponseTypes,2))>=0
        isin_prevTrialResponseTypes=(isin_prevTrialResponseTypes|isnan(trialinfo(:,33)));
    end
    
    
    if sum(ismember(condcfg.imageTypes,0))>0,% Fixation case
        isInCase=(isin_imageTypes&isin_prevBlockImageTypes&isin_prevBlockMemoryTypes);
        isInCaseCurrent=isin_imageTypes;
    else % other cases
        
        isInCase=(isin_imageTypes&isin_memoryTypes&isin_targetTypes&isin_hasResponded&isin_responseTypes&isin_prevBlockImageTypes&isin_prevBlockMemoryTypes&isin_prevTrialTargetTypes&isin_prevTrialHasResponded&isin_prevTrialResponseTypes);
        isInCaseCurrent=isin_imageTypes&isin_memoryTypes&isin_targetTypes&isin_hasResponded&isin_responseTypes;
    end
    
    selection=find(isInCase);
    
    if  isequal(condcfg.prevBlockImageTypes,[1 2 0])&isequal(condcfg.prevBlockMemoryTypes,[1 2 0])&isequal(condcfg.prevTrialTargetTypes,[1 2 3])&isequal(condcfg.prevTrialHasResponded,[0 1])&isequal(condcfg.prevTrialResponseTypes,[0 1 2])
        includeFirstOfAll=1;
    else
        includeFirstOfAll=0;
    end
    
    if includeFirstOfAll
        if (sum(selection==1)==0)&(isInCaseCurrent(1)==1),
            selection=[1; selection];
        end
    end
    
    
    
    %{
condcfg=[];
condcfg.imageTypes=[1 2 ];                %param A % 1. 'face', 2.'tools'
condcfg.memoryTypes=[1 2];                %param B % 1.'0Back' 2.'2Back'
condcfg.targetTypes=[1 2 3];              %param C % 1: target', 2: nontarget', 3: lure ???'
condcfg.hasResponded=[0 1];               %param D % 0. No  1. Yes
condcfg.responseTypes=[0 1 2];            %param E % 0.Mistaken 1. Correct  2.Nan (Non responded or pressed both)
condcfg.prevBlockImageTypes=[1 2 3 4 ];   %param F % 1. 'face', 2. 'place', 3.'body', 4.'tools'
condcfg.prevBlockMemoryTypes=[1 2];       %param G % 1.'0Back' 2.'2Back'
condcfg.prevTrialTargetTypes=[1 2 3];     %param H % 1: target', 2: nontarget', 3: lure ???'
condcfg.prevTrialHasResponded=[0 1];      %param I % 0. No  1. Yes
condcfg.prevTrialResponseTypes=[0 1 2];   %param J % 0.Mistaken 1. Correct  2.Nan (Non responded or pressed both)
    %}
    
    
    
    description=[];
    if ~isequal(condcfg.imageTypes,[1 2])             %param A % 1. 'face', 2.'tools'
        description=['imageTypes: '];
        for iTmp=1:length(condcfg.imageTypes)
            if iTmp==1, tmpSymbol=''; else tmpSymbol='-and-'; end;
            description=[description,descript.imageTypes{condcfg.imageTypes(iTmp)+1}];
        end
        description=[description,'\n'];
    end
    
    if ~isequal(condcfg.memoryTypes,[1 2])                %param B % 1.'0Back' 2.'2Back'
        description=[description, 'memoryTypes: '];
        for iTmp=1:length(condcfg.memoryTypes)
            if iTmp==1, tmpSymbol=''; else tmpSymbol='-and-'; end;
            description=[description,tmpSymbol,descript.memoryTypes{condcfg.memoryTypes(iTmp)}];
        end
        description=[description,'\n'];
    end
    
    if ~isequal(condcfg.targetTypes,[1 2 3])              %param C % 1: target', 2: nontarget', 3: lure ???'
        description=[description, 'targetTypes: '];
        for iTmp=1:length(condcfg.targetTypes)
            if iTmp==1, tmpSymbol=''; else tmpSymbol='-and-'; end;
            description=[description,tmpSymbol,descript.targetTypes{condcfg.targetTypes(iTmp)}];
        end
        description=[description,'\n'];
    end
    
    if ~isequal(condcfg.hasResponded,[0 1])               %param D % 0. No  1. Yes
        description=[description, 'hasResponded: '];
        for iTmp=1:length(condcfg.hasResponded)
            if iTmp==1, tmpSymbol=''; else tmpSymbol='-and-'; end;
            description=[description,tmpSymbol,descript.hasResponded{condcfg.hasResponded(iTmp)}];
        end
        description=[description,'\n'];
    end
    
    if ~isequal(condcfg.responseTypes,[0 1 2])            %param E % 0.Mistaken 1. Correct  2.Nan (Non responded or pressed both)
        description=[description, 'responseTypes: '];
        for iTmp=1:length(condcfg.responseTypes)
            if iTmp==1, tmpSymbol=''; else tmpSymbol='-and-'; end;
            description=[description,tmpSymbol,descript.responseTypes{condcfg.responseTypes(iTmp)}];
        end
        description=[description,'\n'];
    end
    
    if ~isequal(condcfg.prevBlockImageTypes,[1 2 0])   %param F % 1. 'face', 2.'tools'
        description=['prev. Stim. Block - imageTypes: '];
        for iTmp=1:length(condcfg.prevBlockImageTypes)
            if iTmp==1, tmpSymbol=''; else tmpSymbol='-and-'; end;
            description=[description,tmpSymbol,descript.prevBlockImageTypes{condcfg.prevBlockImageTypes(iTmp)}];
        end
        description=[description,'\n'];
    end
    
    if ~isequal(condcfg.prevBlockMemoryTypes,[1 2 0])       %param G % 1.'0Back' 2.'2Back'
        description=['prev. Stim. Block - memoryTypes: '];
        for iTmp=1:length(condcfg.prevBlockMemoryTypes)
            if iTmp==1, tmpSymbol=''; else tmpSymbol='-and-'; end;
            description=[description,tmpSymbol,descript.prevBlockMemoryTypes{condcfg.prevBlockMemoryTypes(iTmp)}];
        end
        description=[description,'\n'];
    end
    
    if ~isequal(condcfg.prevTrialTargetTypes,[1 2 3])     %param H % 1: target', 2: nontarget', 3: lure ???'
        description=['prev. Trial - targetTypes: '];
        for iTmp=1:length(condcfg.prevTrialTargetTypes)
            if iTmp==1, tmpSymbol=''; else tmpSymbol='-and-'; end;
            description=[description,tmpSymbol,descript.prevTrialTargetTypes{condcfg.prevTrialTargetTypes(iTmp)}];
        end
        description=[description,'\n'];
    end
    
    if ~isequal(condcfg.prevTrialHasResponded,[0 1])      %param I % 0. No  1. Yes
        description=['prev. Trial - hasResponded: '];
        for iTmp=1:length(condcfg.prevTrialHasResponded)
            if iTmp==1, tmpSymbol=''; else tmpSymbol='-and-'; end;
            description=[description,tmpSymbol,descript.prevTrialHasResponded{condcfg.prevTrialHasResponded(iTmp)}];
        end
        description=[description,'\n'];
    end
    
    if ~isequal(condcfg.prevTrialResponseTypes,[0 1 2])
        description=['prev. Trial - responseTypes: '];
        for iTmp=1:length(condcfg.prevTrialResponseTypes)
            if iTmp==1, tmpSymbol=''; else tmpSymbol='-and-'; end;
            description=[description,tmpSymbol,descript.prevTrialResponseTypes{condcfg.prevTrialResponseTypes(iTmp)}];
        end
        description=[description,'\n'];
    end
    
    if isempty(description),
        description='allstimuli';
    end
    
    
    
    % === The following are the main general division of analysis
    %freqBands={'D','TH','A','Blow','Bhigh','Glow','Gmid','Ghigh'};
    %pipelines={'eravg','tfavg','srcavglcmv','srcavgdics','conne'};
    %conneCases={'coh','plv','imcoh','psi','powc','orthopowc','xpowc','xpowphase','bigranger'};
    
    %== The following are the main subsets of data to apply analysis
    %datagroups={'TIM','TRESP'};
    %tmpStimCases={'allstim','0B','2B','face','tool','targ','nontarg','lure'};
    
end
    
    function [contr]=newcntrst(varargin)
        % This subfunction creates a default cntrst function
        
        contr =[];
        contr.pipeline    = ft_getopt(varargin,'pipeline',[]);
        %pipelines={'eravg','tfavg','srcavglcmv','srcavgdics','conne'};
        %----------------------------------------------
        contr.lockmode    = ft_getopt(varargin,'lockmode',{[]});
        %for WM datagroups={'TIM','TRESP'};
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
