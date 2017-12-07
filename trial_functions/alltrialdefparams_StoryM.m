function[allTrlCfgs]=alltrialdefparams_StoryM()
%lockedOn = 'TEV': All Events - Onset of Story Sentences, Math number words, Math operand words, Option intro word,Option 1, Option OR, Option 2  (Fixed trial length)
%           'TRESP': Response (Fixed trial length)
%           'BU': Units (Stories or Math Problems including Option Interval) (Variable trial length)

%------------------------
%'TEV': 
countCases=1;
tmpCfg=[];
tmpCfg.badsegmode='remfull'; % This is used by tmeg_preproc to know if trials should be fully removed, or bad segments should be replaced by nans(i.e. 'BU' case)
tmpCfg.trialdef.lockMnem ='TEV';
tmpCfg.trialdef.cutmode = 1;
tmpCfg.trialdef.prestimTime = 1.5;
tmpCfg.trialdef.poststimTime = 4;
allTrlCfgs{countCases}=tmpCfg;
%------------------------
%------------------------
%'TRESP'
countCases=countCases+1;
tmpCfg=[];
tmpCfg.badsegmode='remfull'; % This is used by tmeg_preproc to know if trials should be fully removed, or bad segments should be replaced by nans(i.e. 'BU' case)
tmpCfg.trialdef.lockMnem ='TRESP';
tmpCfg.trialdef.cutmode = 2;
tmpCfg.trialdef.prestimTime = 1.5;
tmpCfg.trialdef.poststimTime = 1.5;
allTrlCfgs{countCases}=tmpCfg;
%------------------------
%------------------------
%'BSENT'
countCases=countCases+1;
tmpCfg=[];
tmpCfg.badsegmode='repnan'; % This is used by tmeg_preproc to know if trials should be fully removed, or bad segments should be replaced by nans(i.e. 'BU' case)
tmpCfg.trialdef.lockMnem ='BSENT';
tmpCfg.trialdef.cutmode = 3;
tmpCfg.trialdef.prestimTime = 1;
tmpCfg.trialdef.poststimTime =1;
allTrlCfgs{countCases}=tmpCfg;
%------------------------
%------------------------
%'BU'
countCases=countCases+1;
tmpCfg=[];
tmpCfg.badsegmode='repnan'; % This is used by tmeg_preproc to know if trials should be fully removed, or bad segments should be replaced by nans(i.e. 'BU' case)
tmpCfg.trialdef.lockMnem ='BU';
tmpCfg.trialdef.cutmode = 4;
tmpCfg.trialdef.prestimTime = 1.5;
tmpCfg.trialdef.poststimTime =1.5;
allTrlCfgs{countCases}=tmpCfg;
%------------------------
