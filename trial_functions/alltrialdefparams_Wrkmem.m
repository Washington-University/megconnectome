function[allTrlCfgs]=alltrialdefparams_Wrkmem()

%lockedOn = 'TIM': TRIALS ON IMAGE ONSET, 'TRESP': TRIALS ON RESPONSE BUTTON ONSET
%------------------------
%'TIM': TRIALS ON IMAGE ONSET
countCases=1;
tmpCfg=[];
tmpCfg.badsegmode='remfull'; % This is used by tmeg_preproc to know if trials should be fully removed, or bad segments should be replaced by nans(i.e. If one trial spans an entire block case)
tmpCfg.trialdef.lockMnem ='TIM';
tmpCfg.trialdef.lockedOn = 'Image';
tmpCfg.trialdef.cutMode = 'trials';
tmpCfg.trialdef.preStimTime = 1.5;
tmpCfg.trialdef.postStimTime = 2.5;
allTrlCfgs{countCases}=tmpCfg;
%------------------------
%------------------------
%'TRESP': TRIALS ON RESPONSE BUTTON ONSET
countCases=countCases+1;
tmpCfg=[];
tmpCfg.badsegmode='remfull'; % This is used by tmeg_preproc to know if trials should be fully removed, or bad segments should be replaced by nans(i.e. If one trial spans an entire block case)
tmpCfg.trialdef.lockMnem ='TRESP';
tmpCfg.trialdef.lockedOn = 'Response';
tmpCfg.trialdef.cutMode = 'trials';
tmpCfg.trialdef.preStimTime = 1.5;
tmpCfg.trialdef.postStimTime = 2.5;
allTrlCfgs{countCases}=tmpCfg;
%------------------------

