function[allTrlCfgs]=alltrialdefparams_Motort()

%lockMnem = 'TFLA': FLASH DRIVER ONSET , 'TEMG': EMG CHANNEL ONSET
allTrlCfgs=[];
%------------------------
%'TEMG': EMG CHANNEL ONSET
countCases=1;
tmpCfg=[];
tmpCfg.badsegmode='remfull'; % This is used by tmeg_preproc to know if trials should be fully removed, or bad segments should be replaced by nans(i.e. If one trial spans an entire block case)
tmpCfg.trialdef.lockMnem ='TEMG';
tmpCfg.trialdef.TrigBasedOnEMG = 'yes';
tmpCfg.trialdef.cutMode = 'trials';
tmpCfg.trialdef.prestimTime = 1.2;
tmpCfg.trialdef.poststimTime = 1.2;
tmpCfg.trialdef.montage=[]; % This is the montage containing the emg channels.To be loaded on the fly.  In .labelnew there should {'EMG_LH','EMG_LF','EMG_RH','EMG_RF'}
allTrlCfgs{countCases}=tmpCfg;
%------------------------
%'TFLA': FLASH DRIVER ONSET
countCases=countCases+1;
tmpCfg=[];
tmpCfg.badsegmode='remfull'; % This is used by tmeg_preproc to know if trials should be fully removed, or bad segments should be replaced by nans(i.e. If one trial spans an entire block case)
tmpCfg.trialdef.lockMnem ='TFLA';
tmpCfg.trialdef.TrigBasedOnEMG = 'no';
tmpCfg.trialdef.cutMode = 'trials';
tmpCfg.trialdef.prestimTime = 1.2;
tmpCfg.trialdef.poststimTime = 1.2;
tmpCfg.trialdef.montage=[]; % This is the montage containing the emg channels.To be loaded on the fly.  In .labelnew there should {'EMG_LH','EMG_LF','EMG_RH','EMG_RF'}
allTrlCfgs{countCases}=tmpCfg;
%------------------------

