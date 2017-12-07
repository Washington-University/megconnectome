function    [trl,trlInfoColDescr,trialSummary,scanStartSamp,scanEndSamp,warninfo]= trialfun_Wrkmem_BaseExtractAll( cfg )
%% This is the Trial Definition Function for Working Memomory experiment.
% It extracts the trial definition for ALL trials within each of different
% datagroup.
% There are 2 different data groups for Working Memory
%
% DATA GROUPS
%----------------------------------
% 1. mnemonic: TIM.     description: Onset of an image that the subject has.to match or not with the target image.
% 2. mnemonic: TRESP. 	description: Onset of button press by the subject.
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
%                          .trialdef.lockedOn = 'Image' or 'Response'; % Defines the 0 reference time of eahc trial
%                                                                      %  according to the desired data group
%                          .trialdef.cutMode = 'trials';               % This representes that data will be cut in trials
%                                                                      %  and not in blocks. The same for both data groups
%                          .trialdef.preStimTime = 1.5;                % Time interval prior to 0 reference point for each trial 
%                          .trialdef.postStimTime = 2.5;               % Time interval after the 0 reference point for each trial 
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
% contained in the trlInfoColDescr output variable. The description is the
% same for both TIM and TRESP data groups.
%
%
% Trialinfo Column Description as presented in trlInfoColDescr output variable
%-------------------------------------------------------------------------------
%========================================================================
% '1. Run Number'
% '2. Block Number within run'
% '3. Nan.  ( This column has been reserved to contain the image ID number which is not encoded in the trigger values. This is not yet implemented.)'
% '4. imgType : 1- Face, 2- Tools  0- Fixation'
% '5. memoryType :  1: 0-Back   2: 2-Back'
% '6. targetType : 1- target, 2- nontarget,  3: lure '
% '7. Trial trigger onset Sample '
% '8. Trial trigger offset Sample'
% '9. Sequence of image in the block'
% '10. isPressed :  0- user did not press any response button, 1- user pressed a response button'
% '11. isPressedLate:  1- If subject responded after the 2 seconds that the image is at the longest displayed and before the next  trial \n    0- If pressed within the presentation time of the image\n ,      NaN: Otherwise'
% '12. isDoubleResponse: 1- user pressed two response buttons in the same trial \n    0: user DID NOT press two response buttons in the same trial'
% '13. pressedCode: Code of the pressed button (If not pressed NaN)'
% '14. isCorrect:  1- If subject has responded that saw a target when a actual target was on or that saw a nontarget when a actual nontarget was on \n      0:  The opposite of the above.\n      NaN:  When subject has not responded or has pressed two  buttons'
% '15. isLureAsCorrect: 1- If subject has responded that saw a target when a lure image of actual target was on \n , 0: In all other cases that a subject has responded\n  NaN:  When subject has not responded or has pressed two buttons'
% '16. respTime: The time from onset of Image to response (sec)'
% '17. respDuration: Duration of button press in seconds'
% '18. isFirstInBlock'
% '19. isLastInBlockk'
% '20. prev. Trial: Run Number'
% '21. prev. Trial: Block Number'
% '22. prev. Trial: Nan.  ( This column has been reserved to contain the image ID number which is not encoded in the trigger values. This is not yet implemented.)'
% '23. prev. Block: imgType : 1- Face, 2-Tools  0- Fixation'
% '24. prev. Block: memoryType :  1: 0-Back   2: 2-Back'
% '25. prev. Trial: targetType : 1- target, 2- nontarget,  3: lure '
% '26. prev. Trial: Trial start Sample'
% '27. prev. Trial: Trial end Sample'
% '28. prev. Trial: Sequence of image in the block'
% '29. prev. Trial: isPressed :  0- user did not press any response button, 1- user pressed a response button'
% '30. prev. Trial: isPressedLate:  1- If subject responded after the 2 seconds that the image is at the longest displayed and before the next  trial \n    0- If pressed within the presentation time of the image\n ,      NaN: Otherwise'
% '31. prev. Trial: isDoubleResponse: 1- user pressed two response buttons in the same trial \n    0: user DID NOT press two response buttons in the same trial'
% '32. prev. Trial: pressedCode: Code of the pressed button (If not pressed NaN)'
% '33. prev. Trial: isCorrect:  1- If subject has responded that saw a target when a actual target was on or that saw a nontarget when a actual nontarget was on \n      0:  The opposite of the above.\n      NaN:  When subject has not responded or has pressed two  buttons'
% '34. prev. Trial: isLureAsCorrect: 1- If subject has responded that saw a target when a lure image of actual target was on \n , 0: In all other cases that a subject has responded\n  NaN:  When subject has not responded or has pressed two buttons'
% '35. prev. Trial: respTime: The time from onset of Image to response (sec)'
% '36. prev. Trial: respDuration: Duration of button press in seconds'
% '37. prev. Trial: isFirstInBlock'
% '38. prev. Trial: isLastInBlock'
% '39. Is button pressed during onset of the stimulus (New field)'};
%==========================================================================

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


% The remaining 39 columns are:
%===== COLUMNS:
trlColDescr={'1. Run Number'
'2. Block Number within run'
'3. Nan.  ( This column has been reserved to contain the image ID number which is not encoded in the trigger values. This is not yet implemented.)'
'4. imgType : 1- Face, 2- Tools  0- Fixation'
'5. memoryType :  1: 0-Back   2: 2-Back'
'6. targetType : 1- target, 2- nontarget,  3: lure '
'7. Trial trigger onset Sample '
'8. Trial trigger offset Sample'
'9. Sequence of image in the block'
'10. isPressed :  0- user did not press any response button, 1- user pressed a response button'
'11. isPressedLate:  1- If subject responded after the 2 seconds that the image is at the longest displayed and before the next  trial \n    0- If pressed within the presentation time of the image\n ,      NaN: Otherwise'
'12. isDoubleResponse: 1- user pressed two response buttons in the same trial \n    0: user DID NOT press two response buttons in the same trial'
'13. pressedCode: Code of the pressed button (If not pressed NaN)'
'14. isCorrect:  1- If subject has responded that saw a target when a actual target was on or that saw a nontarget when a actual nontarget was on \n      0:  The opposite of the above.\n      NaN:  When subject has not responded or has pressed two  buttons'
'15. isLureAsCorrect: 1- If subject has responded that saw a target when a lure image of actual target was on \n , 0: In all other cases that a subject has responded\n  NaN:  When subject has not responded or has pressed two buttons'
'16. respTime: The time from onset of Image to response (sec)'
'17. respDuration: Duration of button press in seconds'
'18. isFirstInBlock'
'19. isLastInBlockk'
'20. prev. Trial: Run Number'
'21. prev. Trial: Block Number'
'22. prev. Trial: Nan.  ( This column has been reserved to contain the image ID number which is not encoded in the trigger values. This is not yet implemented.)'
'23. prev. Block: imgType : 1- Face, 2-Tools  0- Fixation'
'24. prev. Block: memoryType :  1: 0-Back   2: 2-Back'
'25. prev. Trial: targetType : 1- target, 2- nontarget,  3: lure '
'26. prev. Trial: Trial start Sample'
'27. prev. Trial: Trial end Sample'
'28. prev. Trial: Sequence of image in the block'
'29. prev. Trial: isPressed :  0- user did not press any response button, 1- user pressed a response button'
'30. prev. Trial: isPressedLate:  1- If subject responded after the 2 seconds that the image is at the longest displayed and before the next  trial \n    0- If pressed within the presentation time of the image\n ,      NaN: Otherwise'
'31. prev. Trial: isDoubleResponse: 1- user pressed two response buttons in the same trial \n    0: user DID NOT press two response buttons in the same trial'
'32. prev. Trial: pressedCode: Code of the pressed button (If not pressed NaN)'
'33. prev. Trial: isCorrect:  1- If subject has responded that saw a target when a actual target was on or that saw a nontarget when a actual nontarget was on \n      0:  The opposite of the above.\n      NaN:  When subject has not responded or has pressed two  buttons'
'34. prev. Trial:isLureAsCorrect: 1- If subject has responded that saw a target when a lure image of actual target was on \n , 0: In all other cases that a subject has responded\n  NaN:  When subject has not responded or has pressed two buttons'
'35. prev. Trial: respTime: The time from onset of Image to response (sec)'
'36. prev. Trial: respDuration: Duration of button press in seconds'
'37. prev. Trial: isFirstInBlock'
'38. prev. Trial: isLastInBlock'
'39. Is button pressed during onset of the stimulus (New field)'};
%=========================================================================
blockColDescr={''}; % TODO: IF DATA ARE TO BE EXTRACTED IN BLOCKS THEN FORM THE DESCRIPTION OF THE BLOCKS
%=========================================================================
%% INPUTS


datafile=cfg.datafile;
preStimTime=cfg.trialdef.preStimTime;
postStimTime=cfg.trialdef.postStimTime;
cutMode=cfg.trialdef.cutMode;   
lockedOn=cfg.trialdef.lockedOn;

idealISItime=2.5; % This is used to cut Fixation in a similar fashion as the stim blocks
%==========================================================================
%% Read the data file 
hdr = ft_read_header(datafile);
Fsample = hdr.Fs;
trg_resp = ft_read_data(datafile,'chanindx',[1 2],'header',hdr,'eventformat','4d','dataformat','4d');

preStimSamples=floor(Fsample*preStimTime);
postStimSamples=floor(Fsample*postStimTime);
idealISIsamps=floor(Fsample*idealISItime);
%=========================================================================
%% TRIGGER DEFINITIONS

trgMix = trg_resp(1,:);
Nsamples=length(trgMix);

resp = trg_resp(2,:);
% correct possible wrong negative values in response. (Not sure which case this refers to)
ixn = find( resp < 0 );
resp(ixn) = typecast( int16( resp(ixn) ) , 'uint16' );
% Index Finger  (Match to target):      Response Code = 256 ('01 0000 0000')
% Middle Finger (NonMatch to target):   Response Code = 512 ('10 0000 0000')
resp = bitand( resp , bin2dec( '11 0000 0000' ) ); % This ensures that random onset of the most significant bit are removed.

%%
%==========================================
% Separate Photo diode and parallel port triggers
% PhotoDiode only
%--------------------
% Photodiode can take values only 256 ('1 0000 0000') and 0.
trigPhoto = bitand( trgMix , bin2dec('1 0000 0000') );   % This ensures that random onset of the most significant bit are removed.

%Parallel  Port Only
%--------------------
% Parallel port cannot take values outside range 2 (0 0000 0010) and 254 ('0 1111 1110') and 0.
trigPP= bitand(trgMix , bin2dec('0 1111 1110') );   % This ensures that random onset of the most significant bit are removed.

%%
%====================================================
%% Get Events indices in both photodiode and parallel port triggers
difTrigPP=diff(trigPP);
difTrigPhoto=diff(trigPhoto);

indPP_UP=find(difTrigPP>0)+1;
indPP_DOWN=find(difTrigPP<0)+1;
indPP_UPDOWN=find(difTrigPP~=0)+1;

indPhoto_UP=find(difTrigPhoto>0)+1;
indPhoto_DOWN=find(difTrigPhoto<0)+1;
indPhoto_UPDOWN=find(difTrigPhoto~=0)+1;


%----------------------------------------------------------------
%==============================================================
%%
%--------------------------------
% For some scans the photodiode is on at the beginning or/and end of the scan
% So this initial and last interval must not be used in the trial
% extraction

if trigPhoto(1)>0
    trigPhoto(1:indPhoto_DOWN(1)-1)=0;
    difTrigPhoto=diff(trigPhoto);
    [indPhoto_UP]=find((difTrigPhoto>0))+1;
    [indPhoto_DOWN]=find((difTrigPhoto<0))+1;
    [indPhoto_UPDOWN]=find((difTrigPhoto~=0))+1;
end
if trigPhoto(end)>0
    trigPhoto(indPhoto_UP(end):end)=0;
    difTrigPhoto=diff(trigPhoto);
    [indPhoto_UP]=find((difTrigPhoto>0))+1;
    [indPhoto_DOWN]=find((difTrigPhoto<0))+1;
    [indPhoto_UPDOWN]=find((difTrigPhoto~=0))+1;
end

%---------------------------------
%===============================================================
%% First Decode Parallel Port Trigger Events

memTypeTrig=[ 4     % 0-Back
              68    % 2-Back
              2];   % Fixation 
          
imgTypeTrigStep=[ 4   % Face
                 36]; % Tool

targTypeTrigStep=[2    % Non-target
                  4    % Lure
                  6];  % Target  
durTimeCue=2.5;
durTimeImg=2;
durTimeImgFix=0.5;
durTimeFix15=15;
durSampCue=floor(durTimeCue*Fsample);
durSampImg=floor(durTimeImg*Fsample);
durSampImgFix=floor(durTimeImgFix*Fsample);
durSampFix15=floor(durTimeFix15*Fsample);


% The decoding is based on the following trigger protocol used by e-prime
% Cue on trigger  : memTypeTrig
% Img on trigger  : memTypeTrig+imgTypeTrigStep+targTypeTrigStep
% Img off trigger : memTypeTrig+imgTypeTrigStep

    

%{
trlColDescr={

%}    
NeventsPP=length(indPP_UPDOWN);
eventPPInfo=nan(NeventsPP,6);
%{
This structure contains information about all the events in the Parallel
Port
Columns:
'1. event Type: 1. Cue 2. Img On 3. Img Off  4.Fixation 15sec 5. Intro On 6. Intro off
'2. imgType : 1- Face, 2- Place,  3- Body  4- Tools  0- Fixation'
'3. memoryType :  1: 0-Back   2: 2-Back 0-Fixation'
'4. targetType : 1- target, 2- nontarget,  3: lure  0-Cue or Fixation'
'5. Event Sample from PP'
'6. Event Sample from PhotoDiode
%}
    

%%
for iMem=1:length(memTypeTrig), %[ 4     % 0-Back      68    % 2-Back    2];   % Fixation
    
    indIn=find(ismember(trigPP(indPP_UPDOWN),memTypeTrig(iMem)));
    %-----------------------------------------------------
    if isempty(indIn)
        error('No triggers for Working Memory Cue were found');
    else
        if iMem==3
            eventPPInfo(indIn,1)=4;
            eventPPInfo(indIn,2:4)=0;
        else
            eventPPInfo(indIn,1)=1;
            eventPPInfo(indIn,3)=iMem;
            eventPPInfo(indIn,4)=0;
        end
        eventPPInfo(indIn,5)=indPP_UPDOWN(indIn);
    end
    %-----------------------------------------------------
    
    if iMem<3,
        
        tmpImgTypeTrigStep=memTypeTrig(iMem)+imgTypeTrigStep;
        
        for iImg=1:length(tmpImgTypeTrigStep), %memTypeTrig(iMem)+[ 4   % Face   36]; % Tool
            indIn=find(ismember(trigPP(indPP_UPDOWN),tmpImgTypeTrigStep(iImg)));
            %-----------------------------------------------------
            if isempty(indIn)
                error('No triggers for Image off Cue were found');
            else
                eventPPInfo(indIn,1)=3;
                eventPPInfo(indIn,3)=iMem;
                eventPPInfo(indIn,2)=iImg;
                eventPPInfo(indIn,5)=indPP_UPDOWN(indIn);
            end
            %-----------------------------------------------------
            tmpTargTypeTrigStep=tmpImgTypeTrigStep(iImg)+targTypeTrigStep;
            for iTarg=1:length(tmpTargTypeTrigStep), %[2    % Non-target  4    % Lure      6];  % Target
                indIn=find(ismember(trigPP(indPP_UPDOWN),tmpTargTypeTrigStep(iTarg))); 
                %-----------------------------------------------------
                if isempty(indIn)
                    error('No triggers for Image On Cue were found');
                else
                    eventPPInfo(indIn,1)=2;
                    eventPPInfo(indIn,3)=iMem;
                    eventPPInfo(indIn,2)=iImg;
                    if iTarg==1,
                        eventPPInfo(indIn,4)=2;
                    elseif iTarg==2,    
                        eventPPInfo(indIn,4)=3;
                    elseif iTarg==3
                        eventPPInfo(indIn,4)=1;
                    end
                    eventPPInfo(indIn,5)=indPP_UPDOWN(indIn);
                end
                %-----------------------------------------------------
            end
        end
        
    end
end

%======================================================================
% Fill in
% - Image Type for Cue events
% - Target Type for Image Off Events
for iEv=1:NeventsPP-1,
    %-----------------------
    if ((eventPPInfo(iEv,1)==1)&(eventPPInfo(iEv+1,1)==2)),
        eventPPInfo(iEv,2)=eventPPInfo(iEv+1,2);
    end
    %-----------------------
end
for iEv=2:NeventsPP,
    %-----------------------
    if ((eventPPInfo(iEv,1)==3)&(eventPPInfo(iEv-1,1)==2)),
        eventPPInfo(iEv,4)=eventPPInfo(iEv-1,4);
    end
    %-----------------------
end
%======================================================================
%% DO consistency checks on the triggers

warninfo=[];
NeventsPhoto=length(indPhoto_UPDOWN);
rowsPPImg=find(ismember(eventPPInfo(:,1),[2 3]))';
NImgEventsPP=length(rowsPPImg);
indPPImg_UPDOWN=eventPPInfo(rowsPPImg,5)';

if NeventsPhoto~=NImgEventsPP,
    warninfo.dif_ph2pp=[];
    tmpDifP=diff(indPhoto_UPDOWN);
    difThresh=floor((8/60).*Fsample); % 8 monitors refresh frames
    indRand=find(tmpDifP<difThresh);
    if ~isempty(indRand)
        warning(['It seems there are random photodiode events at samples :',num2str(indPhoto_UPDOWN(indRand))]);
        warninfo.dif_ph2pp.samp_extraphotoflips=indPhoto_UPDOWN(indRand);
        NeventsPhoto=NeventsPhoto-2*length(indRand);
        tmpIndPhoto_UPDOWN=indPhoto_UPDOWN;
        tmpIndPhoto_UPDOWN(unique([indRand indRand+1]))=[];
        indPhoto_UPDOWN=tmpIndPhoto_UPDOWN;
   else
        
        if (NeventsPhoto-NImgEventsPP)==2,
            diffTimePhotoUP=(1/Fsample)*diff(indPhoto_UP);
            indDiffLong=find(diffTimePhotoUP>25);
            if length(indDiffLong)>1
                error('There are more than one extra up-down photodiode events');
            end
            if indDiffLong~=1,
                error('There is one extra up-down photodiode event but is not at the beginnin of the scan');
            end
            indPhoto_UPDOWN(1:2)=[];
            indPhoto_UP(1)=[];
            indPhoto_DOWN(1)=[];
            NeventsPhoto=length(indPhoto_UPDOWN);
            
        else
            error('The image onsets-offsets is different in the photodiode and parallel port triggers');
        end
    end
    
end



%------
% Match Parallel port with Photodiode Triggers
difThresh=floor((4/60).*Fsample); % 6 monitors refresh frames
tmpDifPhoto2PP=indPhoto_UPDOWN-indPPImg_UPDOWN;

if  mean(abs(tmpDifPhoto2PP)./Fsample) > (2/60)
    warning('The avergae difference between Photodiode and Parallel Port Triggers is more than 2 monitor refresh rates.'); 
    warninfo.dif_ph2pp.avgdelay=mean(abs(tmpDifPhoto2PP)./Fsample);    
end

indBad=find(abs(tmpDifPhoto2PP)>difThresh);
if ~isempty(indBad)
   warning('For some image onset or offset event the difference between Photodiode and Parallel Port Triggers are more than 6 monitor refresh rates.'); 
   warninfo.dif_ph2pp.samp_longdelay=indPhoto_UPDOWN(indBad);
end
eventPPInfo(rowsPPImg,6)=indPhoto_UPDOWN;
%-----

nanIndx=find(isnan(eventPPInfo(:,1)))';
difNanIndx=diff(nanIndx);

if length(nanIndx)>2
    error('More than 2 PP events remained non assigned');
end
if length(nanIndx)==0
    error('No PP events remained non assigned. There should be one or two at the beginning for the introductory screen and the countdown');
end
if length(find(difNanIndx>1))
    error('The unassigned PP events are not adjacent')
end
if max(nanIndx)>2
       error('Only the first 2 PP events at max can be unassigned')
end

eventPPInfo(nanIndx,5)=indPP_UPDOWN(nanIndx);%Just fill in the onset sample of unassigned events

if length(nanIndx)==1, % This must be a downward trigger event
    if nanIndx~=1,
       error('There is one unassigned event , but is not the first one as expected'); 
    else
       if difTrigPP(eventPPInfo(nanIndx,5)-1)>=0
          error('There is one unassigned event at the beginning but the trigger is not down going') ;
       end
    end
elseif length(nanIndx)==2
    if ~all(nanIndx==[1 2])
         error('There are two unassigned events , but are not the first two as expected'); 
    else
        if ~((difTrigPP(indPP_UPDOWN(nanIndx(1))-1)>0)&(difTrigPP(indPP_UPDOWN(nanIndx(2))-1)<0))
            warning('There are two unassigned events at the beginning but the triggers are not up going and then down going') ;
            warninfo.pp.samp_initTrigNotUpDown=[indPP_UPDOWN(nanIndx(1)) indPP_UPDOWN(nanIndx(2))];
        end
    end  
end

if eventPPInfo(end,1)~=4
   error('The last event should be the start of the last fixation block');
else
    if (eventPPInfo(end,5)+durSampFix15)>Nsamples,
       error('The last fixation block is not as long as expected(15sec). Scan ends earlier.');
    end
end

%======================================
%% If all consistency tests are passed try to read the run number
runNumber=0;
if length(nanIndx)==1, % This must be a downward trigger event
   tmpRunTrig=trigPP(indPP_UPDOWN(nanIndx(1))-1);
elseif length(nanIndx)==2
   tmpRunTrig=trigPP(indPP_UPDOWN(nanIndx(2))-1); 
end
if tmpRunTrig==130,
    runNumber=1;
elseif tmpRunTrig==132,
    runNumber=2;
end
%======================================================================
%% If the code has reached this point with no error then all PP event information must have been decoded,
% and the Image Onset and Offest triggers must be matched between PP and
% Photodiode. No photodiode triggers are available for the cues or for the
% Fixation. So for the decoding, the trigger samples for Cues and Fixation
% blocks will extracted from the PP while the trigger samples for the Image
% onsets will be recovered from the Photodiode triggers.

% FIND BLOCK ON BLOCK OFF - Block starts a the onset of a cue and on sample
% before the onset of the next cue. For the last fixation block the
% end is defined by the default duration of a fixation block.
rowsCueOrFix=find(ismember(eventPPInfo(:,1),[1 4]));
blockPPInfo=eventPPInfo(rowsCueOrFix,:);

indxBlockOn=blockPPInfo(:,5);
indxBlockOff=blockPPInfo(2:end,5)-1;
indxBlockOff(end+1)=indxBlockOn(end)+durSampFix15;
    
Nblocks=length(indxBlockOn);
%======================================================
%% Identify the start and end of scan as the onset of the first block and the end of the last block
% A default temporal padding of 3 second is used in either end
defEndPad=floor(3*Fsample);
scanStartSamp=indxBlockOn(1)-defEndPad;
if scanStartSamp<1 
    scanStartSamp=1;
end
scanEndSamp=indxBlockOff(end)+defEndPad;
if scanEndSamp>Nsamples,
    scanEndSamp=Nsamples;
end
%======================================================
%=======================================================================
 %% GET BASIC TRIALINFO FOR STIM SEQUENCE
  
totalStimInfo=[];
countStims=1;
curBlockNum=0;
for iEv=1:NeventsPP

    if ismember(eventPPInfo(iEv,1),[1 4]),
       curBlockNum=curBlockNum+1;
       curStimSeq=1;
   elseif eventPPInfo(iEv,1)==2,
         if eventPPInfo(iEv+1,1)~=3
             error('Image onset is not followed by an image offset');
         end
       tmpInfo=[runNumber,...               %1.
                 curBlockNum ,...           %2.
                 nan,...                    %3.  % here should be an index refereing to the actual image used
                 eventPPInfo(iEv,2),...     %4.  
                 eventPPInfo(iEv,3),...     %5.
                 eventPPInfo(iEv,4),...     %6.
                 eventPPInfo(iEv,6),...     %7.
                 eventPPInfo(iEv+1,6)-1,... %8.
                 curStimSeq];               %9. %not used. just added for compatibility with glasgow data; Maybe here could be used the sequence of image in the block
             
             totalStimInfo=[totalStimInfo; tmpInfo];             
             curStimSeq=curStimSeq+1;
   end
   
   
end
%==============================================
%% Get Response characteristics
totalStimInfo(:,10:39)=0;

Nstims=size(totalStimInfo,1);
for iStim=1:Nstims,
    indStart=totalStimInfo(iStim,7);
    indEnd=totalStimInfo(iStim,8);
    indTotalEnd=Nsamples;
    tmpResp=resp(indStart:indEnd);
    tmpResp2End=resp(indStart:indTotalEnd);
    tmpDifResp=diff(tmpResp);
    tmpDifResp2End=diff(tmpResp2End);
    [i1,j1]=find(tmpResp==256);
    [k1,l1]=find(tmpResp==512);
    
    if (tmpResp(1)==256)|(tmpResp(1)==512)
        isPressedOnOnset=1;
    else
        isPressedOnOnset=0;
    end
    
    %-----------------------------------------------------
    if (~isempty(i1))&(~isempty(k1)),
        isPressed=1;
        isDoubleResponse=1;
        pressedCode=nan;
        isCorrect=nan;
        isLureAsCorrect=nan;
        respTime=nan;
        respDur=nan; 
        isPressedLate=nan;
        
    elseif (isempty(i1))&(isempty(k1)),
        isPressed=0;
        isDoubleResponse=0;
        pressedCode=nan;
        isCorrect=nan;
        isLureAsCorrect=nan;
        respTime=nan;
        respDur=nan; 
        
        %=======================================
        % Check for late response after the stim is off and before the next
        % stim is on
        if iStim<Nstims,
            indLateStart=totalStimInfo(iStim,8);
            indLateEnd=totalStimInfo(iStim+1,7);
        else
            indLateStart=totalStimInfo(iStim,8);
            indLateEnd=Nsamples;
        end
         tmpLateResp=resp(indLateStart:indLateEnd);
         lat1=find((tmpLateResp==256)|(tmpLateResp==512));
         if ~isempty(lat1)
             isPressedLate=1;
         else
             isPressedLate=nan;
         end
        %=======================================
    elseif (~isempty(i1))&(isempty(k1)),   
        isPressed=1;
        isDoubleResponse=0;
        pressedCode=256;
        if totalStimInfo(iStim,6)==1
             isCorrect=1;
             isLureAsCorrect=0;
        elseif totalStimInfo(iStim,6)==3
             isCorrect=0;
             isLureAsCorrect=1;
        else
             isCorrect=0;
             isLureAsCorrect=0;
        end
        respTime=j1(1)*(1./hdr.Fs);
        [d1,d2]=find(tmpDifResp(j1(1):end)<0,1,'first');
        if isempty(d1),
            [d1,d2]=find(tmpDifResp2End(j1(1):end)<0,1,'first');
        end
        respDur=d2*(1./hdr.Fs);
        isPressedLate=0;
        
    elseif (isempty(i1))&(~isempty(k1)),   
        isPressed=1;
        isDoubleResponse=0;
        pressedCode=512;
        if totalStimInfo(iStim,6)>1
             isCorrect=1;
             isLureAsCorrect=0;
        elseif totalStimInfo(iStim,6)==1
             isCorrect=0;
             isLureAsCorrect=0;
        else
             isCorrect=0;
             isLureAsCorrect=0;
        end
        respTime=l1(1)*(1./hdr.Fs);
        [d1,d2]=find(tmpDifResp(k1(1):end)<0,1,'first');
        if isempty(d1),
            [d1,d2]=find(tmpDifResp2End(k1(1):end)<0,1,'first');
        end
        respDur=d2*(1./hdr.Fs);
        isPressedLate=0;
    end
    %-----------------------------------------------------
            totalStimInfo(iStim,10)=isPressed;
            totalStimInfo(iStim,11)=isPressedLate;
            totalStimInfo(iStim,12)=isDoubleResponse;
            totalStimInfo(iStim,13)=pressedCode;
            totalStimInfo(iStim,14)=isCorrect;
            totalStimInfo(iStim,15)=isLureAsCorrect;
            totalStimInfo(iStim,16)=respTime;
            totalStimInfo(iStim,17)=respDur;
            
            totalStimInfo(iStim,39)=isPressedOnOnset; % New field to check if the button is pressed when the stim comes on; Appended to the end of the original list
    
end



%==========================================
%% Mark the first and the last in Block
for iBlock=1:Nblocks
 [i1,j1]=find(totalStimInfo(:,2)==iBlock);  
 if ~isempty(i1)
    totalStimInfo(i1(1),18)=1;   
    totalStimInfo(i1(end),19)=1;   
 end
end


%===========================================
%% Append trialinfo for previous trial
pastTrialInfo=nan(size(totalStimInfo,1),19);
pastTrialInfo(2:end,:)=totalStimInfo(1:end-1,1:19);
totalStimInfo(:,20:38)=pastTrialInfo;

% --- assign memory and image types of previous block. This is because
% witin each block the memory and image type is always the same. So by
% knowinn the memory and image type of the current , also the ones for the
% previous are known. What can be usefull information is what was the
% memory and image type of the preceding block within the run.
for iBlock=2:Nblocks
 [i1,j1]=find(totalStimInfo(:,2)==iBlock);  
 if ~isempty(i1)
    totalStimInfo(i1,23)=blockPPInfo(iBlock-1,2);   
    totalStimInfo(i1,24)=blockPPInfo(iBlock-1,3);      
 end
end

%====================================================
%====================================================
%====================================================
%====================================================
%====================================================
%====================================================
%============================================================== 
%==============================================================
%==============================================================
%% Fixation  BLocks
% FIXATIONBLOCKINFO
% Columns:
% 1. Previous block ImgType
% 2. Previous block memoryType
% 3. Next block ImgType
% 4. Next block blockType
% 5. Fixation Start Sample
% 6. Fixation End Sample

totalBlockFixInfo=[];
for iBlock=1:Nblocks,

    if blockPPInfo(iBlock,1)==4,
       
        
        if iBlock>1,
           prevInfo= blockPPInfo(iBlock-1,[2:3]);
        else
            prevInfo=[nan nan];
        end
        if iBlock<Nblocks,
            nextInfo=blockPPInfo(iBlock+1,[2:3]);
        else
            nextInfo=[nan nan];
        end
        totalBlockFixInfo=[totalBlockFixInfo;  [prevInfo  nextInfo indxBlockOn(iBlock) indxBlockOff(iBlock)]];
        
    end
end



%%
%==============================================================
%% Stim Blocks
%======================================
% STIMBLOCKINFO
% Columns:
% 1. Nblock
% 2. ImgType
% 3. memoryType
% 4. Previous block ImgType
% 5. Previous block memoryType
% 6. Block Start Sample
% 7. Block End Sample
% 8:17: Images Onset Samples in block
% 18:27: Target type flag for the Images in block
% 28:37: isCorrect binary flag for the subject's responses

rowsStim=find(blockPPInfo(:,1)==1);
blockStimPPInfo=blockPPInfo(rowsStim,:);
NstimBlocks=size(blockStimPPInfo,1);

totalBlockStimInfo=nan(NstimBlocks,37);
countStimBlocks=1;
for iBlock=1:Nblocks,
    if blockPPInfo(iBlock,1)==1
        totalBlockStimInfo(countStimBlocks,1)=iBlock;
        
        [i1,j1]=find(totalStimInfo(:,2)==iBlock);
        
        totalBlockStimInfo(countStimBlocks,2)=totalStimInfo(i1(1),4);
        totalBlockStimInfo(countStimBlocks,3)=totalStimInfo(i1(1),5);
        totalBlockStimInfo(countStimBlocks,4)=totalStimInfo(i1(1),23);
        totalBlockStimInfo(countStimBlocks,5)=totalStimInfo(i1(1),24);
        totalBlockStimInfo(countStimBlocks,6)=indxBlockOn(iBlock);
        totalBlockStimInfo(countStimBlocks,7)=indxBlockOff(iBlock);
        totalBlockStimInfo(countStimBlocks,8:17)=totalStimInfo(i1,7);
        totalBlockStimInfo(countStimBlocks,18:27)=totalStimInfo(i1,6); 
        totalBlockStimInfo(countStimBlocks,28:37)=totalStimInfo(i1,14); % create a binary structure 1 for target , 0 otherwise
        
        countStimBlocks=countStimBlocks+1;
    
    end
    
end
%%
%============================================================
%% ========== Fuse Trials with Fixation blocks in trialinfo =====

NfixBlocks=size(totalBlockFixInfo,1);
fixTrialInfo=nan(NfixBlocks,size(totalStimInfo,2));
fixTrialInfo(:,4)=zeros(NfixBlocks,1);
fixTrialInfo(:,23)=totalBlockFixInfo(:,1);
fixTrialInfo(:,24)=totalBlockFixInfo(:,2);
fixTrialInfo(:,7)=totalBlockFixInfo(:,5);
fixTrialInfo(:,8)=totalBlockFixInfo(:,6);
%===================================================
% FORM FIXATION PSEUDO TRIALS

halfISIsamps=floor(idealISIsamps./2);
psFixTrialInfo=[];
for iFixTrl=1:size(fixTrialInfo,1),
    indxStart   = fixTrialInfo(iFixTrl,7)+halfISIsamps;
        indxEnd     = fixTrialInfo(iFixTrl,8)-halfISIsamps;
        psixtrl=[indxStart:idealISIsamps:indxEnd];
        tmpNtrials=length(psixtrl);
        
        tmpinfo=repmat(fixTrialInfo(iFixTrl,:),tmpNtrials,1);
        tmpinfo=[psixtrl'-preStimSamples psixtrl'+postStimSamples -preStimSamples*ones(tmpNtrials,1) tmpinfo];
    
         psFixTrialInfo=[psFixTrialInfo; tmpinfo];
    
end
 [i1,j1]=find(psFixTrialInfo(:,2)>Nsamples);
 if ~isempty(i1)
     psFixTrialInfo(i1,:)=[];
 end
 [i1,j1]=find(psFixTrialInfo(:,1)<1);
  if ~isempty(i1)
     psFixTrialInfo(i1,:)=[];
 end
%==================================================

%============================================================
%% ------------CREATE OUTPUT
%=========================================================================
trl=[];
if strcmp(cutMode,'trials'),    
outTrialInfo=totalStimInfo;
if strcmp(lockedOn,'Image')
    trlStims=[outTrialInfo(:,7)-preStimSamples outTrialInfo(:,7)+postStimSamples -repmat(preStimSamples,size(outTrialInfo,1),1)  outTrialInfo];
    trl=[trlStims; psFixTrialInfo];
elseif strcmp(lockedOn,'Response')
    outTrialInfo=outTrialInfo((logical(outTrialInfo(:,10)))&(~logical(outTrialInfo(:,12)))&(~logical(outTrialInfo(:,39))),:);        
    trlStims=[outTrialInfo(:,7)+floor(hdr.Fs*outTrialInfo(:,16))-preStimSamples outTrialInfo(:,7)+floor(hdr.Fs*outTrialInfo(:,16))+postStimSamples -repmat(preStimSamples,size(outTrialInfo,1),1)  outTrialInfo];
    trl=[trlStims; psFixTrialInfo];
end

trlInfoColDescr=trlColDescr;
%===================================================================
elseif strcmp(cutMode,'blocks'),
    
   outTrialInfo=totalBlockStimInfo;
   trl=[outTrialInfo(:,6)-preStimSamples outTrialInfo(:,7)+postStimSamples -repmat(preStimSamples,size(outTrialInfo,1),1)  outTrialInfo];

trlInfoColDescr=blockColDescr;

elseif strcmp(cutMode,'fixation'),
    outTrialInfo=totalBlockFixInfo;  
    trl=[outTrialInfo(:,5)-preStimSamples outTrialInfo(:,6)+postStimSamples -repmat(preStimSamples,size(outTrialInfo,1),1)  outTrialInfo];
trlInfoColDescr={''};  
end



%===================================================================

trialSummary=[];
%===========================================
%--- Compute summary of trial information to print is summary file
trialSummary.Nblocks=NstimBlocks+NfixBlocks;
trialSummary.Nblocks_Fixation=NfixBlocks;
trialSummary.Nblocks_0Back=length(find(totalBlockStimInfo(:,3)==1));
trialSummary.Nblocks_2Back=length(find(totalBlockStimInfo(:,3)==2));
trialSummary.Nblocks_faces=length(find(totalBlockStimInfo(:,2)==1));
trialSummary.Nblocks_tools=length(find(totalBlockStimInfo(:,2)==2));
trialSummary.Ntrials=size(totalStimInfo,1);
trialSummary.Ntrials_0Back=length(find(totalStimInfo(:,5)==1));
trialSummary.Ntrials_2Back=length(find(totalStimInfo(:,5)==2));
trialSummary.Ntrials_faces=length(find(totalStimInfo(:,4)==1));
trialSummary.Ntrials_tools=length(find(totalStimInfo(:,4)==2));
trialSummary.NtrialsRESP=length(find((totalStimInfo(:,10)==1)));
trialSummary.NtrialsRESP_0Back=length(find((totalStimInfo(:,5)==1)&(totalStimInfo(:,10)==1)));
trialSummary.NtrialsRESP_2Back=length(find((totalStimInfo(:,5)==2)&(totalStimInfo(:,10)==1)));
trialSummary.NtrialsRESP_faces=length(find((totalStimInfo(:,4)==1)&(totalStimInfo(:,10)==1)));
trialSummary.NtrialsRESP_tools=length(find((totalStimInfo(:,4)==2)&(totalStimInfo(:,10)==1)));
%===========================================
%===================================================================




end
