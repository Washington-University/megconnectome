function    [trl,trlInfoColDescr, trialSummary,scanStartSamp,scanEndSamp,warninfo]  = trialfun_StoryM_BaseExtractAll( cfg )
%% This is the Trial Definition Function for the Story/Math experiment.
% It extracts the trial definition for ALL trials within each of different
% datagroup.
% There are 4 different data groups for Story/Math
%
% DATA GROUPS
%----------------------------------
% 1. mnemonic: TEV.    description:  Onset of any task event during the question
%                                      and option period in stories and math problems.
%                                      All trials have the same fixed duration.
% 2. mnemonic: TRESP.  description:  Onset of button press by the subject. All trials 
%                                      have the same fixed duration.
% 3. mnemonic: BSENT.  description:  Trials contain entire sentences. For stories this
%                                      is a sentence during narration without the option 
%                                      sentence at the end of the story. For Maths this 
%                                      is the sentence of the Math problem excluding the 
%                                      option sentence at the end. Trials
%                                      have variable duration.
% 4. mnemonic: BUN.    description:  Trials containing entire Blocks of stimulus Units. 
%                                      As stimulus unit is defined an entire story or 
%                                      an entire math problem including the option part.
%                                      Trials have variable duration.
%
% 
% INPUT VARIABLE
%----------------------------------
% cfg : This is a structure containing information required for extracting
%       the trials for each of the 4 data groups described above.
%       Fields:
%               .datafile: This is the filename of the raw data file.
%               .trialdef: This is a structure containing the parameters
%                          required to split the data into the trials of each of the data
%                          groups.
%                          Fields:
%                          .trialdef.cutmode = 1 or 2 or 3 or 4;       % Defines the data group cut mode. 1 for TEV, 
%                                                                          2 for TRESP, 3 for BSENT, 4 for BUN.
%                          .trialdef.preStimTime = 1.5;                % Time interval prior to 0 reference point for each trial 
%                          .trialdef.postStimTime = 1.5;               % For fixed length trial data groups TEV and TRESP this is
%                                                                          the time interval after the 0 reference point for each 
%                                                                          trial. For variable length trial data groups BSENT and BU,
%                                                                          this is the time interval after the END of the given event 
%                                                                          to be added to the trial.
%
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
%                  in trials of the main conditions for ALL data groups.
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
%
%
% Trialinfo Column Description as presented in trlInfoColDescr output variable
%-------------------------------------------------------------------------------
%========================================================================
% data group: TEV
%-----------------
%  '1. Block Number within Run'
%  '2. Unit Type : 1.Story 2.Math'
%  '3. Unit Number within Run'
%  '4. Total Number of units (N of Stories or N of Math problems ) in same Run'
%  '5. Unit Number within Block'
%  '6. Total Number of units (N of Stories or N of Math problems ) in same Block'
%  '7. Attribute1: For story this is the story number. For Math is the difficulty level'
%  '8. Unit Narration interval Start Sample - Start with the onset of the  first word trigger or the beginning of the first sentence'
%  '9. Unit Narration interval End Sample'
%  '10. N subunits within Narration interval'
%  '11. Unit Option interval Start Sample'
%  '12. Unit Option interval End Sample'
%  '13. N subunits within Option interval'
%  '14 Correct Option- 1 or 2'
%  '15 Unit Response interval start sample'
%  '16 Unit Response interval end sample'
%  '17 Unit Response sample'
%  '18 is Response Correct'
%  '19 is Response Early'
%  '20 Event Sample'
%  '21 Event Type - 20.Math Narration number word,\n 21.Math Narration operand word,\n 22.Math Option Intro Word,\n 23.Math Option 1 Word\n  24. Math Option OR Word\n 25.Math Option 2 Word\n 10.Story Sentence,\n 12.Story Option Intro\n 13.Math Option 1 Word\n  14. Math Option OR Word\n 15.Math Option 2 Word'
%  '22. Narration Event Number in Narration interval  (Applies only to Sentence or number or operation in Narration interval)'};
%========================================================================
% data group: TRESP or BUN
%-----------------
%  '1. Block Number within Run'
%  '2. Unit Type : 1.Story 2.Math'
%  '3. Unit Number within Run'
%  '4. Total Number of units (N of Stories or N of Math problems ) in same Run'
%  '5. Unit Number within Block'
%  '6. Total Number of units (N of Stories or N of Math problems ) in same Block'
%  '7. Attribute1: For story this is the story number. For Math is the difficulty level'
%  '8. Unit Narration interval Start Sample - Start with the onset of the first word trigger or the beginning of the first sentence'
%  '9.  Unit Narration interval End Sample'
%  '10. N subunits within Narration interval'
%  '11. Unit Option interval Start Sample'
%  '12. Unit Option interval End Sample'
%  '13. N subunits within Option interval'
%  '14. Option Intro Sample - "equals" or "That was about"'
%  '15. Option1 onset sample'
%  '16. OR onset sample'
%  '17. Option2 onset sample'
%  '18. Correct Option- 1 or 2'
%  '19. Unit Response interval start sample'
%  '20. Unit Response interval end sample'
%  '21. Unit Response sample'
%  '22. is Response Correct'
%  '23. is Response Early'};
%========================================================================
% data group: BSENT
%-----------------
%  '1. Block Number within Run'
%  '2. Unit Type : 1.Story 2.Math'
%  '3. Unit Number within Run'
%  '4. Total Number of units (N of Stories or N of Math problems ) in same Run'
%  '5. Unit Number within Block'
%  '6. Total Number of units (N of Stories or N of Math problems ) in same Block'
%  '7. Attribute1: For story this is the story number. For Math is the  difficulty level'
%  '8. N sentences within Narration interval (always one for math as in math one narration interval corresponds to one math sentence)'
%  '9. is Response Correct (For Story this refers to the response at the very end of the sentence)'
%  '10. is Response Early'
%  '11. Narration Sentence Number in Narration interval  (For math always equal to one)'
%  '12. Narration Sentence Start Sample'
%  '13. Narration Sentence End Sample'};
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

if ~isfield(cfg.trialdef, 'plotresults'),  cfg.trialdef.plotresults = 'no';     end

%% =========== Read cfg info
datafile = cfg.datafile;
cutmode=cfg.trialdef.cutmode;
prestimTime = cfg.trialdef.prestimTime;
poststimTime = cfg.trialdef.poststimTime;
if strcmpi(cfg.trialdef.plotresults , 'yes' ), plot_flag = 1; else plot_flag = 0; end

%%
hdr = ft_read_header(datafile);
Fsample = hdr.Fs;
prestimSamples = floor(prestimTime*Fsample);
poststimSamples = floor(poststimTime*Fsample);

detTrig=ft_read_data(datafile,'chanindx',[1],'header',hdr,'eventformat','4d','dataformat','4d');
detTrig=bitand(detTrig,255);
detResp=ft_read_data(datafile,'chanindx',[2],'header',hdr,'eventformat','4d','dataformat','4d');
Nsamps=length(detTrig);



%---------------------------------------

%% TRIGGERS
%==========================================
%-- Trigger definitions from the Eprime
EP=[];
EP.TrigCodStroryMathNumLevel=128;
EP.TrigCodStoryCue=96;
EP.TrigCodMathCue=112;
EP.TrigCodStoryBase=32;
EP.TrigCodStoryTrigStart=24;
EP.TrigCodStoryTrig=26;
EP.TrigCodStoryEndWait=16;
EP.TrigCodStoryThat=14;
EP.TrigCodStoryOptUnCor=12;
EP.TrigCodStoryOptCor=10;
EP.TrigCodStoryOR=8;
EP.TrigCodStoryResp=2;
EP.TrigCodMathBase=64;
EP.TrigCodMathQuestBlock=16;
EP.TrigCodMathQuestTrig=8;
EP.TrigCodMathEqualBlockUnCor=8;
EP.TrigCodMathEqualBlockCor=0;
EP.TrigCodMathEqualTrig=4;
EP.TrigCodMathResp=2;
%==========================================

%==========================================
%-- Form Triggers
trg_headblock=128; % Before each math or story there is header trigger with some superimposed triggers
% For story there are 2 superimposed triggers. The first
% is the story id *2. The second?
% For math there are 3 superimposed triggers. The first
% is the math difficulty level*2. The other two triggers?

trgST_Cue=EP.TrigCodStoryCue; % These 2 cue triggers are sent everytime a change between a story and math
% or vice versa takes place, so they signal
% the big blocks appart from the first one
trgMA_Cue=EP.TrigCodMathCue;
%-------------------------------------------------------------

trgST_start=EP.TrigCodStoryBase+EP.TrigCodStoryTrigStart;
trgST_senton=EP.TrigCodStoryBase+EP.TrigCodStoryTrig;
trgST_endwait=EP.TrigCodStoryBase+EP.TrigCodStoryEndWait;% There is a 500 msec wait after the end of the last sentence? This represents the end of this 500 msec?
trgST_that=EP.TrigCodStoryBase+EP.TrigCodStoryThat; % This is the beginning of the sentence "That story was about"
trgST_OptCor=EP.TrigCodStoryBase+EP.TrigCodStoryOptCor; % This is the correct option;
trgST_OptUnCor=EP.TrigCodStoryBase+EP.TrigCodStoryOptUnCor; % This is the uncorrect option;
trgST_RespBlock=EP.TrigCodStoryBase+EP.TrigCodStoryResp; % This is trigger of the Response Block
trgST_OptOR=EP.TrigCodStoryBase+EP.TrigCodStoryOR;
%--------------------------------------------------------------
trgMA_N_base=EP.TrigCodMathBase+EP.TrigCodMathQuestBlock; %This is the basis trigger in the Math Narration interval upon which the triggers for inidividual words ans numbers are superimposed
trgMA_N_wordon=EP.TrigCodMathBase+EP.TrigCodMathQuestBlock+EP.TrigCodMathQuestTrig; % Math: Narration block: onset of number? There are
% 7 of those and some time after the last one the
% Narration block ends end trigger falls. So it seems
% like the triggers correspond to 4 numbers being
% participating in the addition/subtraction. i.e.
% one plus two plus three minus four

trgMA_O_Cor1=EP.TrigCodMathBase+EP.TrigCodMathEqualBlockCor; % After the narration block is the option block. The base trigger for
% this block
% is different if the 1st or the second option is the correct one.
% The triggers for the  words "equals" option1 "or" option2
% are then  superimposed  on this Base trigger
trgMA_O_Cor2=EP.TrigCodMathBase+EP.TrigCodMathEqualBlockUnCor;
trgMA_O_Cor1wordon=trgMA_O_Cor1+EP.TrigCodMathEqualTrig; % This is the trigger for words "equals" option1 "or" option2
trgMA_O_Cor2wordon=trgMA_O_Cor2+EP.TrigCodMathEqualTrig;
trgMA_RespBlock=EP.TrigCodMathBase+EP.TrigCodMathResp; % This is trigger of the Response Block

durTimeResp=3;
durSampResp=floor(durTimeResp*Fsample);
%---------------------------------------------------------------
% - Left Button should be the response for Option 1
% and Right Button should be the response for Option 2 ???
respLeft=256;
respRight=512;
%-------------------------------------------

%% create indices of up down of trigger and response. Find indices of block changes and header interval onset
difTrig=diff(detTrig);
indxTrigUp=find(difTrig>0)+1;
indxTrigDown=find(difTrig<0)+1;
indxTrigUpDown=find(difTrig~=0)+1;

difResp=diff(detResp);
indxRespUp=find(difResp>0)+1;
indxRespDown=find(difResp<0)+1;

%=====================================================
%=====================================================
%=====================================================
warninfo=[];
%-Identify Story and Math Blocks and make
[i1,j1]=find(ismember(detTrig(2:end),[trgST_Cue trgMA_Cue])&(difTrig>0));
indxBlockChange=j1+1;
if (length(indxBlockChange))~=7,
    disp(['WARNING!!! - It seems that there are NOT 8 blocks in your story math raw data. ']);
end
indxHeadOn=find((detTrig(2:end)==trg_headblock)&(difTrig>0))+1;
Nblocks=length(indxBlockChange)+1;

%-----------------------------------------------------------------

%% Create trlBlockInfo
%-----------------------------
% Create trlBlockInfo
% COLUMNS:
% 1. Block Type : 1.Story 2.Math
% 2. Block Start Sample
% 3. Block End Sample
% 4. Number of units (N of Stories or N of Math problems ) in Block

% Define as the beginning of the block the onset of the first header block
% within the block
% Define the end of each block as the sample before the next block change signal (for the last block is the last sample)
trlBlockInfo=nan(Nblocks,4);
for iBlock=Nblocks-1:-1:1,
    trlBlockInfo(iBlock+1,1)=find(detTrig(indxBlockChange(iBlock))==[trgST_Cue trgMA_Cue]); % 1.Story 2.Math
    trlBlockInfo(iBlock+1,2)=indxHeadOn(find(indxHeadOn>indxBlockChange(iBlock),1,'first'));
    if iBlock==(Nblocks-1),
        
        indRespLast=indxTrigUpDown(find(ismember(detTrig(indxTrigUpDown),[trgMA_RespBlock  trgST_RespBlock]),1,'last'));
        trlBlockInfo(iBlock+1,3)=indRespLast+durSampResp;
    else
        trlBlockInfo(iBlock+1,3)=indxBlockChange(iBlock+1)-1;
    end
end
if trlBlockInfo(2,1)==trgST_Cue
    trlBlockInfo(1,1)=2;
else
    trlBlockInfo(1,1)=1;
end
trlBlockInfo(1,2)=indxHeadOn(find(indxHeadOn>1,1,'first'));
trlBlockInfo(1,3)=indxBlockChange(1)-1;
%-----------------------------
% Find how many units (N of Stories or N of Math problems ) are in each
% block
for iBlock=1:Nblocks,
    trlBlockInfo(iBlock,4)=sum((indxHeadOn>=trlBlockInfo(iBlock,2))&(indxHeadOn<trlBlockInfo(iBlock,3)));
end
%-----------------------------
%======================================================
%% Identify the start and end of scan as the onset of the first block and the end of the last block
% A default temporal padding of 3 second is used in either end
defEndPad=floor(3*Fsample);
scanStartSamp=trlBlockInfo(1,2)-defEndPad;
if scanStartSamp<1 
    scanStartSamp=1;
end
scanEndSamp=trlBlockInfo(end,3)+defEndPad;
if scanEndSamp>Nsamps,
    scanEndSamp=Nsamps;
end
%======================================================



%figure;hold on;
%plot(detTrig);
%plot(trlBlockInfo(:,2),detTrig(trlBlockInfo(:,2)),'.g');
%plot(trlBlockInfo(:,3),detTrig(trlBlockInfo(:,3)),'.r');

%--------------------------------------------------------


%% Create general trlUnitInfo
%=====================================================
%=====================================================
%=====================================================
%=====================================================
% Create general trlUnitInfo
% Columns:
%----------------------------------------------
trlColDescr_UnitInfo={'1. Block Number within Run'
'2. Unit Type : 1.Story 2.Math'
'3. Unit Number within Run'
'4. Total Number of units (N of Stories or N of Math problems ) in same Run'
'5. Unit Number within Block'
'6. Total Number of units (N of Stories or N of Math problems ) in same Block'
'7. Attribute1: For story this is the story number. For Math is the difficulty level'
'8. Unit Narration interval Start Sample - Start with the onset of the first word trigger or the beginning of the first sentence'
'9.  Unit Narration interval End Sample'
'10. N subunits within Narration interval'
'11. Unit Option interval Start Sample'
'12. Unit Option interval End Sample'
'13. N subunits within Option interval'
'14. Option Intro Sample - "equals" or "That was about"'
'15. Option1 onset sample'
'16. OR onset sample'
'17. Option2 onset sample'
'18. Correct Option- 1 or 2'
'19. Unit Response interval start sample'
'20. Unit Response interval end sample'
'21. Unit Response sample'
'22. is Response Correct'
'23. is Response Early'};

Ntotalun=sum(trlBlockInfo(:,4));
trlUnitInfo=nan(Ntotalun,22);
countUn=1;
for iBlock=1:Nblocks,
    tmpNun=trlBlockInfo(iBlock,4);
    
    tmpIndxHeadOn=indxHeadOn(find((indxHeadOn>=trlBlockInfo(iBlock,2))&(indxHeadOn<trlBlockInfo(iBlock,3))));
    
    for iUn=1:tmpNun
        trlUnitInfo(countUn,1)=iBlock;
        trlUnitInfo(countUn,2)=trlBlockInfo(iBlock,1);
        trlUnitInfo(countUn,3)=countUn;
        trlUnitInfo(countUn,4)=Ntotalun;
        trlUnitInfo(countUn,5)=iUn;
        trlUnitInfo(countUn,6)=trlBlockInfo(iBlock,4);
        trlUnitInfo(countUn,7)=0.5*(detTrig(indxTrigUp(find((indxTrigUp>tmpIndxHeadOn(iUn)),1,'first'))) - trg_headblock);
        
        if trlUnitInfo(countUn,2)==1 % Story
            trlUnitInfo(countUn,8)=indxTrigUp(find((indxTrigUp>tmpIndxHeadOn(iUn))&(detTrig(indxTrigUp)==trgST_senton),1,'first'));
            trlUnitInfo(countUn,9)=indxTrigUp(find((indxTrigUp>tmpIndxHeadOn(iUn))&(detTrig(indxTrigUp)==trgST_endwait),1,'first'))-1;
            trlUnitInfo(countUn,10)=length(find((indxTrigUp>=trlUnitInfo(countUn,8))&(indxTrigUp<trlUnitInfo(countUn,9))&(detTrig(indxTrigUp)==trgST_senton)));
            
            trlUnitInfo(countUn,11)=indxTrigUp(find((indxTrigUp>tmpIndxHeadOn(iUn))&(detTrig(indxTrigUp)==trgST_that),1,'first'));
            trlUnitInfo(countUn,12)=indxTrigUp(find((indxTrigUp>tmpIndxHeadOn(iUn))&(detTrig(indxTrigUp)==trgST_RespBlock),1,'first'))-1;
            trlUnitInfo(countUn,13)=length(find((indxTrigUp>=trlUnitInfo(countUn,11))&(indxTrigUp<trlUnitInfo(countUn,12))));
            
            trlUnitInfo(countUn,14)=indxTrigUp(find((indxTrigUp>tmpIndxHeadOn(iUn))&(detTrig(indxTrigUp)==trgST_that),1,'first'));
            tmpCorIndx=indxTrigUp(find((indxTrigUp>tmpIndxHeadOn(iUn))&(detTrig(indxTrigUp)==trgST_OptCor),1,'first'));
            tmpUnCorIndx=indxTrigUp(find((indxTrigUp>tmpIndxHeadOn(iUn))&(detTrig(indxTrigUp)==trgST_OptUnCor),1,'first'));
            if tmpCorIndx>tmpUnCorIndx
                trlUnitInfo(countUn,15)=tmpUnCorIndx;
                trlUnitInfo(countUn,17)=tmpCorIndx;
                trlUnitInfo(countUn,18)=2;
            else
                trlUnitInfo(countUn,15)=tmpCorIndx;
                trlUnitInfo(countUn,17)=tmpUnCorIndx;
                trlUnitInfo(countUn,18)=1;
            end
            trlUnitInfo(countUn,16)=indxTrigUp(find((indxTrigUp>tmpIndxHeadOn(iUn))&(detTrig(indxTrigUp)==trgST_OptOR),1,'first'));
            
            
            
        else % Math
            
            trlUnitInfo(countUn,8)=indxTrigUp(find((indxTrigUp>tmpIndxHeadOn(iUn))&(detTrig(indxTrigUp)==trgMA_N_wordon),1,'first'));
            trlUnitInfo(countUn,9)=indxTrigDown(find((indxTrigDown>trlUnitInfo(countUn,8))&(detTrig(indxTrigDown-1)==trgMA_N_base),1,'first'));
            trlUnitInfo(countUn,10)=length(find((indxTrigUp>=trlUnitInfo(countUn,8))&(indxTrigUp<trlUnitInfo(countUn,9))&(detTrig(indxTrigUp)==trgMA_N_wordon)));
            
            
            if detTrig(trlUnitInfo(countUn,9))==trgMA_O_Cor1
                tmpCorOpt=1;
                tmpOptBase=trgMA_O_Cor1;
                tmp_O_wordon=trgMA_O_Cor1wordon;
            elseif detTrig(trlUnitInfo(countUn,9))==trgMA_O_Cor2
                tmpCorOpt=2;
                tmpOptBase=trgMA_O_Cor2;
                tmp_O_wordon=trgMA_O_Cor2wordon;
            else
                error('In Option INterval of Math Unit , the base trigger was not idenitified');
            end
            
            trlUnitInfo(countUn,11)=indxTrigUp(find((indxTrigUp>tmpIndxHeadOn(iUn))&(detTrig(indxTrigUp)==tmp_O_wordon),1,'first'));
            trlUnitInfo(countUn,12)=indxTrigUpDown(find((indxTrigUpDown>tmpIndxHeadOn(iUn))&(detTrig(indxTrigUpDown)==trgMA_RespBlock),1,'first'))-1;
            tmp_O_WordIndx=indxTrigUp(find((indxTrigUp>=trlUnitInfo(countUn,11))&(indxTrigUp<trlUnitInfo(countUn,12))));
            
            
            if length(tmp_O_WordIndx)~=4,
                error('There should be 4 triggers in option interval of Math')
            end
            trlUnitInfo(countUn,13)=length(tmp_O_WordIndx);
            trlUnitInfo(countUn,14)=tmp_O_WordIndx(1);
            trlUnitInfo(countUn,15)=tmp_O_WordIndx(2);
            trlUnitInfo(countUn,16)=tmp_O_WordIndx(3);
            trlUnitInfo(countUn,17)=tmp_O_WordIndx(4);
            trlUnitInfo(countUn,18)=tmpCorOpt;
            
        end
        
        trlUnitInfo(countUn,19)=trlUnitInfo(countUn,12)+1;
        if countUn<Ntotalun
            trlUnitInfo(countUn,20)=indxTrigDown(find((indxTrigDown>trlUnitInfo(countUn,15))&(detTrig(indxTrigDown)==0),1,'first'))-1;
        else
            indRespLast=indxTrigUpDown(find(ismember(detTrig(indxTrigUpDown),[trgMA_RespBlock  trgST_RespBlock]),1,'last'));
            trlUnitInfo(countUn,20)=indRespLast+durSampResp;
        end
        %-----------------------------
        % Check SUbject's Response
        isRespEarly=0;
        tmpRespIndx=indxRespUp(find((indxRespUp>trlUnitInfo(countUn,19))&(indxRespUp<=trlUnitInfo(countUn,20)),1,'first'));
        if ~isempty(tmpRespIndx),
            trlUnitInfo(countUn,21)=tmpRespIndx;
            if ((detResp(tmpRespIndx)==respLeft)&(trlUnitInfo(countUn,18)==1))|((detResp(tmpRespIndx)==respRight)&(trlUnitInfo(countUn,18)==2))
                trlUnitInfo(countUn,22)=1;
            else
                trlUnitInfo(countUn,22)=0;
            end
        else
            if any(ismember([respLeft,respRight],detResp(trlUnitInfo(countUn,19))))
                isRespEarly=1;
                if ((detResp(trlUnitInfo(countUn,19))==respLeft)&(trlUnitInfo(countUn,18)==1))|((detResp(trlUnitInfo(countUn,19))==respRight)&(trlUnitInfo(countUn,18)==2))
                    trlUnitInfo(countUn,22)=1;
                else
                    trlUnitInfo(countUn,22)=0;
                end
                tmpRespIndx=indxRespUp(find(indxRespUp<=trlUnitInfo(countUn,19),1,'last'));
                trlUnitInfo(countUn,21)=tmpRespIndx;
            end
        end
        
        trlUnitInfo(countUn,23)=isRespEarly;
        countUn=countUn+1;
    end
end
%=====================================================
%=====================================================
%=====================================================
%=====================================================

%% Create detailed trlSubUnitInfo for
%=====================================================
% Create detailed trlSubUnitInfo for
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
% 24. Narration Subunit Number in Narration interval  (Story Sentence or (Math Number or operation))
% 25. Narration Subunit Sample
% 26. Narration Subunit Type: 1:Language (Sentence or math operation word i.e. plus) 2: Math Number

countsubs=1;
trlSubUnitInfo=[];
for iUn=1:Ntotalun,
    
    tmpNsubs=trlUnitInfo(iUn,10);
    
    
    if trlUnitInfo(iUn,2)==1%Story
        tmpIndxSubs=find((indxTrigUp>=trlUnitInfo(iUn,8))&(indxTrigUp<trlUnitInfo(iUn,9))&(detTrig(indxTrigUp)==trgST_senton));
    else
        tmpIndxSubs=find((indxTrigUp>=trlUnitInfo(iUn,8))&(indxTrigUp<trlUnitInfo(iUn,9))&(detTrig(indxTrigUp)==trgMA_N_wordon));
    end
    for iSub=1:tmpNsubs
        trlSubUnitInfo(countsubs,1:23)=trlUnitInfo(iUn,:);
        trlSubUnitInfo(countsubs,24)=iSub;
        trlSubUnitInfo(countsubs,25)=indxTrigUp(tmpIndxSubs(iSub));
        trlSubUnitInfo(countsubs,26)=1;
        if trlUnitInfo(iUn,2)==2%Math
            if mod(iSub,2) % In MAth narration interval the tiggers with odd order 1,3,5,7 correspond to numbers while even to operations i.e. plus
                trlSubUnitInfo(countsubs,26)=2;
            end
        end
        countsubs=countsubs+1;
    end
    
end

%-------------------------------------------------------
%=====================================================
%=====================================================
%=====================================================
%=====================================================
%=====================================================
%=====================================================

%% Create detailed trlAllEventInfo for
%=====================================================
% Create detailed trlAllInfo for
% Columns:
%----------------------------------------------
trlColDescr_AllEvent={'1. Block Number within Run'
'2. Unit Type : 1.Story 2.Math'
'3. Unit Number within Run'
'4. Total Number of units (N of Stories or N of Math problems ) in same Run'
'5. Unit Number within Block'
'6. Total Number of units (N of Stories or N of Math problems ) in same Block'
'7. Attribute1: For story this is the story number. For Math is the difficulty level'
'8. Unit Narration interval Start Sample - Start with the onset of the  first word trigger or the beginning of the first sentence'
'9. Unit Narration interval End Sample'
'10. N subunits within Narration interval'
'11. Unit Option interval Start Sample'
'12. Unit Option interval End Sample'
'13. N subunits within Option interval'
'14 Correct Option- 1 or 2'
'15 Unit Response interval start sample'
'16 Unit Response interval end sample'
'17 Unit Response sample'
'18 is Response Correct'
'19 is Response Early'
'20 Event Sample'
'21 Event Type - 20.Math Narration number word,\n 21.Math Narration operand word,\n 22.Math Option Intro Word,\n 23.Math Option 1 Word\n  24. Math Option OR Word\n 25.Math Option 2 Word\n 10.Story Sentence,\n 12.Story Option Intro\n 13.Math Option 1 Word\n  14. Math Option OR Word\n 15.Math Option 2 Word'
'22. Narration Event Number in Narration interval  (Applies only to Sentence or number or operation in Narration interval)'};

trlAllEventInfo=[];
Nsubs=size(trlSubUnitInfo,1);
countNums=1;
for iSub=1:Nsubs
    
    trlAllEventInfo(countNums,1:13)=trlSubUnitInfo(iSub,1:13);
    trlAllEventInfo(countNums,14:19)=trlSubUnitInfo(iSub,18:23);
    trlAllEventInfo(countNums,20)=trlSubUnitInfo(iSub,25);
    trlAllEventInfo(countNums,22)=trlSubUnitInfo(iSub,24);
    
    tmpNsubs=trlSubUnitInfo(iSub,10);
    
    if trlSubUnitInfo(iSub,2)==1, %Story
        trlAllEventInfo(countNums,21)=10;
        countNums=countNums+1;
        if trlSubUnitInfo(iSub,24)== tmpNsubs
            for iOpt=1:4
                trlAllEventInfo(countNums,1:19)=trlAllEventInfo(countNums-1,1:19);
                trlAllEventInfo(countNums,20)=trlSubUnitInfo(iSub,13+iOpt);
                trlAllEventInfo(countNums,21)=11+iOpt;
                trlAllEventInfo(countNums,22)=nan;
                countNums=countNums+1;
            end
        end
        
    elseif trlSubUnitInfo(iSub,2)==2, %Math
        if trlSubUnitInfo(iSub,26)==1, %operand
            trlAllEventInfo(countNums,21)=21;
        elseif trlSubUnitInfo(iSub,26)==2, %number
            trlAllEventInfo(countNums,21)=20;
        end
        countNums=countNums+1;
        
        if trlSubUnitInfo(iSub,24)== tmpNsubs
            for iOpt=1:4
                trlAllEventInfo(countNums,1:19)=trlAllEventInfo(countNums-1,1:19);
                trlAllEventInfo(countNums,20)=trlSubUnitInfo(iSub,13+iOpt);
                trlAllEventInfo(countNums,21)=21+iOpt;
                trlAllEventInfo(countNums,22)=nan;
                countNums=countNums+1;
            end
        end
        
    end
    
end
%---------------------------------------
% figure;hold on;
% plot(detTrig);
% plot(trlAllEventInfo(:,20),detTrig(trlAllEventInfo(:,20)),'.g');

%=====================================================
%=====================================================

%% Create detailed trlNumberInfo for all numbers
% 1. Block Number within Run
% 2. Unit Sequence Number within Block (unit=story or math problem interval  )
% 3. Total Number of units (N of Stories or N of Math problems ) in same Block
% 4. Attribute1: For story this is the story number. For Math is the
%    difficulty level
% 5. Number Sequence Number in Relative interval  - 1,2,3,4 for
%    Narration numbers and 1,2 for Option Numbers (Use column 7 to distinguish between Narration and Option)
% 6. Number Onset Sample
% 7. Number Type: 0.Number in math narration
%                           11:number in option1 Correct
%                           10:number in option1 UnCorrect
%                           21:number in option2 Correct
%                           20:number in option2 UnCorrect
% 8. is Response Correct % This refers to what the response was at the end
%     of the math operation


trlNumberInfo=[];
Nsubs=size(trlSubUnitInfo,1);
countNums=1;
for iSub=1:Nsubs
    if trlSubUnitInfo(iSub,26)==2
        trlNumberInfo(countNums,1:4)=trlSubUnitInfo(iSub,[1 5 6 7]);
        trlNumberInfo(countNums,5)=floor(trlSubUnitInfo(iSub,26)./2)+1;
        trlNumberInfo(countNums,6)=trlSubUnitInfo(iSub,25);
        trlNumberInfo(countNums,7)=0; %Number in narration;
        trlNumberInfo(countNums,8)=trlSubUnitInfo(iSub,22);
        countNums=countNums+1;
        
        tmpNsubs=trlSubUnitInfo(iSub,10);
        if trlSubUnitInfo(iSub,24)== tmpNsubs;
            for iOpt=1:2
                trlNumberInfo(countNums,1:4)=trlNumberInfo(countNums-1,1:4);
                trlNumberInfo(countNums,5)=iOpt;
                if iOpt==1
                    trlNumberInfo(countNums,6)=trlSubUnitInfo(iSub,15);
                    trlNumberInfo(countNums,7)=10+abs(trlSubUnitInfo(iSub,18)-2);
                    trlNumberInfo(countNums,8)=trlSubUnitInfo(iSub,22);
                else
                    trlNumberInfo(countNums,6)=trlSubUnitInfo(iSub,17);
                    trlNumberInfo(countNums,7)=20+trlSubUnitInfo(iSub,18)-1;
                    trlNumberInfo(countNums,8)=trlSubUnitInfo(iSub,22);
                end
                
                countNums=countNums+1;
            end
        end
        
    end
end
%---------------------------------------
% figure;hold on;
% plot(detTrig);
% plot(trlNumberInfo(:,6),detTrig(trlNumberInfo(:,6)),'.g');
% %plot(trlNumberInfo(:,3),detTrig(trlNumberInfo(:,6)),'.r');
% figure;hold on;
% plot(detTrig);
% plot(trlUnitInfo(:,8),detTrig(trlUnitInfo(:,8)),'.g');
%=====================================================
%=====================================================
%% Create trlSentInfo
% Each trials corresponds to a Narration Sentence , either in story or in
% math (option blocks not included)
% Columns:
%----------------------------------------------
trlColDescr_SentInfo={'1. Block Number within Run'
'2. Unit Type : 1.Story 2.Math'
'3. Unit Number within Run'
'4. Total Number of units (N of Stories or N of Math problems ) in same Run'
'5. Unit Number within Block'
'6. Total Number of units (N of Stories or N of Math problems ) in same Block'
'7. Attribute1: For story this is the story number. For Math is the  difficulty level'
'8. N sentences within Narration interval (always one for math as in math one narration interval corresponds to one math sentence)'
'9. is Response Correct (For Story this refers to the response at the very end of the sentence)'
'10. is Response Early'
'11. Narration Sentence Number in Narration interval  (For math always equal to one)'
'12. Narration Sentence Start Sample'
'13. Narration Sentence End Sample'};
%-----------------------------------------
trlSentInfo=[];
countSent=1;
for iCase=1:size(trlSubUnitInfo,1)
    curUnitRunNum= trlSubUnitInfo(iCase,3);
    
    if trlSubUnitInfo(iCase,2)==1,
        
        trlSentInfo(countSent,1:7)=trlSubUnitInfo(iCase,1:7);
        trlSentInfo(countSent,8)= trlSubUnitInfo(iCase,10);
        trlSentInfo(countSent,9:10)=  trlSubUnitInfo(iCase,22:23);
        trlSentInfo(countSent,11:12)=  trlSubUnitInfo(iCase,24:25);
        if iCase<size(trlSubUnitInfo,1),
            if   trlSubUnitInfo(iCase+1,2)==1,
                if   trlSubUnitInfo(iCase,24)<trlSubUnitInfo(iCase,10),
                    trlSentInfo(countSent,13)=  trlSubUnitInfo(iCase+1,25)-1;
                elseif   trlSubUnitInfo(iCase,24)==trlSubUnitInfo(iCase,10),
                    trlSentInfo(countSent,13)= trlSubUnitInfo(iCase,9);
                end
            else
                trlSentInfo(countSent,13)= trlSubUnitInfo(iCase,9);
            end
        else
            trlSentInfo(countSent,13)= trlSubUnitInfo(iCase,9);
        end
        countSent=countSent+1;
        
        
    elseif trlSubUnitInfo(iCase,2)==2,
        if (curUnitRunNum~=prevUnitRunNum)
            trlSentInfo(countSent,1:7)=trlSubUnitInfo(iCase,1:7);
            trlSentInfo(countSent,8)= 1;
            trlSentInfo(countSent,9:10)=  trlSubUnitInfo(iCase,22:23);
            trlSentInfo(countSent,11)=  1;
            trlSentInfo(countSent,12:13)=trlSubUnitInfo(iCase,8:9);
            countSent=countSent+1;
        end
        
    end
    
    
    prevUnitRunNum=curUnitRunNum;
end

%=====================================================
%%
if cutmode==1       %All EVents
    trialinfo=trlAllEventInfo;
    Ntrls=size(trialinfo,1);
    trl=[trialinfo(:,20)-prestimSamples trialinfo(:,20)+poststimSamples  -repmat(prestimSamples,Ntrls,1) trialinfo];
    trlInfoColDescr=trlColDescr_AllEvent;
elseif cutmode==2   %Response (Fixed trial length)
    
    hasResp=~isnan(trlUnitInfo(:,21));
    
    trialinfo=trlUnitInfo(hasResp,:);
    Ntrls=size(trialinfo,1);
    trl=[trialinfo(:,21)-prestimSamples trialinfo(:,21)+poststimSamples  -repmat(prestimSamples,Ntrls,1) trialinfo];
    trlInfoColDescr=trlColDescr_UnitInfo;
elseif cutmode==3   %Sentences (Story Sentences or Math question sentences(not including option interval)) (Variable trial length)
    
    trialinfo=trlSentInfo;
    Ntrls=size(trialinfo,1);
    trl=[trialinfo(:,12)-prestimSamples trialinfo(:,13)+poststimSamples  -repmat(prestimSamples,Ntrls,1) trialinfo];
    trlInfoColDescr=trlColDescr_SentInfo;
    
elseif cutmode==4   %Units (Stories or Math Problems including Option Interval) (Variable trial length)
    
    trialinfo=trlUnitInfo;
    Ntrls=size(trialinfo,1);
    trl=[trialinfo(:,8)-prestimSamples trialinfo(:,19)-1+poststimSamples  -repmat(prestimSamples,Ntrls,1) trialinfo];
    trlInfoColDescr=trlColDescr_UnitInfo;
end

%----------------------------------
%===========================================================
%% TRIAL SUMMARY INFO
%===========================================================
%===========================================================
%--- Compute summary of trial information to print is summary file
matNnumsInMath=[];
countMathUnits=1;
Nunits=length(unique(trlSubUnitInfo(:,3)));
for iUnit=1:Nunits,
    [indx1]=find(trlSubUnitInfo(:,3)==iUnit,1,'first');
    if trlSubUnitInfo(indx1,2)==2, %Math problem
      [indx2]=find(trlSubUnitInfo(:,3)==iUnit);
      tmpCaseInfo=trlSubUnitInfo(indx2,26);
      matNnumsInMath(countMathUnits)=length(find(tmpCaseInfo==2));
      countMathUnits=countMathUnits+1;
    end
end
%===================================================
trialSummary=[];
trialSummary.Nblocks=size(trlBlockInfo,1);
trialSummary.Nblocks_ST=length(find(trlBlockInfo(:,1)==1));
trialSummary.Nblocks_MA=length(find(trlBlockInfo(:,1)==2));
trialSummary.Nstories=length(find(trlUnitInfo(:,2)==1));
trialSummary.Nstorysent=length(find(trlAllEventInfo(:,21)==10));
trialSummary.Nmathprobs=length(find(trlUnitInfo(:,2)==2));
trialSummary.NnumsInMath=unique(matNnumsInMath); % This show How many number are used in math. It seems the number is fixed. This is to verify this


