opengl software;

% ensure that the time and date of execution are not stored in the provenance information
global ft_default
ft_default.trackcallinfo = 'no';

% allow the user to specify the path where additional data is present, e.g. the channel layout or anatomy files
%if exist('path', 'var')
%    addpath(path)
%end

%%
%-------------------------------------------------------------------------------
% The subject list file should be an ascii file with 5 columns and each row corresponding to a single subject
% COLUMNS:
% 1:   Subject ID
% 2:   0 or 1. flag signifying if this subject has good RESTING STATE data from ALL 3 scans
% 3:   0 or 1. flag signifying if this subject has good MOTOR TASK data from BOTH scans
% 4:   0 or 1. flag signifying if this subject has good WORKING MEMORY data from BOTH scans
% 5:   0 or 1. flag signifying if this subject has good STORY MATH data from BOTH scans
if ~exist( 'subjlistfile','var')
    error('subjlistfile should be specified')
end
subjlistmat=load(subjlistfile);
if size(subjlistmat,2)~=5,
    error('subjlistfile should have 5 columns');
end

%%
if ~exist( 'datamaindir','var')
    error('datamaindir should be specified. This is the main data directory where raw and processed data for all subjects is saved. i.e. /HCP/scratch/meg/intradb/archive1/HCP_Phase2/arc001/');
end

if ~exist( 'avgspecdir','var')
    error('avgspecdir should be specified. This is the directory where the files with the averaged spectra will be saved for all subjects is saved. i.e. /HCP/scratch/meg/frankenstein/workspace/hcp/findpeaks/analysis/');
end
% Files with the averaged spectra are saved with the naming
%[casespectra_',subjID,'_Motort.mat'];
%[casespectra_',subjID,'_Wrkmem.mat'];
%[casespectra_',subjID,'_StoryM.mat'];
%[casespectra_',subjID,'_Restin.mat'];

%%
%--------------------------------------------------------------------


numListMat=subjlistmat;
numList=numListMat(:,1);
subjList=cellfun(@(x) num2str(x) , num2cell(numList),'Uniformoutput',false);

%%
%==========================================================================
% DEFINE SENSOR GROUPS
%subjList={'177746'};
Nsubj=length(subjList);


chanLists=[];
chanLists.Occ_L ={'A104', 'A238', 'A202', 'A220', 'A160', 'A219', 'A186', 'A105', 'A234', 'A215', 'A235', 'A181', 'A137', 'A135', 'A200', 'A102', 'A183', 'A199', 'A203', 'A103', 'A163', 'A236', 'A161', 'A218', 'A201', 'A164', 'A217', 'A136', 'A184', 'A237', 'A182', 'A185', 'A162', 'A133', 'A134', 'A216'}';
chanLists.Occ_R = {'A241', 'A138', 'A207', 'A224', 'A239', 'A165', 'A204', 'A242', 'A168', 'A223', 'A222', 'A188', 'A189', 'A167', 'A240', 'A206', 'A107', 'A139', 'A166', 'A225', 'A205', 'A243', 'A190', 'A208', 'A140', 'A106', 'A187', 'A221'}';
chanLists.Occ_All=[chanLists.Occ_L; chanLists.Occ_R];
chanLists.SeMoAud_L = {'A22', 'A214', 'A71', 'A26', 'A9', 'A8', 'A95', 'A28', 'A233', 'A98', 'A25', 'A70', 'A72', 'A11', 'A47', 'A3', 'A155', 'A10', 'A127', 'A67', 'A12', 'A231', 'A45', 'A180', 'A99', 'A230', 'A24', 'A66', 'A42', 'A96', 'A27', 'A130', 'A100', 'A43', 'A132', 'A21', 'A49', 'A156', 'A128', 'A68', 'A159', 'A4', 'A6', 'A74', 'A232', 'A69', 'A157', 'A97', 'A101', 'A40', 'A179', 'A73', 'A129', 'A131', 'A198', 'A197', 'A46', 'A41', 'A7', 'A23', 'A48', 'A196', 'A158', 'A44'}';
chanLists.SeMoAud_R = {'A114', 'A16', 'A35', 'A170', 'A112', 'A82', 'A13', 'A115', 'A78', 'A31', 'A245', 'A76', 'A50', 'A57', 'A56', 'A80', 'A210', 'A143', 'A113', 'A84', 'A55', 'A32', 'A146', 'A79', 'A54', 'A145', 'A14', 'A15', 'A30', 'A109', 'A172', 'A81', 'A171', 'A173', 'A29', 'A33', 'A147', 'A52', 'A142', 'A211', 'A53', 'A192', 'A226', 'A51', 'A77', 'A83', 'A34', 'A17', 'A18', 'A144', 'A209', 'A110', 'A111', 'A244'}';
chanLists.SeMoAud_All=[chanLists.SeMoAud_L; chanLists.SeMoAud_R];
chanLists.Front_All = {'A93', 'A39', 'A125', 'A175', 'A228', 'A64', 'A177', 'A63', 'A194', 'A176', 'A38', 'A91', 'A86', 'A116', 'A151', 'A120', 'A122', 'A62', 'A60', 'A88', 'A121', 'A61', 'A193', 'A150', 'A227', 'A59', 'A195', 'A124', 'A123', 'A153', 'A178', 'A117', 'A148', 'A87', 'A89', 'A119', 'A92', 'A90', 'A154', 'A149', 'A118', 'A152'}';
%==========================================================================


%((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((
%((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((
%((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((
%((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((

%             START AVERAGING THE SPECTRA IN SENSOR GROUPS

%))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
%))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
%))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
%))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))

%%
%########################################################################################################
%########################################################################################################
%########################################################################################################
% MOTOR TASK

% A variable called casespectra is saved whis contains the speactra for different sensor subsets and different trial cases.
% COLUMN 1:   Member flag ( 1:LH    2:RH    3:LF    4:RF)
% COLUMN 2:   Trial Lock data group ( 1: TMEG    2:TFLA );
% COLUMN 3:   Sensor group   (1: Occ/Posterior  2:SensoryMotor/Auditory Left 3:SensoryMotor/Auditory Right  4:Frontal)
% COLUMN 4:   Spectrum Type (1: Pre-stim(Baseline) 2:Post-stim  3: relative Change in Post-stim vs Pre-stim)
% REMANING COLUMNS :  Spectrum per frequency

columnDescr={
    'COLUMN 1:   Member flag ( 1: LH    2: RH    3: LF     4: RF )'
    'COLUMN 2:   Trial Lock data group ( 1: TMEG    2:TFLA )'
    'COLUMN 3:   Sensor group   (1: Occ/Posterior  2:SensoryMotor/Auditory Left 3:SensoryMotor/Auditory Right  4:Frontal)'
    'COLUMN 4:   Spectrum Type (1: Pre-stim(Baseline) 2:Post-stim  3: relative Change in Post-stim vs Pre-stim)'
    'REMANING COLUMNS :  Spectrum per frequency'};

%----------------------------------------------------------------------------
caseMnems={'LH'
    'RH'
    'LF'
    'RF'};


subjLogFile=[avgspecdir,'SubjectsLog_Motort_tfavg.txt'];
logfid=fopen(subjLogFile,'w+');


outDir=avgspecdir;

for iSubj=1:Nsubj,
    %--------------
    subjID=subjList{iSubj};
    isRESTIN=subjlistmat(iSubj,2);
    isMOTOR=subjlistmat(iSubj,3);
    isWORKMEM=subjlistmat(iSubj,4);
    isSTORYM=subjlistmat(iSubj,5);
    %--------------
    if isMOTOR,
        
        
        
        disp(['===========================================================================================']);
        disp(['MOTOR - ',num2str(iSubj),' : ',subjID]);
        disp(['@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@']);
        %analDir=['/HCP/scratch/meg/intradb/archive1/HCP_Phase2/arc001/',subjID,'_MEG/RESOURCES/analysis/'];
        
        dataDir=[datamaindir,subjID,'_MEG/RESOURCES/analysis/'];
        outputFile=[outDir,'casespectra_',subjID,'_Motort.mat'];
        
        countEntry=1;
        totalSpecCases=[];
        
        isAnyFileFound=0;
        for iCase=1:4,
            
            curMnem=caseMnems{iCase};
            
            datafile1=[dataDir,subjID,'_MEG_Motort_tfavg_[LM-TEMG-',curMnem,']_[MODE-mag].mat'];
            datafile2=[dataDir,subjID,'_MEG_Motort_tfavg_[LM-TFLA-',curMnem,']_[MODE-mag].mat'];
            
            isFile1=exist(datafile1,'file');
            isFile2=exist(datafile2,'file');
            
            
            logChar=[subjID,'   ',num2str(iCase),'   ',num2str(isFile1),'   ',num2str(isFile2)];
            fprintf(logfid,'%s\n',logChar);
            if (isFile1==0)|(isFile2==0),
                logChar=['OOOPS - tfavg files not found for Motor task of subject: ',subjID];
                fprintf(logfid,'%s\n',logChar);
                warning(['OOOPS - tfavg files not found for Motor task of subject: ',subjID]);
                continue;
            else
                isAnyFileFound=1;
            end
            
            load(datafile1);
            
            %===================================================================================================
            tmpdata_post_occ=ft_selectdata(data,'channel',chanLists.Occ_All,'avgoverchan','yes','toilim',[0 0.8],'avgovertime','yes');
            tmpdata_post_semoaud_L=ft_selectdata(data,'channel',chanLists.SeMoAud_L,'avgoverchan','yes','toilim',[0 0.8],'avgovertime','yes');
            tmpdata_post_semoaud_R=ft_selectdata(data,'channel',chanLists.SeMoAud_R,'avgoverchan','yes','toilim',[0 0.8],'avgovertime','yes');
            tmpdata_post_front=ft_selectdata(data,'channel',chanLists.Front_All,'avgoverchan','yes','toilim',[0 0.8],'avgovertime','yes');
            %-----------------------------------
            tmpdata_pre_occ=ft_selectdata(data,'channel',chanLists.Occ_All,'avgoverchan','yes','toilim',[-0.8 0],'avgovertime','yes');
            tmpdata_pre_semoaud_L=ft_selectdata(data,'channel',chanLists.SeMoAud_L,'avgoverchan','yes','toilim',[-0.8 0],'avgovertime','yes');
            tmpdata_pre_semoaud_R=ft_selectdata(data,'channel',chanLists.SeMoAud_R,'avgoverchan','yes','toilim',[-0.8 0],'avgovertime','yes');
            tmpdata_pre_front=ft_selectdata(data,'channel',chanLists.Front_All,'avgoverchan','yes','toilim',[-0.8 0],'avgovertime','yes');
            %===================================================================================================
            totalSpecCases(countEntry,:)=[iCase 1 1 1 tmpdata_pre_occ.powspctrm];countEntry=countEntry+1;
            totalSpecCases(countEntry,:)=[iCase 1 2 1 tmpdata_pre_semoaud_L.powspctrm];countEntry=countEntry+1;
            totalSpecCases(countEntry,:)=[iCase 1 3 1 tmpdata_pre_semoaud_R.powspctrm];countEntry=countEntry+1;
            totalSpecCases(countEntry,:)=[iCase 1 4 1 tmpdata_pre_front.powspctrm];countEntry=countEntry+1;
            
            totalSpecCases(countEntry,:)=[iCase 1 1 2 tmpdata_post_occ.powspctrm];countEntry=countEntry+1;
            totalSpecCases(countEntry,:)=[iCase 1 2 2 tmpdata_post_semoaud_L.powspctrm];countEntry=countEntry+1;
            totalSpecCases(countEntry,:)=[iCase 1 3 2 tmpdata_post_semoaud_R.powspctrm];countEntry=countEntry+1;
            totalSpecCases(countEntry,:)=[iCase 1 4 2 tmpdata_post_front.powspctrm];countEntry=countEntry+1;
            
            totalSpecCases(countEntry,:)=[iCase 1 1 3 -1+tmpdata_post_occ.powspctrm./tmpdata_pre_occ.powspctrm];countEntry=countEntry+1;
            totalSpecCases(countEntry,:)=[iCase 1 2 3 -1+tmpdata_post_semoaud_L.powspctrm./tmpdata_pre_semoaud_L.powspctrm];countEntry=countEntry+1;
            totalSpecCases(countEntry,:)=[iCase 1 3 3 -1+tmpdata_post_semoaud_R.powspctrm./tmpdata_pre_semoaud_R.powspctrm];countEntry=countEntry+1;
            totalSpecCases(countEntry,:)=[iCase 1 4 3 -1+tmpdata_post_front.powspctrm./tmpdata_pre_front.powspctrm];countEntry=countEntry+1;
            
            allFreqs=tmpdata_post_occ.freq;
            %%
            
            load(datafile2);
            
            %===================================================================================================
            tmpdata_post_occ=ft_selectdata(data,'channel',chanLists.Occ_All,'avgoverchan','yes','toilim',[0 0.8],'avgovertime','yes');
            tmpdata_post_semoaud_L=ft_selectdata(data,'channel',chanLists.SeMoAud_L,'avgoverchan','yes','toilim',[0 0.8],'avgovertime','yes');
            tmpdata_post_semoaud_R=ft_selectdata(data,'channel',chanLists.SeMoAud_R,'avgoverchan','yes','toilim',[0 0.8],'avgovertime','yes');
            tmpdata_post_front=ft_selectdata(data,'channel',chanLists.Front_All,'avgoverchan','yes','toilim',[0 0.8],'avgovertime','yes');
            %-----------------------------------
            tmpdata_pre_occ=ft_selectdata(data,'channel',chanLists.Occ_All,'avgoverchan','yes','toilim',[-0.8 0],'avgovertime','yes');
            tmpdata_pre_semoaud_L=ft_selectdata(data,'channel',chanLists.SeMoAud_L,'avgoverchan','yes','toilim',[-0.8 0],'avgovertime','yes');
            tmpdata_pre_semoaud_R=ft_selectdata(data,'channel',chanLists.SeMoAud_R,'avgoverchan','yes','toilim',[-0.8 0],'avgovertime','yes');
            tmpdata_pre_front=ft_selectdata(data,'channel',chanLists.Front_All,'avgoverchan','yes','toilim',[-0.8 0],'avgovertime','yes');
            %===================================================================================================
            totalSpecCases(countEntry,:)=[iCase 2 1 1 tmpdata_pre_occ.powspctrm];countEntry=countEntry+1;
            totalSpecCases(countEntry,:)=[iCase 2 2 1 tmpdata_pre_semoaud_L.powspctrm];countEntry=countEntry+1;
            totalSpecCases(countEntry,:)=[iCase 2 3 1 tmpdata_pre_semoaud_R.powspctrm];countEntry=countEntry+1;
            totalSpecCases(countEntry,:)=[iCase 2 4 1 tmpdata_pre_front.powspctrm];countEntry=countEntry+1;
            
            totalSpecCases(countEntry,:)=[iCase 2 1 2 tmpdata_post_occ.powspctrm];countEntry=countEntry+1;
            totalSpecCases(countEntry,:)=[iCase 2 2 2 tmpdata_post_semoaud_L.powspctrm];countEntry=countEntry+1;
            totalSpecCases(countEntry,:)=[iCase 2 3 2 tmpdata_post_semoaud_R.powspctrm];countEntry=countEntry+1;
            totalSpecCases(countEntry,:)=[iCase 2 4 2 tmpdata_post_front.powspctrm];countEntry=countEntry+1;
            
            totalSpecCases(countEntry,:)=[iCase 2 1 3 -1+tmpdata_post_occ.powspctrm./tmpdata_pre_occ.powspctrm];countEntry=countEntry+1;
            totalSpecCases(countEntry,:)=[iCase 2 2 3 -1+tmpdata_post_semoaud_L.powspctrm./tmpdata_pre_semoaud_L.powspctrm];countEntry=countEntry+1;
            totalSpecCases(countEntry,:)=[iCase 2 3 3 -1+tmpdata_post_semoaud_R.powspctrm./tmpdata_pre_semoaud_R.powspctrm];countEntry=countEntry+1;
            totalSpecCases(countEntry,:)=[iCase 2 4 3 -1+tmpdata_post_front.powspctrm./tmpdata_pre_front.powspctrm];countEntry=countEntry+1;
            
            allFreqs=tmpdata_post_occ.freq;
            %########################################################################################################
        end
        
        if isAnyFileFound
            casespectra=totalSpecCases;
            allfreqs=allFreqs;
            save(outputFile,'casespectra','columnDescr','allfreqs');
        end
        
        
    end
end

fclose(logfid);


%%
%########################################################################################################
%########################################################################################################
%########################################################################################################

% WORKING MEMORY

% A variable called casespectra is saved whis contains the speactra for different sensor subsets and different trial cases.
% 'COLUMN 1:   Member flag ( 1: 0-Back    2: 2-Back )'
% 'COLUMN 2:   Trial Lock data group ( 1: TIM    2:TRESP )'
% COLUMN 3:   Sensor group   (1: Occ/Posterior  2:SensoryMotor/Auditory Left 3:SensoryMotor/Auditory Right  4:Frontal)
% COLUMN 4:   Spectrum Type (1: Pre-stim(Baseline) 2:Post-stim  3: relative Change in Post-stim vs Pre-stim)
% REMANING COLUMNS :  Spectrum per frequency

columnDescr={
    'COLUMN 1:   Member flag ( 1: 0-Back    2: 2-Back )'
    'COLUMN 2:   Trial Lock data group ( 1: TIM    2:TRESP )'
    'COLUMN 3:   Sensor group   (1: Occ/Posterior  2:SensoryMotor/Auditory Left 3:SensoryMotor/Auditory Right  4:Frontal)'
    'COLUMN 4:   Spectrum Type (1: Pre-stim(Baseline) 2:Post-stim  3: relative Change in Post-stim vs Pre-stim)'
    'REMANING COLUMNS :  Spectrum per frequency'};

%----------------------------------------------------------------------------
caseMnems={'0B'
    '2B'};


subjLogFile=[avgspecdir,'SubjectsLog_Wrkmem_tfavg.txt'];
logfid=fopen(subjLogFile,'w+');

outDir=avgspecdir;

for iSubj=1:Nsubj,
    %--------------
    subjID=subjList{iSubj};
    isRESTIN=subjlistmat(iSubj,2);
    isMOTOR=subjlistmat(iSubj,3);
    isWORKMEM=subjlistmat(iSubj,4);
    isSTORYM=subjlistmat(iSubj,5);
    %--------------
    if isWORKMEM,
        
        
        
        disp(['===========================================================================================']);
        disp(['WORKING MEM - ',num2str(iSubj),' : ',subjID]);
        disp(['@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@']);
        %analDir=['/HCP/scratch/meg/intradb/archive1/HCP_Phase2/arc001/',subjID,'_MEG/RESOURCES/analysis/'];
        
        dataDir=[datamaindir,subjID,'_MEG/RESOURCES/analysis/'];
        outputFile=[outDir,'casespectra_',subjID,'_Wrkmem.mat'];
        
        countEntry=1;
        totalSpecCases=[];
        
        isAnyFileFound=0;
        for iCase=1:2,
            
            disp('WM : STAGE  1');
            curMnem=caseMnems{iCase};
            
            datafile1=[dataDir,subjID,'_MEG_Wrkmem_tfavg_[LM-TIM-',curMnem,']_[MODE-mag].mat'];
            datafile2=[dataDir,subjID,'_MEG_Wrkmem_tfavg_[LM-TRESP-',curMnem,']_[MODE-mag].mat'];
            
            isFile1=exist(datafile1,'file');
            isFile2=exist(datafile2,'file');
            
            
            logChar=[subjID,'   ',num2str(iCase),'   ',num2str(isFile1),'   ',num2str(isFile2)];
            fprintf(logfid,'%s\n',logChar);
            if (isFile1==0)|(isFile2==0),
                logChar=['OOOPS - tfavg files not found for Workmem task of subject: ',subjID];
                fprintf(logfid,'%s\n',logChar);
                warning(['OOOPS - tfavg files not found for Workmem task of subject: ',subjID]);
                continue;
            else
                isAnyFileFound=1;
            end
            
            load(datafile1);
            disp('WM : STAGE  2');
            %===================================================================================================
            tmpdata_post_occ=ft_selectdata(data,'channel',chanLists.Occ_All,'avgoverchan','yes','toilim',[0 0.8],'avgovertime','yes');
            tmpdata_post_semoaud_L=ft_selectdata(data,'channel',chanLists.SeMoAud_L,'avgoverchan','yes','toilim',[0 0.8],'avgovertime','yes');
            tmpdata_post_semoaud_R=ft_selectdata(data,'channel',chanLists.SeMoAud_R,'avgoverchan','yes','toilim',[0 0.8],'avgovertime','yes');
            tmpdata_post_front=ft_selectdata(data,'channel',chanLists.Front_All,'avgoverchan','yes','toilim',[0 0.8],'avgovertime','yes');
            %-----------------------------------
            tmpdata_pre_occ=ft_selectdata(data,'channel',chanLists.Occ_All,'avgoverchan','yes','toilim',[-0.8 0],'avgovertime','yes');
            tmpdata_pre_semoaud_L=ft_selectdata(data,'channel',chanLists.SeMoAud_L,'avgoverchan','yes','toilim',[-0.8 0],'avgovertime','yes');
            tmpdata_pre_semoaud_R=ft_selectdata(data,'channel',chanLists.SeMoAud_R,'avgoverchan','yes','toilim',[-0.8 0],'avgovertime','yes');
            tmpdata_pre_front=ft_selectdata(data,'channel',chanLists.Front_All,'avgoverchan','yes','toilim',[-0.8 0],'avgovertime','yes');
            %===================================================================================================
            disp('WM : STAGE  3');
            totalSpecCases(countEntry,:)=[iCase 1 1 1 tmpdata_pre_occ.powspctrm];countEntry=countEntry+1;
            totalSpecCases(countEntry,:)=[iCase 1 2 1 tmpdata_pre_semoaud_L.powspctrm];countEntry=countEntry+1;
            totalSpecCases(countEntry,:)=[iCase 1 3 1 tmpdata_pre_semoaud_R.powspctrm];countEntry=countEntry+1;
            totalSpecCases(countEntry,:)=[iCase 1 4 1 tmpdata_pre_front.powspctrm];countEntry=countEntry+1;
            
            totalSpecCases(countEntry,:)=[iCase 1 1 2 tmpdata_post_occ.powspctrm];countEntry=countEntry+1;
            totalSpecCases(countEntry,:)=[iCase 1 2 2 tmpdata_post_semoaud_L.powspctrm];countEntry=countEntry+1;
            totalSpecCases(countEntry,:)=[iCase 1 3 2 tmpdata_post_semoaud_R.powspctrm];countEntry=countEntry+1;
            totalSpecCases(countEntry,:)=[iCase 1 4 2 tmpdata_post_front.powspctrm];countEntry=countEntry+1;
            
            totalSpecCases(countEntry,:)=[iCase 1 1 3 -1+tmpdata_post_occ.powspctrm./tmpdata_pre_occ.powspctrm];countEntry=countEntry+1;
            totalSpecCases(countEntry,:)=[iCase 1 2 3 -1+tmpdata_post_semoaud_L.powspctrm./tmpdata_pre_semoaud_L.powspctrm];countEntry=countEntry+1;
            totalSpecCases(countEntry,:)=[iCase 1 3 3 -1+tmpdata_post_semoaud_R.powspctrm./tmpdata_pre_semoaud_R.powspctrm];countEntry=countEntry+1;
            totalSpecCases(countEntry,:)=[iCase 1 4 3 -1+tmpdata_post_front.powspctrm./tmpdata_pre_front.powspctrm];countEntry=countEntry+1;
            
            allFreqs=tmpdata_post_occ.freq;
            %%
            disp('WM : STAGE  4');
            load(datafile2);
            disp('WM : STAGE  5');
            %===================================================================================================
            tmpdata_post_occ=ft_selectdata(data,'channel',chanLists.Occ_All,'avgoverchan','yes','toilim',[0 0.8],'avgovertime','yes');
            tmpdata_post_semoaud_L=ft_selectdata(data,'channel',chanLists.SeMoAud_L,'avgoverchan','yes','toilim',[0 0.8],'avgovertime','yes');
            tmpdata_post_semoaud_R=ft_selectdata(data,'channel',chanLists.SeMoAud_R,'avgoverchan','yes','toilim',[0 0.8],'avgovertime','yes');
            tmpdata_post_front=ft_selectdata(data,'channel',chanLists.Front_All,'avgoverchan','yes','toilim',[0 0.8],'avgovertime','yes');
            %-----------------------------------
            tmpdata_pre_occ=ft_selectdata(data,'channel',chanLists.Occ_All,'avgoverchan','yes','toilim',[-0.8 0],'avgovertime','yes');
            tmpdata_pre_semoaud_L=ft_selectdata(data,'channel',chanLists.SeMoAud_L,'avgoverchan','yes','toilim',[-0.8 0],'avgovertime','yes');
            tmpdata_pre_semoaud_R=ft_selectdata(data,'channel',chanLists.SeMoAud_R,'avgoverchan','yes','toilim',[-0.8 0],'avgovertime','yes');
            tmpdata_pre_front=ft_selectdata(data,'channel',chanLists.Front_All,'avgoverchan','yes','toilim',[-0.8 0],'avgovertime','yes');
            %===================================================================================================
            disp('WM : STAGE  6');
            totalSpecCases(countEntry,:)=[iCase 2 1 1 tmpdata_pre_occ.powspctrm];countEntry=countEntry+1;
            totalSpecCases(countEntry,:)=[iCase 2 2 1 tmpdata_pre_semoaud_L.powspctrm];countEntry=countEntry+1;
            totalSpecCases(countEntry,:)=[iCase 2 3 1 tmpdata_pre_semoaud_R.powspctrm];countEntry=countEntry+1;
            totalSpecCases(countEntry,:)=[iCase 2 4 1 tmpdata_pre_front.powspctrm];countEntry=countEntry+1;
            
            totalSpecCases(countEntry,:)=[iCase 2 1 2 tmpdata_post_occ.powspctrm];countEntry=countEntry+1;
            totalSpecCases(countEntry,:)=[iCase 2 2 2 tmpdata_post_semoaud_L.powspctrm];countEntry=countEntry+1;
            totalSpecCases(countEntry,:)=[iCase 2 3 2 tmpdata_post_semoaud_R.powspctrm];countEntry=countEntry+1;
            totalSpecCases(countEntry,:)=[iCase 2 4 2 tmpdata_post_front.powspctrm];countEntry=countEntry+1;
            
            totalSpecCases(countEntry,:)=[iCase 2 1 3 -1+tmpdata_post_occ.powspctrm./tmpdata_pre_occ.powspctrm];countEntry=countEntry+1;
            totalSpecCases(countEntry,:)=[iCase 2 2 3 -1+tmpdata_post_semoaud_L.powspctrm./tmpdata_pre_semoaud_L.powspctrm];countEntry=countEntry+1;
            totalSpecCases(countEntry,:)=[iCase 2 3 3 -1+tmpdata_post_semoaud_R.powspctrm./tmpdata_pre_semoaud_R.powspctrm];countEntry=countEntry+1;
            totalSpecCases(countEntry,:)=[iCase 2 4 3 -1+tmpdata_post_front.powspctrm./tmpdata_pre_front.powspctrm];countEntry=countEntry+1;
            
            allFreqs=tmpdata_post_occ.freq;
            disp('WM : STAGE  7');
            %########################################################################################################
        end
        
        disp('WM : STAGE  8');
        if isAnyFileFound
            casespectra=totalSpecCases;
            allfreqs=allFreqs;
            save(outputFile,'casespectra','columnDescr','allfreqs');
            disp('WM : STAGE  10');
        end
        
        
    end
end

fclose(logfid);

%%

%########################################################################################################
%########################################################################################################
%########################################################################################################
% STORYMATH

% A variable called casespectra is saved whis contains the speactra for different sensor subsets and different trial cases.
% COLUMN 1:   Member flag ( 1: Math Sentence onset    2: Story Sentence onset 3: all (this only for TRESP case))'
% COLUMN 2:   Trial Lock data group ( 1: TEV (For cases 1 and 2 of Columns1 only )    2:TRESP (for case 3 of Column 1 only) )'
% COLUMN 3:   Sensor group   (1: Occ/Posterior  2:SensoryMotor/Auditory Left 3:SensoryMotor/Auditory Right  4:Frontal)
% COLUMN 4:   Spectrum Type (1: Pre-stim(Baseline) 2:Post-stim  3: relative Change in Post-stim vs Pre-stim)
% REMANING COLUMNS :  Spectrum per frequency

columnDescr={
    'COLUMN 1:   Member flag ( 1: Math Sentence onset    2: Story Sentence onset 3: all (this only for TRESP case))'
    'COLUMN 2:   Trial Lock data group ( 1: TEV (For cases 1 and 2 of Columns1 only )    2:TRESP (for case 3 of Column 1 only) )'
    'COLUMN 3:   Sensor group   (1: Occ/Posterior  2:SensoryMotor/Auditory Left 3:SensoryMotor/Auditory Right  4:Frontal)'
    'COLUMN 4:   Spectrum Type (1: Pre-stim(Baseline) 2:Post-stim  3: relative Change in Post-stim vs Pre-stim)'
    'REMANING COLUMNS :  Spectrum per frequency'};

%----------------------------------------------------------------------------
caseMnems={'mathsentnon'
    'storsentnon'
    'all'};


subjLogFile=[avgspecdir,'SubjectsLog_StoryM_tfavg.txt'];
logfid=fopen(subjLogFile,'w+');

outDir=avgspecdir;

for iSubj=1:Nsubj,
    %--------------
    subjID=subjList{iSubj};
    isRESTIN=subjlistmat(iSubj,2);
    isMOTOR=subjlistmat(iSubj,3);
    isWORKMEM=subjlistmat(iSubj,4);
    isSTORYM=subjlistmat(iSubj,5);
    %--------------
    if isSTORYM,
        
        
        
        disp(['===========================================================================================']);
        disp(['STORY MATH - ',num2str(iSubj),' : ',subjID]);
        disp(['@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@']);
        %analDir=['/HCP/scratch/meg/intradb/archive1/HCP_Phase2/arc001/',subjID,'_MEG/RESOURCES/analysis/'];
        
        dataDir=[datamaindir,subjID,'_MEG/RESOURCES/analysis/'];
        outputFile=[outDir,'casespectra_',subjID,'_StoryM.mat'];
        
        countEntry=1;
        totalSpecCases=[];
        
        isAnyFileFound=0;
        for iCase=1:3,
            
            disp('SM : STAGE  1');
            curMnem=caseMnems{iCase};
            
            
            if iCase<3,
                datafile1=[dataDir,subjID,'_MEG_StoryM_tfavg_[LM-TEV-',curMnem,']_[MODE-mag].mat'];
                isFile1=exist(datafile1,'file');
                
                
                
                logChar=[subjID,'   ',num2str(iCase),'   ',num2str(isFile1)];
                fprintf(logfid,'%s\n',logChar);
                if (isFile1==0)
                    logChar=['OOOPS - tfavg files not found for StoryM TEV task of subject: ',subjID];
                    fprintf(logfid,'%s\n',logChar);
                    
                    warning(['OOOPS - tfavg files not found for StoryM TEV task of subject: ',subjID]);
                    disp(datafile1)
                    disp(datafile2);
                    continue;
                else
                    isAnyFileFound=1;
                end
                
                load(datafile1);
                disp('SM : STAGE  2');
                %===================================================================================================
                tmpdata_post_occ=ft_selectdata(data,'channel',chanLists.Occ_All,'avgoverchan','yes','toilim',[0 0.8],'avgovertime','yes');
                tmpdata_post_semoaud_L=ft_selectdata(data,'channel',chanLists.SeMoAud_L,'avgoverchan','yes','toilim',[0 0.8],'avgovertime','yes');
                tmpdata_post_semoaud_R=ft_selectdata(data,'channel',chanLists.SeMoAud_R,'avgoverchan','yes','toilim',[0 0.8],'avgovertime','yes');
                tmpdata_post_front=ft_selectdata(data,'channel',chanLists.Front_All,'avgoverchan','yes','toilim',[0 0.8],'avgovertime','yes');
                %-----------------------------------
                tmpdata_pre_occ=ft_selectdata(data,'channel',chanLists.Occ_All,'avgoverchan','yes','toilim',[-0.8 0],'avgovertime','yes');
                tmpdata_pre_semoaud_L=ft_selectdata(data,'channel',chanLists.SeMoAud_L,'avgoverchan','yes','toilim',[-0.8 0],'avgovertime','yes');
                tmpdata_pre_semoaud_R=ft_selectdata(data,'channel',chanLists.SeMoAud_R,'avgoverchan','yes','toilim',[-0.8 0],'avgovertime','yes');
                tmpdata_pre_front=ft_selectdata(data,'channel',chanLists.Front_All,'avgoverchan','yes','toilim',[-0.8 0],'avgovertime','yes');
                %===================================================================================================
                disp('SM : STAGE  3');
                totalSpecCases(countEntry,:)=[iCase 1 1 1 tmpdata_pre_occ.powspctrm];countEntry=countEntry+1;
                totalSpecCases(countEntry,:)=[iCase 1 2 1 tmpdata_pre_semoaud_L.powspctrm];countEntry=countEntry+1;
                totalSpecCases(countEntry,:)=[iCase 1 3 1 tmpdata_pre_semoaud_R.powspctrm];countEntry=countEntry+1;
                totalSpecCases(countEntry,:)=[iCase 1 4 1 tmpdata_pre_front.powspctrm];countEntry=countEntry+1;
                
                totalSpecCases(countEntry,:)=[iCase 1 1 2 tmpdata_post_occ.powspctrm];countEntry=countEntry+1;
                totalSpecCases(countEntry,:)=[iCase 1 2 2 tmpdata_post_semoaud_L.powspctrm];countEntry=countEntry+1;
                totalSpecCases(countEntry,:)=[iCase 1 3 2 tmpdata_post_semoaud_R.powspctrm];countEntry=countEntry+1;
                totalSpecCases(countEntry,:)=[iCase 1 4 2 tmpdata_post_front.powspctrm];countEntry=countEntry+1;
                
                totalSpecCases(countEntry,:)=[iCase 1 1 3 -1+tmpdata_post_occ.powspctrm./tmpdata_pre_occ.powspctrm];countEntry=countEntry+1;
                totalSpecCases(countEntry,:)=[iCase 1 2 3 -1+tmpdata_post_semoaud_L.powspctrm./tmpdata_pre_semoaud_L.powspctrm];countEntry=countEntry+1;
                totalSpecCases(countEntry,:)=[iCase 1 3 3 -1+tmpdata_post_semoaud_R.powspctrm./tmpdata_pre_semoaud_R.powspctrm];countEntry=countEntry+1;
                totalSpecCases(countEntry,:)=[iCase 1 4 3 -1+tmpdata_post_front.powspctrm./tmpdata_pre_front.powspctrm];countEntry=countEntry+1;
                
                allFreqs=tmpdata_post_occ.freq;
                
            elseif iCase==3
                
                datafile2=[dataDir,subjID,'_MEG_StoryM_tfavg_[LM-TRESP-',curMnem,']_[MODE-mag].mat'];
                
                isFile2=exist(datafile2,'file');
                
                
                logChar=[subjID,'   ',num2str(iCase),'   ',num2str(isFile2)];
                fprintf(logfid,'%s\n',logChar);
                if (isFile2==0),
                    logChar=['OOOPS - tfavg files not found for StoryM TRESP task of subject: ',subjID];
                    fprintf(logfid,'%s\n',logChar);
                    
                    warning(['OOOPS - tfavg files not found for StoryM TRESP task of subject: ',subjID]);
                    disp(datafile1)
                    disp(datafile2);
                    continue;
                else
                    isAnyFileFound=1;
                end
                
                
                
                disp('SM : STAGE  4');
                load(datafile2);
                disp('SM : STAGE  5');
                %===================================================================================================
                tmpdata_post_occ=ft_selectdata(data,'channel',chanLists.Occ_All,'avgoverchan','yes','toilim',[0 0.8],'avgovertime','yes');
                tmpdata_post_semoaud_L=ft_selectdata(data,'channel',chanLists.SeMoAud_L,'avgoverchan','yes','toilim',[0 0.8],'avgovertime','yes');
                tmpdata_post_semoaud_R=ft_selectdata(data,'channel',chanLists.SeMoAud_R,'avgoverchan','yes','toilim',[0 0.8],'avgovertime','yes');
                tmpdata_post_front=ft_selectdata(data,'channel',chanLists.Front_All,'avgoverchan','yes','toilim',[0 0.8],'avgovertime','yes');
                %-----------------------------------
                tmpdata_pre_occ=ft_selectdata(data,'channel',chanLists.Occ_All,'avgoverchan','yes','toilim',[-0.8 0],'avgovertime','yes');
                tmpdata_pre_semoaud_L=ft_selectdata(data,'channel',chanLists.SeMoAud_L,'avgoverchan','yes','toilim',[-0.8 0],'avgovertime','yes');
                tmpdata_pre_semoaud_R=ft_selectdata(data,'channel',chanLists.SeMoAud_R,'avgoverchan','yes','toilim',[-0.8 0],'avgovertime','yes');
                tmpdata_pre_front=ft_selectdata(data,'channel',chanLists.Front_All,'avgoverchan','yes','toilim',[-0.8 0],'avgovertime','yes');
                %===================================================================================================
                disp('SM : STAGE  6');
                totalSpecCases(countEntry,:)=[iCase 2 1 1 tmpdata_pre_occ.powspctrm];countEntry=countEntry+1;
                totalSpecCases(countEntry,:)=[iCase 2 2 1 tmpdata_pre_semoaud_L.powspctrm];countEntry=countEntry+1;
                totalSpecCases(countEntry,:)=[iCase 2 3 1 tmpdata_pre_semoaud_R.powspctrm];countEntry=countEntry+1;
                totalSpecCases(countEntry,:)=[iCase 2 4 1 tmpdata_pre_front.powspctrm];countEntry=countEntry+1;
                
                totalSpecCases(countEntry,:)=[iCase 2 1 2 tmpdata_post_occ.powspctrm];countEntry=countEntry+1;
                totalSpecCases(countEntry,:)=[iCase 2 2 2 tmpdata_post_semoaud_L.powspctrm];countEntry=countEntry+1;
                totalSpecCases(countEntry,:)=[iCase 2 3 2 tmpdata_post_semoaud_R.powspctrm];countEntry=countEntry+1;
                totalSpecCases(countEntry,:)=[iCase 2 4 2 tmpdata_post_front.powspctrm];countEntry=countEntry+1;
                
                totalSpecCases(countEntry,:)=[iCase 2 1 3 -1+tmpdata_post_occ.powspctrm./tmpdata_pre_occ.powspctrm];countEntry=countEntry+1;
                totalSpecCases(countEntry,:)=[iCase 2 2 3 -1+tmpdata_post_semoaud_L.powspctrm./tmpdata_pre_semoaud_L.powspctrm];countEntry=countEntry+1;
                totalSpecCases(countEntry,:)=[iCase 2 3 3 -1+tmpdata_post_semoaud_R.powspctrm./tmpdata_pre_semoaud_R.powspctrm];countEntry=countEntry+1;
                totalSpecCases(countEntry,:)=[iCase 2 4 3 -1+tmpdata_post_front.powspctrm./tmpdata_pre_front.powspctrm];countEntry=countEntry+1;
                allFreqs=tmpdata_post_occ.freq;
            end
            
            
            disp('SM : STAGE  7');
            
            
            %########################################################################################################
        end
        
        disp('SM : STAGE  8');
        if isAnyFileFound
            casespectra=totalSpecCases;
            allfreqs=allFreqs;
            save(outputFile,'casespectra','columnDescr','allfreqs');
            disp('SM : STAGE  10');
        end
        
        
    end
end

fclose(logfid);


%%
%########################################################################################################
%########################################################################################################
%########################################################################################################
% RESTING STATE

% A variable called casespectra is saved whis contains the speactra for different sensor subsets and different trial cases.
% COLUMN 1:   Scan Number ( 1    2  or  3)
% COLUMN 2:   Nan
% COLUMN 3:   Sensor group   (1: Occ/Posterior  2:SensoryMotor/Auditory Left 3:SensoryMotor/Auditory Right  4:Frontal)
% COLUMN 4:   NAN
% REMANING COLUMNS :  Spectrum per frequency

columnDescr={
    'COLUMN 1:   Scan Number ( 1    2  or  3)'
    'COLUMN 2: Nan'
    'COLUMN 3:   Sensor group   (1: Occ/Posterior  2:SensoryMotor/Auditory Left 3:SensoryMotor/Auditory Right  4:Frontal)'
    'COLUMN 4: Nan'
    'REMANING COLUMNS :  Spectrum per frequency'};

%----------------------------------------------------------------------------
caseMnems={'scan1'
    'scan2'
    'scan3'};

defaultFreqs=[1:31 33:2:99];

subjLogFile=[avgspecdir,'SubjectsLog_Restin_powavg.txt'];
logfid=fopen(subjLogFile,'w+');

outDir=avgspecdir;

for iSubj=1:Nsubj,
    %--------------
    subjID=subjList{iSubj};
    isRESTIN=subjlistmat(iSubj,2);
    isMOTOR=subjlistmat(iSubj,3);
    isWORKMEM=subjlistmat(iSubj,4);
    isSTORYM=subjlistmat(iSubj,5);
    %--------------
    if isRESTIN,
        
        
        
        disp(['===========================================================================================']);
        disp(['RESTING STATE - ',num2str(iSubj),' : ',subjID]);
        disp(['@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@']);
        %analDir=['/HCP/scratch/meg/intradb/archive1/HCP_Phase2/arc001/',subjID,'_MEG/RESOURCES/analysis/'];
        
        dataDir=[datamaindir,subjID,'_MEG/RESOURCES/analysis/'];
        outputFile=[outDir,'casespectra_',subjID,'_Restin.mat'];
        
        countEntry=1;
        totalSpecCases=[];
        
        isAnyFileFound=0;
        for iCase=1:3,
            
            disp('SM : STAGE  1');
            curMnem=caseMnems{iCase};
            
            
            
            datafile1=[dataDir,subjID,'_MEG_',num2str(iCase+2),'-Restin_powavg.mat'];
            isFile1=exist(datafile1,'file');
            
            
            
            logChar=[subjID,'   ',num2str(iCase),'   ',num2str(isFile1)];
            fprintf(logfid,'%s\n',logChar);
            if (isFile1==0)
                logChar=['OOOPS - tfavg files not found for Resting state of subject: ',subjID];
                fprintf(logfid,'%s\n',logChar);
                warning(['OOOPS - tfavg files not found for Resting state of subject: ',subjID]);
                disp(datafile1)
                continue;
            else
                isAnyFileFound=1;
            end
            
            load(datafile1,'freq');
            data=freq;
            disp('SM : STAGE  2');
            %===================================================================================================
            tmpdata_rest_occ=ft_selectdata(data,'channel',chanLists.Occ_All,'avgoverchan','yes');
            tmpdata_rest_semoaud_L=ft_selectdata(data,'channel',chanLists.SeMoAud_L,'avgoverchan','yes');
            tmpdata_rest_semoaud_R=ft_selectdata(data,'channel',chanLists.SeMoAud_R,'avgoverchan','yes');
            tmpdata_rest_front=ft_selectdata(data,'channel',chanLists.Front_All,'avgoverchan','yes');
            %===================================================================================================
            
            adjustFreqsIndx=[];
            for iDFs=1:length(defaultFreqs),
                adjustFreqsIndx(iDFs)=nearest(data.freq,defaultFreqs(iDFs));
            end
            
            tmpdata_rest_occ.powspctrm=tmpdata_rest_occ.powspctrm(adjustFreqsIndx);
            tmpdata_rest_occ.freq=defaultFreqs;
            tmpdata_rest_semoaud_L.powspctrm=tmpdata_rest_semoaud_L.powspctrm(adjustFreqsIndx);
            tmpdata_rest_semoaud_L.freq=defaultFreqs;
            tmpdata_rest_semoaud_R.powspctrm=tmpdata_rest_semoaud_R.powspctrm(adjustFreqsIndx);
            tmpdata_rest_semoaud_R.freq=defaultFreqs;
            tmpdata_rest_front.powspctrm=tmpdata_rest_front.powspctrm(adjustFreqsIndx);
            tmpdata_rest_front.freq=defaultFreqs;
            
            
            %===================================================================================================
            disp('SM : STAGE  3');
            totalSpecCases(countEntry,:)=[iCase nan 1 nan tmpdata_rest_occ.powspctrm];countEntry=countEntry+1;
            totalSpecCases(countEntry,:)=[iCase nan 2 nan tmpdata_rest_semoaud_L.powspctrm];countEntry=countEntry+1;
            totalSpecCases(countEntry,:)=[iCase nan 3 nan tmpdata_rest_semoaud_R.powspctrm];countEntry=countEntry+1;
            totalSpecCases(countEntry,:)=[iCase nan 4 nan tmpdata_rest_front.powspctrm];countEntry=countEntry+1;
            
            
            allFreqs=tmpdata_rest_occ.freq;
            %%
            
            
            disp('SM : STAGE  7');
            
            
            %########################################################################################################
        end
        
        disp('SM : STAGE  8');
        if isAnyFileFound
            casespectra=totalSpecCases;
            allfreqs=allFreqs;
            save(outputFile,'casespectra','columnDescr','allfreqs');
            disp('SM : STAGE  10');
        end
        
        
    end
end

fclose(logfid);

%%
%((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((
%((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((
%((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((
%((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((

%             END OF EXTRACTING SPECTRA  IN SENSOR GROUPS

%))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
%))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
%))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
%))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))