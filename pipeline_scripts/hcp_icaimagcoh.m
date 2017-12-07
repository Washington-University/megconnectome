%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2011-2014 by the Human Connectome Project, WU-Minn Consortium (1U54MH091657)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setup the execution environment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opengl software;

% ensure that the time and date of execution are not stored in the provenance information
global ft_default
ft_default.trackcallinfo = 'no';

% allow the user to specify the path where additional data is present, e.g. the channel layout or anatomy files
if exist('path', 'var')
    addpath(path)
end

if ~exist('filename', 'var')
    error('filename should be specified')
end

if ~exist('f', 'var')
    error('frequency should be specified as "f"')
end
freque = f;

% the filename is assumed to be something like
% 'rawdatadir/Phase1MEG/Subjects/CP10018/Experiments/CP10018_MEG/Scans/1-Rnoise_MNN_V1/Resources/4D/c,rfDC'
tok = tokenize(filename, '/');

if ~exist('subjectid', 'var')
    subjectid = tok{end-7};
end

if ~exist('experimentid', 'var')
    experimentid = tok{end-5};
end

if ~exist('scanid', 'var')
    scanid = tok{end-3};
end

if ~exist('pipelinedatadir', 'var')
    pipelinedatadir = hcp_pathdef;
end

% print the matlab and megconnectome version to screen for provenance
ver('megconnectome')

% print the value of all local variables to screen for provenance
w = whos;
w = {w.name};
w = setdiff(w, {'w', 'ans'});
for i=1:length(w)
    fprintf(hcp_printstruct(w{i}, eval(w{i})));
end

% change to the location of the processed data (input and output)
cd(pipelinedatadir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% execute the pipeline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

smodel_type = {'2D';'3D'}; dimindx=1; % 1 for 2D cortical sheet and 2 for 3D gird
griddim = {'4mm';'6mm';'8mm'}; gridindx=1; % $ 1,2,3 for 3D 4mm,6mm and 8mm grid

if(strcmp(smodel_type{dimindx},'2D'))
    sourcemodel_type=smodel_type{dimindx};
elseif(strcmp(smodel_type{dimindx},'3D'))
    sourcemodel_type=[smodel_type{dimindx} griddim{gridindx}];
end

%------------------------
% declare the output file
resultprefix = sprintf('%s_%s', experimentid, scanid);
outputfile = [resultprefix,'_icaimagcoh'];

% ------------------------------------------------------------
% ensure that the output from the previous pipelines is present
% hcp_check_pipelineoutput('baddata', 'subject', subjectid, 'experiment', experimentid, 'scan', scanid);
% hcp_check_pipelineoutput('icaclass', 'subject', subjectid, 'experiment', experimentid, 'scan', scanid);
% hcp_check_pipelineoutput('icamne', 'subject', subjectid, 'experiment', experimentid, 'scan', scanid, 'sourcemodel', sourcemodel_type);

inputfile4 = fullfile([resultprefix '_icaclass_vs.mat']);
hcp_read_matlab(inputfile4);

inputfile5 = fullfile([resultprefix '_icamne.mat']);
hcp_read_matlab(inputfile5,'source');

if(~isfield(comp_class,'trial'))
    
    cfg = [];
    cfg.dataset = filename;
    dataraw=ft_preprocessing(cfg);
    
    badsegments = hcp_read_ascii([resultprefix '_baddata_badsegments.txt']);
    
    disp('bad segments concatenation')
    % collect the results in a structure
    maxsmp2 = max(badsegments.badsegment.ica(:));
    maxsmp3 = max(badsegments.badsegment.manual(:));
    maxsmp  = max([maxsmp2, maxsmp3]);
    badsample = zeros(1,maxsmp);
    if ~isempty(badsegments.badsegment.ica)
        for j=1:size(badsegments.badsegment.ica,1)
            badsample(badsegments.badsegment.ica(j,1):badsegments.badsegment.ica(j,2)) = 1;
        end
    end
    if ~isempty(badsegments.badsegment.manual)
        for j=1:size(badsegments.badsegment.manual,1)
            badsample(badsegments.badsegment.manual(j,1):badsegments.badsegment.manual(j,2)) = 1;
        end
    end
    
    if ~isempty(badsample)
        flank  = diff([0 badsample 0]);
        badsegment_ica=[];
        badsegment_ica(:,1) = find(flank== 1);
        badsegment_ica(:,2) = find(flank==-1) - 1;
        
        badsegments.badsegment.ica       = badsegment_ica;
    end
    
    badsegments = badsegments.badsegment.ica;
    
    badchannels = hcp_read_ascii([resultprefix '_baddata_badchannels.txt']);
    badchannels = badchannels.badchannel.all;
    
    sel_channels={'MEG'};
    grad=dataraw.grad;
    
    
    %------------------------------------------
    
    if ~(isempty([badchannels{:}])),
        for ich=1:size(badchannels,2)
            sel_channels(1,ich+1)={['-' badchannels{1,ich}]};
        end
    end
    
    bandpass = [1.3 150]; % band pass frequency
    if strcmp(subjectid,'CP10128') || strcmp(subjectid,'CP10129')
        bandstop = [49 51 ; 99 101]; % band stop frequency, these were recorded in Glasgow
    else
        bandstop = [59 61 ; 119 121]; % band stop frequency, all others were recorded at SLU
    end
    
    options  = {'dataprefix', resultprefix, 'channels', sel_channels, 'skipped_intervals', badsegments, 'bandpass', bandpass, 'bandstop', bandstop};
    
    data     = hcp_ICA_preprocessing(dataraw, options);
    
    cfg           = [];
    cfg.unmixing  = comp_class.unmixing;
    cfg.topolabel = comp_class.topolabel;
    comp          = ft_componentanalysis(cfg, data);
    comp.class    = comp_class.class;
    comp.topo     = comp_class.topo;
    clear data
else
    comp = comp_class;
end

mixing = comp_class.topo;
if(max(size(source.val))>2)
    for i = 1:size(mixing, 2)
        mixing(:, i) = mixing(:, i)/source.val(i);
        for jic=1:size(comp.trial,2)
            comp.trial{jic}(i,:)=comp.trial{jic}(i,:)*source.val(i);
        end
    end
end
comp.topo = mixing;

% get only the brain ICs
comp_bic=comp;

comp_bic=rmfield(comp_bic,'trial');
for j=1:numel(comp.trial)
    comp_bic.trial{j}(:,:)=comp.trial{j}(comp.class.brain_ic_vs,:);
end
comp_bic.topo=comp.topo(:,comp.class.brain_ic_vs);
comp_bic.unmixing=comp.unmixing(comp.class.brain_ic_vs,:);
comp_bic.label=comp.label(comp.class.brain_ic_vs,:);


% adjust the time axis to avoid memory problems during resampling:
% exact time information is discarded anyway
old_time=comp_bic.time;
for k = 1:numel(comp_bic.time)
    comp_bic.time{k} = comp_bic.time{k} - comp_bic.time{k}(1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the topographies of the components at the level of the
% sources from the mne pipeline

source_bic=source;
source_bic.time=1:comp_bic.class.brain_ic_vs_number;
source_bic.avg.pow=zeros(size(source.pos,1),comp_bic.class.brain_ic_vs_number);
source_bic.avg=rmfield(source_bic.avg,'mom');

for ii=1:numel(source.inside)
    source_bic.avg.mom{source.inside(ii)}=source.avg.mom{source.inside(ii)}(:,comp_bic.class.brain_ic_vs);
    source_bic.avg.pow(source.inside(ii),:)=source.avg.pow(source.inside(ii),comp_bic.class.brain_ic_vs);
end



%%hack to remove balancing
% gradBalanced = comp_bic.grad;
% gradBalanced = ft_apply_montage(gradBalanced, gradBalanced.balance.comp, 'keepunused', 'yes', 'inverse', 'yes');
% gradBalanced = ft_apply_montage(gradBalanced, gradBalanced.balance.Supine, 'keepunused', 'yes', 'inverse', 'yes');
%
% comp_bic.grad=gradBalanced;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   working piece of code to test ft function behaviour on components
cfg = [];
cfg.length     = 1;
cfg.overlap    = 0.5;
cfg.trials    = 'all';
[comp1] = ft_redefinetrial(cfg, comp_bic);

cfg = [];
cfg.output = 'fourier';
cfg.method = 'mtmfft';
cfg.taper  = 'hanning';
[comp2] = ft_freqanalysis(cfg, comp1);
comp2=ft_selectdata(comp2,'foilim',[1:80]);

[ntrials nic nfreq]=size(comp2.fourierspctrm);

cfg = [];
cfg.method = 'csd';
cfg.complex  = 'complex';
[comp3] = ft_connectivityanalysis(cfg, comp2);
%comp_csd dimord=[chan x chan x freq]

nsource = numel(source_bic.inside);
nfreq = length(freque);


nsource = numel(source_bic.inside);
[ndim,nbic] = size(source_bic.avg.mom{source_bic.inside(1)});


rpi=20;

allit = ceil(nsource/rpi);
count = 0;
% tStart = tic;
% dispstat('','init'); % One time only initialization
% dispstat(sprintf('MIM - Beginning ...'),'keepthis','timestamp');

for f=freque
    % allocate memory for single frequency connectomes
    mimf  = zeros(nsource, nsource,'single');
   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    WlfAic = cat(1,source_bic.avg.mom{source_bic.inside(:)'});
    WlfAic = reshape( permute( reshape(WlfAic,ndim,nsource,[]),[2 1 3]),nsource*ndim,[]) ;
    ndum = size(WlfAic,1);
    
    CS = comp3.crsspctrm(:,:,f);
    
    cs_bb = hcp_cs2csvv(CS,WlfAic,ndim);
    pro_prcsbb = zeros(ndim,ndum);
    for iv=1:nsource
        pro_prcsbb(:,(iv-1)*ndim + [1:ndim]) =  pinv(real( full(cs_bb((iv-1)*ndim + [1:ndim],(iv-1)*ndim + [1:ndim]))));
    end
    pro_prcsbb = repmat(pro_prcsbb,1,rpi);
    
    iBD = reshape(padarray(reshape(1:ndum*rpi,ndim,[]),ndim,'circular'),[],1);
    jBD = reshape(padarray(1:ndum*rpi,1,'replicate'),[],1);
    prcsbb = sparse(iBD,jBD,reshape(pro_prcsbb,[],1),ndum*rpi,ndum*rpi);
    
    iv_set = reshape( padarray([1:nsource]',(rpi*ceil(nsource/rpi))-nsource,'post'),rpi,[])' ;
    for iset = 1:size(iv_set,1)
        
        iv = nonzeros(iv_set(iset,:))';
        jv  = iv(1):nsource;
        
        niv = length(iv);
        njv = length(jv);
        
        pro_prcsaa = zeros(ndim,ndim,niv);
        for i = 1:niv
            pro_prcsaa(:,:,i) = prcsbb((iv(i)-1)*ndim + [1:ndim],(iv(i)-1)*ndim + [1:ndim]);
        end
        pro_prcsaa = repmat(pro_prcsaa,[1 njv 1]);
        prcsaa = sparse(iBD(1:ndim*ndim*njv*niv),jBD(1:ndim*ndim*njv*niv),reshape(pro_prcsaa,[],1),ndim*njv*niv,ndim*njv*niv);
        
        indA = reshape(repmat(iv',1,ndim) + repmat([0 nsource 2*nsource],niv,1),1,[]);
        indB = reshape(repmat(jv',1,ndim) + repmat([0 nsource 2*nsource],njv,1),1,[]);
        
        CSso_iv = WlfAic(indA,:)* CS * (WlfAic(indB,:)');
        
        nvox_iv = size(CSso_iv,1)/ndim;
        nvox_jv = size(CSso_iv,2)/ndim;
        
        CSso_iv = reshape(CSso_iv,niv,ndim,ndim*njv);
        CSso_iv = permute(CSso_iv,[1 3 2]);
        CSso_iv = reshape(CSso_iv,niv,njv,ndim,ndim);
        CSso_iv = permute(CSso_iv,[1 2 4 3]);
        
        cs_ab = sparse(iBD(1:ndim*ndim*njv*niv),jBD(1:ndim*ndim*njv*niv),reshape(permute(CSso_iv,[3 4 2 1]),[],1),ndim*njv*niv,ndim*njv*niv);
        
        sub_ind = zeros(niv,njv*3);
        for i=1:niv
            sub_ind(i,:) = (i-1)*nsource*ndim + reshape( ( repmat((jv-1)'*ndim,1,ndim) + repmat([1:ndim],njv,1) )',1,[]);
        end
        sub_ind = reshape(sub_ind',1,[]);
        
        pro_MIM = full( diag(prcsaa * imag(cs_ab) * prcsbb(sub_ind,sub_ind) * imag(cs_ab)')) ;
        
        mimf(iv,jv) = cast( squeeze(permute(sum( reshape(pro_MIM,ndim,[],niv),1),[3 2 1])),'single');
        
        count = count+1;
        %         dispstat(sprintf('Progress
        %         %d%%',round(100*count/allit)),'timestamp');
        
    end
    mimf = ( triu(mimf) + triu(mimf,1)' );
    mimf = mimf - diag(0.5*diag(mimf));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    imagcoh=[];
    imagcoh.freq   = comp2.freq(f);
    imagcoh.dimord = 'pos_pos';
    imagcoh.mimspctrm = mimf;
    imagcoh.pos = source.pos(source.inside,:);
    
    hcp_write_matlab([outputfile,'_',num2str(round(comp3.freq(f))),'Hz'], 'imagcoh');
    hcp_check_pipelineoutput('icaimagcoh', 'subject', subjectid, 'experiment', experimentid, 'scan', scanid, 'freq', num2str(round(comp3.freq(f))));
    
end % nfreq

% dispstat('Finished.','keepprev');
% tElapsed = toc(tStart);
% disp(['Elapsed time = ',num2str(floor(tElapsed/60)),' min and ',num2str(round(tElapsed - floor(tElapsed/60)*60)),' sec.']);
