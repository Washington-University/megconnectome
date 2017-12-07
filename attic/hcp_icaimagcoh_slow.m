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

if ~exist('freque', 'var')
    error('frequency of interest should be specified')
end

resultprefix = sprintf('%s_%s', experimentid, scanid);

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
outputfile = [resultprefix,'_icaimagcoh_' sourcemodel_type];


%------------------------------------------------------------
% ensure that the output from the previous pipelines is present
% hcp_check_pipelineoutput('baddata', 'subject', subjectid, 'experiment', experimentid, 'scan', scanid);
% hcp_check_pipelineoutput('icaclass', 'subject', subjectid, 'experiment', experimentid, 'scan', scanid);
%     hcp_check_pipelineoutput('icamne', 'subject', subjectid, 'experiment', experimentid, 'scan', scanid, 'sourcemodel', sourcemodel_type);


inputfile4 = fullfile([resultprefix,'_icaclass_vs.mat']);
hcp_read_matlab(inputfile4);

inputfile5 = fullfile([resultprefix,'_icamne_' sourcemodel_type]);
hcp_read_matlab(inputfile5,'source');

if(~isfield(comp_class,'trial'))
    
    cfg = [];
    cfg.dataset = filename;
    dataraw=ft_preprocessing(cfg);
    
    badsegments = hcp_read_ascii([resultprefix '_baddata_badsegments.txt']);
    %     badsegments = badsegments.badsegment.all;
    
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
        bandstop = [49 51 ; 99 101]; % band stop frequency
    else
        bandstop = [59 61 ; 119 121]; % band stop frequency
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
    comp=comp_class;
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
comp.topo=mixing;



% get only the brain ICs
comp_bic=comp;
if(~isfield(comp_bic.class,'brain_ic_vs')) comp_bic.class.brain_ic_vs=comp_bic.class.brain_ic; end

comp_bic=rmfield(comp_bic,'trial');
for j=1:numel(comp.trial)
    comp_bic.trial{j}(:,:)=comp.trial{j}(comp.class.brain_ic_vs,:).*1e15;
end
comp_bic.topo=comp.topo(:,comp.class.brain_ic_vs);
comp_bic.unmixing=comp.unmixing(comp.class.brain_ic_vs,:);
comp_bic.label=comp.label(comp.class.brain_ic_vs,:);

% channel-data is not needed any more
clear data

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

comp2=ft_selectdata(comp2,'foilim',freque);
%   comp2=ft_selectdata(comp2,'rpt',[1:10]);

[ntrials nic nfreq]=size(comp2.fourierspctrm);

cfg = [];
cfg.method = 'csd';
cfg.complex  = 'complex';
[comp3] = ft_connectivityanalysis(cfg, comp2);
%comp_csd dimord=[chan x chan x freq]

nsource = numel(source_bic.inside);
%   nfreq = length(comp3.freq);
%nfreq = 20; % FIXME why is this hardcoded?
estimate = nsource * nsource * nfreq * 8; % 2 because complex, 8 because double precision
fprintf('estimated memory requirement = %f GB\n', estimate/(1024^3));

nsource = numel(source_bic.inside);
[ndim,nbic] = size(source_bic.avg.mom{source_bic.inside(1)});

% freque=10;
% disp(freque)



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

comp2=ft_selectdata(comp2,'foilim',[1:30]);
%   comp2=ft_selectdata(comp2,'rpt',[1:10]);

[ntrials nic nfreq]=size(comp2.fourierspctrm);

cfg = [];
cfg.method = 'csd';
cfg.complex  = 'complex';
[comp3] = ft_connectivityanalysis(cfg, comp2);
%comp_csd dimord=[chan x chan x freq]

nsource = numel(source_bic.inside);
%   nfreq = length(comp3.freq);
%nfreq = 20; % FIXME why is this hardcoded?
estimate = nsource * nsource * nfreq * 8; % 2 because complex, 8 because double precision
fprintf('estimated memory requirement = %f GB\n', estimate/(1024^3));

nsource = numel(source_bic.inside);


for f=freque
    % allocate memory for single frequency connectomes
    mimf        = zeros(nsource, nsource,'single');
    %     powenvcorrf = zeros(nsource, nsource);
    
    % show the progress
    % concatenate all fourier coefficients into a single matrix
    allfourier = cat(1,source_bic.avg.mom{source_bic.inside(:)'})*comp2.fourierspctrm(:,:,f).';
    
    % concatenate all powerenvelope values (demeaned) into a single matrix
    %     allpowenv  = abs(allfourier).^2;
    %     allpowenv  = allpowenv - mean(allpowenv,2)*ones(1,size(allpowenv,2));
    
    % create a cell-array containing the within source csd and one
    % cell-array with the within source power envelope covariance
    % this needs to be computed only once for all sources
    autocsd    = cell(1,nsource);
    
    %     autopowcov = cell(1,nsource);
    gim        = zeros(1,nsource,'single');
    n          = zeros(1,nsource);
    pow        = zeros(1,nsource);
    %     powpow     = zeros(1,nsource);
    indx = 0;
    for i=1:nsource
        n(i)       = size(source_bic.avg.mom{source_bic.inside(i)},1);
        indx       = indx(end)+(1:n(i));
        autocsd{i} = allfourier(indx,:)*allfourier(indx,:)';
    end
    
    %     for j=1:nsource
    %         [eigvec{j} eigval{j}]=eig(autocsd{j});
    %         [junk indxeig{j}]=sort(real(diag(eigval{j})),'descend');
    %         autocsd_new{j}=eigvec{j}(indxeig{j}(1:2),indxeig{j}(1:2))*eigval{j}(indxeig{j}(1:2),indxeig{j}(1:2))*inv(eigvec{j}(indxeig{j}(1:2),indxeig{j}(1:2)));
    %     n(j)=2;
    %     n2(j)=3;
    %     end
%     clear autocsd;
%     autocsd=autocsd_new;
%     clear eigvec eigval indxeig
    
    for j=1:numel(indx)
        autocsd{i}(j,j)=real(autocsd{i}(j,j));
    end
    for i=1:nsource
        pow(i)     = trace((autocsd{i}));
        
        %         autopowcov{i} = allpowenv(indx,:)*allpowenv(indx,:)';
        %         powpow(i)  = trace(autopowcov{i});
        gim(i)     = cast(hcp_connectivity_mim([autocsd{i} autocsd{i};autocsd{i}' autocsd{i}], 'indices', [ones(1,n(i)) ones(1,n(i))*2]),'single');
    end
    
    indx = 0;
    for i=1:nsource
        %         disp(i);
        % treat source i as the reference source and compute the
        % metrics of interest against all other sources. compute the
        % cross-terms only once
        indx = indx(end)+(1:n(i));
        crosscsd    = allfourier*allfourier(indx,:)';
        
        %         crosspowcov = allpowenv*allpowenv(indx,:)';
        
        indx2 = cumsum(n(1:i));
        for j=(i+1):nsource
            % now compute the metric per voxel pair
            indx2 = indx2(end)+(1:n(j));
            
            %         [eigvec eigval]=eig(crosscsd(indx2,:));
            %         [junk indxeig]=sort(real(diag(eigval)),'descend');
            %         crosscsd_new=eigvec(indxeig(1:2),indxeig(1:2))*eigval(indxeig(1:2),indxeig(1:2))*inv(eigvec(indxeig(1:2),indxeig(1:2)));
            %
            %
            % combine the auto and cross-voxel csd matrices into a
            % single matrix, for the mim computation
            %               avecsd = [autocsd{i} crosscsd';crosscsd autocsd{j}];
            avecsd = [autocsd{i} crosscsd(indx2,:)';crosscsd(indx2,:) autocsd{j}];
            mim    = hcp_connectivity_mim(avecsd, 'indices', [ones(1,n(i)) ones(1,n(j))*2]);
            
            % power envelope correlation
            %             powenvcorr = trace(crosspowcov(indx2,:))./sqrt(powpow(i)*powpow(j)); %FIXME the trace-operator does not work with non-square matrices
            
            % fill the connectivity matrix
            mimf(i,j)=cast(mim,'single');
            mimf(j,i)=cast(mim,'single');
            
            %             powenvcorrf(i,j) = powenvcorr;
            %             powenvcorrf(j,i) = powenvcorr;
            
        end % nsource
    end % nsource
    mimf        = mimf + cast(diag(gim),'single');
    
    source2=source_bic;
    source2=rmfield(source2,'avg');
    source2=rmfield(source2,'time');
    %         source2        = [];
    %         source2.pos    = source.pos(source.inside,:);
    %         source2.inside = 1:nsource;
    source2.freq   = comp2.freq(f);
    source2.dimord = 'pos_pos_freq';
    source2.mimspctrm = mimf;
    %     source2.powenvcorrspctrm = powenvcorrf;
    
    imagcoh=[];
    imagcoh.freq   = comp2.freq(f);
    imagcoh.dimord = 'pos_pos_freq';
    imagcoh.mimspctrm = mimf;
%     
%     save([outputfile,'_freq',num2str(f)], 'imagcoh','-v7.3');
%     
%     imagcoh.autocsd = autocsd;
%     save([outputfile,'_full_freq',num2str(f)], 'imagcoh','-v7.3');
    
    %     imagcoh.powenvcorrspctrm = powenvcorrf;

    
        hcp_write_matlab([outputfile,'_freq',num2str(f)], 'imagcoh');
    
    
%     clear source2 imagcoh mimf allfourier allpowenv autocsd gim n pow
%     clear crosscsd crosspowcov avecsd mim
     hcp_check_pipelineoutput('icaimagcoh', 'subject', subjectid, 'experiment', experimentid, 'scan', scanid, 'sourcemodel', sourcemodel_type, 'freq', num2str(f));
    
    
    clear source2 imagcoh mimf allfourier allpowenv autocsd gim n pow
    clear crosscsd crosspowcov avecsd mim
    % hcp_check_pipelineoutput('icaimagcoh', 'subject', subjectid, 'experiment', experimentid, 'scan', scanid, 'sourcemodel', sourcemodel_type, 'freq', num2str(f));
    
end % nfreq

