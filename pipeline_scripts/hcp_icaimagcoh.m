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

% This is not default anymore, as of October 28
%if ~exist('freque', 'var')
%  freque=[1:80];
%end

if exist('freque', 'var')
    blp_bands = freque(:);
    aband=1:size(blp_bands);
    for i = 1:numel(freque)
        band_prefix{i,1} = sprintf('%02dHz',freque(i));
    end
else
    blp_bands   = [1 4;4 8;8 15;15 26;26 35;35 50;50 76;76 120];
    band_prefix={ 'delta'; 'theta' ; 'alpha' ; 'betalow' ; 'betahigh' ; 'gammalow' ; 'gammamid' ;  'gammahigh'};
    
    if ~exist('aband', 'var')
        aband=[1 2 3 4 5 6 7 8];
    end
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

%------------------------
% declare the output file
resultprefix = sprintf('%s_%s', experimentid, scanid);
outputfile = [resultprefix,'_icaimagcoh'];

% ------------------------------------------------------------
% ensure that the output from the previous pipelines is present
hcp_check_pipelineoutput('baddata', 'subject', subjectid, 'experiment', experimentid, 'scan', scanid);
hcp_check_pipelineoutput('icaclass', 'subject', subjectid, 'experiment', experimentid, 'scan', scanid);
hcp_check_pipelineoutput('icamne', 'subject', subjectid, 'experiment', experimentid, 'scan', scanid);

inputfile1 = fullfile([resultprefix '_icaclass_vs.mat']);
hcp_read_matlab(inputfile1);

inputfile2 = fullfile([resultprefix '_icamne.mat']);
hcp_read_matlab(inputfile2,'source');

% rename for convenience
comp = comp_class;
clear comp_class;

% weight the topographies and the component time series (in opposite
% directions, so the mixing stays correct)
mixing   = comp.topo;
unmixing = comp.unmixing;
if(max(size(source.val))>2)
    for i = 1:size(mixing, 2)
        unmixing(i, :) = unmixing(i, :)*source.val(i);
        mixing(:, i)   = mixing(:, i)/source.val(i);
        for jic = 1:size(comp.trial,2)
            comp.trial{jic}(i,:) = comp.trial{jic}(i,:)*source.val(i);
        end
    end
end
comp.topo     = mixing;
comp.unmixing = unmixing;
clear mixing unmixing

% get only the brain ICs
comp_bic = rmfield(comp, 'trial');

for k = 1:numel(comp.trial)
    comp_bic.trial{k}(:,:) = comp.trial{k}(comp.class.brain_ic_vs,:);
end
comp_bic.topo     = comp.topo(:,comp.class.brain_ic_vs);
comp_bic.unmixing = comp.unmixing(comp.class.brain_ic_vs,:);
comp_bic.label    = comp.label(comp.class.brain_ic_vs,:);
clear comp

% adjust the time axis to avoid memory problems during resampling:
% exact time information is discarded anyway
old_time = comp_bic.time;
for k = 1:numel(comp_bic.time)
    comp_bic.time{k} = comp_bic.time{k} - comp_bic.time{k}(1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the topographies of the components at the level of the
% sources from the mne pipeline

source_bic         = source;
source_bic.time    = 1:comp_bic.class.brain_ic_vs_number;
source_bic.avg.pow = zeros(size(source.pos,1),comp_bic.class.brain_ic_vs_number);
source_bic.avg     = rmfield(source_bic.avg,'mom');

inside_indices = find(source.inside(:))';

for k = inside_indices(:)'
    source_bic.avg.mom{k}   = source.avg.mom{k}(:,comp_bic.class.brain_ic_vs);
    source_bic.avg.pow(k,:) = source.avg.pow(k,comp_bic.class.brain_ic_vs);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create component-level cross-spectral density matrices

% create 1-second snippets of data with 50% overlap
cfg         = [];
cfg.length  = 1;
cfg.overlap = 0.5;
cfg.trials  = 'all';
comp1       = ft_redefinetrial(cfg, comp_bic);
clear comp_bic;

% spectral analysis
cfg        = [];
cfg.output = 'fourier';
cfg.method = 'mtmfft';
cfg.taper  = 'hanning';
cfg.foilim = [1 max(blp_bands(:))];
comp2      = ft_freqanalysis(cfg, comp1);
clear comp1;


% get some sizes
[ntrials, nic, dummy] = size(comp2.fourierspctrm);
nfreq                 = numel(comp2.freq);
nsource               = numel(find(source_bic.inside));
[ndim, nbic]          = size(source_bic.avg.mom{source_bic.inside(1)});


% compute csd-matrices
cfg          = [];
cfg.method   = 'csd';
cfg.complex  = 'complex';
comp3        = ft_connectivityanalysis(cfg, comp2);
%comp_csd dimord=[chan x chan x freq]
clear comp2;


% allocate some other stuff
rpi   = 20;
allit = ceil(nsource/rpi);
count = 0;

for fband=aband
    
    % allocate memory for single frequency-band connectome
    mimf_band  = zeros(nsource, nsource,'single');
    
    fbeg = blp_bands(fband,1)
    fend = blp_bands(fband,end)
    
    % decision: if the lower edge is in between two frequency bins, round to the ceiling
    % decision: if the upper edge is in between two frequency bins, round to the floor
    % rationale: a frequency bin captures energy at frequencies +/- 0.5*Rayleigh frequency
    
    % decision: if a single frequency bin is requested which is in between the sampled frequencies,
    % take the nearest
    if fbeg==fend
        fbegindx = nearest(comp3.freq, fbeg);
        findx    = fbegindx;
    else
        fbegindx = find(comp3.freq-fbeg>=0, 1, 'first');
        fendindx = find(comp3.freq-fend<=0, 1, 'last');
        findx    = fbegindx:fendindx;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for f=findx
        disp(['f ' num2str(f)])
        mimf = zeros(nsource, nsource, 'single');
        
        WlfAic = cat(1,source_bic.avg.mom{inside_indices(:)});
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
            
        end
        mimf = ( triu(mimf) + triu(mimf,1)' );
        mimf = mimf - diag(0.5*diag(mimf));
        
        mimf_band = mimf_band + mimf;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % normalize with the number of frequency bins
    mimf_band = mimf_band./numel(findx);
    
    % get the data in 'full' representation, i.e. also allocate nans to the
    % outside vertices
    tmp = nan+zeros(size(source.pos,1));
    tmp(inside_indices, inside_indices) = mimf_band;
    
    % create the output structure
    imagcoh           = [];
    imagcoh.freq      = mean(comp3.freq(findx));
    imagcoh.dimord    = 'pos_pos';
    imagcoh.mimspctrm = tmp;
    imagcoh.pos       = source.pos;
    if isfield(source, 'tri'),     imagcoh.tri     = source.tri;     end
    if isfield(source, 'inside'),  imagcoh.inside  = source.inside;  end
    if isfield(source, 'outside'), imagcoh.outside = source.outside; end
    if isfield(source, 'brainstructure'), imagcoh.brainstructure = source.brainstructure; end
    if isfield(source, 'brainstructure'), imagcoh.brainstructurelabel = source.brainstructurelabel; end
    
    % write as a cifti
    hcp_write_cifti([outputfile,'_',band_prefix{fband}], imagcoh, 'parameter', 'mimspctrm', 'type', 'dconn');
    
    % also write as a mat file
    hcp_write_matlab([outputfile,'_',band_prefix{fband}], 'imagcoh');
    
    % check the output
    hcp_check_pipelineoutput('icaimagcoh', 'subject', subjectid, 'experiment', experimentid, 'scan', scanid, 'band', band_prefix{fband});
    
end % fband
