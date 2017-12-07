function [ source ] = hcp_beamforming_blp(gridAllLF, datain, vol, grad, options_blp)

% HCP_ICA_BLP allows the calculation of Band Limited Power.

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



subject   = ft_getopt(options_blp, 'dataprefix');
band_prefix = ft_getopt(options_blp, 'band_prefix');
dofiltcheck = ft_getopt(options_blp, 'dofiltcheck', 'no');
blp_band     = ft_getopt(options_blp, 'blp_band', [6.3 16.5]);
bf_lambda = ft_getopt(options_blp, 'bf_lambda', '75%');

if strcmp(band_prefix,'betalow')
    N_hp=9; N_lp=10 ;
elseif strcmp(band_prefix,'betahigh')
    N_hp=9; N_lp=12 ;
elseif strcmp(band_prefix,'alpha')
    N_hp=7; N_lp=8 ;
elseif strcmp(band_prefix,'theta')
    N_hp=6; N_lp=7 ;
elseif strcmp(band_prefix,'delta')
    N_hp=4; N_lp=7 ;
elseif strcmp(band_prefix,'gammalow')
    N_hp=9; N_lp=13 ;
elseif strcmp(band_prefix,'gammamid')
    N_hp=11; N_lp=14 ;
elseif strcmp(band_prefix,'gammahigh')
    N_hp=13; N_lp=14 ;
elseif strcmp(band_prefix,'whole')
    N_hp=4; N_lp=14 ;
end

step   = ft_getopt(options_blp, 'blp_step',    20);
window = ft_getopt(options_blp, 'blp_window', 400);
band_str=[num2str(blp_band(1,1)) '-' num2str(blp_band(1,2))];

disp(['blp band -> ' band_str ' Hz'])
disp(['power calculation step -> ' num2str(step) ' ms'])
disp(['power calculation window -> ' num2str(window) ' ms'])

% filtering
cfgin          = [];
cfgin.hpfilter = 'yes';
cfgin.hpfiltord=N_hp;
cfgin.hpfreq   = blp_band(1,1);
junk     = ft_preprocessing(cfgin,datain);

cfgin          = [];
cfgin.lpfilter = 'yes';
cfgin.lpfiltord=N_lp;
cfgin.lpfreq   = blp_band(1,2);
databand     = ft_preprocessing(cfgin,junk);
% % % % junk=databand;
dummydata=databand;
dummydata.trial={cell2mat(databand.trial)};
dummydata.time = {[1:size(dummydata.trial{1},2)]/databand.fsample};
dummydata.sampleinfo=[1 size(dummydata.trial{1},2)];


clear junk
% % cfg=[];
% % cfg.trials    = 'all'
% %  cfg.length    = 25;
% %     cfg.overlap = 0;
% % junk=ft_redefinetrial(cfg, dummydata);

ncountseg=dummydata.fsample*25;

cfg=[];

ntrial=ceil(dummydata.sampleinfo(1,2)/ncountseg);
begsample=(ones(1,ntrial)+([0:(ntrial-1)]*ncountseg))';
endsample=(([1:(ntrial)]*ncountseg))';
endsample(end,1)=dummydata.sampleinfo(1,2);
cfg.trl=[begsample endsample zeros(numel(endsample),1)]
% cfg.trials    = 'all'
% cfg.length    = 25;
%     cfg.overlap = 0;
junk=ft_redefinetrial(cfg, dummydata);

cfg=[];
cfg.covariance='yes';
cfg.keeptrials='no';
cfg.removemean='yes';
cfg.vartrllength=2;
datavg=ft_timelockanalysis(cfg,junk);


%----- Set initial source localization settings
srcCfg=[];
srcCfg.vol=vol;
srcCfg.method='lcmv';
srcCfg.fixedori='no';
srcCfg.feedback='text';
srcCfg.lambda = bf_lambda;
srcCfg.projectnoise='no';
srcCfg.keepfilter='yes';
srcCfg.keepleadfield='yes';
srcCfg.keepmom='no';


%--- Select only channels of current data set in already computed leadfields -----------
[indA,indB]=match_str(datavg.label,gridAllLF.label);
Nsrc=length(gridAllLF.inside);
gridCase=gridAllLF;
gridCase.leadfield(gridCase.inside)=cellfun(@(x,y) x(y,:),gridAllLF.leadfield(gridAllLF.inside),repmat({indB},1,Nsrc),'UniformOutput',false);

%----- Add grid to source localization settings
srcCfg.grid=gridCase;
srcCfg.grid.label=gridCase.label(indB);
%----- Compute Inverse Solution from Covariance matrix ------------

%  this provides the spatial filters that will project the data for the
%  time points and time periods of interest
datavg.grad=grad;
disp('beamforming starting...');
datavg.grad=grad;
source=ft_sourceanalysis(srcCfg,datavg);
disp('...beamforming ended');
% source.inside=[1:numel(source.inside)];

%===================================================================
%%

TC=cell2mat(databand.trial);
TC=TC.*10^15;

% TC=TC(:,1:size(datavg.time,2));
nsource = numel(find(source.inside));

pIC=size(TC,2);
source_sig=zeros(3,pIC);
sigt=zeros(1,pIC);

% difines step and window in points;
step_pnt=round(datain.fsample*step/1000);
window_pnt=round(datain.fsample*window/1000);

nwin=fix((size(sigt,2)-window_pnt)/step_pnt);

for k=1:nwin-1
    time_power(k)=(1/datain.fsample)*mean((k-1)*step_pnt+1:(k-1)*step_pnt+window_pnt);       % in seconds
end

power=zeros(nsource,nwin-1);
ft_progress('init', 'text',  'Please wait...');
str=['evaluating power for ' subject ' band ' band_prefix];
disp(str)

% create a sparse matrix that is essentially an averaging operator
% and prune the IC time course matrix
ix = zeros(window_pnt*(nwin-1),1);
iy = zeros(window_pnt*(nwin-1),1);
for k = 1:(nwin-1)
    indx     = (k-1)*window_pnt+(1:window_pnt);
    ix(indx) = (k-1)*step_pnt+(1:window_pnt);
    iy(indx) = k;
end
iz = ones(numel(iy),1)./window_pnt;
P  = sparse(ix,iy,iz);
TC = TC(:,1:size(P,1)); % the last bunch of samples are not used anyway
source_sig= source_sig(:,1:size(P,1));

inside_indices = find(source.inside(:))';

for is=1:numel(inside_indices)
    ft_progress(is/nsource, 'evaluating power voxel %d from %d\n', is, nsource);
    srcFilter=source.avg.filter{inside_indices(is)};
    source_sig(:,:)=srcFilter*TC;
    sigt=sum(source_sig.^2,1);
    power(is,:) = sigt*P;
end

% remove some unwanted stuff from the data structure that was output of hcp_icamne.m
if isfield(source,'avg'); source=rmfield(source,'avg'); end
if isfield(source,'time'); source=rmfield(source,'time'); end
if isfield(source,'val'); source=rmfield(source,'val'); end
if isfield(source,'snr'); source=rmfield(source,'snr'); end

% add some relevant stuff
source.power=power;
source.blp_band=blp_band;
source.step=step;
source.window=window;
source.step_pnt=step_pnt;
source.window_pnt=window_pnt;
source.time=time_power;
source.time_power=time_power;
end

