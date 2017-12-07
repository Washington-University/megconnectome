function [ source ] = hcp_ica_blp(source, comp, options_blp)

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

if strcmp(band_prefix,'betalow')
    N_hp=9; N_lp=10 ;
elseif strcmp(band_prefix,'betahigh')
    N_hp=9; N_lp=12 ;
elseif strcmp(band_prefix,'alpha')
    N_hp=7; N_lp=8 ;
elseif strcmp(band_prefix,'theta')
    N_hp=6; N_lp=7 ;
elseif strcmp(band_prefix,'delta')
    N_hp=6; N_lp=7 ;
elseif strcmp(band_prefix,'gammalow')
    N_hp=9; N_lp=13 ;
elseif strcmp(band_prefix,'gammamid')
    N_hp=11; N_lp=14 ;
elseif strcmp(band_prefix,'gammahigh')
    N_hp=13; N_lp=14 ;
end

step   = ft_getopt(options_blp, 'blp_step',    20);
window = ft_getopt(options_blp, 'blp_window', 400);
band_str=[num2str(blp_band(1,1)) '-' num2str(blp_band(1,2))];

disp(['blp band -> ' band_str ' Hz'])
disp(['power calculation step -> ' num2str(step) ' ms'])
disp(['power calculation window -> ' num2str(window) ' ms'])

if(~isfield(comp.class,'brain_ic_vs')) 
    comp.class.brain_ic_vs=comp.class.brain_ic;
end

brainic_indx=comp.class.brain_ic_vs;
nsource = numel(source.inside);


% filtering

if(strcmp(band_prefix,'delta'))
    cfgin          = [];
    cfgin.lpfilter = 'yes';
    cfgin.lpfreq   = blp_band(1,2);
    cfgin.lpfiltord=N_lp;
    comp_blp     = ft_preprocessing(cfgin,comp);
elseif(strcmp(band_prefix,'whole'))
    comp_blp=comp;
else
    cfgin          = [];
    cfgin.hpfilter = 'yes';
    cfgin.hpfiltord=N_hp;
    cfgin.hpfreq   = blp_band(1,1);
    junk     = ft_preprocessing(cfgin,comp);
    
    cfgin          = [];
    cfgin.lpfilter = 'yes';
    cfgin.lpfiltord=N_lp;
    cfgin.lpfreq   = blp_band(1,2);
    comp_blp     = ft_preprocessing(cfgin,junk);
end
clear junk

if strcmp(dofiltcheck,'yes')
    compappo=comp;
    compappo.trial=comp_blp.trial;
    options   = {'doplot', 'no', 'grad', comp.grad, 'plottype', 'components'};
    disp('STARTING ft_ICA_freq');
    comp_freq = hcp_ICA_freq(compappo, options);
    disp('DONE ft_ICA_freq');
    
    options   = {'doplot', 'no', 'grad', comp.grad, 'plottype', 'components'};
    disp('STARTING ft_ICA_freq');
    comp_freq2 = hcp_ICA_freq(comp, options);
    disp('DONE ft_ICA_freq');
    
    imgname=[subject '_blp_iccheck_' band_prefix];
    options={'plottype','components','component',9,'saveres','no','grad', comp.grad,'modality','MEG','saveformat','png','fileout',imgname,'visible','on'};
    hcp_ICA_plot(comp_freq2,options) % summary plots of the IC
    
    mspec     = sqrt(comp_freq.freq_comp.powspctrm(9,:));
    F         = comp_freq.freq_comp.freq;
    
    subplot(2,2,[3 4])
    hold on
    plot(F,mspec,'r')
end

junk=cell2mat(comp_blp.trial);

IC=junk(brainic_indx,:);
pIC=size(IC,2);
source_sig=zeros(3,pIC);
sigt=zeros(1,pIC);


% difines step and window in points;
step_pnt=round(comp.fsample*step/1000);
window_pnt=round(comp.fsample*window/1000);

nwin=fix((size(sigt,2)-window_pnt)/step_pnt);

for k=1:nwin-1
    time_power(k)=(1/comp.fsample)*mean((k-1)*step_pnt+1:(k-1)*step_pnt+window_pnt);       % in seconds
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
IC = IC(:,1:size(P,1)); % the last bunch of samples are not used anyway
source_sig= source_sig(:,1:size(P,1));

for is=1:nsource
    ft_progress(is/nsource, 'evaluating power voxel %d from %d\n', is, nsource);
    source_sig(:,:)=source.avg.mom{source.inside(is)'}(:,brainic_indx)*IC;
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

