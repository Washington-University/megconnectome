function hcp_ica_plotreport(comp,options_ICA,datain)

% HCP_ICA_PLOTREPORT allows the plots of Independent
% Components automatically classified by HCP_ICA_RMEG_CLASSIFICATION function. This
% function is able to classify Independent Components into Brain Activity
% and ECG/EOG related artifacts provided that ECG and EOR reference
% channels are used.
%
% Use as
%   hcp_ica_plotreport(comp,options_ICA,datain)
%
%     where the input comp structure should be obtained from
%     HCP_ICA_RMEG_CLASSIFICATION and the datain structure is a fieldtrip data
%     structure contining the channel time courses used as imput of the
%     FT_COMPONENTANALYSIS .
%
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

resultprefix =ft_getopt(options_ICA, 'resultprefix');
subject      = ft_getopt(options_ICA, 'subject');
bandpass     = ft_getopt(options_ICA, 'bandpass');
bandstop     = ft_getopt(options_ICA, 'bandstop');
tsize        = ft_getopt(options_ICA, 'textsize',22)
modality     = ft_getopt(options_ICA, 'modality', 'MEG'); % GIORGOS
grad=ft_getopt(options_ICA, 'grad');

if strcmp(modality,'MEG')
    layout='4D248.mat';
elseif strcmp(modality,'EEG')
    layout='EEG1010.lay';
end

% if(~isfield(comp,'trial'))
%     for i=1:size(datain.trial,2)
%         comp.trial{i} = comp.unmixing*datain.trial{i};
%     end
% end
% 
% comp.time={cell2mat(datain.time)};


if(~isfield(comp,'pow_data'))
    doplot='no'; 
    options = {'doplot',doplot, 'grad', grad,'plottype','components'};
    
    [comp_freq]=hcp_ICA_freq(comp, options, datain);
    comp.pow_data=comp_freq.pow_data;
    comp.freq_comp=comp_freq.freq_comp;
end
n_IC=size(comp.unmixing,1);

frange=bandpass; 

cfgin =[];
cfgin.grad=grad;
cfgin.zlim='maxabs';
cfgin.component=[];
cfgin.parameter = 'topo';
cfgin.comment = 'no';
cfgin.colorbar = 'yes';


mspec = sqrt(comp.freq_comp.powspctrm);
F = comp.freq_comp.freq;
pow_data=comp.pow_data;

pIC = cell2mat(pow_data.trial).^2;
time_pIC=cell2mat(pow_data.time);


leg={'OneOverF (>0.5)  ' ; 'Specflat (<2.0)  ' ; 'Kurtosis (>15)   ' ; 'elecSig (>0.1)  '; 'elecPow (>0.25) '; 'elecSpe (>0.95) '};
thrs=[0.5;2.0;15;0.1;0.25;0.95];
class=comp.class;

for ix=1:n_IC
    str(1,:)={num2str(abs(comp.class.spectrum_one_over_f(1,ix)),'%5.3f')};
    str(2,:)={num2str(abs(comp.class.spectrum_flat(1,ix)),'%5.1f')};
    str(3,:)={num2str(abs(comp.class.time_kurtosis(1,ix)),'%5.1f')};
    str(4,:)={num2str(abs(comp.class.elc_signal_correlation(1,ix)),'%5.3f')};
    str(5,:)={num2str(abs(comp.class.elc_power_correlation(1,ix)),'%5.3f')};
    str(6,:)={num2str(abs(comp.class.elc_spectrum_correlation(1,ix)),'%5.3f')};
    strcolor=['k' ;'k'; 'k'; 'k'; 'k'; 'k'];
    for i=[1 3 4 5 6]
        if(str2num(str{i})>thrs(i)) strcolor(i)='b';end
    end
    if(str2num(str{2})<thrs(2)) strcolor(2)='b';end
    
    
    str_ic = 'ARTIFACT';
    
    if(~isempty(find(class.brain_ic==ix)) && ~isempty(find(class.brain_ic_vs==ix)))
        str_ic='BRAIN';
    elseif(~isempty(find(class.physio==ix)))
        str_ic='PHYSIOLOGICAL ARTIFACT CORRECTED';
    elseif(isempty(find(class.brain_ic==ix)) && ~isempty(find(class.ecg_eog_ic==ix)) && isempty(find(class.physio==ix)))
        str_ic='EOG/ECG RELATED ARTIFACT';
    elseif(~isempty(find(class.brain_ic==ix)) && isempty(find(class.brain_ic_vs==ix)) && isempty(find(class.ecg_eog_ic==ix)))
        str_ic='ARTIFACT CORRECTED';
    elseif(isempty(find(class.brain_ic==ix)) && ~isempty(find(class.brain_ic_vs==ix)))
        str_ic='BRAIN CORRECTED';
    end
   
    
    figure('Menubar','none');
    %
    set(gca, 'XTickLabel',[], 'YTickLabel',[], ...
        'Units','normalized', 'Position',[0 0 1 1])

    set(gcf, 'visible', 'off')
    set(gcf, 'PaperUnits','inches')

    set(gcf, 'paperposition', [1 1 20 12]);
    
    subplot(2,3,1)
    cfgin.component=ix;
    cfgin.layout=layout;
    cfgin.colorbar='yes';
    ft_topoplotIC(cfgin, comp);
    title('IC sensor map','FontSize',12);
    
    subplot(2,3,4)
    plot(F,mspec(ix,:));
    title(['IC power spectral density']);
    Prms = mspec(ix,:);
    lags = find(F>1 & F<frange(1,2));
    max_val = max(Prms(lags))+0.1*(max(Prms(lags))-min(Prms(lags)));
    axis([0 75 0 max_val]);
    ylabel(['power']);
    xlabel(['frequency (Hz)'],'color','k');
    
    str_icnum = {['ic ' num2str(ix) ' = ' str_ic]};
    
 
    text(1.2 ,0.65,leg{1},'Units','normalized','FontSize',tsize,'color',strcolor(1))
    text(2,0.65, str{1},'Units','normalized','FontSize',tsize,'color',strcolor(1))
    text(1.2,0.59,leg{2},'Units','normalized','FontSize',tsize,'color',strcolor(2))
    text(2,0.59,str{2},'Units','normalized','FontSize',tsize,'color',strcolor(2))
    text(1.2,0.53,leg{3},'Units','normalized','FontSize',tsize,'color',strcolor(3))
    text(2,0.53,str{3},'Units','normalized','FontSize',tsize,'color',strcolor(3))
    text(1.2,0.465,leg{4},'Units','normalized','FontSize',tsize,'color',strcolor(4))
    text(2,0.465,str{4},'Units','normalized','FontSize',tsize,'color',strcolor(4))
    text(1.2,0.40,leg{5},'Units','normalized','FontSize',tsize,'color',strcolor(5))
    text(2,0.40,str{5},'Units','normalized','FontSize',tsize,'color',strcolor(5))
    text(1.2,0.335,leg{6},'Units','normalized','FontSize',tsize,'color',strcolor(6))
    text(2,0.335,str{6},'Units','normalized','FontSize',tsize,'color',strcolor(6))
    
    text(1.2,0.8,str_icnum,'Units','normalized','FontSize',tsize,'color','b')
    
    subplot(2,3,[2 3])
    timeic=cell2mat(comp.time);
    trialic=cell2mat(comp.trial);
    plot(timeic(1,:),trialic(ix,:));
    sig = trialic(ix,:);
    base = timeic;
    axis([base(1) base(end) min(sig)-0.1*(max(sig)-min(sig)) max(sig)+0.1*(max(sig)-min(sig))]);
    ylabel(['amplitude']);
    xlabel(['time (Sec)']);
    title('IC time course','FontSize',12);
    
    subplot(2,3,6)
    plot(time_pIC,pIC(ix,:));
    axis([time_pIC(1) time_pIC(end) 0 2*max(pIC(ix,:))]);
    ylabel(['power']);
    xlabel(['time (sec)']);
    title('IC power timecourse','FontSize',12);
    
    
    imgname=[resultprefix  '_icaclass_vs_' num2str(ix) '.png'];
    hcp_write_figure(imgname, gcf, 'resolution', 300)
    
    close
end

end
