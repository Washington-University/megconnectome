function [comp_best] = hcp_ICA_RMEG_classification(ref_data,options,iteration,datain)

% FUNCTION hcp_ICA_RMEG_classification
%
% hcp_ICA_RMEG_classification performs the classification of the output
% of an Independent Component Analyis arranged in an "iteration" structure.
% This structure must contains "comp" fields coming from an iterative procedure
% using ft_componentanalysis (like hcp_doICARMEG).  In general the iteration structure
% contains a number of "comp" fields as many as the number of ICA iteration performed.
%
% Use as
%   ICA = hcp_ICA_RMEG_classification(options, iteration, datain)
%
% where datain is data structure used as input in the ft_componentanalysis
% and options is a cell-array specifying the behaviour of the
% algorithm. Cell-arrays need to be organized as sets of key-value pairs,
% i.e. {'key1', 'value1', ...}.
%
% Options needs to contain the following key:
%   dataset: dataset containing the reference signals
%   ref_channels: reference channel labels (default hcp-labels {'E31', 'E32'})
%
% Plotting and storing options
%   showbestit: summary plots of the selected best iteration (default='yes')
%   plottype:   plot layout (summary or components, default = 'sumamry' [see hcp_plotICA])
%   saveres:    save the summary plot (default = 'no') using the 'saveformat' file extention (defualt = 'fig')
%   store_bestcomp:   store the comp structure for the best iteration with the time courses and
%                the time vector (default = 'yes')
%
% Otions by default retrieved from the comp structure to be modified with care if required in particular case
%   bandpass: pass-band interval ([f1 f2])
%   bandstop: stop-band intervals ([fs1 fs2])
%   grad: gradiometer structure
%
% Example use:
% options = {'dataset', '0' , 'ref_channels', {'E31', 'E32'}, 'saveres', 'yes', 'saveformat', 'jpg'};
% ICA = hcp_ICA_RMEG_classification(options,iteration,megdata);

% set the defaults
% cfgin.dataset      = ft_getopt(options, 'dataset');
% cfgin.dataset=fname;
subject   = ft_getopt(options, 'subject');

cfgin.ref_channels      = ft_getopt(options, 'ref_channels');
cfgin.bpfreq        = ft_getopt(options, 'bandpass');
cfgin.bsfreq         = ft_getopt(options, 'bandstop');
cfgin.skipped_intervals = ft_getopt(options, 'skipped_intervals');
cfgin.grad         = ft_getopt(options, 'grad');

cfgin.showbestit  = ft_getopt(options, 'showbestit');
cfgin.plottype = ft_getopt(options, 'plottype');

cfgin.store_bestcomp     = ft_getopt(options, 'store_bestcomp');
cfgin.saveres  = ft_getopt(options, 'saveres');
cfgin.saveformat     = ft_getopt(options, 'saveformat');

if isempty(cfgin.ref_channels)  , cfgin.ref_channels = {'E31', 'E32'}; end
if isempty(cfgin.bpfreq) , cfgin.bpfreq = iteration(1,1).comp(1,1).bandpass ; end
if isempty(cfgin.bsfreq) , cfgin.bsfreq = iteration(1,1).comp(1,1).bandstop ; end
if isempty(cfgin.showbestit), cfgin.showbestit = 'no'; end
if isempty(cfgin.plottype)  , cfgin.plottype = 'summary'; end
if isempty(cfgin.store_bestcomp)    , cfgin.store_bestcomp = 'yes';             end
if isempty(cfgin.saveres)    , cfgin.saveres = 'no';             end
if isempty(cfgin.saveformat) , cfgin.saveformat = 'fig';         end
cfgin.zlim='maxabs';

nref=0;

if isfield(cfgin,'ref_channels')
    
    nref= size(ref_data.trial{1},1);
    ntrl= size(ref_data.trial,2)
    
    %   ----- Normalize ECG and EOG data
    for i = 1:nref
        for k=1:ntrl
            sig=ref_data.trial{k}(i,:);
            if std(sig) > 0
                sig=(sig-mean(sig))/std(sig);
                ref_data.trial{k}(i,:)=sig;
            end
        end
    end
    elec=cell2mat(ref_data.trial);
    
    h=figure;
    set(h, 'visible', 'off','paperposition', [1 1 10 7]);
    for i = 1:nref
        subplot(ceil(nref/2),2,i); plot(cell2mat(ref_data.time),elec(i,:)); axis tight;
        ylabel(ref_data.label(i));
        xlabel('seconds');
        title('Reference data timecourses','FontSize',12);
    end
    imgname=[subject '_icaclass_refch'];
    hcp_write_figure(imgname, h);
    clear h;
    
    %   ----- Time Course of Power of Electric Channels -----------------------
    
    cfg = [];
    
    cfg.bpfilter      = 'yes';
    cfg.bpfreq        = cfgin.bpfreq;
    cfg.hilbert='abs';
    
    pow_data                 = ft_preprocessing(cfg,ref_data);
    
    pelec = cell2mat(pow_data.trial).^2;
    
    %   ----- Power Spectral Estimation of Electric Channels ------------------
    
    win=2048/(ref_data.fsample*2^(round(log2(1024/ref_data.fsample))));
    
    cfg = [];
    cfg.length               = win;
    ref_data_seg                 = ft_redefinetrial(cfg, ref_data);
    
    cfg=[];
    cfg.method = 'mtmfft';
    cfg.output = 'pow';
    cfg.trials = 'all';
    cfg.taper = 'hamming';
    
    [freq_ref_data] = ft_freqanalysis(cfg,ref_data_seg);
    
    selec = sqrt(freq_ref_data.powspctrm);
end


n_iter=size(iteration,2);
for j=1:n_iter
    clear pow_data A W IC pIC pow_f_IC freq_comp_tr good_ic mspettro G G2 G3 H H2 H3 K K2 K3R R2 R3 ave T  T1 T2 ave_on ave_off clas potenza order spec_on spec_off
    
    comp_iter = iteration(1,j).comp;
    comp_iter.time=datain.time;
    for it=1:ntrl
        comp_iter.trial{it}=comp_iter.unmixing*datain.trial{it};
    end
    Nc = size(comp_iter.topo,2);
    
    
    
    %------ POWER SPECTRAL DENSITY ESTIMATION  -----
    cfgold=comp_iter.cfg;
    comp_iter= hcp_ICA_freq(comp_iter, options, datain);
    comp_iter.cfgold=cfgold;
    freq_comp =comp_iter.freq_comp;
    temp_struc(j).freq_comp=freq_comp;
    
    mspettro = sqrt(freq_comp.powspctrm);
    F = freq_comp.freq;
    
    %--------------------------------
    %----> IC POWER TIME COURSE <----
    %--------------------------------
    
    pow_data=comp_iter.pow_data;
    temp_struc(j).pow_data=pow_data;
    pIC = cell2mat(pow_data.trial).^2;
    
    %-------------------------------------
    %--------> IC CLASSIFICATION <--------
    %-------------------------------------
    
    if Nc==0
        good_ic=[]; art_con=1; G=[]; G2=[]; G3=[];
        H=[]; H2=[]; H3=[];
    else
        
        %   -----Correlation between ICs and electric channels <---
        
        IC=cell2mat(comp_iter.trial);
        if isfield(cfgin,'ref_channels')
            if nref > 0
                
                %   ----- Correlation between the ICs and electric channels
                R = corrcoef([elec',IC']);
                K = R(nref+1:nref+size(comp_iter.topo,2),1:nref);
                H = max(abs(K'),[],1);
                
                %   ----- Correlation between the ICs' and electric channels' power in window of 400 ms duration (step = 20 ms)
                R2 = corrcoef([pelec',pIC']);
                K2 = R2(nref+1:nref+size(comp_iter.topo,2),1:nref);
                H2 = max(abs(K2'),[],1);
                
                %   ----- Correlation between ICs' electric channels'power spectra between 3 and 40 Hz
                
                vspec = find(F > 3 & F < 40);
                R3 = corrcoef([selec(:,vspec)',mspettro(:,vspec)']);
                K3 = R3(nref+1:nref+size(comp_iter.topo,2),1:nref);
                H3 = max(abs(K3'),[],1);
            else
                H = zeros(1,size(comp_iter.topo,2)); H2=zeros(1,size(comp_iter.topo,2)); H3=zeros(1,size(comp_iter.topo,2));
            end
        else
            H=zeros(1,Nc);
            H2=zeros(1,Nc);
            H3=zeros(1,Nc);
        end
        
        
        %   ----- Fit with 1/f spectrum, flat spectrum "quantification", IC time kurtosis
        
        frq = [cfgin.bpfreq(1,1) 13];
        frq2 = [2 78];
        vec = find(F > frq(1) & F < frq(2));
        vec2 = find(F > frq2(1) & F < frq2(2));
        for ix = 1:Nc
            sspettro = mspettro(ix,:);
            sspettro_fit=sspettro./sspettro(vec(1,1));
            [fitx,goodness] = hcp_fitspectrum(F(vec)',sspettro_fit(vec)');
            if(fitx.a>0)
                G(ix) = goodness.rsquare;
            else
                G(ix)=0.11;
            end
            ss = sspettro(vec2);
            G2(ix)= prctile(ss,95)/(median(ss));
            G3(ix)= kurtosis(IC(ix,:));
        end
        
        %   ---- Set thresholds for ICs' classification
        thres_a = 0.10; % Elc-IC Signal Correlation
        thres_b = 0.25; % Elc-IC Power Correlation
        thres_c = 0.95; % Elc-IC Spectrum Correlation
        thres_d = 0.5;  % PSD 1/f
        thres_e = 2.02; % Spectrum Flat
        thres_f = 15;   % IC Time Kurtosis
        
        
        %----------------------------------------------------------
        %       classification: 0=artifact   1=brain signal
        %----------------------------------------------------------
        
        clas = zeros(1,Nc);
        
        clas_ecg_eog = zeros(1,Nc);
        
        %------ Find correlation with ECG and EOG
        clas_ecg_eog(find(H > thres_a | H2 > thres_b | H3 > thres_c)) = 1;  %---> Elc_signal_correlation - Elc_power_correlation - Elc_spectrum_correlation
        %------
        ecg_eog_ic = find(clas_ecg_eog==1);
        
        
        clas(find(G2 > thres_e & G3 < thres_f & H3 < thres_c)) = 1; %---> Spectrum_flat - Time_kurtosis - Elc_spectrum_correlation
        
        clas(find(G > thres_d | H > thres_a | H2 > thres_b)) = 0;   %---> Spectrum_one_over_f - Elc_signal_correlation - Elc_power_correlation
        
        good_ic = find(clas==1);
        
        if nref > 0
            art_con = mean(H2(good_ic));
        else
            art_con = 0;
        end
        
    end
    
    
    
    %    -------- Save Classification Parameters --------------
    
    
    iteration(j).comp.class.total_ic_number=size(comp_iter.topo,2);
    iteration(j).comp.class.brain_ic_number=length(good_ic);
    iteration(j).comp.class.brain_ic=good_ic;
    iteration(j).comp.class.ecg_eog_ic=ecg_eog_ic;
    iteration(j).comp.class.artifact_contamination=art_con;
    iteration(j).comp.class.elc_signal_correlation=H;
    iteration(j).comp.class.elc_power_correlation=H2;
    iteration(j).comp.class.elc_spectrum_correlation=H3;
    iteration(j).comp.class.spectrum_one_over_f=G;
    iteration(j).comp.class.spectrum_flat=G2;
    iteration(j).comp.class.time_kurtosis=G3;
    
    goodness_ic(j,1) = length(good_ic);
    goodness_ic(j,2) = 1-art_con;
    
    ica_nonlinearity=comp_iter.cfgold.fastica.g;
    disp([' run ' num2str(j) ' \ ' num2str(n_iter)]);
    disp(['total number of ICs : ' num2str(size(comp_iter.topo,2)) ' -  number of brain ICs : '  num2str(length(good_ic)) ' -  artifact contamination : ' num2str(art_con)]);
    
end

ICA.iteration = iteration;

clear comp_iter A W IC pIC good_ic mspettro G G2 G3 H H2 H3 K K2 K3 R R2 R3 ave ave_on ave_off T T1 T2 clas  spec_on spec_off;


% ------- Save Best Iteration Data --------------

[tmp,pos_new]=max(goodness_ic(:,1).*goodness_ic(:,2));

W = ICA.iteration(pos_new).comp.unmixing;
good_ic = ICA.iteration(pos_new).comp.class.brain_ic;
comp=ICA.iteration(pos_new).comp;

comp_best=ICA.iteration(pos_new).comp;

comp.time=datain.time;
for iic=1:size(datain.trial,2)
    comp.trial{iic}=W*datain.trial{iic};
end
comp.freq_comp=temp_struc(pos_new).freq_comp;
comp.pow_data=temp_struc(pos_new).pow_data;
comp.grad=cfgin.grad;

% ------- Display Best Iteration --------------
if(strcmp(cfgin.showbestit,'yes'))
    hcp_ICA_plot(comp,options);
end

if ~isempty(W)
    ICA.best_iter.index = pos_new;
    ICA.best_iter.nonlinearity = ica_nonlinearity;
    ICA.best_iter.brain_ic_number = length(good_ic);
    ICA.best_iter.brain_ic = good_ic;
    if strcmp(cfgin.store_bestcomp,'yes')
        ICA.best_iter.comp=comp;
    end
else
    disp('!!!ATTENTION:ICA processing completed, but no brain IC!');
end

end
