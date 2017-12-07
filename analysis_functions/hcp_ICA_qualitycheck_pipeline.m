function [ bad_segments, bad_channels , iter_res] = hcp_ICA_qualitycheck_pipeline(dataraw,options_ICA)

% HCP_ICA_QUALITYCHECK_PIPELINE allows the rejection of channels and intervals
% harmful for an Independent Component Analysis. It is babsed on an
% iterative procedure updating the list of bad segments and channels each
% iteration ending when no channels and intervals are found.
%
% Output data are vectors containing the bad channels and intervals
% selected.
%
% Use as
% [skint, bch, iter_res] = hcp_ICA_qualitycheck_pipeline(filename, options_ICA)
%
% where filename is a string that points to a raw data file in the database
% and options_ICA is a cell-array specifying the behaviour of the
% algorithm. Cell-arrays need to be organized as sets of key-value pairs,
% i.e. {'key1', 'value1', ...}.
%
% Options needs to contain the following key:
%   channels: channel selection , it should be a row cell
%   knownbadchannels:  bad channels already known to be bad, it should be a row cell -GIORGOS
%   skipped_intervals: Nx2 matrix of skipped intervals previusoly selected (t11 t12 ; t21 t22; ... ; tn1 tn2)
%   bandpass (optional): band pass intervals ([f1 f2] default [1 150])
%   bandstop (optional): band stop intervals ([fs1 fs2] default [59 61 ; 119 121])
%
% The following steps are performed:
% -reading in all data from disk
% -band pass and if required band stop filtering of the selected channels
% -FastICA analysis
% -Find Bad Channels according to the weight of the Mixing matrix
% -Find large artifact in the IC time courses according to a local variance
% -Selects the intervals to be cutted
%
% Example use:
% fname='0';
% options_ICA = {'channels', 'MEG', 'skipped_intervals', [], 'bandpass', [1 150], 'bandstop', [59 61 ; 119 121]};
% [skipped_intervals, bad_channels ] = hcp_ICA_qualitycheck_pipeline(fname,options_ICA)

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

resultprefix   = ft_getopt(options_ICA, 'resultprefix');
skintpt     = ft_getopt(options_ICA, 'skipped_intervals');        % skipped intervals previously identified
skint= skintpt/dataraw.fsample;
plotics   = ft_getopt(options_ICA, 'plotics',       'no');      % 'yes' or 'no' to plot the ICA result each iteration (default='no')
ch_method = ft_getopt(options_ICA, 'bad_ch_method', 'std_thr'); % 'std_thr' or 'max_weight' (default='std_thr')
selmode   = ft_getopt(options_ICA, 'mode',          'auto');    % skipped interval selection 'auto' or 'user' (default='auto')
saveit    = ft_getopt(options_ICA, 'iter_results');             % save the results of each iteration
modality   = ft_getopt(options_ICA, 'modality','MEG'); %Francesco
resamplefs     = ft_getopt(options_ICA, 'resamplefs');

% Add some options to the options_ICA cell-array that will be passed to
% lower level functions. Should these be configurable?
options_ICA = ft_setopt(options_ICA, 'save_comp',   'yes');
options_ICA = ft_setopt(options_ICA, 'ica_iterations', 1);

decfactor=1;
if(~isempty(resamplefs))
    decfactor=dataraw.fsample/resamplefs; % decimation factor (double)
end

%------------------------
% GIORGOS
if strcmp(modality,'MEG')
    layout='4D248.mat';
    threshold_std=10;
elseif strcmp(modality,'EEG')
    layout='EEG1010.lay';
    threshold_std=15; % The original threshold for MEG was 10 but many channels were identified as bad. Increased the threshold to be adapted to the EEG signal characteristics
end
%-----------------------------

% These are hard coded thresholds. Should we consider to make them
% configurable?
num_sig     = 12; % number of standard deviation out of the mean from the artifact detection
sic_thr_par = 5;  % factor of multiplication of the threshold for the determination of the intervals to be skipped


% size_s   = get(0,'ScreenSize');
size_s = [1,1 1200 700];
flag     = 1; % flag that checks convergence
it       = 0;  % iteration counter
bch2     = []; % some variable
segcount = 0;  % some variable
while flag
    clear datain
    disp('doing ICA preprocessing')
    datain = hcp_ICA_preprocessing(dataraw, options_ICA);
    disp('STARTING hcp_ICA_unmix');
    iteration = hcp_ICA_unmix(datain, options_ICA);
    disp('DONE hcp_ICA_unmix');
    comp   = iteration.comp;
    
    % normalize the time courses to unit standard deviation (and zero mean)
    % -> this shouldn't matter because the component time courses are by
    % definition zero mean
    disp('STARTING ft_channelnormalise');
    tmp        = ft_channelnormalise([], comp);
    disp('DONE ft_channelnormalise');
    comp.trial = tmp.trial;
    
    % perform spectral analysis on ICs
    options   = {'doplot', 'no', 'grad', dataraw.grad, 'plottype', 'summary'}; % perform frequency analysis
    disp('STARTING ft_ICA_freq');
    comp_freq = hcp_ICA_freq(comp, options);
    disp('DONE ft_ICA_freq');
    
    disp(['plotics:',plotics])
    if strcmp(plotics, 'yes')
        options = {'plottype', 'components'};
        disp('STARTING ft_ICA_plot');
        hcp_ICA_plot(comp_freq, options) % summary plots of the IC
        disp('DONE ft_ICA_plot');
    end
    
    % HACK JM work it back into the 1 trial case in order for the rest of the
    % code to work: I suggest however to clean this up at some point
    disp('STAGE POST A');
    fs = 0;
    nsmp = 0;
    for k = 1:numel(comp.trial)
        fs   = fs   + sum(diff(comp.time{k}));
        nsmp = nsmp + numel(comp.time{k}) - 1;
    end
    fs = nsmp./fs;
    
    %   timelim  = [comp.time{1}(1) comp.time{end}(end)];
    %   timeaxis = timelim(1):1/fs:timelim(2);
    if(isempty(resamplefs))
        nSamples = dataraw.sampleinfo(1,2);
    else
        nSamples = round(dataraw.sampleinfo(1,2)*(datain.fsample/dataraw.fsample));
    end
    disp('STAGE POST B');
    timelim  = [1 nSamples];
    timeaxis = [0:timelim(2)]/fs;
    
    Nc    = size(comp.topo,2);
    ntim  = numel(timeaxis);
    
    % Reconstruct the whole IC signal and power adding zeros in the skipped intervals
    disp('STAGE POST C');
    IC  = zeros(Nc, nSamples);
    pIC = zeros(Nc, nSamples);
    for i=1:numel(comp.trial)
        begindx = nearest(timeaxis, comp.time{i}(1));
        endindx = nearest(timeaxis, comp.time{i}(end));
        IC(:,begindx:endindx)  = comp_freq.trial{i};
        pIC(:,begindx:endindx) = comp_freq.pIC{i};
    end
    disp('STAGE POST D');
    % Smoothing of the IC Power Time Course (used in the peak analysis)
    sIC = zeros(size(pIC));
    for i=1:size(IC,1)
        % this relies on the curve fitting toolbox
        %     sIC(i,:)=(smooth(pIC(i,:),101))';
        sIC(i,:)=ft_preproc_smooth(pIC(i,:), 101);
    end
    
    
    %%%%%%%%%% Unworking channel identification %%%%%%%%%%%%%%%
    disp('STAGE POST E');
    disp('unworking channel identification...');
    
    % declare some variables
    ch_ic_list    = [];
    bad_channels  = [];
    bad_channels2 = [];
    size_fx  = size_s(1,3)/3;
    size_fy  = size_s(1,4)-(size_s(1,4)/5);
    border_x = (size_s(1,4)/20);
    A        = comp_freq.topo;
    for i=1:Nc
        clear chan_lab;
        flag_bad=0;
        [vt,order] = sort(abs(A(:,i))); % Sort the channels according to the weight on each IC
        vt = flipud(vt);
        order = flipud(order);
        
        if (strcmp(ch_method,'max_weight'))
            
            if vt(1)/vt(2) > 10 % One channel is discarded if the IC weight is 10 times grater than the others
                flag_bad=1;
                mth='val';
            end
        else
            %         [histo_ch bins_ch]=hist(abs(A(:,i)),400);
            mean_weight=prctile(vt,50,1);
            std_weight=std(vt(2:end));
            if vt(1)> mean_weight+threshold_std*std_weight % One channels is discarded if the IC weight is 10 times std out of the mean
                h_fig=figure;
                %         set(h_fig, 'Position', [(size_s(1,3)-(border_x+size_fx)) border_x size_fx size_fy])
                set(h_fig, 'visible', 'off','paperposition', [1 1 10 7]);
                [histo_ch,bin_ch]=hist(vt,400);
                hist(vt,400);
                xlabel('weight value ','color','k');
                title('Mixing Matrix Weight distribution','FontSize',12);
                line([mean_weight+10*std_weight mean_weight+10*std_weight], [0 max(histo_ch)],'color','r')
                pause(3);
                flag_bad=1;
                mth='std';
            end
        end
        if(flag_bad==1)
            posi=[border_x border_x size_fx size_fy];
            %       options={'plottype','components','component',i,'position',posi};
            %       hcp_ICA_plot(comp_freq,options) % summary plots of the IC
            if(strcmp(mth,'std'))
                n_chan=size(find(vt>mean_weight+10*std_weight),1);
                chan_lab(1,1:n_chan) = iteration(1,1).comp.topolabel(order(1:n_chan));
            else
                chan_lab = iteration(1,1).comp.topolabel(order(1));
            end
            str_chan=char(chan_lab(1,1));
            if(n_chan>1)
                for is=2:n_chan
                    str_chan=[str_chan, '-' ,char(chan_lab(1,is))];
                end
            end
            
            if(strcmp(selmode,'user'))
                string_1=['method ', mth , ': Do you want elminiate the channel/s  ', str_chan ];
                ButtonName2=questdlg(string_1, ...
                    'Question', ...
                    'Yes','No','Yes');
            else
                ButtonName2='Yes';
            end
            if strcmp(ButtonName2,'Yes')
                imgname=[resultprefix '_icaqc_badchannel_' str_chan];
                options={'plottype','components','component',i,'saveres','yes','grad', dataraw.grad,'modality',modality,'saveformat','png','fileout',imgname,'visible','off'};
                hcp_ICA_plot(comp_freq,options) % summary plots of the IC
                ch_ic_list=[ch_ic_list i];
                chan=strcat('-',chan_lab);
            else
                chan_lab=[];
            end
        else
            chan_lab=[];
        end
        if(~isempty(chan_lab))
            bad_channels=[bad_channels,chan];
            bad_channels2=[bad_channels2, chan_lab];
        end
        close all
    end
    
    if(isempty(bad_channels))
        disp('no channel identified by ICA');
    else
        disp('bad channels identified ');
        disp(bad_channels2);
    end
    
    if(~isempty(ch_ic_list))
        IC(ch_ic_list,:)=[]; sIC(ch_ic_list,:)=[];
    end
    
    %%%%%%% Find IC intervals with large artifacts according to a local variance %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% The variance of the time course is evaluated in different segments %%%%%%%%%%%%%%%
    disp('evaluating IC variances');
    
    ic_ind=[];
    [Nc data_length] = size(IC);
    bin=round(comp.fsample/decfactor); % width of the segment
    nstep=floor(data_length/bin); % number of segments
    clear std_IC mean_sIC max_std max_val max_ind
    
    for i=1:Nc
        for j=1:nstep
            std_IC(i,j)  =std(IC(i,(bin*(j-1)+1):(bin*j))); % standard deviation of the i-th segment
            mean_sIC(i,j)=mean(sIC(i,(bin*(j-1)+1):(bin*j))); % mean fo the smoothed pIC
        end
        max_std(i)=max(std_IC(i,:)); % find the segment with largest standard deviation
        [max_val(i),max_ind(i)] = max(abs(IC(i,:)));
    end
    
    for i=1:Nc
        [junk c0]=find(std_IC(i,:)==0);
        std_IC(:,c0)=[]; % delete sements with zero std
        [histo bins]=hist(std_IC(i,:),400); % create an histogram for the std
        mean_std=mean(std_IC(i,:));
        std_std=std(std_IC(i,:));
        
        if(max_std(i)>(mean_std+(num_sig*std_std/decfactor)))
            ic_ind=i;
            break
        else
            ic_ind=0;
        end
    end
    
    
    %%%%%%%%%%%%% Graphically or automatically select the time intervals to be excluded %%%%%%%%%
    int_vect=[];
    while (ic_ind)
        disp('performing bad segments selection ');
        
        ButtonName='No';
        while strcmp(ButtonName,'No')
            
            % Choosing a window around the selected artifact
            offset=10000;
            start_point=max_ind(ic_ind)-offset;
            end_point=max_ind(ic_ind)+offset;
            if(start_point<=0)
                start_point=1;
            end
            if(end_point>=size(IC,2))
                end_point=size(IC,2);
            end
            IC2=abs(IC(ic_ind,start_point:end_point));
            sIC2=sIC(ic_ind,start_point:end_point);
            [junk,max_ind2] = max(abs(IC2(1,:)));
            
            % Defining the interval to be cutted around the artifact peak by requiring that the
            % smooted power is over a given threshold
            flag_inf=0;
            flag_sup=0;
            clear inf_thr sup_thr
            ref_sIC = prctile(mean_sIC(ic_ind,:),50,2); % median of the smoothed IC power
            for i=1:20
                infseg=max_ind2-i*500; if(infseg<1), infseg=1; end
                supseg=max_ind2+i*500; if(supseg>size(IC2,2)), supseg=size(IC2,2); end
                %         [junk over_thr]=find(sIC2(infseg+1:supseg)>sic_thr_par*ref_sIC);
                [junk over_thr]=find(sIC2(infseg:supseg)>sic_thr_par*ref_sIC);
                
                if(~isempty(over_thr))
                    if(flag_inf==0)
                        inf_thr(i)=min(over_thr)+infseg;
                        if(i>1)
                            if(inf_thr(i)==inf_thr(i-1)), flag_inf=1; end
                        end
                    end
                    if(flag_sup==0)
                        sup_thr(i)=max(over_thr)+infseg;
                        if(i>1)
                            if(sup_thr(i)==sup_thr(i-1)), flag_sup=1;end
                        end
                    end
                else
                    inf_thr(i)=infseg;
                    flag_inf=1;
                    sup_thr(i)=supseg;
                    flag_sup=1;
                end
                if(i>1)
                    if(flag_inf==1 && flag_sup==1)
                        break
                    end
                end
            end
            
            % Plots the selected artifact interval
            disp('plotting bad segments ');
            
            ymax=max(IC2);
            points_all=(1:size(IC,2));
            
            cfg=[];
            cfg.component=ic_ind;
            cfg.layout=layout;
            
            h1_fig=figure;
            %       set(h1_fig, 'Position', [75 75 size_s(1,3)-150 size_s(1,4)-150])
            set(h1_fig, 'visible', 'off','paperposition', [1 1 10 7]);
            
            subplot(2,3,3) ; plot(points_all,IC(ic_ind,:))
            hold on
            plot(points_all(1,start_point:end_point),IC(ic_ind,start_point:end_point),'k')
            legend('IC time course','interval around the peak','Location','NorthOutside')
            subplot(2,3,4)
            ft_topoplotIC(cfg, comp_freq)
            subplot(2,3,5); plot(bins,histo)
            line([mean_std+10*std_std mean_std+10*std_std], [0 max(histo)],'Color','r')
            legend('std distribution','mean and peak thr','Location','NorthOutside')
            line([mean_std mean_std], [0 max(histo)],'Color','r')
            %             subplot(2,3,[1 2]) ; plot(IC2,'k')
            subplot(2,3,6) ; plot(sIC2,'k')
            line([1 size(IC2,2)], [sic_thr_par*ref_sIC sic_thr_par*ref_sIC],'Color','r')
            axis([1 length(sIC2) 0 max(sIC2)])
            subplot(2,3,[1 2]) ; plot(IC2,'k')
            hold on
            plot(min(inf_thr):min(max(sup_thr),size(IC2,2)),IC2(1,min(inf_thr):min(max(sup_thr),size(IC2,2))),'r')
            legend('abs of the IC','suggested cut')
            axis([1 length(IC2) 0 ymax])
            title('Select time interval to be excluded...');
            
            
            %             time_points=[1:cut_interval(1) cut_interval(2):length(IC)];
            if(strcmp(selmode,'user'))
                
                [X,Y] = ginput(2);
                cut_interval=round(X)+start_point;
                
                subplot(2,3,3) ; plot(points_all,IC(ic_ind,points_all))
                hold on
                plot(points_all(1,start_point:end_point),IC(ic_ind,start_point:end_point),'k')
                plot(points_all(1,cut_interval(1,1):cut_interval(2,1)),IC(ic_ind,cut_interval(1,1):cut_interval(2,1)),'g')
                legend('IC time course','interval around the peak','selected interval','Location','NorthOutside')
                
                ButtonName=questdlg('Are you satisfied with the interval selection?', ...
                    'Question', ...
                    'Yes','No','skip','Yes');
            else
                segcount=segcount+1;
                ButtonName='Yes';
                cut_interval=[min(inf_thr)+start_point;max(sup_thr)+start_point];
                imgname=[resultprefix '_icaqc_badsegment_' num2str(segcount)];
                hcp_write_figure(imgname, h1_fig);
            end
            close(h1_fig);
            
        end
        close all
        disp('Done bad segments selection ');
        
        if(cut_interval(2)-cut_interval(1)<100)
            cut_interval(1)=cut_interval(1)-100;
            cut_interval(2)=cut_interval(2)+100;
        end
        if(cut_interval(1)<=0)
            cut_interval(1)=1;
        end
        if(cut_interval(2)>size(IC,2))
            cut_interval(2)=size(IC,2);
        end
        IC(:,cut_interval(1):cut_interval(2))=0;
        sIC(:,cut_interval(1):cut_interval(2))=0;
        if(strcmp(ButtonName,'Yes'))
            int_vect=[int_vect; cut_interval(1) cut_interval(2)];
        end
        
        clear  std_IC mean_sIC max_std max_val max_ind
        for i=1:Nc
            for j=1:nstep
                std_IC(i,j)=std(IC(i,(bin*(j-1)+1):(bin*j)));
                mean_sIC(i,j)=mean(sIC(i,(bin*(j-1)+1):(bin*j)));
            end
            max_std(i)=max(std_IC(i,:));
            [max_val(i),max_ind(i)] = max(abs(IC(i,:)));
        end
        
        for i=1:Nc
            [junk c0]=find(std_IC(i,:)==0);
            std_IC(:,c0)=[];
            [histo bins]=hist(std_IC(i,:),400);
            mean_std=mean(std_IC(i,:));
            std_std=std(std_IC(i,:));
            
            if(max_std(i)>(mean_std+(num_sig*std_std)))
                ic_ind=i;
                break
            else
                ic_ind=0;
            end
        end
        
    end
    disp('Verify the iteration results ');
    
    %%%%%%%%%%%%% Verify the iteration results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    flag_int=0;
    if(size(int_vect))
        int_vect_time=int_vect/comp_freq.fsample;
        skint=vertcat(skint,int_vect_time);
        [junk ord]=sort(skint(:,1));
        skint=skint(ord,:);
        i=0;
        flag_sk=0;
        while i<size(skint,1)-1+flag_sk
            flag_sk=0;
            for i=1:(size(skint,1)-1)
                if (skint(i,2)>skint(i+1,1) && skint(i,2)<skint(i+1,2))
                    skint(i,2)=skint(i+1,2);
                    skint(i+1,:)=[];
                    flag_sk=1;
                    break
                end
                if (skint(i,2)>skint(i+1,1) && skint(i,2)>skint(i+1,2))
                    skint(i+1,:)=[];
                    flag_sk=1;
                    break
                end
            end
        end
        flag_sk=0;
        i=0;
        while i<size(skint,1)-1+flag_sk
            flag_sk=0;
            for i=1:(size(skint,1)-1)
                if (skint(i+1,1)-skint(i,2) < 1)
                    skint(i,2)=skint(i+1,2);
                    skint(i+1,:)=[];
                    flag_sk=1;
                    break
                end
            end
        end
        
        %modified in order to properly manage bad segments in samples. Should
        %be done better
        options_ICA = ft_setopt(options_ICA, 'skipped_intervals', skint*dataraw.fsample);
        flag_int=1;
    end
    disp(' arranging temporary results ');
    
    flag_bch=0;
    if(~isempty(bad_channels))
        flag_bch=1;
        pv_bad_ch=ft_getopt(options_ICA, 'knownbadchannels'); % GIORGOS
        bad_channels=horzcat(pv_bad_ch,bad_channels);
        options_ICA = ft_setopt(options_ICA, 'knownbadchannels', bad_channels); % GIORGOS
        bch2=horzcat(bch2,bad_channels2);
    end
    flag=(flag_bch+flag_int)-(flag_bch*flag_int);
    pause(3);
    it=it+1;
    iter_res(1,it).bch=ft_getopt(options_ICA, 'knownbadchannels'); % GIORGOS
    iter_res(1,it).bch2=bch2;
    iter_res(1,it).skint=ft_getopt(options_ICA, 'skipped_intervals');
end
disp(' Finishing ICAqc ');
bch= ft_getopt(options_ICA, 'knownbadchannels'); % GIORGOS

bad_segments=round(skint*dataraw.fsample);
bad_channels=unique(bch2);
if isempty(bad_channels)
    bad_channels=[]; % This is because in Matlab 2013a unique of an empty matrix return a 0x1 empty matrix
end
imgname=[resultprefix '_icaqc_results' ];
options={'plottype','components','saveres','yes','grad', dataraw.grad,'modality',modality,'saveformat','png','fileout',imgname,'visible','off'};
hcp_ICA_plot(comp_freq,options) % summary plots of the IC
options={'plottype','components','saveres','yes','grad', dataraw.grad,'modality',modality,'saveformat','fig','fileout',imgname,'visible','off'};
hcp_ICA_plot(comp_freq,options) % summary plots of the IC
disp(' ICAqc ended ');
end
