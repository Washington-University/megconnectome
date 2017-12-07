function [connect_mcw] = hcp_icaemcw_estimate(source, sourcemodel, options)

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

subject= ft_getopt(options, 'subjectid');
outputfile = ft_getopt(options, 'outputfile');
label=ft_getopt(options, 'label');
step = ft_getopt(options, 'mcw_step', 200);
window     = ft_getopt(options, 'mcw_window', 10000);
dofig = ft_getopt(options, 'dofig', 'yes');
ac=ft_getopt(options, 'ac', 0.1);
thresh_tn=ft_getopt(options, 'thresh', 0.6);
extnode=ft_getopt(options, 'extnode');
% outputfile=[outputfile '_test']

disp('reading mri template')
mri=ft_read_mri('mni_icbm152_t1_tal_nlin_sym_09a.nii');
mri=ft_convert_units(mri, 'mm');

sourcemodel=ft_convert_units(sourcemodel, 'mm');
sourcemodel.inside=source.inside;
sourcemodel.outside=source.outside;

source=ft_convert_units(source, 'mm');
% oldpos=source.pos;

source.pos=sourcemodel.pos;
source.coordsys = 'spm';

X = 45;                  %
Y = 20;                  %
xMargin = 1;               %
yMargin = 1;               %
xSize = X - 2*xMargin;     %
ySize = Y - 2*yMargin;     %

%________________ MCW sampling details _______
Fs = 50; % hard coded... to be retrievied from the blp structure
step2=round(Fs*step/1000);
window2=round(Fs*window/1000); % in points;
%______________________________________________


% % % % % % label1=label(1,:)';
% % % % % % label2=label(2,:)';
% % % % % % label3=label(3,:)';
% % % % % % label4=label(4,:)';
% % % % % % % %____________ EXT NODE
% % % % % % label5=label(5,:)';


%_____________ SETTINGS (THESE WILL BE MANAGED FROM SOMEWHERE ELSE) _______

% % % % % % % % % % % % xdim = source.dim(1);
% % % % % % % % % % % % ydim = source.dim(2);
% % % % % % % % % % % % zdim = source.dim(3);


%__________________________________________________________________________
%_______ all power is mean removed and divided for std once and for all
%_______ here, not needed later
disp('removing mean from power')
source.power=source.power.*1e15;
source.power = (source.power-...
    repmat(squeeze(mean(source.power,2)),1,size(source.power,2)))./...
    (repmat(std(source.power,1,2),1,size(source.power,2)));

%________________________________________________________________________
position = source.pos(source.inside,:);
if ~isempty(extnode)
% label(5,:)=[-35, -88, 4];
label(5,:)=extnode;
end
for inode=1:5
    ref_p=label(inode,:)';
    
    disp('finding seed voxels')
    abs2=sqrt((position(:,1)-ref_p(1)*ones(length(position(:,1)),1)).^2+...
        (position(:,2)-ref_p(2)*ones(length(position(:,2)),1)).^2+...
        (position(:,3)-ref_p(3)*ones(length(position(:,3)),1)).^2);
    
    [x, iref(inode)]=min(abs2);
    
    ref(inode,:) = source.power(iref(inode),:);
    seed_pos(inode,:)=position(iref(inode),:);
end

ntp=size(ref,2);
nwin=fix((ntp-window2)/step2);


disp('estimating seed correlations')
for k=1:nwin
    vect1=[(k-1)*step2+1:(k-1)*step2+window2];
    corr_wins(k,:)=[(k-1)*step2+1 (k-1)*step2+window2];
    time_corr(k)=(1/Fs)*mean(vect1);       % in seconds
    
    tmp=find(vect1 >= 1 & vect1 <= ntp);
    
    r=corrcoef(ref(1,vect1(tmp)),ref(2,vect1(tmp))); c_p1(1,k)=r(1,2);
    r=corrcoef(ref(1,vect1(tmp)),ref(3,vect1(tmp))); c_p1(2,k)=r(1,2);
    r=corrcoef(ref(1,vect1(tmp)),ref(4,vect1(tmp))); c_p1(3,k)=r(1,2);
    r=corrcoef(ref(1,vect1(tmp)),ref(5,vect1(tmp))); c_p1(4,k)=r(1,2);
    
    r=corrcoef(ref(2,vect1(tmp)),ref(1,vect1(tmp))); c_p2(1,k)=r(1,2);
    r=corrcoef(ref(2,vect1(tmp)),ref(3,vect1(tmp))); c_p2(2,k)=r(1,2);
    r=corrcoef(ref(2,vect1(tmp)),ref(4,vect1(tmp))); c_p2(3,k)=r(1,2);
    r=corrcoef(ref(2,vect1(tmp)),ref(5,vect1(tmp))); c_p2(4,k)=r(1,2);
    
    r=corrcoef(ref(3,vect1(tmp)),ref(1,vect1(tmp))); c_p3(1,k)=r(1,2);
    r=corrcoef(ref(3,vect1(tmp)),ref(2,vect1(tmp))); c_p3(2,k)=r(1,2);
    r=corrcoef(ref(3,vect1(tmp)),ref(4,vect1(tmp))); c_p3(3,k)=r(1,2);
    r=corrcoef(ref(3,vect1(tmp)),ref(5,vect1(tmp))); c_p3(4,k)=r(1,2);
    
    r=corrcoef(ref(4,vect1(tmp)),ref(1,vect1(tmp))); c_p4(1,k)=r(1,2);
    r=corrcoef(ref(4,vect1(tmp)),ref(2,vect1(tmp))); c_p4(2,k)=r(1,2);
    r=corrcoef(ref(4,vect1(tmp)),ref(3,vect1(tmp))); c_p4(3,k)=r(1,2);
    r=corrcoef(ref(4,vect1(tmp)),ref(5,vect1(tmp))); c_p4(4,k)=r(1,2);
end
c_sn(1,:)=min(abs(c_p1(1:3,:)));
c_sn(2,:)=min(abs(c_p2(1:3,:)));
c_sn(3,:)=min(abs(c_p3(1:3,:)));
c_sn(4,:)=min(abs(c_p4(1:3,:)));
c_out=[c_p1(4,:) ; c_p2(4,:) ; c_p3(4,:) ; c_p4(4,:)];

disp('finding mcw')
iii=[]; 

junk=[];

for inode=1:1
    t_n=max(abs(c_sn(inode,:)));
    while(t_n>thresh_tn)
        t_e=t_n-0.9;
        while(t_e<(thresh_tn-ac))
            junk=[find(abs(c_sn(inode,:))>t_n & abs(c_out(inode,:))<t_e)];
            iii=horzcat(iii,junk);
            t_e=t_e+0.01;
        end
        t_n=t_n-0.01;
    end
end
iii=unique(iii);

tmpindx=[]; cindx=1; iw=[];
% % % % % if isempty(iii)==0
% % % % %     %     if isempty(find(diff(iii)>5))==0
% % % % %     tmpindx=iii(1); iw(1)=iii(1);
% % % % %     if(size(iii,2)>1)
% % % % %         for ip=2:size(iii,2)
% % % % %             if(iii(ip)>(tmpindx+5))
% % % % %                 cindx=cindx+1;
% % % % %                 iw(cindx)=iii(ip);
% % % % %                 tmpindx=iii(ip);
% % % % %             end
% % % % %         end
% % % % %     end
% % % % %     %     else
% % % % %     %         iw=iii(1);
% % % % %     %     end
% % % % % end


    if isempty(iii)==0
        if isempty(find(diff(iii)>5))==0
            iw=[iii(1) iii(find(diff(iii)>5)+ones(1,length(find(diff(iii)>5))))];

        else
            iw=iii(1);
        end
    end
    
    time_corr_win=(1/Fs)*window2/2;

if(~isempty(iw))
    % if(strcmp(dofig,'yes'))
    h=figure;         set(gcf, 'Units','centimeters', 'Position',[1 1 xSize ySize],'Color',[0.8 0.8 0.8])
    set(h, 'paperposition', [1 1 24 6]);
    plot(time_corr,c_p1(1,:),'b'); hold on;
    plot(time_corr,c_p1(2,:),'r')
    plot(time_corr,c_p1(3,:),'g')
    plot(time_corr,c_p1(4,:),'k')
    % plot(time_corr,c_mean,'k:')
    axis([0 time_corr(1,end) 0 1])
%     plot(time_corr,repmat(0.6,1,size(time_corr,2)),'-')
 line([time_corr(1) time_corr(end)], [0.6 0.6],'color','k')
  line([time_corr(1) time_corr(end)], [0.5 0.5],'color','k')
    for il=1:size(iii,2)
        line([time_corr(iii(il)) time_corr(iii(il))], [0 c_p1(4,iii(il))],'color','y')
    end
%        for il=1:size(iw,2)
%           max_corr= max([c_p1(1,iw(il)) c_p1(2,iw(il)),c_p1(3,iw(il))])
%         line([time_corr(iw(il)) time_corr(iw(il))], [0 1],'color','k')
%     end
    % grid on
    xlabel('time (s)');
    ylabel('corr over windows');
    % legend('rpips','lfef','rfef','out','mean')
    legend('rpips','lfef','rfef','out','Location','NorthEastOutside')
    plot(time_corr(iw),c_p1(1,iw),'*');
    plot(time_corr(iw),c_p1(2,iw),'r*')
    plot(time_corr(iw),c_p1(3,iw),'g*')
    plot(time_corr(iw),c_p1(4,iw),'m*')
        plot(time_corr(iw),c_p1(1,iw),'o');
    plot(time_corr(iw),c_p1(2,iw),'ro')
    plot(time_corr(iw),c_p1(3,iw),'go')
    plot(time_corr(iw),c_p1(4,iw),'mo')
    % if(~isempty(tmpindx))
    % plot(time_corr(tmpindx),c_p1(tmpindx),'ro')
    % end
    set(gca,'Color',[0.5 0.5 0.5]);
    imgname = [outputfile '_seed_correlations.png'];
    hcp_write_figure(imgname, h);
    close(h)
    % end
end

if isempty(iii)==0
    xp=time_corr(iw);
    clear extr extr2
    
    %____ extremes in seconds _________________________________________
    extr=[xp'-window/2000*ones(size(xp))' xp'+window/2000*ones(size(xp))'];
    %____ extremes in points _____________________
    extr2=round(Fs*extr);
    if extr2(1)==0
        extr2(1)=1;
    end
    
    ref1 = ref(1,:);
    %_________________________________________________________________
    
    %_____________________ CONNECTIVITY MAP COMPUTATION ______________
    
    dm1_pre = zeros(size(extr2,1),size(source.power,1));
    z_fisher_dm1 = zeros(size(extr2,1),size(source.power,1));
    
    disp('Connectivity map Computation')
    for w=1:size(extr2,1)
        disp(num2str(w))
                
        
        dm1_pre(w,:) = corr(ref1(1,extr2(w,1):extr2(w,2))',source.power(:,extr2(w,1):extr2(w,2))');
        
        
        z_fisher_dm1(w,:)=0.5.*log((1+dm1_pre(w,:))./(1-dm1_pre(w,:)));
        
%                     if(strcmp(dofig,'yes'))
%         
%                         h=figure;         set(gcf, 'Units','centimeters', 'Position',[5 5 xSize ySize],'Color',[0.8 0.8 0.8])
%                         plot(dm1_pre(w,:),'m') ; hold on
%                         plot(iref,dm1_pre(w,iref),'k*')
%                         plot(iref2,dm1_pre(w,iref2),'*')
%                         plot(iref3,dm1_pre(w,iref3),'r*')
%                         plot(iref4,dm1_pre(w,iref4),'g*')
%                         plot(iref5,dm1_pre(w,iref5),'y*')
%                         xlabel('voxel');
%                         ylabel('correlation with seed');
%                         legend('seed', 'ref1','ref2','ref3','ref4','out')
%                         set(gca,'Color',[0.5 0.5 0.5]);
%                         line([1 size(dm1_pre,2)], [0.6 0.6])
%                         imgname = [outputfile '_voxel_correlations_' num2str(w) '.png'];
%                         hcp_write_figure(imgname, h);
%                         close(h)
%                     end
    end
    
%                 for izf=1:size(z_fisher_dm1,1)
%                 [junk seed_indx]=max(z_fisher_dm1(izf,:));
%                 z_fisher_dm1(izf,seed_indx)=0;
%                 z_fisher_dm1(izf,seed_indx)=max(z_fisher_dm1(izf,:));
%             end
    
    connect_mcw=source;
    connect_mcw=rmfield(connect_mcw,'power')
    connect_mcw=rmfield(connect_mcw,'time_power')
    if(isfield(connect_mcw,'avg')) connect_mcw=rmfield(connect_mcw,'avg'); end
    
    connect_mcw.time=[1:size(dm1_pre,1)];
    connect_mcw.mcw=dm1_pre;
    connect_mcw.mcwextr=extr2;
    
    connect_mcw.seeds_indx=iref;
    
    disp('Evaluating Z-Score MAP')
    
%         mu_base=mean(mean(z_fisher_dm1,2),1);
%         mu_s=mean(z_fisher_dm1,1);
%         sigma_s=std(z_fisher_dm1,1,1);
%         z_s=(mu_s-mu_base).*sqrt(size(z_fisher_dm1,1))./sigma_s;
%         
        mu_base=mean(mean(dm1_pre,2),1);
        mu_s=mean(dm1_pre,1);
        sigma_s=std(dm1_pre,1,1);
        z_s=(mu_s-mu_base).*sqrt(size(dm1_pre,1))./sigma_s;
        
        for iz=1:5
        [junk seed_indx(iz)]=max(z_s);
        z_s(seed_indx(iz))=0;
        end
        z_s(seed_indx)=max(z_s);
        connect_mcw.z_s=z_s;
       

    
    if(strcmp(dofig,'yes'))
        h=figure;         set(gcf, 'Units','centimeters', 'Position',[5 5 xSize ySize],'Color',[0.8 0.8 0.8])
        
        plot(z_s,'m') ; hold on
        plot(iref(1,1),z_s(1,iref(1,1)),'k*')
        plot(iref(1,2),z_s(1,iref(1,2)),'*')
        plot(iref(1,3),z_s(1,iref(1,3)),'r*')
        plot(iref(1,4),z_s(1,iref(1,4)),'g*')
        plot(iref(1,5),z_s(1,iref(1,5)),'y*')
        xlabel('voxel');
        ylabel('z-score');
        legend('seed', 'ref1','ref2','ref3','ref4','out')
        set(gca,'Color',[0.5 0.5 0.5]);
        %             line([1 size(dm1_pre,2)], [0.6 0.6])
        imgname = [outputfile '_voxel_z_score.png'];
        hcp_write_figure(imgname, h);
        close(h)
    end
    
    if(strcmp(dofig,'yes'))
        
        intcfg=[];
        intcfg.parameter = 'avg.pow';
        tmp = connect_mcw;
        %         tmp.pos=oldpos;
        tmp.avg.pow=zeros(size(source.pos,1),1);
        tmp.avg.pow(source.inside) = z_s;
        s_interp = ft_sourceinterpolate(intcfg, tmp, mri);
        [junk thr_ind junk2 thr_ind_dep]= hcp_mcw_fdr(z_s,0.05,size(connect_mcw.mcw,1)); % [z_thr thr3 z_thr_dep thr3_dep]
        
        
%         indxmask=find(s_interp.avg.pow>thr_ind_dep);
        indxmask=find(s_interp.avg.pow>0);

        s_interp.avg.mask=zeros(size(s_interp.avg.pow));
        s_interp.avg.mask(indxmask)=1;
        %         s_interp.avg.mask=zeros(size(s_interp.avg.pow));
        %         maxabs=max(max(max(s_interp.avg.pow)));
        %         indxmask=find(s_interp.avg.pow>1.5);
        %         s_interp.avg.mask(indxmask)=1;
        %
        intcfg=[];
        intcfg.parameter = 'avg.pow';
        %         intcfg.downsample = 1;
        
        tmp = connect_mcw;
        
        for k = 1:numel(connect_mcw.time)
            tmp.avg.pow=zeros(size(source.pos,1),1);
            tmp.avg.pow(connect_mcw.inside) = connect_mcw.mcw(k,:);
            s_interp2(k) = ft_sourceinterpolate(intcfg, tmp, mri);
            s_interp2(k).avg.mask=zeros(size(s_interp2(k).avg.pow));
            maxabs=max(max(max(s_interp2(k).avg.pow)));
            indxmask=find(s_interp2(k).avg.pow>0.4*maxabs);
            s_interp2(k).avg.mask(indxmask)=1;
        end
        
        
        tmp = connect_mcw;
        tmp.time=1;
        
        tmp.avg.pow=zeros(size(tmp.pos,1),1);
        tmp.avg.pow(:,:)=0;
        
        tmp.avg.pow(tmp.inside(iref(1,1))) = -10;
        tmp.avg.pow(tmp.inside(iref(1,2))) = 10;
        tmp.avg.pow(tmp.inside(iref(1,3))) = 10;
        tmp.avg.pow(tmp.inside(iref(1,4))) = 10;
        tmp.avg.pow(tmp.inside(iref(1,5))) = 10;
        
        s_interp_seed = ft_sourceinterpolate(intcfg, tmp, mri);
        s_interp_seed.avg.mask=zeros(size(s_interp_seed.avg.pow));
        maxabs=max(max(max(s_interp_seed.avg.pow)));
        indxmask=find(abs(s_interp_seed.avg.pow)>0.4*maxabs);
        s_interp_seed.avg.mask(indxmask)=1;

        %         netstring=net_seed{seed_indx};
        %         tmpindx=strfind(netstring, '_')
        %         netstring(tmpindx)='-';
        %         imgname = ['mcw_netseed_' net_seed{seed_indx}];
        
        
        
        cfg = [];
        cfg.method = 'slice';
        cfg.funcolormap='jet';
        cfg.nslices=36;
        %   cfg.method = 'surface';
        %   cfg.anaparameter = 'anatomy';
        %   cfg.coordsys = 'mni';
        %   cfg.anaparameter = 'mri';
        cfg.funparameter = 'avg.pow';
        %   cfg.method = 'ortho';
        cfg.projmethod = 'nearest';
        cfg.interactive = 'no';
% % % % %         cfg.funcolorlim=[-10 10];        
                        cfg.maskparameter='avg.mask';
% % % % %         ft_sourceplot(cfg, s_interp_seed);
% % % % %         
% % % % %         ax1 = gca;
        
        cfg.funcolorlim=[0 10];
        % % % % %             cfg.maskparameter='avg.mask';
        
        imgname = [outputfile '_z_score'];
        ft_sourceplot(cfg, s_interp);
%         hcp_write_figure(imgname, gcf, 'format', 'fig')
        disp('step1')
% % % % %         ax2 = gca;
        h3=gcf;
% % % % %         h3 = figure; %create new figure
% % % % %         set(gca, 'XTickLabel',[], 'YTickLabel',[], ...
% % % % %             'Units','normalized', 'Position',[0 0 1 1])
% % % % %         
% % % % %         %# figure size on screen (50% scaled, but same aspect ratio)
% % % % %         set(gcf, 'Units','centimeters', 'Position',[5 5 xSize ySize])
% % % % %         
% % % % %         %# figure size printed on paper
% % % % %         set(gcf, 'visible', 'on')
% % % % %         set(gcf, 'PaperUnits','centimeters')
% % % % %         set(gcf, 'PaperSize',[X Y])
% % % % %         set(gcf, 'PaperPosition',[xMargin yMargin xSize ySize])
% % % % %         set(gcf, 'PaperOrientation','portrait')
        
        disp('step2')

% % % % %         s1 = subplot(1,2,1); %create and get handle to the subplot axes
% % % % %         fig1 = get(ax1,'children'); %get handle to all the children in the figure
% % % % %         copyobj(fig1,s1); %copy children to new parent axes i.e. the subplot axes
% % % % %         caxis([-10,10])
% % % % %         colorbar
       
        disp('step3')

% % % % %         s2 = subplot(1,2,2);
% % % % %         fig2 = get(ax2,'children');
% % % % %         copyobj(fig2,s2);
% % % % %         colorbar
% % % % %         caxis(cfg.funcolorlim)
        
        disp('step4')

        imgname = [outputfile '_z_score.png'];

        hcp_write_figure(imgname, h3);
        hcp_write_figure(imgname, h3, 'format' ,'fig');

        
        disp('step5')

        close(h3)
                
        disp('step6')

%         close all
        
        for i = 1:size(s_interp2, 2)
            imgname = [outputfile '_' num2str(i) '.png'];
            %                 maxVal=max(s_interp2(i).avg.pow(:));
            %                 minVal=min(s_interp2(i).avg.pow(:));
            %                 maxAbs=max([maxVal;minVal])
            %                 limSign=unique(sign([minVal maxVal]));
            %                 if (length(limSign)==1)&limSign>0,
            %                     funLim=maxAbs*[0 1];
            %                 elseif (length(limSign)==1)&limSign<0,
            %                     funLim=maxAbs*[-1 0];
            %                 else
            %                     funLim=maxAbs*[-1 1];
            %                 end
            funLim=[0 1];
            cfg.funcolorlim=funLim;
            
            ft_sourceplot(cfg, s_interp2(i));
            h1=gcf;
%             set(h1, 'paperposition', [1 1 10 7]);
            hcp_write_figure(imgname, h1);
            close(h1)
        end
    end
    
    
else
    connect_mcw=[];
end %_over k ref points
close all
clear time_corr c_p1 c_p2 c_p3 c_p4
end



