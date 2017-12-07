function [connect_mcw] = hcp_icamcw_estimate(source, sourcemodel, options)

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
net_seeds=ft_getopt(options, 'net_seeds');
step = ft_getopt(options, 'mcw_step', 200);
window     = ft_getopt(options, 'mcw_window', 10000);
dofig = ft_getopt(options, 'dofig', 'yes');
ac=ft_getopt(options, 'ac', 0.1);
thresh_tn=ft_getopt(options, 'thresh', 0.6);
extnode=ft_getopt(options, 'extnode');
% outputfile=[outputfile '_test']

mni_pos=net_seeds.pos;

X = 45; Y = 20;
xMargin = 1;  yMargin = 1;               
xSize = X - 2*xMargin;  ySize = Y - 2*yMargin; 

%________________ MCW sampling details _______
Fs = 50; % hard coded... to be retrievied from the blp structure
step2=round(Fs*step/1000);
window2=round(Fs*window/1000); % in points;
%______________________________________________

%__________________________________________________________________________
%_______ all power is mean removed and divided for std once and for all
%_______ here, not needed later
disp('removing mean from power')
% source.power=source.power.*1e15;
source.power = (source.power-...
    repmat(squeeze(mean(source.power,2)),1,size(source.power,2)))./...
    (repmat(std(source.power,1,2),1,size(source.power,2)));

%________________________________________________________________________
position = sourcemodel.pnt(source.inside,:);
if ~isempty(extnode)
    mni_pos(5,:)=extnode;
end

inode=5
ref_p=mni_pos(inode,:)';
disp('finding seed out')
abs2=sqrt((position(:,1)-ref_p(1)*ones(length(position(:,1)),1)).^2+...
    (position(:,2)-ref_p(2)*ones(length(position(:,2)),1)).^2+...
    (position(:,3)-ref_p(3)*ones(length(position(:,3)),1)).^2);

[x, iref(inode)]=min(abs2);

iref(1:4)=net_seeds.cortex_index
ref=source.power(iref,:);
seed_pos=position(iref,:);

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
if isempty(iii)==0
    %     if isempty(find(diff(iii)>5))==0
    tmpindx=iii(1); iw(1)=iii(1);
    if(size(iii,2)>1)
        for ip=2:size(iii,2)
            if(iii(ip)>(tmpindx+9))
                cindx=cindx+1;
                iw(cindx)=iii(ip);
                tmpindx=iii(ip);
            end
        end
    end
    %     else
    %         iw=iii(1);
    %     end
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
    legend(net_seeds.seeds_string{2:5},'Location','NorthEastOutside')
    
    %     legend('rpips','lfef','rfef','out','Location','NorthEastOutside')
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

load('parcellations_VGD11b_4k.mat')
mask1=find(strcmp(atlas.parcellation2label,'L_MEDIAL.WALL'));
mask2=find(strcmp(atlas.parcellation2label,'R_MEDIAL.WALL'));
mask_wall=[find(atlas.parcellation2==mask1); find(atlas.parcellation2==mask2)];

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
    
    mcw_corrmap = zeros(size(extr2,1),size(source.power,1));
    
    disp('Connectivity map Computation')
    for w=1:size(extr2,1)
        disp(num2str(w))
        
        mcw_corrmap(w,:) = corr(ref1(1,extr2(w,1):extr2(w,2))',source.power(:,extr2(w,1):extr2(w,2))');
    end
    
    connect_mcw=source;
    connect_mcw=rmfield(connect_mcw,'power')
    connect_mcw=rmfield(connect_mcw,'time_power')
    if(isfield(connect_mcw,'avg')) connect_mcw=rmfield(connect_mcw,'avg'); end
    
    connect_mcw.time=[1:size(mcw_corrmap,1)];
    connect_mcw.mcw=mcw_corrmap;
    connect_mcw.mcwextr=extr2;
    connect_mcw.seeds_indx=iref;
    
    disp('Evaluating Z-Score MAP')
    
    tmp_mcwcorr=mcw_corrmap;
    tmp_mcwcorr(:,mask_wall)=NaN;
    mu_base=nanmean(nanmean(tmp_mcwcorr,2),1);
    mu_s=mean(mcw_corrmap,1);
    sigma_s=std(mcw_corrmap,1,1);
    z_s=(mu_s-mu_base).*sqrt(size(mcw_corrmap,1))./sigma_s;
    
    for iz=1:5
        [junk seed_indx(iz)]=max(z_s);
        z_s(seed_indx(iz))=0;
    end
    z_s(seed_indx)=max(z_s);
    connect_mcw.z_s=z_s;
    
    if(strcmp(dofig,'yes'))
        h=figure;
        set(gcf, 'Units','centimeters', 'Position',[5 5 xSize ySize],'Color',[0.8 0.8 0.8])
        
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
        imgname = [outputfile '_voxel_z_score.png'];
        hcp_write_figure(imgname, h);
        close(h)
    end
    
    if(strcmp(dofig,'yes'))
        
        imgname = [outputfile '_z_score_view.png'];
        color_extr= [0 5];
        options_pl={'outputfile',imgname, 'mask','yes','color_extr',color_extr,'color_map','jet'};
        
        hcp_icaplotcortex(z_s, subject, options_pl)
                
    end
    
    
else
    connect_mcw=[];
end 
close all
end



