function []=hcp_icaplotconnectome(connectome, sourcemodel2d, options)

outputfile=ft_getopt(options, 'outputfile');
sortindx=ft_getopt(options, 'sorting','yes');
rm_wall=ft_getopt(options, 'mask', 'yes');
rm_dist=ft_getopt(options, 'mask_edist','yes');
edist_radius=ft_getopt(options, 'edist_radius', 3.5);
plot_parc=ft_getopt(options, 'plot_parc','yes');
color_extr=ft_getopt(options, 'color_extr');
fig_type=ft_getopt(options, 'fig_type','png');
yeo_mode=ft_getopt(options, 'yeo_mode',[]);

if strcmp(sortindx,'yes')
    
    % reading Yeo parcellation for vertex ordering
    atlas_yeo=ft_read_cifti('Yeo2011_17Networks.LR.min50sqmm.4k_fs_LR.dlabel.nii');
    
    if ~isempty(yeo_mode) & strcmp(yeo_mode,'compact')
    % ordering of the Yeo 17Networks as will be displayed (DMN, MOT, FP, DAN ...)
    ordering_n=[17 18 16 4 5 15 9 12 13 14 7 6 8 10 11 2 3 1];
    
    NetLabels={ 'Medial_Wall'
        '17Networks_1_VIS_1'
        '17Networks_2_VIS_2'
        '17Networks_3_MOT_1'
        '17Networks_4_MOT_2'
        '17Networks_5_DAN_2'
        '17Networks_6_DAN_1'
        '17Networks_7_VAN_1'
        '17Networks_8_FP_1'
        '17Networks_9_LIM_1'
        '17Networks_10_LIM_2'
        '17Networks_11_FP_2'
        '17Networks_12_FP_3'
        '17Networks_13_FP_4'
        '17Networks_14_MOT_3'
        '17Networks_15_DMN_3'
        '17Networks_16_DMN_1'
        '17Networks_17_DMN_2'};
    
    % assign to each parcel (patch) the corresponding original 17Network parcel
    parcel_net= [0 2 16 16 16 16 0 0 6 6 6 6 12 12 12 12 12 12 10 5 13 13 13 13 13 13 13 ...
        17 17 17 17 17 17 15 15 15 1 4 9 11 11 3 8 8 8 8 8 8 7 7 7 7 14 14 ...
        2 16 16 16 16 16 16 0 0 6 6 6 6 12 12 12 12 12 12 10 6 5 13 13 13 13 13 ...
        17 17 17 17 17 15 15 15 1 4 9 11 11 11 3 8 8 8 8 8 8 7 7 7 7 14 14];
    
    atlasindex=zeros(numel(atlas_yeo.x1),1);
    
    for i=0:17
        
        indextmp=find(parcel_net==i)-1 ;
        
        index_parcels=[];
        
        for ipar=1:numel(indextmp)
            junk=find(atlas_yeo.x1==indextmp(ipar));
            index_parcels=[index_parcels' junk']';
            
            atlasindex(index_parcels,1)=repmat(i,numel(index_parcels),1);
            
            
        end
        
    end
    else
        ordering_n=[1:numel(atlas_yeo.x1label)+1];
        NetLabels{1}='Medial_Wall';
        for ij=1:numel(atlas_yeo.x1label) ; NetLabels{ij+1}=atlas_yeo.x1label{ij}; end
        atlasindex=atlas_yeo.x1;
    end
    
    for it=1:numel(NetLabels) ; atlaslabel{it}=NetLabels{ordering_n(it)} ; end
    
    indexvoxels=[];
    jn=0;
    for in=1:numel(atlaslabel)
        jn=jn+1;
        networkl{jn}=atlaslabel{in};
        nindx=ordering_n(in)-1;
        if ~isempty(nindx)
            indexvoxels=[indexvoxels ; find(atlasindex==nindx)];
            indxvoxels_net{jn}=find(atlasindex==nindx);
        end
        label_net_n{jn}=atlaslabel{in};
        if (jn==1)
            net_max(jn,:)=[1 numel(indexvoxels)];
        else
            net_max(jn,:)=[(net_max(jn-1,2)+1) numel(indexvoxels)];
        end
    end
    
else
    indexvoxels=[1:size(connectome,1)];
end

if strcmp(rm_dist,'yes')
    eudist=squareform(pdist(sourcemodel2d.pos));
    mask_c=find(eudist<edist_radius);
    connectome(mask_c)=NaN;
end

if strcmp(rm_wall,'yes')
    mask1_wall=find(strcmp(NetLabels,'Medial_Wall'))-1;
    mask_wall=[find( atlasindex==mask1_wall)];
    connectome(mask_wall,:)=NaN;
    connectome(:,mask_wall)=NaN;
end

if isempty(color_extr)
    color_extr=[nanmin(nanmin(connectome)) nanmax(nanmax(connectome))];
end
indexvoxels_p=setdiff(indexvoxels,mask_wall,'stable');
connectome_ordered=connectome(indexvoxels_p,indexvoxels_p);
figure
h=imagesc(connectome_ordered);
set(h,'alphadata',~isnan(connectome_ordered))
caxis(color_extr)
colorbar
imgname = [outputfile '.' fig_type];
hcp_write_figure(imgname, gcf)

close all

connectome_ordered=connectome(indexvoxels,indexvoxels);
if strcmp(plot_parc,'yes')
    connect_parc=zeros(numel(networkl));
    for ip=1:numel(networkl)
        for jp=1:numel(networkl)
            connect_parc(ip,jp)=nanmean(nanmean(connectome_ordered(net_max(ip,1):net_max(ip,2),net_max(jp,1):net_max(jp,2))));
        end
    end
    
    saveparc='no';
    if(strcmp(saveparc,'yes'))
        fname = [outputfile '_parc'];
        save(fname,'connect_parc')
    end
    
    connect_parc=connect_parc(1:end-1,1:end-1);
    
    figure
    h=imagesc(connect_parc);
    set(h,'alphadata',~isnan(connect_parc))
    caxis(color_extr)
    colorbar
    imgname = [outputfile '_parc' '.' fig_type];
    if(numel(networkl)<20)
    set(gca,'ytick',[1:numel(networkl)])
    set(gca,'yticklabel',networkl)
    end
    hcp_write_figure(imgname, gcf)
    
    close all
    
end
end