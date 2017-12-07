function []=hcp_icaplotconnectome(connectome, sourcemodel2d, options)

outputfile=ft_getopt(options, 'outputfile');
sortindx=ft_getopt(options, 'sorting','yes');
atlas_label=ft_getopt(options, 'parcel_type','RSN');
rm_wall=ft_getopt(options, 'mask', 'yes');
rm_dist=ft_getopt(options, 'mask_edist','yes');
edist_radius=ft_getopt(options, 'edist_radius', 3.5);
plot_parc=ft_getopt(options, 'plot_parc','yes');
color_extr=ft_getopt(options, 'color_extr');
fig_type=ft_getopt(options, 'fig_type','png');

if strcmp(sortindx,'yes')
    if strcmp(atlas_label,'RSN')
        atlas_rsn_l=ft_read_atlas('RSN-networks.L.4k_fs_LR.label.gii');
        atlas_rsn_r=ft_read_atlas('RSN-networks.R.4k_fs_LR.label.gii');
        ordering_n=[4 7 8 15 25 5 17 16 6 21 12 14 3 10 20 11 18 26 13 22 23 27 24 9 19 2 1];
        atlaslabel_rtemp=atlas_rsn_r.parcellation4label; atlasindex_r=atlas_rsn_r.parcellation4;
        for it=1:numel(atlaslabel_rtemp) ; atlaslabel_r{it}=atlaslabel_rtemp{ordering_n(it)} ; end 
        atlaslabel_ltemp=atlas_rsn_l.parcellation4label; atlasindex_l=atlas_rsn_l.parcellation4;
        for it=1:numel(atlaslabel_ltemp) ; atlaslabel_l{it}=atlaslabel_ltemp{ordering_n(it)} ; end 

    elseif strcmp(atlas_label,'VGD11b')
        atlas_rsn_l=ft_read_atlas('parcellations_VGD11b.L.4k_fs_LR.label.gii');
        atlas_rsn_r=ft_read_atlas('parcellations_VGD11b.R.4k_fs_LR.label.gii');
        atlaslabel_r=atlas_rsn_r.parcellation2label; atlasindex_r=atlas_rsn_r.parcellation2;
        atlaslabel_l=atlas_rsn_l.parcellation2label; atlasindex_l=atlas_rsn_l.parcellation2;
        ordering_n=[1:numel(atlaslabel_r)];
    end
    
    indexvoxels=[];
    jn=0;
    for in=1:numel(atlaslabel_r)
        jn=jn+1;
        networkl{jn}=atlaslabel_r{in};
%         nindxl=find(strcmp(atlaslabel_l,atlaslabel_r{in}));
        nindxl=ordering_n(in);         nindxr=ordering_n(in); 
        if ~isempty(nindxl)
            indexvoxels=[indexvoxels ; find(atlasindex_l==nindxl)];
            indxvoxels_net{jn}=find(atlasindex_l==nindxl);
        end
        indexvoxels=[indexvoxels ; (find(atlasindex_r==nindxr)+numel(atlasindex_l))];
        indxvoxels_net{jn}=[indxvoxels_net{jn} ; (find(atlasindex_r==nindxr)+numel(atlasindex_l))];
        label_net_n{jn}=atlaslabel_r{in};
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
    load('parcellations_VGD11b_4k.mat')
    mask1_wall=find(strcmp(atlas.parcellation2label,'L_MEDIAL.WALL'));
    mask2_wall=find(strcmp(atlas.parcellation2label,'R_MEDIAL.WALL'));
    mask_wall=[find(atlas.parcellation2==mask1_wall); find(atlas.parcellation2==mask2_wall)];
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
    figure
    h=imagesc(connect_parc);
    set(h,'alphadata',~isnan(connect_parc))
    caxis(color_extr)
    colorbar
    imgname = [outputfile '_parc' '.' fig_type];
    set(gca,'ytick',[1:numel(networkl)])
    set(gca,'yticklabel',networkl)
    hcp_write_figure(imgname, gcf)
    
    close all
    
end
end