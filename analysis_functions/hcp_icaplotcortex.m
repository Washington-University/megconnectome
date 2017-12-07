function []=hcp_icaplotcortex(map_val, subject, options)


imgname=ft_getopt(options, 'outputfile');
rm_wall=ft_getopt(options, 'mask', 'yes');
color_extr=ft_getopt(options, 'color_extr');
c_map= ft_getopt(options, 'color_map', 'jet');
savefig= ft_getopt(options, 'savefig', 'yes');
browsefig= ft_getopt(options, 'browsefig', 'no');
cortex = ft_getopt(options, 'cortex', 'midt');

X = 45; Y = 20;
xMargin = 1;  yMargin = 1;
xSize = X - 2*xMargin;  ySize = Y - 2*yMargin;

if isempty(color_extr)
    color_extr=[nanmin(nanmin(map_val)) nanmax(nanmax(map_val))];
end

if strcmp(rm_wall,'yes')
    atlas_yeo=ft_read_cifti('Yeo2011_17Networks.LR.min50sqmm.4k_fs_LR.dlabel.nii');
    
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
    
    atlasindex=zeros( numel(atlas_yeo.x1),1);
    
    for i=0:17
        
        indextmp=find(parcel_net==i)-1 ;
        
        index_parcels=[];
        
        for ipar=1:numel(indextmp)
            junk=find(atlas_yeo.x1==indextmp(ipar));
            index_parcels=[index_parcels' junk']';
            
            atlasindex(index_parcels,1)=repmat(i,numel(index_parcels),1);
            
            
        end
        
    end
    
    mask1_wall=find(strcmp(NetLabels,'Medial_Wall'))-1;
    mask_wall=find(atlasindex==mask1_wall);
    map_val(1,mask_wall)=-(max(abs(map_val))*10);
    
end

if(strcmp(cortex,'midt'))
head=ft_read_headshape(['Q1-Q6_R440.L.midthickness.4k_fs_LR.surf.gii']);
elseif(strcmp(cortex,'infl'))
head=ft_read_headshape(['Q1-Q6_R440.L.inflated.4k_fs_LR.surf.gii']);
else
head=ft_read_headshape([subject '.L.midthickness.4k_fs_LR.surf.gii']);
end

tmp2=[];
tmp2.tri=head.tri;
tmp2.pos=head.pnt;
tmp2.avg.pow=map_val(1:size(tmp2.pos,1))';

if(strcmp(cortex,'midt'))
head2=ft_read_headshape(['Q1-Q6_R440.R.midthickness.4k_fs_LR.surf.gii']);
elseif(strcmp(cortex,'infl'))
head2=ft_read_headshape(['Q1-Q6_R440.R.inflated.4k_fs_LR.surf.gii']);
else
head2=ft_read_headshape([subject '.R.midthickness.4k_fs_LR.surf.gii']);
end

tmp3=[];
tmp3.tri=head2.tri;
tmp3.pos=head2.pnt;
tmp3.avg.pow=map_val((size(tmp2.pos,1))+1:end)';

figure
colorbar_extr=color_extr;
h1=gcf;
set(gca, 'XTickLabel',[], 'YTickLabel',[], ...
    'Units','normalized', 'Position',[0 0 1 1])

%# figure size on screen (50% scaled, but same aspect ratio)
set(gcf, 'Units','centimeters', 'Position',[2 2 xSize+5 ySize+5])

%# figure size printed on paper
set(gcf, 'visible', 'on')
set(gcf, 'PaperUnits','centimeters')
set(gcf, 'PaperSize',[X Y])
set(gcf, 'PaperPosition',[xMargin yMargin xSize ySize])
set(gcf, 'PaperOrientation','portrait')

subplot(2,4,2)
ft_plot_mesh(tmp2,'vertexcolor',tmp2.avg.pow)
caxis(colorbar_extr) ; colormap(c_map)
view(0,90)

subplot(2,4,6)
ft_plot_mesh(tmp2,'vertexcolor',tmp2.avg.pow)
caxis(colorbar_extr) ; colormap(c_map)
view(0,0)
colorbar('location','NorthOutside')

subplot(2,4,3)
ft_plot_mesh(tmp3,'vertexcolor',tmp3.avg.pow)
caxis(colorbar_extr) ; colormap(c_map)
view(0,90)

subplot(2,4,7)
ft_plot_mesh(tmp3,'vertexcolor',tmp3.avg.pow)
caxis(colorbar_extr) ; colormap(c_map)
view(0,0)
colorbar('location','NorthOutside')

subplot(2,4,1)
ft_plot_mesh(tmp2,'vertexcolor',tmp2.avg.pow)
caxis(colorbar_extr) ; colormap(c_map)
view(90,0)

subplot(2,4,5)
ft_plot_mesh(tmp2,'vertexcolor',tmp2.avg.pow)
caxis(colorbar_extr) ; colormap(c_map)
view(-90,0)

subplot(2,4,4)
ft_plot_mesh(tmp3,'vertexcolor',tmp3.avg.pow)
caxis(colorbar_extr) ; colormap(c_map)
view(-90,0)

subplot(2,4,8)
ft_plot_mesh(tmp3,'vertexcolor',tmp3.avg.pow, 'facealpha', mask_wall)
caxis(colorbar_extr) ; colormap(c_map)
view(90,0)

if strcmp(savefig,'yes')
hcp_write_figure(imgname, h1)
end
if strcmp(browsefig,'no')
close(h1)
end
end
