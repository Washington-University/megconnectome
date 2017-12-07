function []=hcp_icaplotcortex(map_val, subject, options)


imgname=ft_getopt(options, 'outputfile');
rm_wall=ft_getopt(options, 'mask', 'yes');
color_extr=ft_getopt(options, 'color_extr');
c_map= ft_getopt(options, 'color_map', 'jet');

X = 45; Y = 20;
xMargin = 1;  yMargin = 1;
xSize = X - 2*xMargin;  ySize = Y - 2*yMargin;

if isempty(color_extr)
    color_extr=[nanmin(nanmin(map_val)) nanmax(nanmax(map_val))];
end

load('parcellations_VGD11b_4k.mat')
mask1_wall=find(strcmp(atlas.parcellation2label,'L_MEDIAL.WALL'));
mask2_wall=find(strcmp(atlas.parcellation2label,'R_MEDIAL.WALL'));
mask_wall=[find(atlas.parcellation2==mask1_wall); find(atlas.parcellation2==mask2_wall)];

map_val(1,mask_wall)=-(max(abs(map_val))*10);

head=ft_read_headshape([subject '.L.midthickness.4k_fs_LR.surf.gii']);

tmp2=[];
tmp2.tri=head.tri;
tmp2.pos=head.pnt;
tmp2.avg.pow=map_val(1:size(tmp2.pos,1))';

head2=ft_read_headshape([subject '.R.midthickness.4k_fs_LR.surf.gii']);

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
ft_plot_mesh(tmp3,'vertexcolor',tmp3.avg.pow)
caxis(colorbar_extr) ; colormap(c_map)
view(90,0)

hcp_write_figure(imgname, h1)
close(h1)
end
