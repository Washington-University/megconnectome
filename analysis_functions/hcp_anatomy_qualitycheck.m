function hcp_anatomy_qualitycheck(filename, fiducials, landmarks, transform, mri, headshape, headshapemri, headmodel, sourcemodel2d, sourcemodel3d)

% HCP_ANATOMY_QUALITYCHECK saves a bunch of figures from the output of
% the anatomy pipeline, which need to be inspected visually in order to
% check the quality of the output
%
% Use as
%   hcp_anatomy_qualitycheck(filename)
%
% For example
%   hcp_anatomy_qualitycheck('CP10018_anatomy.mat')
%
% See also HCP_ANATOMY_PIPELINE, HCP_ENSURE_UNITS, HCP_ENSURE_COORDSYS

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

close all

% FIXME this should be assigned in hcp_anatomy_pipeline
fiducials.coordsys = 'vox';
landmarks.coordsys = 'vox';
headshape.coordsys = 'bti';
headshapemri.coordsys = 'bti';

% ensure that they all are represented in a consistent coordinate system
mri           = hcp_ensure_coordsys(mri,           transform, 'bti');
fiducials     = hcp_ensure_coordsys(fiducials,     transform, 'bti');
headmodel     = hcp_ensure_coordsys(headmodel,     transform, 'bti');
sourcemodel2d = hcp_ensure_coordsys(sourcemodel2d, transform, 'bti');
sourcemodel3d = hcp_ensure_coordsys(sourcemodel3d, transform, 'bti');
headshape     = hcp_ensure_coordsys(headshape,     transform, 'bti');
headshapemri  = hcp_ensure_coordsys(headshapemri,  transform, 'bti');

% ensure that they all have consistent units
mri           = hcp_ensure_units(mri,           'mm');
%fiducials     = hcp_ensure_units(fiducials,     'mm');
headmodel     = hcp_ensure_units(headmodel,     'mm');
sourcemodel2d = hcp_ensure_units(sourcemodel2d, 'mm');
sourcemodel3d = hcp_ensure_units(sourcemodel3d, 'mm');
headshape     = hcp_ensure_units(headshape,     'mm');
headshapemri  = hcp_ensure_units(headshapemri,  'mm');

% define some helper functions
viewtop    = @() view(  0,  90);
viewbottom = @() view(180, -90);
viewleft   = @() view(180,   0);
viewright  = @() view(  0,   0);
viewfront  = @() view( 90,   0); % FIXME this might also be back

figure;
subplot(2,2,1); hold on; viewbottom(); ft_plot_vol(headmodel); ft_plot_fiducials(fiducials);
subplot(2,2,2); hold on; viewtop();    ft_plot_vol(headmodel); ft_plot_fiducials(fiducials);
subplot(2,2,3); hold on; viewleft();   ft_plot_vol(headmodel); ft_plot_fiducials(fiducials);
subplot(2,2,4); hold on; viewright();  ft_plot_vol(headmodel); ft_plot_fiducials(fiducials);
axis on; grid on;
hcp_write_figure([filename,'_headmodel.png']);

figure;
subplot(2,2,1); hold on; viewbottom(); ft_plot_mesh(sourcemodel3d.pos(sourcemodel3d.inside,:)); ft_plot_fiducials(fiducials);
subplot(2,2,2); hold on; viewtop();    ft_plot_mesh(sourcemodel3d.pos(sourcemodel3d.inside,:)); ft_plot_fiducials(fiducials);
subplot(2,2,3); hold on; viewfront();  ft_plot_mesh(sourcemodel3d.pos(sourcemodel3d.inside,:)); ft_plot_fiducials(fiducials);
subplot(2,2,4); hold on; viewright();  ft_plot_mesh(sourcemodel3d.pos(sourcemodel3d.inside,:)); ft_plot_fiducials(fiducials);
axis on; grid on;
hcp_write_figure([filename,'_sourcemodel_3d.png']);

figure;
options = {'transform',mri.transform,'intersectmesh',{sourcemodel2d headmodel.bnd}};
subplot(2,2,1); hold on; ft_plot_slice(mri.anatomy, 'location', [0  0 60], 'orientation', [0 0 1], options{:}); viewtop();
subplot(2,2,2); hold on; ft_plot_slice(mri.anatomy, 'location', [0  0 20], 'orientation', [0 0 1], options{:}); viewtop();
subplot(2,2,3); hold on; ft_plot_slice(mri.anatomy, 'location', [0 20  0], 'orientation', [1 0 0], options{:}); viewfront();
subplot(2,2,4); hold on; ft_plot_slice(mri.anatomy, 'location', [0 20  0], 'orientation', [0 1 0], options{:}); viewright();
set(gcf, 'Renderer', 'zbuffer')
hcp_write_figure([filename,'_sourcemodel_2d.png']);

options = {'transform', mri.transform,'nslice',16,'intersectmesh',{sourcemodel2d headmodel.bnd},...
           'intersectlinewidth',1,'slicesize',[300 300]};

figure;
ft_plot_montage(mri.anatomy, 'location', [0 0 0], 'orientation', [0 0 1], 'slicerange', [-20 120], options{:});
set(gcf, 'Renderer', 'zbuffer'); 
hcp_write_figure([filename,'_slice1.png'], 'resolution', 500);

figure;
ft_plot_montage(mri.anatomy, 'location', [0 0 0], 'orientation', [0 1 0], 'slicerange', [-60 60], options{:});
set(gcf, 'Renderer', 'zbuffer');
hcp_write_figure([filename,'_slice2.png'], 'resolution', 500);

figure;
ft_plot_montage(mri.anatomy, 'location', [0 0 0], 'orientation', [1 0 0], 'slicerange', [-70 110], options{:});
set(gcf, 'Renderer', 'zbuffer');
hcp_write_figure([filename,'_slice3.png'], 'resolution', 500);

v = headshapemri.bnd.pnt;
f = headshapemri.bnd.tri;
[f,v]=reducepatch(f,v, 0.2);
headshapemri.bnd.pnt = v;
headshapemri.bnd.tri = f;
figure
subplot(2,2,1);hold on; ft_plot_mesh(headshapemri.bnd,'edgecolor','none','facecolor','b','fidcolor','y'); 
ft_plot_headshape(headshape); viewbottom();
subplot(2,2,2);hold on; ft_plot_mesh(headshapemri.bnd,'edgecolor','none','facecolor','b','fidcolor','y'); 
ft_plot_headshape(headshape); viewtop();
subplot(2,2,3);hold on; ft_plot_mesh(headshapemri.bnd,'edgecolor','none','facecolor','b','fidcolor','y'); 
ft_plot_headshape(headshape); viewfront();
subplot(2,2,4);hold on; ft_plot_mesh(headshapemri.bnd,'edgecolor','none','facecolor','b','fidcolor','y'); 
ft_plot_headshape(headshape); viewright();
axis on;
grid on;
hcp_write_figure([filename,'_headshape.png'], 'resolution', 500);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function writetextfile(filename, str)
fid = fopen(filename, 'wt');
fwrite(fid, str);
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ft_plot_fiducials(fiducials)
x = fiducials.nas(1);
y = fiducials.nas(2);
z = fiducials.nas(3);
plot3(x, y, z, 'b*');
x = fiducials.lpa(1);
y = fiducials.lpa(2);
z = fiducials.lpa(3);
plot3(x, y, z, 'b*');
x = fiducials.rpa(1);
y = fiducials.rpa(2);
z = fiducials.rpa(3);
plot3(x, y, z, 'b*');


