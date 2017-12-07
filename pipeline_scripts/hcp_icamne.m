%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setup the execution environment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opengl software;

% ensure that the time and date of execution are not stored in the provenance information
global ft_default
ft_default.trackcallinfo = 'no';

% allow the user to specify the path where additional data is present, e.g. the channel layout or anatomy files
if exist('path', 'var')
    addpath(path)
end

if ~exist('filename', 'var')
    error('filename should be specified')
end

% the filename is assumed to be something like
% 'rawdatadir/Phase1MEG/Subjects/CP10018/Experiments/CP10018_MEG/Scans/1-Rnoise_MNN_V1/Resources/4D/c,rfDC'
tok = tokenize(filename, '/');

if ~exist('subjectid', 'var')
    subjectid = tok{end-7};
end

if ~exist('experimentid', 'var')
    experimentid = tok{end-5};
end

if ~exist('scanid', 'var')
    scanid = tok{end-3};
end

if ~exist('pipelinedatadir', 'var')
    pipelinedatadir = hcp_pathdef;
end

resultprefix = sprintf('%s_%s', experimentid, scanid);

% change to the location of the processed data (input and output)
cd(pipelinedatadir)

% hcp_check_pipelineoutput('baddata', 'subject', subjectid, 'experiment', experimentid, 'scan', scanid);
% hcp_check_pipelineoutput('icaclass', 'subject', subjectid, 'experiment', experimentid, 'scan', scanid);
% hcp_check_pipelineoutput('anatomy', 'subject', subjectid);

% print the matlab and megconnectome version to screen for provenance
ver('megconnectome')

% print the value of all local variables to screen for provenance
w = whos;
w = {w.name};
w = setdiff(w, {'w', 'ans'});
for i=1:length(w)
    fprintf(hcp_printstruct(w{i}, eval(w{i})));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% execute the pipeline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hcp_read_matlab([resultprefix '_icaclass_vs.mat'])
% comp_class.grad=rmfield(comp_class.grad,'balance');

% read the source and volume conduction model from current dir with
% outputs of previous pipelines

hcp_read_matlab([subjectid '_MEG_anatomy_sourcemodel_2d']);
sourcemodel2d=ft_convert_units(sourcemodel2d, 'cm');
sourcemodel2d.inside = 1:size(sourcemodel2d.pos,1);
sourcemodel2d.outside = [];
sourcemodelsubj = sourcemodel2d;

hcp_read_matlab(sprintf('%s.mat', [subjectid '_MEG_anatomy_headmodel']));
headmodel = ft_convert_units(headmodel, 'cm');

mri = ft_read_mri([subjectid '_MEG_anatomy_anatomical.nii']);
mri=ft_convert_units(mri, 'cm');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ensure the correct geometrical units
% sourcemodel = ft_convert_units(sourcemodelsubj, 'mm');
grad = ft_read_sens(filename);
gradBalanced = grad;
gradBalanced = ft_apply_montage(gradBalanced, gradBalanced.balance.Supine, 'keepunused', 'yes', 'inverse', 'yes');
grad=gradBalanced;
grad = ft_convert_units(grad, 'cm');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deal with the component data in order for ft_sourceanalysis to be able to
% swallow it
mixing = comp_class.topo;
channels = comp_class.topolabel;
% normalisation of the topographies
for i = 1:size(mixing, 2)
    val(i) = 0.01*max(abs(mixing(:, i)));
    mixing(:, i) = mixing(:, i)/val(i);
end


% create a 'timelock' structure
tlck = [];
tlck.label = channels;
tlck.cov = eye(numel(tlck.label)); % perhaps this one should be scaled
tlck.time=1;
tlck.grad = grad;
tlck.dimord = 'chan_time';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare the leadfields: this needs to be done only once:
% if I understand the method well, each component will be reconstructed by
% using an MNE with different regularisation

% keep parcel info if present
if isfield(sourcemodelsubj, 'label')
    %     parcellation.label = sourcemodelsubj.label;
    %     parcellation.labelindx = sourcemodelsubj.labelindx;
    %     parcellation.annotation = sourcemodelsubj.annotation;
    %     parcellation.ctable = sourcemodelsubj.ctable;
    parcellation.label = sourcemodelsubj.aparclabel;
    parcellation.labelindx = sourcemodelsubj.aparc;
    %     parcellation.annotation = sourcemodelsubj.annotation;
    %     parcellation.ctable = sourcemodelsubj.ctable;
end


cfg = [];
cfg.vol = headmodel;
cfg.grid = sourcemodelsubj;
cfg.grad = grad;
cfg.channel = channels;
cfg.normalize = 'yes';
cfg.reducerank = 2;
gridLF = ft_prepare_leadfield(cfg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do the source reconstruction. In this case all components are
% reconstructed using the same amount of regularisation

noise_level=8;

cfg = [];
cfg.method = 'mne';
cfg.grid = gridLF;
cfg.vol = headmodel;
cfg.channel = channels;
cfg.mne.prewhiten = 'yes';
cfg.mne.noisecov=eye(numel(channels))*noise_level;

nic=size(mixing,2);
for i=1:nic
    tlck.avg = mixing(:,i);
    
    mgfp(i)=sqrt(mean((mixing(:,i)-mean(mixing(:,i))).^2));
    cfg.mne.snr=mgfp(i)/noise_level;
    
    noisevec(i)=cfg.mne.snr;
    
    tmp = ft_sourceanalysis(cfg, tlck);
    if i==1
        source=tmp;
    else
        % concatenate the mom here
        for k = 1:numel(tmp.inside)
            source.avg.mom{source.inside(k)} = cat(2,source.avg.mom{source.inside(k)}, tmp.avg.mom{source.inside(k)});
        end
        source.avg.pow = horzcat(source.avg.pow,tmp.avg.pow);
    end
end
source.val=val;
source.time=1:size(mixing,2);
source.noise=noisevec;


if(isfield(sourcemodelsubj,'tri')) source.tri=sourcemodelsubj.tri; end

if exist('parcellation', 'var')
    % store them back in the source structure
    source.label = parcellation.label;
    source.labelindx = parcellation.labelindx;
    %     source.annotation = parcellation.annotation;
end


cfgtopo =[];
cfgtopo.grad=grad;
cfgtopo.zlim='maxabs';
cfgtopo.component=[];
cfgtopo.parameter = 'topo';
cfgtopo.comment = 'no';
cfgtopo.colorbar = 'yes';
cfgtopo.layout='4D248.mat';
tmpclass=comp_class;
tmpclass.trial{1}(1,1)=0;
tmpclass.time{1}(1,1)=0;
tmpclass.topo=mixing;

X = 29.7;                  %# A3 paper size
Y = 14.35;                  %# A3 paper size
xMargin = 1;               %# left/right margins from page borders
yMargin = 1;               %# bottom/top margins from page borders
xSize = X - 2*xMargin;     %# figure size on paper (widht & hieght)
ySize = Y - 2*yMargin;     %# figure size on paper (widht & hieght)
%     figure('Menubar','none');


cfg = [];
cfg.funparameter = 'avg.pow'; 
cfg.method = 'surface';
cfg.funcolormap='jet';
cfg.interactive = 'no';

tmp = source;
for k = 1:numel(source.time)
    tmp.avg.pow = source.avg.pow(:, k);
    tmp.time=source.time(k);
    %         tmp.avg.mask=zeros(size(tmp.avg.pow));
    maxabs=max(max(max(tmp.avg.pow)));
    %         indxmask=find(tmp.avg.pow>0.4*maxabs);
    %         tmp.avg.mask(indxmask)=1;
    imgname = [resultprefix '_icamne_' num2str(k) '.png'];
    
    ft_plot_mesh(tmp,'vertexcolor',tmp.avg.pow)
    view(-90,90)
    
    h1=gcf;
    set(gca, 'XTickLabel',[], 'YTickLabel',[], ...
        'Units','normalized', 'Position',[0 0 1 1])
    
    %# figure size on screen (50% scaled, but same aspect ratio)
    set(gcf, 'Units','centimeters', 'Position',[5 5 xSize ySize])
    
    %# figure size printed on paper
    set(gcf, 'visible', 'on')
    set(gcf, 'PaperUnits','centimeters')
    set(gcf, 'PaperSize',[X Y])
    set(gcf, 'PaperPosition',[xMargin yMargin xSize ySize])
    set(gcf, 'PaperOrientation','portrait')
    
    fax(1)=gca;
    set(fax(1),'position',[0.25 0.1 0.95 0.8]);
    fax(2)= axes('position',[0 0.1 0.5 0.8]);
    
    cfgtopo.component=k;
    ft_topoplotIC(cfgtopo, tmpclass);
    
    hcp_write_figure(imgname, h1)
    close(h1)
end
hcp_write_matlab([resultprefix,'_icamne'],'source');



% ensure that the expected output files were created
% hcp_check_pipelineoutput('icamne', 'subject', subjectid, 'experiment', experimentid, 'scan', scanid);