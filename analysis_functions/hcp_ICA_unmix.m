function [iteration, data_meg] = hcp_ICA_unmix(data_meg, options)

% hcp_ICA_unmix performs ICA on a dataset.
% Output data is a structure containing the ft_componentanalysis result of
% each iteration.
%
% Use as
% [iteration] = hcp_ICA_unmix(data_meg, options)
%
% where data_meg is a data-structure containing the data to be unmixed
% and options is a cell-array specifying the behaviour of the
% algorithm. Cell-arrays need to be organized as sets of key-value pairs,
% i.e. {'key1', 'value1', ...}.
%
% Options needs to contain the following key:
%   channel: channel selection
%   ica_iterations: number of ICA iterations to be performed
%
% The following steps are performed:
% -FastICA analysis
%
%
% Example use:
% options = {'channel', 'MEG', 'ica_iterations', 1};
% iteration = hcp_ICA_unmix(filename, options);

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

% call the old-fashioned way (i.e. also do preprocessing) when first input
% is a string
if ischar(data_meg)
    data_meg = hcp_ICA_preprocessing(data_meg, options);
end

%%%%%%%%% deal with the input options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the options
ica_iterations = ft_getopt(options, 'ica_iterations', 1);
saveIC         = ft_getopt(options, 'save_comp', 'yes');
selchan        = ft_getopt(options, 'channels',  'MEG');

%%%%%%%%% Performing ICA  using the FastICA algorithm%%%%%%%%%%%%%%%%%%%%%%

disp('performing ICA...');

cfg                  = [];
cfg.method           = 'fastica';
cfg.fastica.approach = 'defl';
cfg.fastica.g        ='tanh';
cfg.fastica.maxNumIterations = 200;
cfg.fastica.numOfIC     =  size(data_meg.label,1);
cfg.fastica.verbose     = 'on';
cfg.fastica.displayMode = 'off';
cfg.demean              = 'no';

for j = 1:ica_iterations
    
    pack;
    randn('state',j);
    cfg.fastica.initGuess = randn(numel(data_meg.label));
    
    comp = hcp_componentanalysis(cfg, data_meg);
    A    = comp.topo;
    W    = comp.unmixing;
    
    [chan, Nc] = size(A);
    
    %%% Ordering ICs according to the Power
    if Nc > 0
        power          = mean(A.^2);
        [~, order] = sort(power, 'descend');
        
        A = A(:,order);
        W = W(order,:);
    end
    
    % copy the stuff back into the output structure
    comp.order    = order;
    comp.topo     = A;
    comp.unmixing = W;
    if strcmp(saveIC,'yes')
        for iic=1:size(comp.sampleinfo,1), comp.trial{iic} = comp.trial{iic}(order,:); end
    else
        comp = rmfield(comp, 'trial');
        comp = rmfield(comp, 'time');
    end
    
    % put in the output structure
    iteration(j).comp = comp; % structure containing all the iterations
    
    clear comp A W IC power order
end
