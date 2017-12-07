function [net_seeds] = hcp_mcw_netdef(net,options)

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

% network= ft_getopt(options, 'network', 'DMN');
if isempty(net) net='DMN'; end
custom_seeds = ft_getopt(options, 'custom_seeds');

net_seeds.pos=zeros(5,3);
%____________ EXT NODE
net_seeds.pos(5,:)=[9 42 53];

load('hcp_cortex_seeds')


if isempty(custom_seeds)
if(strcmp(net,'DMN'))
%'DMN'   
label{1}='preCunPC' ;  hemi{1}='L';
label{2}='mPFC2';      hemi{2}='L';    
label{3}='AG';         hemi{3}='R';
label{4}='AG';         hemi{4}='L';
net_name={'a3_Default_mode'};
par_name={'L_G_pariet_inf-Angular'};

elseif(strcmp(net,'DAN'))
%'DMN'   
label{1}='pIPS-SPL' ;   hemi{1}='L';
label{2}='pIPS-SPL';    hemi{2}='R';    
label{3}='FEF';         hemi{3}='L';
label{4}='FEF';         hemi{4}='R';
net_name={'a15_Dorsal_attention'};


elseif(strcmp(net,'MN'))
%'DMN'   
label{1}='CS' ;   hemi{1}='L';
label{2}='CS';    hemi{2}='R';    
label{3}='S2';    hemi{3}='L';
label{4}='vCS';   hemi{4}='R'; % this is not exactly the same as the previous implementation
% label{4}='vPoCe';   hemi{4}='R'; % this is not exactly the same as the previous implementation
% label{4}='mI2';   hemi{4}='R'; % this is not exactly the same as the previous implementation
net_name={'a4_"Hand"_somatosensory-motor' 'a16_"Mouth"_somatosensory-motor'};

elseif(strcmp(net,'VIS'))
%'DMN'   
label{1}='CS' ;   hemi{1}='L';
label{2}='CS';    hemi{2}='R';    
label{3}='S2';    hemi{3}='L';
label{4}='vCS';   hemi{4}='R'; % this is not exactly the same as the previous implementation
% label{4}='vPoCe';   hemi{4}='R'; % this is not exactly the same as the previous implementation
% label{4}='mI2';   hemi{4}='R'; % this is not exactly the same as the previous implementation
net_name={'a5_Visual'};

elseif(strcmp(net,'AUD'))
%'DMN'   
label{1}='CS' ;   hemi{1}='L';
label{2}='CS';    hemi{2}='R';    
label{3}='S2';    hemi{3}='L';
label{4}='vCS';   hemi{4}='R'; % this is not exactly the same as the previous implementation
% label{4}='vPoCe';   hemi{4}='R'; % this is not exactly the same as the previous implementation
% label{4}='mI2';   hemi{4}='R'; % this is not exactly the same as the previous implementation
net_name={'a24_Auditory'};
end

% 
% %___ DMN
% 
% if(strcmp(net_seeds,'dmn_pcc'))
% %%_________ DMN NET
% %__ pcc
% % label(1,:)=[-3 -54 31];
% % DMN preCunPC   L	   
% label(1,:)=[-7.78  -49.59   30.99];
% label_name{1}= 'lPCPC';
% 
% %__ lfp
% %label(2,:)=[-2 50 2];
% %DMN mPFC2      L	   
% label(2,:)=[-1.27   46.05   23.16];
% label_name{2}='LmPFC'; %This is a little bit hihger than the one we used before
% %DMN mPFC1      R	    
% % label(2,:)=[3.71   49.05   11.26]; % This is actually RmPFC. Perhaps we can use 
% % label_name(2,:)='RmPFC';
% 
% %__ rag
% label(3,:)=[51 -64 32];
% % DMN AG         R	   
% label(3,:)=[47.97  -64.56   42.25];
% label_name{3}= 'rAG';
% 
% %__ lag
% % label(4,:)=[-43 -76 35];
% % DMN AG         L	  
% label(4,:)=[-44.42  -62.52   35.72];
% label_name{4}= 'lAG';
% 
% 
% %___ DAN
% 
% elseif(strcmp(net_seeds,'dan_lpips'))
% %_________ DAN NET
% %__ lpips
% label(1,:)=[-25 -67 48];
% %DAN pIPS-SPL   L	  
% label(1,:)=[-15.73 -64.01 53.13];
% label_name{1}= 'lpIPS';
% 
% %__ rpips
% label(2,:)=[23 -69 49];
% %DAN pIPS-SPL   R	   
% label(2,:)=[21.85 -65.83 46.07];
% label_name{2}= 'rpIPS';
% 
% %__ lfef
% label(3,:)=[-26 -12 53];
% %DAN FEF        L	  
% label(3,:)=[-18.52 -4.02 58.16];
% label_name{3}= 'lFEF';
% 
% %__ rfef
% label(4,:)=[30 -13 53];
% %DAN FEF        R	   
% label(4,:)=[25.44 -3.59 55.05];
% label_name{4}= 'rFEF';
% 
% %___ DMN 
% elseif(strcmp(net_seeds,'dmn_lmpfc'))
% %%_________ DMN NET
% %__ lag
% label(2,:)=[-43 -76 35];
% %__ pcc
% label(4,:)=[-3 -54 31];
% %__ rag
% label(3,:)=[51 -64 32];
% %__ lfp
% label(1,:)=[-2 50 2];
%  
% 
% 
% %___ VIS
% 
% elseif(strcmp(net_seeds,'vis_lv1'))
% % % %__________ VIS RSN NODES
% % 
% % % % % _______ LV1
% %   label(1,:)=[-2.4	-99.2	-3.5];
% % VPN V1v        L	   
% label(1,:)=[-8.73  -89.12   -4.73];
% label_name{1}= 'lV1V';
% %   VPN V1d-V2d    L	   -2.10  -91.93   11.14
% 
% % % % %__________ RV1 ______________
% %   label(2,:)=[10.3 -91 3.5];
% %  VPN V1         R	   
% label(2,:)=[14.68  -86.50   10.45];
% label_name{2}= 'rV1';
% %   VPN V1d        R	    6.39  -77.08   13.25
% 
% %  % % %_________RV7 _____________
% %   label(3,:)=[30.4	-78.7	22.7];
% %   VPN POSd       R	   
%   label(3,:)=[19.32  -75.62   28.73]; % This is not exactly RV7
%   label_name{3}= 'rPOSd(v7)';
% %   VPN V3-V3A     R	   12.10  -85.04   35.78
% %   VPN V4v        R	   23.72  -74.13   -9.49
% 
% % % % % _______ LV7 dorsal
% %   label(4,:)=[-22.4	-77.9	22.5];
% %   VPN V7         L	  
%   label(4,:)=[-23.05  -71.18    8.28];
%   label_name{4}= 'lV7)';
% 
% 
% % % VPN V3-V3A     L	  -13.40  -93.97   22.07
% % % VPN V3A        L	  -14.16  -93.76   34.75
% % % VPN POSd       L	  -24.24  -55.58   -1.06
% % % VPN VP         L	  -11.81  -74.43   -4.46
% % % VPN V7-POSd    L	   -6.67  -86.22   40.49
% % % VPN V3-V3A     R	   21.18  -92.55   21.22
% % % VPN V1v-RV2v   R	   10.21  -77.33   -3.72
% % % VPN VP         R	   17.58  -63.60   -4.68
% % % VPN POSv       R	   20.81  -48.82   -3.50
% 
% 
% elseif(strcmp(net_seeds,'vis_rv1'))
% 
% %___ VIS RSN NODES
% % 
% % % % % _______ LV1
%   label(2,:)=[-2.4	-99.2	-3.5];
% % % % %__________ RV1 ______________
%   label(1,:)=[10.3 -91 3.5];
% %  % % %_________RV7 _____________
%   label(3,:)=[30.4	-78.7	22.7];
% % % % % _______ LV7 dorsal
%   label(4,:)=[-22.4	-77.9	22.5];
% 
% 
% 
% elseif(strcmp(net_seeds,'vis_rv7'))
% 
% %___ VIS RSN NODES
% % 
% % % % % _______ LV1
%   label(3,:)=[-2.4	-99.2	-3.5];
% % % % %__________ RV1 ______________
%   label(2,:)=[10.3 -91 3.5];
% %  % % %_________RV7 _____________
%   label(1,:)=[30.4	-78.7	22.7];
% % % % % _______ LV7 dorsal
%   label(4,:)=[-22.4	-77.9	22.5];   
% 
% 
% 
% elseif(strcmp(net_seeds,'vis_lv7'))
%     
% %___ VIS RSN NODES
% % 
% % % % % _______ LV1
%   label(4,:)=[-2.4	-99.2	-3.5];
% % % % %__________ RV1 ______________
%   label(2,:)=[10.3 -91 3.5];
% %  % % %_________RV7 _____________
%   label(3,:)=[30.4	-78.7	22.7];
% % % % % _______ LV7 dorsal
%   label(1,:)=[-22.4	-77.9	22.5];
% 
%   
% 
% %___ VAN
% 
% elseif(strcmp(net_seeds,'van_rstg'))
% % % %__________ VAN RSN NODES
% % 
% % % %__________ RSTG ______________
% label(1,:)=[57.8 -48.2 10.4];
%  % % %_________RMFG _____________
% label(2,:)=[41.8 17.2 31.3];
% % % % _______ RPCS
% label(3,:)=[41.1 1.8  50.2];
%  % _______ RVFG
% label(4,:)=[39.9 20.8 -3.8];
% 
% 
% 
% elseif(strcmp(net_seeds,'van_rmfg'))
% % % %__________ VAN RSN NODES
% % 
% % % %__________ RSTG ______________
% label(3,:)=[57.8 -48.2 10.4];
%  % % %_________RMFG _____________
% label(1,:)=[41.8 17.2 31.3];
% % % % _______ RPCS
% label(2,:)=[41.1 1.8  50.2];
%  % _______ RVFG
% label(4,:)=[39.9 20.8 -3.8];
% 
% 
% 
% elseif(strcmp(net_seeds,'van_rpcs'))
% % % %__________ VAN RSN NODES
% % 
% % % %__________ RSTG ______________
% label(4,:)=[57.8 -48.2 10.4];
%  % % %_________RMFG _____________
% label(2,:)=[41.8 17.2 31.3];
% % % % _______ RPCS
% label(1,:)=[41.1 1.8  50.2];
%  % _______ RVFG
% label(3,:)=[39.9 20.8 -3.8];
% 
% 
% 
% elseif(strcmp(net_seeds,'van_rvfg'))
% % % %__________ VAN RSN NODES
% % 
% % % %__________ RSTG ______________
% label(4,:)=[57.8 -48.2 10.4];
%  % % %_________RMFG _____________
% label(2,:)=[41.8 17.2 31.3];
% % % % _______ RPCS
% label(3,:)=[41.1 1.8  50.2];
%  % _______ RVFG
% label(1,:)=[39.9 20.8 -3.8];
% 
% 
% 
% %___ DMN
% 
% elseif(strcmp(net_seeds,'dmn_lag'))
% %%_________ DMN NET
% %__ lag
% label(1,:)=[-43 -76 35];
% %__ pcc
% label(4,:)=[-3 -54 31];
% %__ rag
% label(3,:)=[51 -64 32];
% %__ lfp
% label(2,:)=[-2 50 2];
% 
% 
% 
% elseif(strcmp(net_seeds,'dmn_rag'))
% %%_________ DMN NET
% %__ lag
% label(4,:)=[-43 -76 35];
% %__ pcc
% label(3,:)=[-3 -54 31];
% %__ rag
% label(1,:)=[51 -64 32];
% %__ lfp
% label(2,:)=[-2 50 2];
% 
% 
% 
% elseif(strcmp(net_seeds,'dmn_lfp'))
% %%_________ DMN NET
% %__ lag
% label(4,:)=[-43 -76 35];
% %__ pcc
% label(3,:)=[-3 -54 31];
% %__ rag
% label(2,:)=[51 -64 32];
% %__ lfp
% label(1,:)=[-2 50 2];
% 
% 
% 
% %___ MOT
% 
% elseif(strcmp(net_seeds,'mot_ls2'))
% %__________ MOT RSN NODES
% %______________ LS2
% label(1,:)=[-39.4	-26.7	18.2];
% % % _______ RS2
% label(2,:)=[36.4	-22.7	20.6];
% % % %__________ lcs_stef ______________
% label(3,:)=[-32 -25 55];
% % %__________ rcs_stef ______________
% label(4,:)=[35 -26 55];
% 
% 
% 
% elseif(strcmp(net_seeds,'mot_rs2'))
% %__________ MOT RSN NODES
% %______________ LS2
% label(2,:)=[-39.4	-26.7	18.2];
% % % _______ RS2
% label(1,:)=[36.4	-22.7	20.6];
% % % %__________ lcs_stef ______________
% label(3,:)=[-32 -25 55];
% % %__________ rcs_stef ______________
% label(4,:)=[35 -26 55];
% 
% 
% 
% elseif(strcmp(net_seeds,'mot_lcs'))
% %__________ MOT RSN NODES
% %______________ LS2
% label(2,:)=[-39.4	-26.7	18.2];
% % % _______ RS2
% label(3,:)=[36.4	-22.7	20.6];
% % % %__________ lcs_stef ______________
% label(1,:)=[-32 -25 55];
% % %__________ rcs_stef ______________
% label(4,:)=[35 -26 55];
% 
% 
% 
% elseif(strcmp(net_seeds,'mot_rcs'))
% %__________ MOT RSN NODES
% %______________ LS2
% label(2,:)=[-39.4	-26.7	18.2];
% % % _______ RS2
% label(3,:)=[36.4	-22.7	20.6];
% % % %__________ lcs_stef ______________
% label(4,:)=[-32 -25 55];
% % %__________ rcs_stef ______________
% label(1,:)=[35 -26 55];
% 
% 
% 
% %___ DAN
% 
% elseif(strcmp(net_seeds,'dan_lpips'))
% %_________ DAN NET
% %__ lpips
% label(1,:)=[-25 -67 48];
% %__ rpips
% label(2,:)=[23 -69 49];
% %__ lfef
% label(3,:)=[-26 -12 53];
% %__ rfef
% label(4,:)=[30 -13 53];
% 
% 
% 
% elseif(strcmp(net_seeds,'dan_rpips'))
% %_________ DAN NET
% %__ lpips
% label(2,:)=[-25 -67 48];
% %__ rpips
% label(1,:)=[23 -69 49];
% %__ lfef
% label(3,:)=[-26 -12 53];
% %__ rfef
% label(4,:)=[30 -13 53];
% 
% 
% 
% elseif(strcmp(net_seeds,'dan_lfef'))
% %_________ DAN NET
% %__ lpips
% label(2,:)=[-25 -67 48];
% %__ rpips
% label(3,:)=[23 -69 49];
% %__ lfef
% label(1,:)=[-26 -12 53];
% %__ rfef
% label(4,:)=[30 -13 53];
% 
% 
% 
% elseif(strcmp(net_seeds,'dan_rfef'))
% %_________ DAN NET
% %__ lpips
% label(2,:)=[-25 -67 48];
% %__ rpips
% label(3,:)=[23 -69 49];
% %__ lfef
% label(4,:)=[-26 -12 53];
% %__ rfef
% label(1,:)=[30 -13 53];
%  
% 
end
for i=1:4
    indx_hemi = find(strcmp(hcp_seeds.hemisphere,hemi{i}));
    indx_lab=find(strcmp(hcp_seeds.label,label{i}));
    
    [junk temp_indx]=ismember(indx_lab,indx_hemi);
    [junk junk2 temp_indx]=find(temp_indx);
    indx_seed=indx_hemi(temp_indx);
    
    net_seeds.pos(i,:)=hcp_seeds.pos(indx_seed,:);
    net_seeds.seeds_string{i}=[hemi{i} '-' label{i}];
    net_seeds.cortex_pos(i,:)=hcp_seeds.pos_cortex(indx_seed,:);
    net_seeds.cortex_index(i,:)=hcp_seeds.cortex_indx_both(indx_seed,1);
    net_seeds.cortex_index_hemi(i,:)=hcp_seeds.cortex_indx_lf(indx_seed,1);
end
net_seeds.seeds_string{5}=['out' '-' 'seed'];

net_seeds.hemi=hemi;
net_seeds.label=label;
net_seeds.net_name=net_name;
% net_seeds.par_name=par_name;
end


