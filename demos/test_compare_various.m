% EXPERIMENT 3: COMPARISON OF METRICS WITH STEREO %
% Copyright (c) 2018 Clemence Prevost, Konstantin Usevich
% https://github.com/cprevost4/HSR_Tucker
% Contact: clemence.prevost@univ-lorraine.fr
addpath ../utils
addpath ../methods
addpath ../metrics


%% Indian Pines (removed regions of water absorption):
% 1. load data
SRI = cell2mat(struct2cell(load('Indian_pines_corrected.mat')));
SRI = SRI(2:end,2:end,:);  % crop image to 144x144
% 2. spectral degradation
Pm = spectral_deg(SRI,"LANDSAT");
d1 = 4; d2 = 4;
[P1,P2] = spatial_deg(SRI, 9, d1, d2);
MSI = tmprod(SRI,Pm,3);
HSI = tmprod(tmprod(SRI,P1,1),P2,2);
% 3. run metrics
methods = {'STEREO' 'stereo3' '50'; ...
           'STEREO' 'stereo3' '100'; ...
           'SCOTT' 'scott' '[40,40,6]'; ...
           'SCOTT' 'scott' '[30,30,16]'; ...
           'SCOTT' 'scott' '[24,24,25]';...
           'B-SCOTT' 'bscott_b1_adaptor' '[40,40,6]'; ...                 
           'B-SCOTT' 'bscott_b1_adaptor' '[60,60,6]'; ...                 
           'B-SCOTT' 'bscott_b1_adaptor' '[100,100,6]'};      
DegMat = struct('Pm', Pm, 'P1', P1, 'P2', P2);         
res = compare_methods(SRI, HSI, MSI, DegMat, [d1 d2], methods);           

%% Salinas-A scene


%% Pavia dataset
SRI = cell2mat(struct2cell(load('PaviaU.mat')));
SRI(1:2,:,:) = []; SRI(:,1:4,:) = [];
%SRI = SRI(:,1:320,:);
%SRI(1,:,:) = []; SRI(:,1,:) = [];
% 
Pm = spectral_deg(SRI,"Quickbird");
MSI = tmprod(SRI,Pm,3);
d1 = 4; d2 = 4; q = 9;
[P1,P2] = spatial_deg(SRI, q, d1, d2);
HSI = tmprod(tmprod(SRI,P1,1),P2,2);
        
methods2 = {'SCUBA'  'scuba_adaptor' '[120,3]' [4,4]; ...
           'B-SCOTT' 'bscott_adaptor' '[100,100,4]' [1,1]; ...                 
           'B-SCOTT' 'bscott_adaptor' '[60,60,3]' [4,4]; ...                 
           'B-SCOTT' 'bscott_adaptor' '[120,60,4]' [4,4]};      
res = compare_methods(SRI, HSI, MSI, struct('Pm', Pm), [d1 d2], methods2);           


