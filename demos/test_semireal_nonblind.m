% TABLES FOR SEMIREAL DATA %
% Copyright (c) 2018 Clemence Prevost, Konstantin Usevich, Pierre Comon, David Brie
% https://github.com/cprevost4/HSR_Tucker
% Contact: clemence.prevost@univ-lorraine.fr

%% Table for Indian Pines

SRI = cell2mat(struct2cell(load('Indian_pines.mat')));
SRI(:,:,[104:108 150:163 220]) = []; %Regions of water absorption
SRI(1,:,:) = []; SRI(:,1,:) = [];
Pm = spectral_deg(SRI,"LANDSAT");
MSI = tmprod(SRI,Pm,3);
d1 = 4; d2 = 4; q = 9;
[P1,P2] = spatial_deg(SRI, q, d1, d2);
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
T = cell2mat(res(:,2:end)); save('exp3_table2.txt','T','-ASCII')


%% Table for Salinas A-scene

% 1. load data
SRI = cell2mat(struct2cell(load('Salinas.mat')));
SRI = crop(SRI,[80,84,size(SRI,3)]);
SRI(:,:,[108:112 154:167 224]) = []; %Regions of water absorption (Salinas)
% 2. degradation
Pm = spectral_deg(SRI,"LANDSAT");
d1 = 4; d2 = 4; q = 9;
[P1,P2] = spatial_deg(SRI, q, d1, d2);
MSI = tmprod(SRI,Pm,3); HSI = tmprod(tmprod(SRI,P1,1),P2,2);

% 3. run metrics
methods = {'STEREO' 'stereo3' '50'; ...
           'STEREO' 'stereo3' '100'; ...
           'SCOTT' 'scott' '[40,40,6]'; ...
           'SCOTT' 'scott' '[14,14,15]'; ...
           'SCOTT' 'scott' '[10,15,25]';...
           'SCOTT' 'scott' '[30,30,6]';...                
           'SCOTT' 'scott' '[58,58,6]'};      
DegMat = struct('Pm', Pm, 'P1', P1, 'P2', P2);         
res = compare_methods(SRI, HSI, MSI, DegMat, [d1 d2], methods); 
T = cell2mat(res(:,2:end)); save('exp3_table2_sal.txt','T','-ASCII')