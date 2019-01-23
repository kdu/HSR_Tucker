% EXPERIMENT 1: HOSVD FOR VARIOUS RANKS %
% Copyright (c) 2018 Clemence Prevost, Konstantin Usevich, Pierre Comon, David Brie
% https://github.com/cprevost4/HSR_Tucker
% Contact: clemence.prevost@univ-lorraine.fr

%% Load Pavia University

% 1. load data
SRI = cell2mat(struct2cell(load('PaviaU.mat')));
SRI(1:2,:,:) = []; SRI(:,1:4,:) = [];
% 2. degradation
Pm = spectral_deg(SRI,"Quickbird");
MSI = tmprod(SRI,Pm,3);
d1 = 4; d2 = 4; q = 7;
[P1,P2] = spatial_deg(SRI, q, d1, d2);
HSI = tmprod(tmprod(SRI,P1,1),P2,2);
% 3. add noise
for k=1:size(HSI,3)
    HSI(:,:,k) = awgn(HSI(:,:,k),15);
end
for k=1:size(MSI,3)
    MSI(:,:,k) = awgn(MSI(:,:,k),25);
end

%% Load Cuprite

load Cuprite.mat
SRI = double(X); clear X; %Convert groundtruth to double
SRI(:,:,[1:2 104:113 148:167 221:224]) = []; %Regions of water absorption
d1 = 4; d2 = 4; %Spatial downsampling ratios
[P1,P2] = spatial_deg(SRI,7,d1,d2); Pm = spectral_deg(SRI,"LANDSAT");
HSI = tmprod(SRI,{P1,P2},[1,2]); MSI = tmprod(SRI,Pm,3);
for k=1:size(HSI,3)
    HSI(:,:,k) = awgn(HSI(:,:,k),15);
end
for k=1:size(MSI,3)
    MSI(:,:,k) = awgn(MSI(:,:,k),25);
end

%% Table for Pavia University

methods2 = {'SCUBA'  'scuba_adaptor' '[120,3]' [4,4]; ...
           'B-SCOTT' 'bscott_adaptor' '[size(HSI,1),size(HSI,2),3]' [4,4]; ...                 
           'SCUBA'  'scuba_adaptor' '[300,9]' [1,1]; ... 
           'B-SCOTT' 'bscott_adaptor' '[size(HSI,1),size(HSI,2),size(MSI,3)]' [1,1]; ...
           'SCUBA'  'scuba_adaptor' '[120,3]' [8,8]; ... 
           'B-SCOTT' 'bscott_adaptor' '[size(MSI,1)/8,size(MSI,2)/8,3]' [8,8];};      
res = compare_methods(SRI, HSI, MSI, struct('Pm', Pm), [d1 d2], methods2);
T = cell2mat(res(:,2:end)); save('table_BPavia.txt','T','-ASCII')

%% Table for Cuprite

methods2 = {'SCUBA'  'scuba_adaptor' '[45,3]' [4,4]; ...
           'B-SCOTT' 'bscott_adaptor' '[45,45,3]' [4,4]; ...                 
           'SCUBA'  'scuba_adaptor' '[150,10]' [1,1]; ... 
           'B-SCOTT' 'bscott_adaptor' '[150,150,6]' [1,1]; ...
           'SCUBA'  'scuba_adaptor' '[45,3]' [8,8]; ... 
           'B-SCOTT' 'bscott_adaptor' '[45,45,3]' [8,8];};      
res = compare_methods(SRI, HSI, MSI, struct('Pm', Pm), [d1 d2], methods2);
T = cell2mat(res(:,2:end)); save('table_BCuprite.txt','T','-ASCII')



