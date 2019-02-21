% VISUAL FOR SEMIREAL DATA %
% Copyright (c) 2018 Clemence Prevost, Konstantin Usevich, Pierre Comon, David Brie
% https://github.com/cprevost4/HSR_Tucker
% Contact: clemence.prevost@univ-lorraine.fr

%% Simulations for Indian Pines

SRI = cell2mat(struct2cell(load('Indian_pines.mat')));
SRI(:,:,[104:108 150:163 220]) = []; %Regions of water absorption
SRI(1,:,:) = []; SRI(:,1,:) = [];
Pm = spectral_deg(SRI,"LANDSAT");
MSI = tmprod(SRI,Pm,3);
d1 = 4; d2 = 4; q = 9;
[P1,P2] = spatial_deg(SRI, q, d1, d2);
HSI = tmprod(tmprod(SRI,P1,1),P2,2);

methods = {'STEREO' 'stereo3' '50' []; ...
           'SCOTT' 'scott' '[40,40,6]' []; ...
           'HySure','hysure_b1_adaptor','[]',16};      
DegMat = struct('Pm', Pm, 'P1', P1, 'P2', P2);         
[res, est] = compare_methods(SRI, HSI, MSI, DegMat, [d1 d2], methods); 

%% Simulations for Salinas A-scene

% 1. load data
SRI = cell2mat(struct2cell(load('SalinasA.mat')));
%SRI = crop(SRI,[80,84,size(SRI,3)]);
SRI(:,:,[108:112 154:167 224]) = []; %Regions of water absorption (Salinas)
% 2. degradation
Pm = spectral_deg(SRI,"LANDSAT");
d1 = 4; d2 = 4; q = 9;
[P1,P2] = spatial_deg(SRI, q, d1, d2);
MSI = tmprod(SRI,Pm,3); HSI = tmprod(tmprod(SRI,P1,1),P2,2);

methods = {'STEREO' 'stereo3' '100' []; ...         
           'SCOTT' 'scott' '[58,58,6]' [];...
           'HySure','hysure_b1_adaptor','[]',6};      
DegMat = struct('Pm', Pm, 'P1', P1, 'P2', P2);         
[~,est] = compare_methods(SRI, HSI, MSI, DegMat, [d1 d2], methods); 

%% Make figure

figure
subplot(2,2,1); imagesc(SRI(:,:,44)); title('Groundtruth SRI'); colorbar; lim = caxis; axis off
subplot(2,2,2); imagesc(real(est{1}(:,:,44))); title(sprintf('%s, F=%s',methods{1,1},methods{1,3})); colorbar; caxis(lim); axis off
subplot(2,2,3); imagesc(real(est{2}(:,:,44))); title(sprintf('%s, R=%s',methods{2,1},methods{2,3})); colorbar; caxis(lim); axis off
subplot(2,2,4); imagesc(real(est{3}(:,:,44))); title(sprintf('%s, p=%d',methods{3,1},methods{3,4})); colorbar; caxis(lim); axis off
