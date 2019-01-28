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

%% Images for Pavia

methods = {'SCUBA'  'scuba_adaptor' '[120,3]' [4,4]; ...
           'B-SCOTT' 'bscott_adaptor' '[size(HSI,1),size(HSI,2),3]' [4,4]; ...                 
           'SCUBA'  'scuba_adaptor' '[300,9]' [1,1]; ... 
           'B-SCOTT' 'bscott_adaptor' '[size(HSI,1),size(HSI,2),size(MSI,3)]' [1,1]; ...
           'SCUBA'  'scuba_adaptor' '[120,3]' [8,8]; ... 
           'B-SCOTT' 'bscott_adaptor' '[size(MSI,1)/8,size(MSI,2)/8,3]' [8,8];};      
[~, est] = compare_methods(SRI, HSI, MSI, struct('Pm', Pm), [d1 d2], methods);

figure
subplot(2,3,1); imagesc(real(est{1}(:,:,44))); ...
    title(sprintf('%s, %s, [%d %d]',methods{1,1},methods{1,3},methods{1,4}(1),methods{1,4}(2)))
subplot(2,3,4); imagesc(real(est{2}(:,:,44)));...
    title(sprintf('%s, %s, [%d %d]',methods{2,1},methods{2,3},methods{2,4}(1),methods{2,4}(2)))
subplot(2,3,2); imagesc(real(est{3}(:,:,44)));...
    title(sprintf('%s, %s, [%d %d]',methods{3,1},methods{3,3},methods{3,4}(1),methods{3,4}(2)))
subplot(2,3,5); imagesc(real(est{4}(:,:,44)));...
    title(sprintf('%s, %s, [%d %d]',methods{4,1},methods{4,3},methods{4,4}(1),methods{4,4}(2)))
subplot(2,3,3); imagesc(real(est{5}(:,:,44)));...
    title(sprintf('%s, %s, [%d %d]',methods{5,1},methods{5,3},methods{5,4}(1),methods{5,4}(2)))
subplot(2,3,6); imagesc(real(est{6}(:,:,44)));...
    title(sprintf('%s, %s, [%d %d]',methods{6,1},methods{6,3},methods{6,4}(1),methods{6,4}(2)))

%% Images for Cuprite

methods = {'SCUBA'  'scuba_adaptor' '[45,3]' [4,4]; ...
           'B-SCOTT' 'bscott_adaptor' '[45,45,3]' [4,4]; ...                 
           'SCUBA'  'scuba_adaptor' '[150,10]' [1,1]; ... 
           'B-SCOTT' 'bscott_adaptor' '[150,150,6]' [1,1]; ...
           'SCUBA'  'scuba_adaptor' '[45,3]' [8,8]; ... 
           'B-SCOTT' 'bscott_adaptor' '[45,45,3]' [8,8];};      
[~,est] = compare_methods(SRI, HSI, MSI, struct('Pm', Pm), [d1 d2], methods);

figure
subplot(2,3,1); imagesc(real(est{1}(:,:,44))); ...
    title(sprintf('%s, %s, [%d %d]',methods{1,1},methods{1,3},methods{1,4}(1),methods{1,4}(2)))
subplot(2,3,4); imagesc(real(est{2}(:,:,44)));...
    title(sprintf('%s, %s, [%d %d]',methods{2,1},methods{2,3},methods{2,4}(1),methods{2,4}(2)))
subplot(2,3,2); imagesc(real(est{3}(:,:,44)));...
    title(sprintf('%s, %s, [%d %d]',methods{3,1},methods{3,3},methods{3,4}(1),methods{3,4}(2)))
subplot(2,3,5); imagesc(real(est{4}(:,:,44)));...
    title(sprintf('%s, %s, [%d %d]',methods{4,1},methods{4,3},methods{4,4}(1),methods{4,4}(2)))
subplot(2,3,3); imagesc(real(est{5}(:,:,44)));...
    title(sprintf('%s, %s, [%d %d]',methods{5,1},methods{5,3},methods{5,4}(1),methods{5,4}(2)))
subplot(2,3,6); imagesc(real(est{6}(:,:,44)));...
    title(sprintf('%s, %s, [%d %d]',methods{6,1},methods{6,3},methods{6,4}(1),methods{6,4}(2)))

