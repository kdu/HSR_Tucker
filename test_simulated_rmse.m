% EXPERIMENT ON SIMULATED DATA %
% Copyright (c) 2018 Clemence Prevost, Konstantin Usevich
% https://github.com/cprevost4/HSR_Tucker
% Contact: clemence.prevost@univ-lorraine.fr

%% Load 6-material dataset

% 1. abundance maps
I = 240; J = 240; %spatial dimensions
A1 = zeros(I,J); A2 = A1; A3 = A1; A4 = A1; 
A5 = A1; A6 = A1; %initialize abundance maps
g1 = fspecial('gaussian',40, 15);
A1(1:40,1:40) = g1; 
A2(41:80,1:40) = g1; A2(1:40,41:80) = g1;
A3(81:120,1:40) = g1; A3(41:80,41:80) = g1; A3(1:40,81:120) = g1; 
A4(121:160,1:40) = g1; A4(1:40,121:160) = g1; 
A5(161:200,1:40) = g1; A5(121:160,41:80) = g1; A5(41:80,121:160) = g1; A5(1:40,161:200) = g1;
A6(201:240,1:40) = g1; A6(161:200,41:80) = g1; A6(121:160, 81:120) = g1;
    A6(81:120,121:160) = g1; A6(41:80,161:200) = g1; A6(1:40,201:240) = g1;
% 2. spectra
K = 200; %spectral dimension
SRI = cell2mat(struct2cell(load('Indian_pines.mat')));
SRI(:,:,[104:108 150:163 220]) = []; %Regions of water absorption
SRI(1,:,:) = []; SRI(:,1,:) = [];
load('Indian_pines_gt.mat'); indian_pines_gt(:,1) = []; indian_pines_gt(1,:) = [];
mat = tens2mat(SRI,3,[]);
ind = find(reshape(indian_pines_gt,1,[]) == 3);
s1 = real(mean(mat(:,ind),2)); 
ind = find(reshape(indian_pines_gt,1,[]) == 5);
s2 = real(mean(mat(:,ind),2)); 
ind = find(reshape(indian_pines_gt,1,[]) == 11);
s3 = real(mean(mat(:,ind),2)); 
ind = find(reshape(indian_pines_gt,1,[]) == 13);
s4 = real(mean(mat(:,ind),2)); 
ind = find(reshape(indian_pines_gt,1,[]) == 14);
s5 = real(mean(mat(:,ind),2)); 
ind = find(reshape(indian_pines_gt,1,[]) == 15);
s6 = real(mean(mat(:,ind),2)); 
clear SRI; clear mat; clear ind;
% 3; generate images
SRI = outprod(A1,s1) + outprod(A2,s2) + outprod(A3,s3) + outprod(A4,s4) + outprod(A5,s5) + outprod(A6,s6);
Pm = spectral_deg(SRI,"Quickbird");
MSI_gt = tmprod(SRI,Pm,3);
d1 = 4; d2 = 4; q = 9;
[P1,P2] = spatial_deg(SRI, q, d1, d2);
HSI_gt = tmprod(tmprod(SRI,P1,1),P2,2); 

%% Make figures - input SNR 20dB

for k=1:size(MSI_gt,3)
    MSI(:,:,k) = awgn(MSI_gt(:,:,k),20);
end
for k=1:size(HSI_gt,3)
    HSI(:,:,k) = awgn(HSI_gt(:,:,k),20);
end

methods = {'SCUBA'  'scuba_adaptor' '[6,4]' [1,1]; ...
           'SCUBA'  'scuba_adaptor' '[1,1]' [6,6]; ...
           'B-SCOTT' 'bscott_adaptor' '[6,6,4]' [1,1]; ...                 
           'B-SCOTT' 'bscott_adaptor' '[1,1,1]' [6,6];};               
[~,est] = compare_methods(SRI, HSI, MSI, struct('Pm', Pm), [4,4], methods);

err1 = zeros(1,size(SRI,3)); err2 = err1; err3 = err1; err4 = err1;
for k=1:size(SRI,3)
    err1(1,k) = norm(SRI(:,:,k)-est{1}(:,:,k),'fro')/norm(SRI(:,:,k),'fro');
    err2(1,k) = norm(SRI(:,:,k)-est{2}(:,:,k),'fro')/norm(SRI(:,:,k),'fro');
    err3(1,k) = norm(SRI(:,:,k)-est{3}(:,:,k),'fro')/norm(SRI(:,:,k),'fro');
    err4(1,k) = norm(SRI(:,:,k)-est{4}(:,:,k),'fro')/norm(SRI(:,:,k),'fro');
end

figure
subplot(1,2,1); plot(err1); hold on; plot(err3); xlabel('Spectral bands'); ylabel('RMSE'); title('(1,1) blocks'); legend('SCUBA','BSCOTT')
subplot(1,2,2); plot(err2); hold on; plot(err4); xlabel('Spectral bands'); ylabel('RMSE'); title('(6,6) blocks'); legend('SCUBA','BSCOTT')

%% Make figures - input SNR 60dB

for k=1:size(MSI_gt,3)
    MSI(:,:,k) = awgn(MSI_gt(:,:,k),60);
end
for k=1:size(HSI_gt,3)
    HSI(:,:,k) = awgn(HSI_gt(:,:,k),60);
end

methods = {'SCUBA'  'scuba_adaptor' '[6,4]' [1,1]; ...
           'SCUBA'  'scuba_adaptor' '[3,1]' [6,6]; ...
           'B-SCOTT' 'bscott_adaptor' '[6,6,4]' [1,1]; ...                 
           'B-SCOTT' 'bscott_adaptor' '[2,2,2]' [6,6]; ...                 
           'B-STEREO' 'stereo_blind' '3' [1,1]};               
[~,est] = compare_methods(SRI, HSI, MSI, struct('Pm', Pm), [4,4], methods);

err1 = zeros(1,size(SRI,3)); err2 = err1; err3 = err1; err4 = err1;
for k=1:size(SRI,3)
    err1(1,k) = norm(SRI(:,:,k)-est{1}(:,:,k),'fro')/norm(SRI(:,:,k),'fro');
    err2(1,k) = norm(SRI(:,:,k)-est{2}(:,:,k),'fro')/norm(SRI(:,:,k),'fro');
    err3(1,k) = norm(SRI(:,:,k)-est{3}(:,:,k),'fro')/norm(SRI(:,:,k),'fro');
    err4(1,k) = norm(SRI(:,:,k)-est{4}(:,:,k),'fro')/norm(SRI(:,:,k),'fro');
end

figure
subplot(1,2,1); plot(err1); hold on; plot(err3); xlabel('Spectral bands'); ylabel('RMSE'); title('(1,1) blocks'); legend('SCUBA','BSCOTT')
subplot(1,2,2); plot(err2); hold on; plot(err4); xlabel('Spectral bands'); ylabel('RMSE'); title('(6,6) blocks'); legend('SCUBA','BSCOTT')
