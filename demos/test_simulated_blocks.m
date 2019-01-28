% EXPERIMENT 1: HOSVD FOR VARIOUS RANKS %
% Copyright (c) 2018 Clemence Prevost, Konstantin Usevich, Pierre Comon, David Brie
% https://github.com/cprevost4/HSR_Tucker
% Contact: clemence.prevost@univ-lorraine.fr

%% Load 2-material dataset

% 1. abundance maps
I = 120; J = 120; %spatial dimensions
A1 = zeros(I,J); A2 = A1; %initialize abundance maps
g1 = fspecial('gaussian',60, 20);
A1(1:60,61:120) = g1; A1(61:120,1:60) = g1; A2(1:60,1:60) = g1;
% 2. spectra
K = 200; %spectral dimension
SRI = cell2mat(struct2cell(load('Indian_pines.mat')));
SRI(:,:,[104:108 150:163 220]) = []; %Regions of water absorption
SRI(1,:,:) = []; SRI(:,1,:) = [];
load('Indian_pines_gt.mat'); indian_pines_gt(:,1) = []; indian_pines_gt(1,:) = [];
mat = tens2mat(SRI,3,[]); 
ind = find(reshape(indian_pines_gt,1,[]) == 3); s1 = real(mean(mat(:,ind),2)); 
ind = find(reshape(indian_pines_gt,1,[]) == 5); s2 = real(mean(mat(:,ind),2)); 
clear SRI; clear mat; clear ind;
% 3. generate images
SRI = outprod(A1,s1) + outprod(A2,s2);
Pm = spectral_deg(SRI,"LANDSAT"); MSI_gt = tmprod(SRI,Pm,3);
d1 = 4; d2 = 4; q = 9;
[P1,P2] = spatial_deg(SRI, q, d1, d2);
HSI_gt = tmprod(tmprod(SRI,P1,1),P2,2);

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

%% R-SNR for Nblocks - 2 materials

lvl = ["20","35","60"]; ranks1 = [3,2]; ranks2 = [2,2,2];
for j=1:length(lvl)
    
    snr1 = []; snr2 = [];
    for k=1:size(MSI_gt,3)
        eval(sprintf('MSI(:,:,k) = awgn(MSI_gt(:,:,k),%s);',lvl(j)))
    end
    for k=1:size(HSI_gt,3)
        eval(sprintf('HSI(:,:,k) = awgn(HSI_gt(:,:,k),%s);',lvl(j)))
    end
    
    for i=1:12

            [i,i]
            opts.Nblocks = [i,i];
            [SRI_hat1,~] = scuba(MSI, HSI, Pm,ranks1,opts);
            snr1(i,1) = r_snr(SRI,SRI_hat1);
            [SRI_hat2,~] = bscott(MSI, HSI, Pm,ranks2,opts);
            snr2(i,1) = r_snr(SRI,SRI_hat2);

    end
    
    figure(1)
    subplot(1,3,j)
    plot(snr1); hold on; plot(snr2); legend('SCUBA','BSCOTT')
    
end

%% R-SNR for Nblocks - 6 materials

lvl = ["20","35","60"]; ranks1 = [7,4]; ranks2 = [6,6,4];
for j=1:length(lvl)
    
    snr1 = []; snr2 = [];
    for k=1:size(MSI_gt,3)
        eval(sprintf('MSI(:,:,k) = awgn(MSI_gt(:,:,k),%s);',lvl(j)))
    end
    for k=1:size(HSI_gt,3)
        eval(sprintf('HSI(:,:,k) = awgn(HSI_gt(:,:,k),%s);',lvl(j)))
    end
    
    for i=1:12

            [i,i]
            opts.Nblocks = [i,i];
            [SRI_hat1,~] = scuba(MSI, HSI, Pm,ranks1,opts);
            snr1(i,1) = r_snr(SRI,SRI_hat1);
            [SRI_hat2,~] = bscott(MSI, HSI, Pm,ranks2,opts);
            snr2(i,1) = r_snr(SRI,SRI_hat2);

    end
    
    figure(1)
    subplot(1,3,j)
    plot(snr1); hold on; plot(snr2); legend('SCUBA','BSCOTT')
    
end