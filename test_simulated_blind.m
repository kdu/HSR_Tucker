% EXPERIMENT 1: HOSVD FOR VARIOUS RANKS %
% Copyright (c) 2018 Clemence Prevost, Konstantin Usevich, Pierre Comon, David Brie
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

%% Experiments - 6 materials

lvl = ["20","60"]; opts.Nblocks = [1,1]; opts2.Nblocks = [6,6];
snr_scuba = []; snr_scuba2 = []; snr_bscott = []; snr_bscott2 = [];

for j=1:length(lvl)

    for k=1:size(MSI_gt,3)
        eval(sprintf('MSI(:,:,k) = awgn(MSI_gt(:,:,k),%s);',lvl(j)))
    end
    for k=1:size(HSI_gt,3)
        eval(sprintf('HSI(:,:,k) = awgn(HSI_gt(:,:,k),%s);',lvl(j)))
    end
    
     for F=1:25
        for r3=1:size(MSI,3)
            [F r3]
            [SRI_hat, ~] = scuba(MSI, HSI,Pm, [F,r3], opts); 
            [SRI_hat1, ~] = scuba(MSI, HSI,Pm, [F,r3], opts2);
            snr_scuba(F,r3) = r_snr(SRI,SRI_hat); snr_scuba2(F,r3) = r_snr(SRI,SRI_hat1);
        end
     end
     
    figure
    surf(1:size(MSI,3),1:25,snr_scuba,'FaceAlpha',0.6); colormap summer; freezeColors; hold on
    surf(1:size(MSI,3),1:25,snr_scuba2,'FaceAlpha',0.6); colormap spring
    xlabel('R_3'); ylabel('F'); zlabel('R-SNR (dB)'); legend('(1,1) blocks','(6,6) blocks')
    title(sprintf('SCUBA, 6 sources, Km=4, %s dB noise',lvl(j)))
    
     for r1=1:40
        for r3=1:size(MSI,3)
            R = [r1,r1,r3]
            [SRI_hat2,~] = bscott(MSI, HSI, Pm,R,opts);
            [SRI_hat3,~] = bscott(MSI, HSI, Pm,R,opts2);
            snr_bscott(r1,r3) = r_snr(SRI,SRI_hat2);
            snr_bscott2(r1,r3) = r_snr(SRI,SRI_hat3);
            if not((r3<=size(MSI,3) || r1<=size(HSI,1)) && (r1<=min(r3,size(MSI,3))*r1 && r3<=min(r1,size(HSI,1))^2))
                snr_bscott(r1,r3) = NaN; snr_bscott2(r1,r3) = NaN;
            end  
        end
     end
    
        figure
        surf(1:size(MSI,3),1:40,snr_bscott, 'FaceAlpha', 0.6); colormap summer; freezeColors; hold on
        surf(1:size(MSI,3),1:40,snr_bscott2, 'FaceAlpha',0.6); colormap spring
        xlabel('R_3'); ylabel('R_1 = R_2'); zlabel('SNR (dB)'); legend('(1,1) blocks','(6,6) blocks')
        title(sprintf('BSCOTT, 6 sources, Km=4, %sdB noise',lvl(j)))
    
end