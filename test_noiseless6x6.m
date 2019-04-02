% EXPERIMENT 1: R-SNR FOR NOISELESS SIMULATED SRI 6x6x4 %
% Copyright (c) 2018 Clemence Prevost, Konstantin Usevich, Pierre Comon, David Brie
% https://github.com/cprevost4/HSR_Tucker
% Contact: clemence.prevost@univ-lorraine.fr

%% Generate data

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

X = outprod(A1,s1) + outprod(A2,s2) + outprod(A3,s3) + outprod(A4,s4) + outprod(A5,s5) + outprod(A6,s6);
Pm = spectral_deg(X,"Quickbird");
%Pm = spectral_deg(X,"LANDSAT");
MSI_true = tmprod(X,Pm,3); MSI=MSI_true;
d1 = 4; d2 = 4; q = 9;
[P1,P2] = spatial_deg(X, q, d1, d2);
HSI_true = tmprod(tmprod(X,P1,1),P2,2); HSI=HSI_true;


%     for k=1:4
%         MSI(:,:,k) = awgn(MSI_true(:,:,k),35);
%     end
%     for k=1:size(HSI_true,3)
%         HSI(:,:,k) = awgn(HSI_true(:,:,k),35);
%     end

%% Run simulations

snr_stereo = []; snr_scott = [];
 for F=1:40%70
    F
    try
        [SRI_hat2, ~] = stereo3(HSI, MSI, P1,P2,Pm, F);
        SRI_hat2 = real(SRI_hat2); snr_stereo(F,1) = r_snr(X,SRI_hat2);
    catch
        snr_stereo(F,1) = NaN;
        continue
    end
 end

 for r1=1:40%size(HSI,1)+10
    for r3=1:15
        R = [r1,r1,r3]

        if ((r3<=size(MSI,3) || r1<=size(HSI,1)) && r3<=(min(r1,size(HSI,1)))^2)%identifiability
            [SRI_hat2,~] = scott(HSI, MSI, P1, P2, Pm, R);
            snr_scott(r1,r3) = r_snr(X,SRI_hat2);
        else
             snr_scott(r1,r3) = NaN;
        end
    end
end

%% Make figure

figure
subplot(1,2,1); plot(1:40,snr_stereo) ;xlabel('F'); ylabel('SNR (dB)');...
    xlim([1 40]); title('STEREO')
subplot(1,2,2); surf(1:15,1:40,snr_scott); xlabel('R_3');...
    ylabel('R_1 = R_2'); title('SCOTT')

%% 

% figure
% subplot(1,2,1); imagesc(X(:,:,44)); title('Spectral band 44')
% subplot(1,2,2); imagesc(X(:,:,160)); title('Spectral band 160')

%% ADD CURVE TO THE STEREO PLOT
% 
% r3=6;
%  for r1=1:size(HSI,1)+10
%         R = [r1,r1,r3]
%         %try
%             [SRI_hat2,~] = scott2(HSI, MSI, P1, P2, Pm, R);
%             snr_scott(r1,r3) = r_snr(X,SRI_hat2);
%         %catch
%             %snr_scott(r1,r3) = NaN;
%             %continue
%         %end
% 
%         if not((r3<=size(MSI,3) || r1<=size(HSI,1)) && r3<=(min(r1,size(HSI,1)))^2)%identifiability
%              snr_scott(r1,r3) = NaN;
%         end
%  end

