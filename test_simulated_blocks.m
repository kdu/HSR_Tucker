%% EXPERIMENT ON SIMULATED DATA %%
% CASE WHERE THE CP MODEL IS IDENTIFIABLE %

% Copyright (c) 2018 Clemence Prevost, Konstantin Usevich, Pierre Comon, David Brie
% https://github.com/cprevost4/HSR_Tucker
% Contact: clemence.prevost@univ-lorraine.fr

%% GENERATE SIMULATED DATA
 
% 1. abundance maps
I = 120; J = 120; g1 = fspecial('gaussian',10, 4); %spatial dimensions
A1 = zeros(I,J); A2 = A1; A3 = A1; A4 = A1; A5 = A1; A6 = A1; %init. maps

A1(1:10,1:10) = g1; A1(11:20,11:20) = g1; 
A2(21:30,21:30) = g1; A2(31:40,31:40) = g1;
A3(41:50,41:50) = g1; A3(51:60,51:60) = g1;
A4(61:70,61:70) = g1; A4(71:80,71:80) = g1; 
A5(81:90,81:90) = g1; A5(91:100,91:100) = g1; 
A6(101:110,101:110) = g1; A6(111:120,111:120) = g1;

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
Pm = spectral_deg(SRI,"Quickbird"); MSI = tmprod(SRI,Pm,3);
d1 = 4; d2 = 4; q = 9; [P1,P2] = spatial_deg(SRI, q, d1, d2);
HSI = tmprod(tmprod(SRI,P1,1),P2,2); 


%% RUN EXPERIMENTS

snr_stereo = []; snr_scott = [];
 for F=1:40
    F
    try
        [SRI_hat2, ~] = stereo3(HSI, MSI, P1,P2,Pm, F);
        SRI_hat2 = real(SRI_hat2); snr_stereo(F,1) = r_snr(SRI,SRI_hat2);
    catch
        snr_stereo(F,1) = NaN;
        continue
    end
 end

 for r1=1:size(HSI,1)+10
    for r3=1:15
        R = [r1,r1,r3]
        %try
            [SRI_hat2,~] = scott2(HSI, MSI, P1, P2, Pm, R);
            snr_scott(r1,r3) = r_snr(SRI,SRI_hat2);
        %catch
            %snr_scott(r1,r3) = NaN;
            %continue
        %end

        if not((r3<=size(MSI,3) || r1<=size(HSI,1)) && r3<=(min(r1,size(HSI,1)))^2)%identifiability
             snr_scott(r1,r3) = NaN;
        end
    end
 end
 
 %% MAKE FIGURES
 
 figure
subplot(1,2,1); plot(1:40,snr_stereo) ;xlabel('F'); ylabel('SNR (dB)');...
    xlim([1 40]); title('STEREO')
subplot(1,2,2); surf(1:15,1:size(HSI)+10,snr_scott); xlabel('R_3');...
    ylabel('R_1 = R_2'); title('SCOTT')

%% 

figure
subplot(1,2,1); imagesc(SRI(:,:,44)); title('Spectral band 44')
subplot(1,2,2); imagesc(SRI(:,:,160)); title('Spectral band 160')

%% TRY OUT SOME SHIT

X3 = tens2mat(SRI,3,[]);




 