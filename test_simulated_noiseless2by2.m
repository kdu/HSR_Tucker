% EXPERIMENT 1: R-SNR FOR NOISELESS SIMULATED SRI 2x2x2 %
% Copyright (c) 2018 Clemence Prevost, Konstantin Usevich, Pierre Comon, David Brie
% https://github.com/cprevost4/HSR_Tucker
% Contact: clemence.prevost@univ-lorraine.fr

%% Generate data

I = 120; J = 120; %spatial dimensions
A1 = zeros(I,J); A2 = A1; %initialize abundance maps

g1 = fspecial('gaussian',60, 20);
A1(1:60,61:120) = g1; A1(61:120,1:60) = g1;
A2(1:60,1:60) = g1;

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

clear SRI; clear mat; clear ind;

X = outprod(A1,s1) + outprod(A2,s2);

Pm = spectral_deg(X,"LANDSAT");
MSI = tmprod(X,Pm,3);
d1 = 4; d2 = 4; q = 9;
[P1,P2] = spatial_deg(X, q, d1, d2);
HSI = tmprod(tmprod(X,P1,1),P2,2);

%% Run simulations

snr_stereo = []; snr_scott = [];
 for F=1:40
    F
    try
        [SRI_hat2, ~] = stereo3(HSI, MSI, P1,P2,Pm, F);
        SRI_hat2 = real(SRI_hat2); snr_stereo(F,1) = r_snr(X,SRI_hat2);
    catch
        snr_stereo(F,1) = NaN;
        continue
    end
 end

 for r1=1:size(HSI,1)+10
    for r3=1:10
        R = [r1,r1,r3]
        %try
            [SRI_hat2,~] = scott2(HSI, MSI, P1, P2, Pm, R);
            snr_scott(r1,r3) = r_snr(X,SRI_hat2);
        %catch
            %snr_scott(r1,r3) = NaN;
            %continue
        %end

        if not((r3<=size(MSI,3) || r1<=size(HSI,1)) && r3<=(min(r1,size(HSI,1)))^2)%identifiability
             snr_scott(r1,r3) = NaN;
        end
    end
end

%% Make figure

figure
plot(1:40,snr_stereo); hold on; plot(1:40, snr_scott(:,2))
legend('STEREO','SCOTT')
xlabel('F'); ylabel('SNR (dB)');...
    xlim([1 40]); title('STEREO and SCOTT (R_3=N)')
figure
surf(1:10,1:size(HSI)+10,snr_scott); xlabel('R_3');...
    ylabel('R_1 = R_2'); title('SCOTT')

%% 

% figure
% subplot(1,2,1); imagesc(X(:,:,44)); title('Spectral band 44')
% subplot(1,2,2); imagesc(X(:,:,160)); title('Spectral band 160')

%% ADD CURVE TO THE STEREO PLOT

r3=2;
 for r1=1:size(HSI,1)+10
        R = [r1,r1,r3]
        %try
            [SRI_hat2,~] = scott2(HSI, MSI, P1, P2, Pm, R);
            snr_scott(r1,r3) = r_snr(X,SRI_hat2);
        %catch
            %snr_scott(r1,r3) = NaN;
            %continue
        %end

        if not((r3<=size(MSI,3) || r1<=size(HSI,1)) && r3<=(min(r1,size(HSI,1)))^2)%identifiability
             snr_scott(r1,r3) = NaN;
        end
 end




