%% 1st EXAMPLE ON SIMULATED DATA %

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
MSI_true = tmprod(X,Pm,3);
d1 = 4; d2 = 4; q = 9;
[P1,P2] = spatial_deg(X, q, d1, d2);
HSI_true = tmprod(tmprod(X,P1,1),P2,2);

%% 2nd EXAMPLE ON SIMULATED DATA

I = 180; J = 180; %spatial dimensions
A1 = zeros(I,J); A2 = A1; A3 = A1; %initialize abundance maps

g1 = fspecial('gaussian',60, 20);
A1(1:60,1:60) = g1; 
A2(61:120,1:60) = g1; A2(1:60,61:120) = g1;
A3(121:180,1:60) = g1; A3(61:120,61:120) = g1; A3(1:60,121:180) = g1; 


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
ind = find(reshape(indian_pines_gt,1,[]) == 12);
s3 = real(mean(mat(:,ind),2)); 

clear SRI; clear mat; clear ind;

X = outprod(A1,s1) + outprod(A2,s2) + outprod(A3,s3);

Pm = zeros(2,size(X,3));
Pm(1,1:100) = (1/100)*ones(1,100); Pm(2,101:200) = (1/100)*ones(1,100);
MSI_true = tmprod(X,Pm,3);
d1 = 4; d2 = 4; q = 9;
[P1,P2] = spatial_deg(X, q, d1, d2);
HSI_true = tmprod(tmprod(X,P1,1),P2,2);


%% 3rd EXAMPLE ON SIMULATED DATA

I = 240; J = 240; %spatial dimensions
A1 = zeros(I,J); A2 = A1; A3 = A1; A4 = A1; %initialize abundance maps

g1 = fspecial('gaussian',60, 20);
A1(1:60,1:60) = g1; 
A2(61:120,1:60) = g1; A2(1:60,61:120) = g1;
A3(121:180,1:60) = g1; A3(1:60,121:180) = g1; 
A4(181:240,1:60) = g1; A4(121:180,61:120) = g1; A4(61:120,121:180) = g1; A4(1:60,181:240) = g1;  


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
ind = find(reshape(indian_pines_gt,1,[]) == 12);
s3 = real(mean(mat(:,ind),2)); 
ind = find(reshape(indian_pines_gt,1,[]) == 15);
s4 = real(mean(mat(:,ind),2)); 

clear SRI; clear mat; clear ind;

X = outprod(A1,s1) + outprod(A2,s2) + outprod(A3,s3) + outprod(A4,s4);

Pm = zeros(2,size(X,3));
Pm(1,1:100) = (1/100)*ones(1,100); Pm(2,101:200) = (1/100)*ones(1,100);
MSI_true = tmprod(X,Pm,3);
d1 = 4; d2 = 4; q = 9;
[P1,P2] = spatial_deg(X, q, d1, d2);
HSI_true = tmprod(tmprod(X,P1,1),P2,2);


%% 4rd EXAMPLE ON SIMULATED DATA

I = 360; J = 360; %spatial dimensions
A1 = zeros(I,J); A2 = A1; A3 = A1; A4 = A1; 
A5 = A1; A6 = A1; %initialize abundance maps

g1 = fspecial('gaussian',60, 20);
A1(1:60,1:60) = g1; 
A2(61:120,1:60) = g1; A2(1:60,61:120) = g1;
A3(121:180,1:60) = g1; A3(61:120,61:120) = g1; A3(1:60,121:180) = g1; 
A4(181:240,1:60) = g1; A4(1:60,181:240) = g1; 
A5(241:300,1:60) = g1; A5(181:240,61:120) = g1; A5(61:120,181:240) = g1; A5(1:60,241:300) = g1;
A6(301:360,1:60) = g1; A6(241:300,61:120) = g1; A6(181:240, 121:180) = g1;
    A6(121:180,181:240) = g1; A6(61:120,241:300) = g1; A6(1:60,301:360) = g1;
    

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
% Pm = zeros(2,size(X,3));
% Pm(1,1:100) = (1/100)*ones(1,100); Pm(2,101:200) = (1/100)*ones(1,100);
% Pm = zeros(4,size(X,3));
% Pm(1,1:50) = (1/50)*ones(1,50); Pm(2,51:100) = (1/50)*ones(1,50);
% Pm(3,101:150) = (1/50)*ones(1,50); Pm(4,151:200) = (1/50)*ones(1,50);
Pm = spectral_deg(X,"LANDSAT");
MSI_true = tmprod(X,Pm,3);
d1 = 4; d2 = 4; q = 9;
[P1,P2] = spatial_deg(X, q, d1, d2);
HSI_true = tmprod(tmprod(X,P1,1),P2,2);


%%
snr_stereo = []; %MSI = MSI_true; HSI = HSI_true;
for n=1:1
    
    for k=1:size(MSI_true,3)
        MSI(:,:,k) = awgn(MSI_true(:,:,k),60);
    end
    for k=1:size(HSI_true,3)
        HSI(:,:,k) = awgn(HSI_true(:,:,k),60);
    end
    
    Niter = 10; lambda = 1; 
    for F=1:50
        F
        [SRI_hat, ~] = stereo3( HSI, MSI, P1,P2,Pm, F);
        SRI_hat = real(SRI_hat); snr_stereo(F,n) = r_snr(X,SRI_hat);
    end
    
    for i=1:70
        for j=1:10
            R = [i,i,j]
            [SRI_hat,~] = scott(HSI, MSI, P1, P2, Pm, R);
            snr_scott(i,j) = r_snr(X,SRI_hat);
        end
    end

end

%%

snr_stereo = []; %MSI = MSI_true; HSI = HSI_true;
for n=1:10
    
    for k=1:size(MSI_true,3)
        MSI(:,:,k) = awgn(MSI_true(:,:,k),Inf);
    end
    for k=1:size(HSI_true,3)
        HSI(:,:,k) = awgn(HSI_true(:,:,k),Inf);
    end
    
    Niter = 10; lambda = 1; 
    for F=1:50
        F
        [SRI_hat, ~] = stereo3( HSI, MSI, P1,P2,Pm, F);
        SRI_hat = real(SRI_hat); snr_stereo(F,n) = r_snr(X,SRI_hat);
    end
    
    for i=1:100
        if i<=size(HSI,1)
            R = [i,i,size(MSI,3)]
            [SRI_hat,~] = scott(HSI, MSI, P1, P2, Pm, R);
            snr_scott(i,n) = r_snr(X,SRI_hat);
        else
            snr_scott(i,n) = NaN;
        end
    end

end


 %% MAKE PLOTS
 
% figure(1)
% for i=1:size(snr_stereo,2)
%     plot(snr_stereo(:,i))
%     hold on
% end
% xlabel('F'); ylabel('SNR (dB)')
% figure(2)
% surfc(1:10,1:68,snr_scott)
% xlabel('R_3'); ylabel('R_1 = R_2'); zlabel('SNR (dB)')
figure(1)
boxplot(snr_stereo')
figure(2)
boxplot(snr_scott(1:size(HSI,1),:)')


%% RUN SCOTT AUTOMATICALLY

lvl = ["20","35","60"];
%param = [2,6; 3,2; 4,2; 6,4];
param = [3,2; 4,2; 6,4];
Niter = 10; lambda = 1; 

for i=1:4
    for j=1:3
        snr_stereo = []; snr_scott = [];
        eval(sprintf('load(''data_%dS%dB%sdB.mat'')',param(i,1),param(i,2),lvl(j)))
        
        for F=1:50
            F
            [SRI_hat, ~] = stereo3(HSI, MSI, P1,P2,Pm, F);
            SRI_hat = real(SRI_hat); snr_stereo(F,1) = r_snr(X,SRI_hat);
        end
        
        figure
        plot(1:50,snr_stereo)
        xlabel('F'); ylabel('SNR (dB)')
        title(sprintf('STEREO, %d sources, Km=%d, %sdB noise',param(i,1),param(i,2),lvl(j)))
        savefig(gcf, sprintf('SNR_CP%dM%sdB',param(i,1),lvl(j)))
    
        for r1=1:size(HSI,1)
            for r3=1:10
                R = [r1,r1,r3]
                [SRI_hat,~] = scott(HSI, MSI, P1, P2, Pm, R);
                snr_scott(r1,r3) = r_snr(X,SRI_hat);
            end
        end
        
        figure
        surfc(1:10,1:79,snr_scott)
        xlabel('R_3'); ylabel('R_1 = R_2'); zlabel('SNR (dB)')
        title(sprintf('SCOTT, %d sources, Km=%d, %sdB noise',param(i,1),param(i,2),lvl(j)))
        savefig(gcf, sprintf('SNR_SVD%dM%sdB',param(i,1),lvl(j)))
        
           
    end
end


%% CASE WHERE SNR = Inf

param = [2,6; 3,2; 4,2; 6,4];

for i=1:4
    snr_stereo = []; snr_scott = [];
    eval(sprintf('load(''data_%dS%dBInfdB.mat'')',param(i,1),param(i,2)))
    

        for F=1:50
            F
            try
                [SRI_hat, ~] = stereo3(HSI, MSI, P1,P2,Pm, F);
                SRI_hat = real(SRI_hat); snr_stereo(F,1) = r_snr(X,SRI_hat);
            catch
                snr_stereo(F,1) = NaN;
                continue
            end
        end
        
        figure
        plot(1:50,snr_stereo)
        xlabel('F'); ylabel('SNR (dB)'); xlim([1 50])
        title(sprintf('STEREO, %d sources, Km=%d, Inf dB noise',param(i,1),param(i,2)))
        savefig(gcf, sprintf('SNR_CP%dMInfdB',param(i,1)))
        
        for r1=1:size(HSI,1)
            for r3=1:10
                R = [r1,r1,r3]
                try
                    [SRI_hat,~] = scott(HSI, MSI, P1, P2, Pm, R);
                    snr_scott(r1,r3) = r_snr(X,SRI_hat);
                catch
                    snr_scott(r1,r3) = NaN;
                    continue
                end
            end
        end
        
        figure
        surfc(1:10,1:size(HSI),snr_scott)
        xlabel('R_3'); ylabel('R_1 = R_2'); zlabel('SNR (dB)')
        title(sprintf('SCOTT, %d sources, Km=%d, Inf dB noise',param(i,1),param(i,2)))
        savefig(gcf, sprintf('SNR_SVD%dMInfdB',param(i,1)))



end

