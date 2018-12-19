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
MSI_true = tmprod(X,Pm,3); MSI=MSI_true;
d1 = 4; d2 = 4; q = 9;
[P1,P2] = spatial_deg(X, q, d1, d2);
HSI_true = tmprod(tmprod(X,P1,1),P2,2); HSI=HSI_true;

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
MSI_true = tmprod(X,Pm,3); MSI=MSI_true;
d1 = 4; d2 = 4; q = 9;
[P1,P2] = spatial_deg(X, q, d1, d2);
HSI_true = tmprod(tmprod(X,P1,1),P2,2); HSI=HSI_true;

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
MSI_true = tmprod(X,Pm,3); MSI=MSI_true;
d1 = 4; d2 = 4; q = 9;
[P1,P2] = spatial_deg(X, q, d1, d2);
HSI_true = tmprod(tmprod(X,P1,1),P2,2); HSI=HSI_true;

%% 4rd EXAMPLE ON SIMULATED DATA

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
%         MSI(:,:,k) = awgn(MSI_true(:,:,k),20);
%     end
%     for k=1:size(HSI_true,3)
%         HSI(:,:,k) = awgn(HSI_true(:,:,k),20);
%     end

%% RUN SCOTT AUTOMATICALLY

lvl = ["20","35","60"];
%param = [2,6; 3,2; 4,2; 6,4];
param = [6,4; 6,6];
Niter = 10; lambda = 1; 

for i=1:4
    for j=1:3
        snr_stereo = []; snr_scott = [];
        eval(sprintf('load(''data_%dS%dB%sdB.mat'')',param(i,1),param(i,2),lvl(j)))
        
        for F=1:25
            F
            [SRI_hat2, ~] = stereo3(HSI, MSI, P1,P2,Pm, F);
            SRI_hat2 = real(SRI_hat2); snr_stereo(F,1) = r_snr(X,SRI_hat2);
        end
        
        
        figure
        plot(1:25,snr_stereo)
        xlabel('F'); ylabel('SNR (dB)')
        title(sprintf('STEREO, %d sources, Km=%d, %sdB SNR',param(i,1),param(i,2),lvl(j)))
        savefig(gcf, sprintf('SNR_CP%dM%sdB',param(i,1),lvl(j)))
    
        for r1=1:size(HSI,1)+10
            for r3=1:10
                R = [r1,r1,r3]
                [SRI_hat2,~] = scott(HSI, MSI, P1, P2, Pm, R);
                snr_scott(r1,r3) = r_snr(X,SRI_hat2);
                
                if (r3>size(MSI,3) && r1>size(HSI,1))%identifiability
                     snr_scott(r1,r3) = NaN;
                end
            end
        end
        
        figure
        surfc(1:10,1:size(HSI,1)+10,snr_scott)
        xlabel('R_3'); ylabel('R_1 = R_2'); zlabel('SNR (dB)')
        title(sprintf('SCOTT, %d sources, Km=%d, %sdB SNR',param(i,1),param(i,2),lvl(j)))
        savefig(gcf, sprintf('SNR_SVD%dM%sdB',param(i,1),lvl(j)))
        
           
    end
end

%% CASE WHERE SNR = Inf

param = [2,6; 3,2; 4,2; 6,4];

for i=1:4
    snr_stereo = []; snr_scott = [];
    eval(sprintf('load(''data_%dS%dBInfdB.mat'')',param(i,1),param(i,2)))
    

        for F=1:25
            F
            try
                [SRI_hat2, ~] = stereo3(HSI, MSI, P1,P2,Pm, F);
                SRI_hat2 = real(SRI_hat2); snr_stereo(F,1) = r_snr(X,SRI_hat2);
            catch
                snr_stereo(F,1) = NaN;
                continue
            end
        end
        
        figure
        plot(1:25,snr_stereo)
        xlabel('F'); ylabel('SNR (dB)'); xlim([1 50])
        title(sprintf('STEREO, %d sources, Km=%d, infinite SNR',param(i,1),param(i,2)))
        savefig(gcf, sprintf('SNR_CP%dMInfdB',param(i,1)))
        
        for r1=1:size(HSI,1)+10
            for r3=1:10
                R = [r1,r1,r3]
                try
                    [SRI_hat2,~] = scott(HSI, MSI, P1, P2, Pm, R);
                    snr_scott(r1,r3) = r_snr(X,SRI_hat2);
                catch
                    snr_scott(r1,r3) = NaN;
                    continue
                end
                
                if (r3>size(MSI,3) && r1>size(HSI,1))%identifiability
                     snr_scott(r1,r3) = NaN;
                end
            end
        end
        
        figure
        surfc(1:10,1:size(HSI)+10,snr_scott)
        xlabel('R_3'); ylabel('R_1 = R_2'); zlabel('SNR (dB)')
        title(sprintf('SCOTT, %d sources, Km=%d, infinite SNR',param(i,1),param(i,2)))
        savefig(gcf, sprintf('SNR_SVD%dMInfdB',param(i,1)))



end

%% PLOT LOG OF SINGULAR VALUES

lvl = ["20","35","60","Inf"];
param = [2,6; 3,2; 4,2; 6,4];

for i=1:4
    figure(i)
    for j=1:4
        eval(sprintf('load(''data_%dS%dB%sdB.mat'')',param(i,1),param(i,2),lvl(j)))
        
        
        subplot(4,3,3*j-2)
            x = log(svd(tens2mat(MSI,1,[])));
            ylim auto
            plot(x(1:30),'-s','MarkerFaceColor','red','MarkerSize',2,'MarkerEdgeColor','red')
            title(sprintf('1st mode MSI, %d mat., Km=%d, %sdB SNR',param(i,1),param(i,2),lvl(j)))
        subplot(4,3,3*j-1)
            x = log(svd(tens2mat(MSI,2,[])));
            ylim auto
            plot(x(1:30),'-s','MarkerFaceColor','red','MarkerSize',2,'MarkerEdgeColor','red')
            title(sprintf('2nd mode MSI, %d mat., Km=%d, %sdB SNR',param(i,1),param(i,2),lvl(j)))
        subplot(4,3,3*j)
            x = log(svd(tens2mat(HSI,3,[])));
            ylim auto
            plot(x(1:30),'-s','MarkerFaceColor','red','MarkerSize',2,'MarkerEdgeColor','red')
            title(sprintf('3rd mode HSI, %d mat., Km=%d, %sdB SNR',param(i,1),param(i,2),lvl(j)))
            
            %sgtitle(sprintf('Log of singular values, %d materials',param(i,1)))
        
      
    end
end

%% PLOT LOG OF SINGULAR VALUES 2

lvl = ["20","35","60","Inf"];
param = [2,6; 3,2; 4,2; 6,4];
leg = {'20dB SNR', '35dB SNR','60dB SNR','Inf. SNR'};
newPosition = [0.4 0.4 0.2 0.2];

for i=1:4
    figure(i)
    sgtitle(sprintf('Log of singular values, %d materials, Km=%d',param(i,1),param(i,2)))
    for j=1:4
        eval(sprintf('load(''data_%dS%dB%sdB.mat'')',param(i,1),param(i,2),lvl(j)))
        
        
       subplot(1,3,1)
        x = log(svd(tens2mat(MSI,1,[])));
        ylim auto
        plot(x(1:10),'-s','MarkerSize',3,'Linewidth',1.1)
        title('1st unfolding MSI')
        legend(leg(1:j),'Location', 'southoutside')
        hold on
       subplot(1,3,2)
        x = log(svd(tens2mat(MSI,2,[])));
        ylim auto
        plot(x(1:10),'-s','MarkerSize',3,'Linewidth',1.1)
        title('2nd unfolding MSI')
        legend(leg(1:j),'Location', 'southoutside')
        hold on
       subplot(1,3,3)
        x = log(svd(tens2mat(HSI,3,[])));
        ylim auto
        plot(x(1:10),'-s','MarkerSize',3,'Linewidth',1.1)
        title('3rd unfolding HSI')
        legend(leg(1:j),'Location', 'southoutside')
        hold on
    end
end

%% BLIND ALGORITHMS

lvl = ["20","35","60"];
param = [2,6; 3,2; 4,2; 6,4];
Niter = 10; lambda = 1; opts.Nblocks = [4,4];

for i=1:4
    for j=1:3
        snr_stereo_blind = []; snr_scott_blind = [];
        eval(sprintf('load(''data_%dS%dB%sdB.mat'')',param(i,1),param(i,2),lvl(j)))
        
        for F=1:30
            F
            [SRI_hat2, ~] = stereo_blind(HSI, MSI,Pm, F);
            SRI_hat2 = real(SRI_hat2); snr_stereo_blind(F,1) = r_snr(X,SRI_hat2);
        end
        
        figure
        plot(1:30,snr_stereo_blind)
        xlabel('F'); ylabel('SNR (dB)')
        title(sprintf('BLIND STEREO, %d sources, Km=%d, %sdB noise',param(i,1),param(i,2),lvl(j)))
        savefig(gcf, sprintf('SNR_BCP%dM%sdB',param(i,1),lvl(j)))
    
        for r1=1:size(HSI,1)
            for r3=1:size(MSI,3)
                R = [r1,r1,r3]
                [SRI_hat2,~] = bscott(MSI, HSI, Pm,R,opts);
                snr_scott_blind(r1,r3) = r_snr(X,SRI_hat2);
            end
        end
        
        figure
        surfc(1:size(MSI,3),1:size(HSI,1),snr_scott_blind)
        xlabel('R_3'); ylabel('R_1 = R_2'); zlabel('SNR (dB)')
        title(sprintf('BSCOTT, %d sources, Km=%d, %sdB noise',param(i,1),param(i,2),lvl(j)))
        savefig(gcf, sprintf('SNR_BSVD%dM%sdB',param(i,1),lvl(j)))
        
        
    end
end

%% BLIND CASE WHERE SNR = Inf

param = [2,6; 3,2; 4,2; 6,4]; opts.Nblocks = [4,4];

for i=1:4
    snr_stereo_blind = []; snr_scott_blind = [];
    eval(sprintf('load(''data_%dS%dBInfdB.mat'')',param(i,1),param(i,2)))
    

        for F=1:30
            F
            try
            [SRI_hat2, ~] = stereo_blind(HSI, MSI,Pm, F);
            SRI_hat2 = real(SRI_hat2); snr_stereo_blind(F,1) = r_snr(X,SRI_hat2);
            catch
                snr_stereo_blind(F,1) = NaN;
                continue
            end
        end
        
        figure
        plot(1:30,snr_stereo_blind)
        xlabel('F'); ylabel('SNR (dB)'); xlim([1 30])
        title(sprintf('BLIND STEREO, %d sources, Km=%d, Inf dB noise',param(i,1),param(i,2)))
        savefig(gcf, sprintf('SNR_BCP%dMInfdB',param(i,1)))
        
        for r1=1:size(HSI,1)
            for r3=1:size(MSI,3)
                R = [r1,r1,r3]
                try
                [SRI_hat2,~] = bscott(MSI, HSI, Pm,R,opts);
                snr_scott_blind(r1,r3) = r_snr(X,SRI_hat2);
                catch
                    snr_scott_blind(r1,r3) = NaN;
                    continue
                end
            end
        end
        
        figure
        surfc(1:size(MSI,3),1:size(HSI),snr_scott_blind)
        xlabel('R_3'); ylabel('R_1 = R_2'); zlabel('SNR (dB)')
        title(sprintf('BSCOTT, %d sources, Km=%d, Inf dB noise',param(i,1),param(i,2)))
        savefig(gcf, sprintf('SNR_BSVD%dMInfdB',param(i,1)))



end

%% MAKE DOUBLE FIGURES

% Load saved figures
a=hgload('SNR_CP6MInfdB.fig');
b=hgload('SNR_SVD6MInfdB.fig');
% c=hgload('SNR_CP2M35dB.fig');
% d=hgload('SNR_SVD2M35dB.fig');
% e=hgload('SNR_CP2M60dB.fig');
% f=hgload('SNR_SVD2M60dB.fig');
% Prepare subplots
figure
%sgtitle('SNR between groundtruth and estimate, 2 mat., Km=6, SNR=20dB')
h(1)=subplot(1,2,1); title('STEREO'); xlabel('F')
h(2)=subplot(1,2,2); title('SCOTT'); xlabel('R_3'); ylabel('R_1=R_2')
% h(3)=subplot(3,2,3);
% h(4)=subplot(3,2,4);
% h(5)=subplot(3,2,5);
% h(6)=subplot(3,2,6);
% Paste figures on the subplots
copyobj(allchild(get(a,'CurrentAxes')),h(1));
copyobj(allchild(get(b,'CurrentAxes')),h(2));
% copyobj(allchild(get(c,'CurrentAxes')),h(3));
% copyobj(allchild(get(d,'CurrentAxes')),h(4));
% copyobj(allchild(get(e,'CurrentAxes')),h(5));
% copyobj(allchild(get(f,'CurrentAxes')),h(6));

%% MAKE BOXPLOTS

lvl = ["20","35","60"];
%param = [2,6; 3,2; 4,2; 6,4];
param = [6,4;];
Niter = 10; lambda = 1; Nreal = 10;

for i=1:3
    for j=1:3
        
        eval(sprintf('load(''rawdata_%dS%dB.mat'')',param(i,1),param(i,2)))
        MSI_true = MSI; HSI_true = HSI;
        R3 = size(MSI_true,3);%R3 = [param(i,1) param(i,2)];
        
        for n = 1:Nreal
            
            % Generate one realization of noise
            for k=1:size(MSI_true,3)
                eval(sprintf('MSI(:,:,k) = awgn(MSI_true(:,:,k),%s);',lvl(j)))
            end
            for k=1:size(HSI_true,3)
                eval(sprintf('HSI(:,:,k) = awgn(HSI_true(:,:,k),%s);',lvl(j)))
            end
            
            % Run STEREO
            for F=1:25
                F
                [SRI_hat2, ~] = stereo3(HSI, MSI, P1,P2,Pm, F);
                SRI_hat2 = real(SRI_hat2); snr_stereo(n,F) = r_snr(X,SRI_hat2);
            end
            
            
            for r1=1:size(HSI,1)+10
                R = [r1,r1,R3]
                [SRI_hat2,~] = scott(HSI, MSI, P1, P2, Pm, R);
                snr_scott(n,r1) = r_snr(X,SRI_hat2);

                if (R3>size(MSI,3) && r1>size(HSI,1))%identifiability
                     snr_scott(n,r1) = NaN;
                end
            end          
        end
        
        figure
        boxplot(snr_stereo)
        xlabel('F');
        title(sprintf('STEREO, %d sources, Km=%d, %sdB SNR',param(i,1),param(i,2),lvl(j)))
        savefig(gcf, sprintf('BP_CP%dM%sdB',param(i,1),lvl(j)))
        
        figure
        boxplot(snr_scott)
        xlabel('R_1 = R_2');
        title(sprintf('SCOTT, %d sources, Km=%d, %sdB SNR, R_3 = %d',param(i,1),param(i,2),lvl(j),R3))
        savefig(gcf, sprintf('BP_SVD%dM%sdB%dRs',param(i,1),lvl(j),R3))

           
    end
end

%% FIGURES OF SPECTRAL BANDS 

% FIND ORIGINAL SPECTRUM FOR MATERIAL 1
SRI = cell2mat(struct2cell(load('Indian_pines.mat')));
SRI(:,:,[104:108 150:163 220]) = []; %Regions of water absorption
SRI(1,:,:) = []; SRI(:,1,:) = [];
load('Indian_pines_gt.mat'); indian_pines_gt(:,1) = []; indian_pines_gt(1,:) = [];
mat = tens2mat(SRI,3,[]);
ind = find(reshape(indian_pines_gt,1,[]) == 3);
s1 = real(mean(mat(:,ind),2)); 

% RUN ALGORITHMS
load('data_6S4BInfdB.mat')
F = 25; R = [6,6,6];
[SRI_hat1, ~] = stereo3(HSI, MSI, P1,P2,Pm, F);
SRI_hat1 = real(SRI_hat1);
r_snr(X,SRI_hat1)
[SRI_hat2, ~] = scott(HSI, MSI, P1, P2, Pm, R);
r_snr(X,SRI_hat2)

err1 = zeros(1,200); err2 = err1;
for k=1:size(X,3)
    err1(1,k) = norm(X(1:40,1:40,k)-SRI_hat1(1:40,1:40,k),'fro');
    err2(1,k) = norm(X(1:40,1:40,k)-SRI_hat2(1:40,1:40,k),'fro');
end
    
ind = find(err1 == max(err1));
% find(err2 == max(err2))

% MAKE PLOTS
figure(1)
subplot(1,3,1)
    imagesc(X(1:40,1:40,ind))
    title('Groundtruth SRI')
subplot(1,3,2)
    imagesc(SRI_hat1(1:40,1:40,ind))
    title('STEREO')
subplot(1,3,3)
    imagesc(SRI_hat2(1:40,1:40,ind))
    title('SCOTT')

figure(2)
plot(err1)
hold on
plot(err2)
% hold on 
% plot(s1'/1000)
    









