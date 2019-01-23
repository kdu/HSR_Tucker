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
param = [3,2; 4,2; 6,4];
Niter = 10; lambda = 1; 

for i=1:4
    for j=1:3
        cost1 = [NaN*ones(param(i,1)-1,1)]; cost2 = [NaN*ones(size(HSI,1)+10,param(i,1)-1)];% snr_stereo = []; snr_scott = [];
        eval(sprintf('load(''data_%dS%dB%sdB.mat'')',param(i,1),param(i,2),lvl(j)))
        
%         for F=param(i,1):25
%             F
%             [~, info] = stereo3(HSI, MSI, P1,P2,Pm, F);
%             %SRI_hat2 = real(SRI_hat2); snr_stereo(F,1) = r_snr(X,SRI_hat2);
%             cost1(F,1) = cost_stereo(info.factors, HSI, MSI, lambda, P1,P2,Pm);
%             
%         end
%         
%         
%         
%         figure
% %         plot(1:25,snr_stereo)
% %         xlabel('F'); ylabel('SNR (dB)')
% %         title(sprintf('STEREO, %d sources, Km=%d, %sdB SNR',param(i,1),param(i,2),lvl(j)))
% %         savefig(gcf, sprintf('SNR_CP%dM%sdB',param(i,1),lvl(j)))
%         plot(1:25,cost1)
%         xlabel('F'); ylabel('Value of cost function')
%         title(sprintf('STEREO, %d sources, Km=%d, %sdB SNR',param(i,1),param(i,2),lvl(j)))
%         savefig(gcf, sprintf('cost_CP%dM%sdB',param(i,1),lvl(j)))
    
        for r1=1:size(HSI,1)+10
            for r3=param(i,1):10
                R = [r1,r1,r3]
                [~,info] = scott(HSI, MSI, P1, P2, Pm, R);
                %snr_scott(r1,r3) = r_snr(X,SRI_hat2);
                cost2(r1,r3) = cost_scott(info.factors,info.core, HSI, MSI, lambda, P1,P2,Pm);
                
                if (r3>size(MSI,3) && r1>size(HSI,1))%identifiability
                     cost2(r1,r3) = NaN; %snr_scott(r1,r3) = NaN;
                end
            end
        end
        
        cost2(1:param(i,1)-1,:) = NaN;
        figure
        surfc(1:10,1:size(HSI,1)+10,cost2)
        xlabel('R_3'); ylabel('R_1 = R_2'); zlabel('Value of cost function')
        title(sprintf('SCOTT, %d sources, Km=%d, %sdB SNR',param(i,1),param(i,2),lvl(j)))
        savefig(gcf, sprintf('cost_SVD%dM%sdB',param(i,1),lvl(j)))
        
           
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

lvl = "Inf"; %lvl = ["20","35","60"];
%param = [2,6; 3,2; 4,2; 6,4]; 
param = [6,4]; %lvl = "20";
Niter = 10; lambda = 1; opts.Nblocks = [1,1]; opts2.Nblocks = [6,6];

for i=1:length(param)
    %for j=1:length(lvl)
        snr_scuba = []; snr_scuba2 = []; snr_bscott = []; snr_bscott2 = [];
        eval(sprintf('load(''data_%dS%dBInfdB.mat'')',param(i,1),param(i,2)))
        
        for F=1:25
            for r3=1:size(MSI,3)
                [F r3]
                %try
                [SRI_hat, ~] = scuba(MSI, HSI,Pm, [F,r3], opts); 
                [SRI_hat1, ~] = scuba(MSI, HSI,Pm, [F,r3], opts2);
                snr_scuba(F,r3) = r_snr(X,SRI_hat); snr_scuba2(F,r3) = r_snr(X,SRI_hat1);
                %catch
                    %snr_scuba(F,r3) = NaN; snr_scuba2(F,r3) = NaN;
                    %continue
                %end
            end
        end
%         
        figure
        surf(1:size(MSI,3),1:25,snr_scuba,'FaceAlpha',0.6); colormap summer; freezeColors; hold on
        surf(1:size(MSI,3),1:25,snr_scuba2,'FaceAlpha',0.6); colormap spring
        xlabel('R_3'); ylabel('F'); zlabel('R-SNR (dB)'); legend('(1,1) blocks','(4,4) blocks')
        title(sprintf('SCUBA, %d sources, Km=%d, Inf dB noise',param(i,1),param(i,2)))
        %savefig(gcf, sprintf('SNR_BCP%dM%sdB',param(i,1),lvl(j)))
%     
%         for r1=1:40%+10
%             for r3=1:size(MSI,3)
%                 R = [r1,r1,r3]
%                 [SRI_hat2,~] = bscott(MSI, HSI, Pm,R,opts);
%                 [SRI_hat3,~] = bscott(MSI, HSI, Pm,R,opts2);
%                 snr_bscott(r1,r3) = r_snr(X,SRI_hat2);
%                 snr_bscott2(r1,r3) = r_snr(X,SRI_hat3);
                if not((r3<=size(MSI,3) || r1<=size(HSI,1)) && (r1<=min(r3,size(MSI,3))*r1 && r3<=min(r1,size(HSI,1))^2))
                    snr_bscott(r1,r3) = NaN; snr_bscott2(r1,r3) = NaN;
                end  
%             end
%         end
%         
%         figure
%         surf(1:size(MSI,3),1:40,snr_bscott, 'FaceAlpha', 0.6); colormap summer; freezeColors; hold on
%         surf(1:size(MSI,3),1:40,snr_bscott2, 'FaceAlpha',0.6); colormap spring
%         xlabel('R_3'); ylabel('R_1 = R_2'); zlabel('SNR (dB)'); legend('(1,1) blocks','(4,4) blocks')
%         title(sprintf('BSCOTT, %d sources, Km=%d, %sdB noise',param(i,1),param(i,2),lvl(j)))
%         savefig(gcf, sprintf('SNR_BSVD%dM%sdB',param(i,1),lvl(j)))
      
        
    %end
end

%% VISUALIZE FIGURES - BLIND

lvl = ["20","35","60"]; 
F = 2; R = [6,6,4]; opts.Nblocks = [2,2];

for i=1:length(lvl)
    
    
    eval(sprintf('load(''data_2S6B%sdB.mat'')',lvl(i)))
    
    [SRI_hat1, ~] = scuba(MSI, HSI, Pm, [F,2], opts);
    SRI_hat1 = real(SRI_hat1);
    [SRI_hat2, ~] = bscott(MSI, HSI, Pm, R, opts);
    
    [nmse(X,SRI_hat1) nmse(X,SRI_hat2)]

end

%% TEST ON NUMBER OF BLOCKS

lvl = ["20","35","60"];
for j=1:length(lvl)
    
    eval(sprintf('load(''data_6S6B%sdB.mat'')',lvl(j)))
    snr1 = []; snr2 = [];
    for i=1:40
        %for j=1:size(HSI,2)

            [i,i]
            opts.Nblocks = [i,i]; ranks1 = [7,6]; ranks2 = [6,6,6];
            [SRI_hat1,~] = scuba(MSI, HSI, Pm,ranks1,opts);
            snr1(i,1) = r_snr(X,SRI_hat1);
            [SRI_hat2,~] = bscott(MSI, HSI, Pm,ranks2,opts);
            snr2(i,1) = r_snr(X,SRI_hat2);

        %end
    end
    
    figure
    subplot(1,2,1)
    plot(snr1)
    subplot(1,2,2)
    plot(snr2)
    savefig(gcf, sprintf('Nblocks_6M6B%sdB',lvl(j)))
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
a=hgload('SNR_BCP2M20dB.fig');
b=hgload('SNR_BSVD2M20dB.fig');
% c=hgload('Nblocks_6M4B60dB.fig');
% d=hgload('SNR_SVD2M35dB.fig');
% e=hgload('SNR_CP2M60dB.fig');
% f=hgload('SNR_SVD2M60dB.fig');
% Prepare subplots
figure
%sgtitle('SNR between groundtruth and estimate, 2 mat., Km=6, SNR=20dB')
h(1)=subplot(1,2,1); xlabel('R_3'); ylabel('F'); zlabel('R-SNR (dB)'); legend('(1,1) blocks','(4,4) blocks')
h(2)=subplot(1,2,2); xlabel('R_3'); ylabel('R_1 = R_2'); zlabel('R-SNR (dB)'); legend('(1,1) blocks','(4,4) blocks')
%h(3)=subplot(1,3,3); xlabel('Number of blocks'); ylabel('R-SNR (dB)'); legend('SCUBA','BSCOTT')
% h(3)=subplot(3,2,3);
% h(4)=subplot(3,2,4);
% h(5)=subplot(3,2,5);
% h(6)=subplot(3,2,6);
% Paste figures on the subplots
copyobj(allchild(get(a,'CurrentAxes')),h(1));
copyobj(allchild(get(b,'CurrentAxes')),h(2));
%copyobj(allchild(get(c,'CurrentAxes')),h(3));
% copyobj(allchild(get(d,'CurrentAxes')),h(4));
% copyobj(allchild(get(e,'CurrentAxes')),h(5));
% copyobj(allchild(get(f,'CurrentAxes')),h(6));


% 
% sp1 = findobj(subplot(1,2,1), 'Type', 'line');
% sp2 = findobj(subplot(1,2,2), 'Type', 'line');
% figure
% plot(sp1.YData,'Linewidth',1); hold on; plot(sp2.YData,'Linewidth',1)
% xlabel('Number of blocks in both dimensions'); ylabel('R-SNR (dB)'); legend('SCUBA','BSCOTT')

%% MAKE BOXPLOTS

lvl = "35"; %lvl = ["20","35","60"];
%param = [2,6; 3,2; 4,2; 6,4]; PENSER A FAIRE 4M
param = [6,4;];
Niter = 10; lambda = 1; Nreal = 10;

for i=1:length(param)
    for j=1:length(lvl)
        
        eval(sprintf('load(''rawdata_%dS%dB.mat'')',param(i,1),param(i,2)))
        MSI_true = MSI; HSI_true = HSI;
        R3 = param(i,1);%R3 = [param(i,1) param(i,2)];
        snr_scott = zeros(Nreal,size(HSI,1)+10);
        
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
            
            
%             parfor r1=1:size(HSI,1)+10
%                 R = [r1,r1,R3]
%                 [SRI_hat2,~] = scott(HSI, MSI, P1, P2, Pm, R);
%                 snr_scott(n,r1) = r_snr(X,SRI_hat2);
% 
%                 if (R3>size(MSI,3) && r1>size(HSI,1))%identifiability
%                      snr_scott(n,r1) = NaN;
%                 end
%             end          
        end
        
        figure
        boxplot(snr_stereo)
        xlabel('F');
        title(sprintf('STEREO, %d sources, Km=%d, %sdB SNR',param(i,1),param(i,2),lvl(j)))
        savefig(gcf, sprintf('BP_CP%dM%sdB',param(i,1),lvl(j)))
%         
%         figure
%         boxplot(snr_scott)
%         xlabel('R_1 = R_2');
%         title(sprintf('SCOTT, %d sources, Km=%d, %sdB SNR, R_3 = %d',param(i,1),param(i,2),lvl(j),R3))
%         savefig(gcf, sprintf('BP_SVD%dM%sdB%dRs',param(i,1),lvl(j),R3))

           
    end
end

%% TABLES OF METRICS 

lvl = ["20","35","60","Inf"]; F = [2,2,2,3]; R = [2,2,2];

for i=1:length(lvl)
    
    eval(sprintf('load(''data_2S6B%sdB.mat'')',lvl(i)))
    
    tic; [SRI_hat1, ~] = stereo3(HSI, MSI, P1,P2,Pm, F(i)); t1=toc;
    SRI_hat1 = real(SRI_hat1);
    err1 = compute_metrics(X,SRI_hat1,4,4); err1 = [cell2mat(err1) t1];
    tic; [SRI_hat2, ~] = scott(HSI, MSI, P1, P2, Pm, R); t2 = toc;
    err2 = compute_metrics(X,SRI_hat2,4,4);  err2 = real([cell2mat(err2) t2]);
    err1(3) = []; err2(3) = [];
    T = [err1; err2]; save(sprintf('table_2S6B%sdB.txt',lvl(i)),'T','-ASCII')
    
end



% C = info1.factors{3}; W = info2.factors{3}; 
% theta = subspace(C,W)
% 
% err1 = zeros(1,200); err2 = err1;
% for k=1:size(X,3)
%     err1(1,k) = norm(X(:,:,k)-SRI_hat1(:,:,k),'fro');
%     err2(1,k) = norm(X(:,:,k)-SRI_hat2(:,:,k),'fro');
% end
      
%% VISUALIZE FIGURES

lvl = ["60"]; 
%F = [2,2,2,3]; R = [2,2,2];
%F = [6,10,13,25]; R = [6,6,4];
opts.Nblocks = [1,1]; opts2.Nblocks = [6,6];

for i=1:length(lvl)
    
    eval(sprintf('load(''data_6S4B%sdB.mat'')',lvl(i)))
    
    tic; [SRI_hat1, ~] = scuba(MSI, HSI,Pm, [1,4], opts); t1 = toc;
    %SRI_hat1 = real(SRI_hat1);
    err1 = compute_metrics(X,SRI_hat1,4,4); err1 = [cell2mat(err1) t1];
    
    tic; [SRI_hat2, ~] = scuba(MSI, HSI,Pm, [1,2], opts2); t2 = toc;
    err2 = compute_metrics(X,SRI_hat2,4,4); err2 = [cell2mat(err2) t2];
    
    tic; [SRI_hat3,~] = bscott(MSI, HSI, Pm,[6,6,4],opts); t3 = toc;
    err3 = compute_metrics(X,SRI_hat3,4,4); err3 = [cell2mat(err3) t3];
    
    tic; [SRI_hat4,~] = bscott(MSI, HSI, Pm,[2,2,2],opts2); t4 = toc;
    err4 = compute_metrics(X,SRI_hat4,4,4); err4 = [cell2mat(err4) t4];
    
    tic; [SRI_hat5,~] = stereo_blind(HSI, MSI, Pm, 3); t5 = toc;
    err5 = compute_metrics(X,SRI_hat5,4,4); err5 = [cell2mat(err5) t5];
    
    err1(3) = []; err2(3) = []; err3(3) = []; err4(3) = []; err5(3) = [];
    T = [err1; err2; err3; err4; err5]; save(sprintf('table_B6S4B%sdB.txt',lvl(i)),'T','-ASCII')
    
    
    
%     figure(1)
%         %for k=1:size(X,3)
%     subplot(2,3,1); imagesc(SRI_hat1(:,:,44)); title('SCUBA (1,1)')
%     subplot(2,3,2); imagesc(SRI_hat2(:,:,44)); title('SCUBA (6,6)')
%     subplot(2,3,4); imagesc(SRI_hat3(:,:,44)); title('BSCOTT (1,1)')
%     subplot(2,3,5); imagesc(SRI_hat4(:,:,44)); title('BSCOTT (6,6)')
%     subplot(2,3,3); imagesc(X(:,:,44)); title('Groundtruth')
%     subplot(2,3,6); imagesc(SRI_hat5(:,:,44)); title('Blind-STEREO F=3')
%     %pause(0.4)
%         %end
    
    
    
%     err1 = zeros(1,size(X,3)); err2 = err1; err3 = err1; err4 = err1;
%     for k=1:size(X,3)
%         err1(1,k) = norm(X(:,:,k)-SRI_hat1(:,:,k),'fro')/norm(X(:,:,k),'fro');
%         err2(1,k) = norm(X(:,:,k)-SRI_hat2(:,:,k),'fro')/norm(X(:,:,k),'fro');
%         err3(1,k) = norm(X(:,:,k)-SRI_hat3(:,:,k),'fro')/norm(X(:,:,k),'fro');
%         err4(1,k) = norm(X(:,:,k)-SRI_hat4(:,:,k),'fro')/norm(X(:,:,k),'fro');
%     end
%     
%     figure(i)
%     subplot(1,2,1); plot(err1); hold on; plot(err3); xlabel('Spectral bands'); ylabel('RMSE'); title('(1,1) blocks'); legend('SCUBA','BSCOTT')
%     subplot(1,2,2); plot(err2); hold on; plot(err4); xlabel('Spectral bands'); ylabel('RMSE'); title('(6,6) blocks'); legend('SCUBA','BSCOTT')
%     savefig(gcf, sprintf('bandwise_error6M4B%sdB',lvl(i)))
    

end

%% TEST FOR SVD - SEMI-REAL DATA

% %%%%%%%%%%
% %LOAD INDIAN PINES DATASET % 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SRI = cell2mat(struct2cell(load('Indian_pines.mat')));
% SRI(:,:,[104:108 150:163 220]) = []; %Regions of water absorption
% SRI(1,:,:) = []; SRI(:,1,:) = [];
% Pm = spectral_deg(SRI,"LANDSAT");
% MSI = tmprod(SRI,Pm,3);
% d1 = 4; d2 = 4; q = 9;
% [P1,P2] = spatial_deg(SRI, q, d1, d2);
% HSI = tmprod(tmprod(SRI,P1,1),P2,2);

% %%%%%%%%%%%
% % LOAD SALINAS DATASET % 
% %%%%%%%%%%%%%%%%%%%%%%%%
% SRI = cell2mat(struct2cell(load('Salinas.mat')));
% SRI = crop(SRI,[80,84,size(SRI,3)]);
% SRI(:,:,[108:112 154:167 224]) = []; %Regions of water absorption (Salinas)
% Pm = spectral_deg(SRI,"LANDSAT");
% MSI = tmprod(SRI,Pm,3);
% d1 = 4; d2 = 4; q = 9;
% [P1,P2] = spatial_deg(SRI, q, d1, d2);
% HSI = tmprod(tmprod(SRI,P1,1),P2,2);

% %%%%%%%%%%%
% % LOAD PAVIA DATASET % 
% %%%%%%%%%%%%%%%%%%%%%%
% SRI = cell2mat(struct2cell(load('PaviaU.mat')));
% Pm = spectral_deg(SRI,"Quickbird");
% MSI = tmprod(SRI,Pm,3);
% d1 = 4; d2 = 4; q = 9;
% [P1,P2] = spatial_deg(SRI, q, d1, d2);
% HSI = tmprod(tmprod(SRI,P1,1),P2,2);

%%%%%%%%%%%
% LOAD CUPRITE DATASET % 
%%%%%%%%%%%%%%%%%%%%%%%%
SRI = double(cell2mat(struct2cell(load('Cuprite.mat'))));
Pm = spectral_deg(SRI,"LANDSAT");
MSI = tmprod(SRI,Pm,3);
d1 = 4; d2 = 4; q = 9;
[P1,P2] = spatial_deg(SRI, q, d1, d2);
HSI = tmprod(tmprod(SRI,P1,1),P2,2);

% for k=1:size(MSI,3)
%     MSI(:,:,k) = awgn(MSI(:,:,k),20);
% end
% for k=1:size(HSI,3)
%     HSI(:,:,k) = awgn(HSI(:,:,k),20);
% end

figure
subplot(1,3,1)
x = svd(tens2mat(MSI,1,[]));
semilogy(x)
subplot(1,3,2)
x = svd(tens2mat(MSI,2,[]));
semilogy(x)
subplot(1,3,3)
x = svd(tens2mat(HSI,3,[]));
semilogy(x)

%% TEST ON CURVATURE

lambda = 1;
for r1=1:size(HSI,1)+10
    for r3=1:10
        R = [r1,r1,r3]
        [SRI_hat2,info] = scott(HSI, MSI, P1, P2, Pm, R);
        cost2(r1,r3) = cost_scott(info.factors,info.core, HSI, MSI, lambda, P1,P2,Pm);
        snr_scott(r1,r3) = r_snr(X,SRI_hat2);
        if not((r3<=size(MSI,3) || r1<=size(HSI,1)) && (r1<=min(r3,size(MSI,3))*r1 && r3<=min(r1,size(HSI,1))^2))
             cost2(r1,r3) = NaN; snr_scott(r1,r3) = NaN;
        end
    end
end

[X,Y] = meshgrid(1:10, 1:size(HSI,1)+10);
[K,H] = surfature(X,Y,cost2); 
figure
surf(X,Y,cost2,H,'facecolor','interp'); 
set(gca,'clim',[-1,1])
figure
surf(X,Y,snr_scott)


