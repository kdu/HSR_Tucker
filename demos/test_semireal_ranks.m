% EXPERIMENT 1: HOSVD FOR VARIOUS RANKS %
% Copyright (c) 2018 Clemence Prevost, Konstantin Usevich, Pierre Comon, David Brie
% https://github.com/cprevost4/HSR_Tucker
% Contact: clemence.prevost@univ-lorraine.fr


%% Load Indian pines 

SRI = cell2mat(struct2cell(load('Indian_pines.mat')));
SRI(:,:,[104:108 150:163 220]) = []; %Regions of water absorption
SRI(1,:,:) = []; SRI(:,1,:) = [];
Pm = spectral_deg(SRI,"LANDSAT");
MSI = tmprod(SRI,Pm,3);
d1 = 4; d2 = 4; q = 9;
[P1,P2] = spatial_deg(SRI, q, d1, d2);
HSI = tmprod(tmprod(SRI,P1,1),P2,2);

%% Simulations for Indian Pines 
lambda = 1;

%R1=R2, R3 varies
R1 = 10:50; R3 = 2:25;
for i=1:length(R1)
    for j=1:length(R3)
        R = [R1(i),R1(i),R3(j)];
        filename = sprintf('data_exp1_%d_%d_%d_IP',R1(i),R1(i),R3(j));
        if not((R3(j)<=size(MSI,3) || R1(i)<=size(HSI,1)) && (R1(i)<=min(R3(j),size(MSI,3))*R1(i) && R3(j)<=min(R1(i),size(HSI,1))^2))
            snr = NaN; cost = NaN;
        else  
            [SRI_hat,info] = scott(HSI, MSI, P1, P2, Pm, R);
            snr = r_snr(SRI,SRI_hat); 
            cost = cost_scott(info.factors,info.core, HSI, MSI, lambda, P1,P2,Pm);
        end
        save(filename,'SRI_hat','cost','snr');
    end
end

% R3 fixed to 6
R1 = 10:40; R2 = 10:40; R3 = 6;
lambda = 1;
for i=1:length(R1)
    for j=1:length(R2)
        R = [R1(i), R2(j), R3];
        filename = sprintf('data_exp1_%d_%d_%d_IP',R1(i),R2(j),R3);
        if not((R3<=size(MSI,3) || (i<=size(HSI,1) && j<=size(HSI,2))) && (i<=min(R3,size(MSI,3))*j && j<=min(R3,size(MSI,3))*i && R3<=min(i,size(HSI,1))*min(j,size(HSI,2))))
            snr = NaN; cost = NaN;
        else  
            [SRI_hat,info] = scott(HSI, MSI, P1, P2, Pm, R); 
            cost = cost_scott(info.factors,info.core, HSI, MSI, lambda, P1,P2,Pm);
        end
        save(filename,'SRI_hat','cost','snr');
    end
end

% R3 fixed to 16
R1 = 10:40; R2 = 10:40; R3 = 16;
lambda = 1;
for i=1:length(R1)
    for j=1:length(R2)
        R = [R1(i), R2(j), R3];
        filename = sprintf('data_exp1_%d_%d_%d_IP',R1(i),R2(j),R3);
        if not((R3<=size(MSI,3) || (i<=size(HSI,1) && j<=size(HSI,2))) && (i<=min(R3,size(MSI,3))*j && j<=min(R3,size(MSI,3))*i && R3<=min(i,size(HSI,1))*min(j,size(HSI,2))))
            snr = NaN; cost = NaN;
        else  
            [SRI_hat,info] = scott(HSI, MSI, P1, P2, Pm, R); 
            cost = cost_scott(info.factors,info.core, HSI, MSI, lambda, P1,P2,Pm);
        end
        save(filename,'SRI_hat','cost','snr');
    end
end

%% Figure, R1=R2

snr1 = []; cost1 = [];
for i=10:50
    for j=2:25
        eval(sprintf('load(''data_exp1_%d_%d_%d_IP'')',i,i,j));
        snr1(i-9,j-1) = snr; cost1(i-9,j-1) = cost;
    end
end

R3 = 2:25; R1 = 10:50;
figure(1)
surfc(R3,R1,snr1,'FaceColor','interp')
ylabel('R1=R2'); xlabel('R3'); zlabel('SNR(dB)');
title('SNR between SRI and estimate for R1=R2 and R3'); 
%saveas(gcf,'fig_exp1_snr_R2f_IP','fig')
figure(2)
surfc(R3,R1,cost1)
ylabel('R1=R2'); xlabel('R3'); zlabel('Value of cost function');
title('Cost function value between SRI and estimate for R1=R2 and R3')
%saveas(gcf,'fig_exp1_cost_R2f_IP','fig')

%% Figure, R3 = 6

snr2 = []; cost2 = [];
for i=10:40
    for j=10:40
        eval(sprintf('load(''data_exp1_%d_%d_%d_IP'')',i,j,6));
        snr2(i-9,j-9) = snr; cost2(i-9,j-9) = cost;
    end
end

R1 = 10:40; R2 = 10:40;
figure(3)
surfc(R2,R1,snr2)
xlabel('R2'); ylabel('R1'); zlabel('SNR(dB)');
title('SNR between SRI and estimate for R1 and R2 (R3=6)')
saveas(gcf,'fig_exp1_snr_R3f6','fig')
figure(4)
surfc(R1,R2,cost2)
xlabel('R1'); ylabel('R2'); zlabel('Value of cost function');
title('Cost function value between SRI and estimate for R1 and R2 (R3=6)')
saveas(gcf,'fig_exp1_cost_R3f6','fig')

%% Figure, R3 = 16

snr3 = []; cost3 = [];
for i=10:40
    for j=10:40
        eval(sprintf('load(''data_exp1_%d_%d_%d_IP'')',i,j,16));
        snr3(i-9,j-9) = snr; cost3(i-9,j-9) = cost;
    end
end

R1 = 10:40; R2 = 10:40;
figure(5)
surfc(R1,R2,snr3)
xlabel('R1'); ylabel('R2'); zlabel('SNR(dB)');
title('SNR between SRI and estimate for R1 and R2 (R3=16)')
saveas(gcf,'fig_exp1_snr_R3f16','fig')
figure(6)
surfc(R1,R2,cost3)
xlabel('R1'); ylabel('R2'); zlabel('Value of cost function');
title('Cost function value between SRI and estimate for R1 and R2 (R3=16)')
saveas(gcf,'fig_exp1_cost_R3f16','fig')

