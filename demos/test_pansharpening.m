% EXPERIMENTS ON WOOD DATA %
% Copyright (c) 2018 Clemence Prevost, Konstantin Usevich, Pierre Comon, David Brie
% https://github.com/cprevost4/HSR_Tucker
% Contact: clemence.prevost@univ-lorraine.fr

%% Load Indian Pines

SRI = cell2mat(struct2cell(load('Indian_pines.mat')));
SRI(:,:,[104:108 150:163 220]) = []; %Regions of water absorption (Indian pines)
SRI(1,:,:) = []; SRI(:,1,:) = [];
% 2. degradation
d1 = 4; d2 = 4; q = 9;
[P1,P2] = spatial_deg(SRI, q, d1, d2);
HSI = tmprod(tmprod(SRI,P1,1),P2,2);
Km = 1; Pm = zeros(Km,size(SRI,3));
spec_range = linspace(400,2500,size(SRI,3));
Pm(1,:) = 1/length(spec_range);
MSI = tmprod(SRI,Pm,3); 

%% Make table

methods = {'SCOTT' 'scott' '[24,24,25]' [];...
           'SCOTT' 'scott' '[30,30,16]' [];...
           'SCOTT' 'scott' '[35,35,6]' [];
           'HySure','hysure_b1_adaptor','[]',16};      
DegMat = struct('Pm', Pm, 'P1', P1, 'P2', P2);         
res = compare_methods(SRI, HSI, MSI, DegMat, [d1 d2], methods); 
T = cell2mat(res(:,2:end)); save('exp4_table2_ip.txt','T','-ASCII')

%% Run simulations

R1 = 1:36; R3 = 1:25; lambda = 1;
for i=1:length(R1)
    for j=1:length(R3)
        R = [R1(i), R1(i), R3(j)];
        filename = sprintf('HSPan_exp4_%d_%d_%d_IP',R1(i),R1(i),R3(j));
        if not((j<=size(MSI,3) || i<=size(HSI,1)) && (i<=min(j,size(MSI,3))*i && j<=min(i,size(HSI,1))^2))
            snr = NaN; cost = NaN;
        else  
            [SRI_hat,info] = scott(HSI, MSI, P1, P2, Pm, R);
            snr = r_snr(SRI,SRI_hat); 
            cost = cost_scott(info.factors,info.core, HSI, MSI, lambda, P1,P2,Pm);
        end
        %save(filename,'SRI_hat','cost','snr');
    end
end

%% Make figure

snr1 = []; cost1 = [];
for i=1:36
    for j=1:25
        eval(sprintf('load(''HSPan_exp4_%d_%d_%d_IP'')',i,i,j));
        snr1(i,j) = snr; cost1(i,j) = cost;
    end
end

R3 = 1:25; R1 = 1:36;
figure(1)
surfc(R3,R1,snr1)
ylabel('R1=R2'); xlabel('R3'); zlabel('SNR(dB)');
title('SNR between SRI and estimate for R1=R2 and R3')
saveas(gcf,'fig_exp4_snr_R2f_IP','fig')
figure(2)
surfc(R3,R1,cost1)
ylabel('R1=R2'); xlabel('R3'); zlabel('Value of cost function');
title('Cost function value between SRI and estimate for R1=R2 and R3')
saveas(gcf,'fig_exp4_cost_R2f_IP','fig')
