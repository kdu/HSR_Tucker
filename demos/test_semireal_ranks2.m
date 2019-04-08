% EXPERIMENT 1: HOSVD FOR VARIOUS RANKS %
% Copyright (c) 2018 Clemence Prevost, Konstantin Usevich, Pierre Comon, David Brie
% https://github.com/cprevost4/HSR_Tucker
% Contact: clemence.prevost@univ-lorraine.fr

%% Load Salinas-A

% 1. load data
SRI = cell2mat(struct2cell(load('SalinasA.mat')));
SRI(:,:,[108:112 154:167 224]) = []; %Regions of water absorption (Salinas)
SRI = crop(SRI,[80,84,size(SRI,3)]);
% 2. degradation
Pm = spectral_deg(SRI,"LANDSAT");
d1 = 4; d2 = 4; q = 9;
[P1,P2] = spatial_deg(SRI, q, d1, d2);
MSI = tmprod(SRI,Pm,3); HSI = tmprod(tmprod(SRI,P1,1),P2,2);

%% Simulations

lambda = 1; 

%R1=R2, R3 varies
R1 = 10:50; R3 = 2:25; snr = []; cost = [];
for i=1:length(R1)
    for j=1:length(R3)
        R = [R1(i),R1(i),R3(j)]
        %filename = sprintf('data_exp1_%d_%d_%d_Sal',R1(i),R1(i),R3(j));
        if not((R3(j)<=size(MSI,3) || R1(i)<=size(HSI,1)) && (R1(i)<=min(R3(j),size(MSI,3))*R1(i) && R3(j)<=min(R1(i),size(HSI,1))^2))
            snr(i+9,j+1) = NaN; cost(i+9,j+1) = NaN;
        else  
            [SRI_hat,info] = scott(HSI, MSI, P1, P2, Pm, R);
            snr(i+9,j+1) = r_snr(SRI,SRI_hat); 
            cost(i+9,j+1) = cost_scott(info.factors,info.core, HSI, MSI, lambda, P1,P2,Pm);
        end
        %save(filename,'SRI_hat','cost','snr');
    end
end

%% Make figure

R3 = 2:25; R1 = 10:50;

snr(find(snr==0)) = 24;
cost(find(cost==0)) = 0.2e9;



figure(1)
shading interp
pcolor(snr)
xlabel("R_3");ylabel("R_1 = R_2");
rectangle('Position',[7 21 50 50],'FaceColor',[1,1,1])
xlim([2 25]); ylim([10 50]);


 
%title('SNR between SRI and estimate for R1=R2 and R3'); 
%saveas(gcf,'fig_exp1_snr_R2f_IP','fig')
% figure(2)
% surfc(R3,R1,cost1)
% ylabel('R1=R2'); xlabel('R3'); zlabel('Value of cost function');colormap
%title('Cost function value between SRI and estimate for R1=R2 and R3'); 
%saveas(gcf,'fig_exp1_cost_R2f_IP','fig')
figure(2)
pcolor(cost)
xlabel("R_3");ylabel("R_1 = R_2");
rectangle('Position',[7 21 50 50],'FaceColor',[1,1,1])
xlim([2 25]); ylim([10 50]);
