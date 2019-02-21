% EXPERIMENT 1: HOSVD OF UNFOLDINGS %
% Copyright (c) 2018 Clemence Prevost, Konstantin Usevich, Pierre Comon, David Brie
% https://github.com/cprevost4/HSR_Tucker
% Contact: clemence.prevost@univ-lorraine.fr

%% Load Pavia University

% 1. load data
SRI = cell2mat(struct2cell(load('PaviaU.mat')));
SRI(1:2,:,:) = []; SRI(:,1:4,:) = [];
% 2. degradation
Pm = spectral_deg(SRI,"Quickbird");
MSI = tmprod(SRI,Pm,3);
d1 = 4; d2 = 4; q = 9;
[P1,P2] = spatial_deg(SRI, q, d1, d2);
HSI = tmprod(tmprod(SRI,P1,1),P2,2);

%% SVD

figure
        
   subplot(1,3,1)
    x = log(svd(tens2mat(MSI,1,[])));
    ylim auto
    loglog(x(1:10),'-s','MarkerSize',3,'Linewidth',1.1)
    title('1st unfolding MSI')
   subplot(1,3,2)
    x = log(svd(tens2mat(MSI,2,[])));
    ylim auto
    loglog(x(1:10),'-s','MarkerSize',3,'Linewidth',1.1)
    title('2nd unfolding MSI')
   subplot(1,3,3)
    x = log(svd(tens2mat(HSI,3,[])));
    ylim auto
    loglog(x(1:10),'-s','MarkerSize',3,'Linewidth',1.1)
    title('3rd unfolding HSI')