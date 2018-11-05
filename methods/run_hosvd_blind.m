function [SRI_hat,cost, err] = run_hosvd_blind(SRI,MSI,HSI,R,P1,P2,Pm, alpha)

% RUN_HOSVD runs the HOSVD algorithm for specified rank R
% [SRI_hat,cost, err] = RUN_HOSVD(SRI,MSI,HSI,R,P1,P2,Pm, alpha) returns 
% estimation of SRI, value of cost function and metrics in the cell array err
% 
% INPUT ARGUMENTS:
%     SRI, MSI, HSI: input datasets (resp. groundtruth SRI, MSI and HSI
%     R: specified multilinear rank
%     P1,P2,Pm: spatial and spectral degratation matrices
%     alpha: possible regularization (usually set to zero when fct is employed
% OUTPUT ARGUMENTS:
%     SRI_hat: estimated SRI
%     cost: value of the cost function
%     err: cell array of metrics
% Copyright (c) 2018 Clémence Prévost, Konstantin Usevich, Pierre Comon, David Brie
% https://github.com/cprevost4/HSR_Tucker

tic,
[fact,S] = mlsvd(MSI, R); % hosvd
U = cell2mat(fact(1)); V = cell2mat(fact(2)); W_tilde = cell2mat(fact(3));    
[W_H, ~, ~] = svds(tens2mat(HSI,3,[]), R(3));
T = (Pm * W_H) \ W_tilde;

lam = 1;

W = W_H * T;
SRI_hat = lmlragen({U,V,W},S); 
time = toc

cost = frob(HSI - lmlragen({P1*U,P2*V,W}, S),'squared') + frob(MSI - lmlragen({U,V,Pm*W},S),'squared');
err = {cost nmse(SRI,SRI_hat), cc(SRI,SRI_hat), sam(SRI,SRI_hat), ergas(SRI,SRI_hat,1/4), r_snr(SRI,SRI_hat), time};
%err = r_snr(SRI,SRI_hat);

end

