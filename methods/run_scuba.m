function [SRI_hat,cost, err] = run_scuba(SRI,MSI,HSI,F,R,P1,P2,Pm, Nblocks)

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
% Copyright (c) 2018 Clemence Prevost, Konstantin Usevich, Pierre Comon, David Brie
% https://github.com/cprevost4/HSR_Tucker

step_MSI = size(MSI); step_MSI = step_MSI(1:2) ./ Nblocks 
step_HSI = size(HSI); step_HSI = step_HSI(1:2) ./ Nblocks 
SRI_hat = zeros(size(MSI,1), size(MSI,2), size(HSI,3));

tic,
for i1=1:Nblocks(1)
  for i2=1:Nblocks(2) 
    options.Display = true;  
    options.MaxIter = 25; %options.Display = true;
    U = cpd(MSI((1:step_MSI(1)) + (i1-1)*step_MSI(1), ...
                         (1:step_MSI(2)) + (i2-1)*step_MSI(2), :), F,options);
    A = cell2mat(U(1)); B = cell2mat(U(2)); C_tilde = cell2mat(U(3));

    [W_H, ~, ~] = svds(tens2mat(HSI((1:step_HSI(1)) + (i1-1)*step_HSI(1), ...
                                (1:step_HSI(2)) + (i2-1)*step_HSI(2), :),3,[]), R);
    T = (Pm * W_H) \ C_tilde;

    lam = 1;
    C = W_H * T;
    SRI_hat((1:step_MSI(1)) + (i1-1)*step_MSI(1), ...
            (1:step_MSI(2)) + (i2-1)*step_MSI(2), :) = cpdgen({A,B,C}); 
  end
end 
time = toc
    
    
cost = 0;
err = {cost nmse(SRI,SRI_hat), cc(SRI,SRI_hat), sam(SRI,SRI_hat), ergas(SRI,SRI_hat,1/4), r_snr(SRI,SRI_hat),time};

end

