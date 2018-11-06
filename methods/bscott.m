function [SRI_hat,info] = bscott(MSI,HSI,R,Pm,opts)

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
% Contact: clemence.prevost@univ-lorraine.fr

if nargin==4
    Nblocks = [1,1]; 
elseif nargin==5
    Nblocks = opts.Nblocks; 
end

step_MSI = size(MSI); step_MSI = step_MSI(1:2) ./ Nblocks; 
step_HSI = size(HSI); step_HSI = step_HSI(1:2) ./ Nblocks ;
SRI_hat = zeros(size(MSI,1), size(MSI,2), size(HSI,3));

for i1=1:Nblocks(1)
  for i2=1:Nblocks(2) 
    [fact,S] = mlsvd(MSI((1:step_MSI(1)) + (i1-1)*step_MSI(1), ...
                         (1:step_MSI(2)) + (i2-1)*step_MSI(2), :), R); % hosvd
    U = cell2mat(fact(1)); V = cell2mat(fact(2)); W_tilde = cell2mat(fact(3));    
    [W_H, ~, ~] = svds(tens2mat(HSI((1:step_HSI(1)) + (i1-1)*step_HSI(1), ...
                                (1:step_HSI(2)) + (i2-1)*step_HSI(2), :),3,[]), R(3));
    T = (Pm * W_H) \ W_tilde;

    W = W_H * T;
    SRI_hat((1:step_MSI(1)) + (i1-1)*step_MSI(1), ...
            (1:step_MSI(2)) + (i2-1)*step_MSI(2), :) = lmlragen({U,V,W},S); 
  end
end 

info.factors = {'U','V','W'};
info.core = {'S'};
info.rank = {'R'};

end

