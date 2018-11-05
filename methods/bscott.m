function [SRI_hat,info] = bscott(HSI, MSI, ranks, Pm, opts)

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

if nargin==6
    Display = 'false';
elseif nargin==7
    Display = opts.Display;
end

[fact,S] = mlsvd(MSI, ranks); % hosvd
U = cell2mat(fact(1)); V = cell2mat(fact(2)); W_tilde = cell2mat(fact(3));    
[W_H, ~, ~] = svds(tens2mat(HSI,3,[]), ranks(3));
T = (Pm * W_H) \ W_tilde;

W = W_H * T;
SRI_hat = lmlragen({U,V,W},S); 
info.factors = {'U','V','W'};
info.core = {'S'};
info.rank = {'ranks'};

end

