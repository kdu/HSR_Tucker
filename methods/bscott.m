function [SRI_hat,info] = bscott(MSI,HSI,Pm,R,opts)

% BSCOTT runs the BSCOTT algorithm for specified rank R
% [SRI_hat,info] = BSCOTT(MSI,HSI,R,Pm,opts) returns 
% estimation of SRI and info structure
% 
% INPUT ARGUMENTS:
%     MSI, HSI: input datasets (resp. MSI and HSI)
%     R: specified multilinear rank
%     Pm: spectral degradation matrix
%     opts: options structure
% OUTPUT ARGUMENTS:
%     SRI_hat: estimated SRI
%     info: informative structure
%
% Copyright (c) 2018 Clemence Prevost, Konstantin Usevich, Pierre Comon, David Brie
% https://github.com/cprevost4/HSR_Tucker
% Contact: clemence.prevost@univ-lorraine.fr

if ~exist('opts','var')
    opts = struct();
end
if ~isfield(opts,'Nblocks') || isempty(opts.Nblocks)
    opts.Nblocks = [1,1];
end

step_MSI = size(MSI); step_HSI = size(HSI);
step_MSI = step_MSI(1:2) ./ opts.Nblocks; step_HSI = step_HSI(1:2) ./ opts.Nblocks;
step_MSI = floor(step_MSI); step_HSI = floor(step_HSI);
SRI_hat = zeros(size(MSI,1), size(MSI,2), size(HSI,3));

for i1=1:opts.Nblocks(1)-1
  for i2=1:opts.Nblocks(2)-1
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

for i2 = 1:opts.Nblocks(2)-1
    [fact,S] = mlsvd(MSI(1+(opts.Nblocks(1)-1)*step_MSI(1):size(MSI,1), ...
                         (1:step_MSI(2)) + (i2-1)*step_MSI(2), :), R); % hosvd
    U = cell2mat(fact(1)); V = cell2mat(fact(2)); W_tilde = cell2mat(fact(3));    
    [W_H, ~, ~] = svds(tens2mat(HSI(1+(opts.Nblocks(1)-1)*step_HSI(1):size(HSI,1), ...
                                (1:step_HSI(2)) + (i2-1)*step_HSI(2), :),3,[]), R(3));
    T = (Pm * W_H) \ W_tilde;

    W = W_H * T;
    SRI_hat(1+(opts.Nblocks(1)-1)*step_MSI(1):size(MSI,1), ...
            (1:step_MSI(2)) + (i2-1)*step_MSI(2), :) = lmlragen({U,V,W},S); 
end

for i1 = 1:opts.Nblocks(1)-1
    [fact,S] = mlsvd(MSI((1:step_MSI(1)) + (i1-1)*step_MSI(1), ...
                         1+(opts.Nblocks(2)-1)*step_MSI(2):size(MSI,2), :), R); % hosvd
    U = cell2mat(fact(1)); V = cell2mat(fact(2)); W_tilde = cell2mat(fact(3));    
    [W_H, ~, ~] = svds(tens2mat(HSI((1:step_HSI(1)) + (i1-1)*step_HSI(1), ...
                                1+(opts.Nblocks(2)-1)*step_HSI(2):size(HSI,2), :),3,[]), R(3));
    T = (Pm * W_H) \ W_tilde;

    W = W_H * T;
    SRI_hat((1:step_MSI(1)) + (i1-1)*step_MSI(1), ...
            1+(opts.Nblocks(2)-1)*step_MSI(2):size(MSI,2), :) = lmlragen({U,V,W},S); 
end


    [fact,S] = mlsvd(MSI(1+(opts.Nblocks(1)-1)*step_MSI(1):size(MSI,1), ...
                         1+(opts.Nblocks(2)-1)*step_MSI(2):size(MSI,2), :), R); % hosvd
    U = cell2mat(fact(1)); V = cell2mat(fact(2)); W_tilde = cell2mat(fact(3));    
    [W_H, ~, ~] = svds(tens2mat(HSI(1+(opts.Nblocks(1)-1)*step_HSI(1):size(HSI,1), ...
                                1+(opts.Nblocks(2)-1)*step_HSI(2):size(HSI,2), :),3,[]), R(3));
    T = (Pm * W_H) \ W_tilde;

    W = W_H * T;
    SRI_hat(1+(opts.Nblocks(1)-1)*step_MSI(1):size(MSI,1), ...
           1+(opts.Nblocks(2)-1)*step_MSI(2):size(MSI,2), :) = lmlragen({U,V,W},S); 

info.factors = {U,V,W};
info.core = {S};
info.rank = {R};

end

