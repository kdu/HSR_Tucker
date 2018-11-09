function [SRI_hat,info] = scuba(MSI,HSI,Pm, ranks, opts)

% SCUBA runs the SCUBA algorithm for specified rank R
% [SRI_hat,info] = SCUBA(MSI,HSI,ranks,Pm, opts) returns 
% estimation of SRI and informative structure
% 
% INPUT ARGUMENTS:
%     MSI, HSI: input datasets (resp. MSI and HSI)
%     ranks: field of CP rank + spectral rank
%     Pm: spectral degratation matrix
%     opts: options structure
% OUTPUT ARGUMENTS:
%     SRI_hat: estimated SRI
%     info: informative structure
%     
% Copyright (c) 2018 Clemence Prevost, Konstantin Usevich, Pierre Comon, David Brie
% https://github.com/cprevost4/HSR_Tucker
% Contact: clemence.prevost@univ-lorraine.fr

F = ranks(1); R = ranks(2);

if ~exist('opts','var')
    opts = struct();
end
if ~isfield(opts,'Nblocks') || isempty(opts.Nblocks)
    opts.Nblocks = [1,1];
end

step_MSI = size(MSI); step_MSI = step_MSI(1:2) ./ opts.Nblocks ;
step_HSI = size(HSI); step_HSI = step_HSI(1:2) ./ opts.Nblocks ;
SRI_hat = zeros(size(MSI,1), size(MSI,2), size(HSI,3));


for i1=1:opts.Nblocks(1)
  for i2=1:opts.Nblocks(2) 
    %options.Display = true;  
    options.MaxIter = 25; %options.Display = true;
    U = cpd(MSI((1:step_MSI(1)) + (i1-1)*step_MSI(1), ...
                         (1:step_MSI(2)) + (i2-1)*step_MSI(2), :), F,options);
    A = cell2mat(U(1)); B = cell2mat(U(2)); C_tilde = cell2mat(U(3));

    [W_H, ~, ~] = svds(tens2mat(HSI((1:step_HSI(1)) + (i1-1)*step_HSI(1), ...
                                (1:step_HSI(2)) + (i2-1)*step_HSI(2), :),3,[]), R);
    T = (Pm * W_H) \ C_tilde;

    C = W_H * T;
    SRI_hat((1:step_MSI(1)) + (i1-1)*step_MSI(1), ...
            (1:step_MSI(2)) + (i2-1)*step_MSI(2), :) = cpdgen({A,B,C}); 
  end
end 

info.steps = {step_MSI,step_HSI};
info.factors = {A,B,C};
info.rank = {R, F};


end

