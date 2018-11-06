function [SRI_hat, info] = scott(HSI, MSI, P1, P2, Pm, R, opts)

% SCOTT runs the SCOTT algorithm for specified rank R
% [SRI_hat,info] = SCOTT(HSI, MSI, P1, P2, Pm, ranks, opts) returns 
% estimation of SRI and informative structure
% 
% INPUT ARGUMENTS:
%     MSI, HSI: input datasets (resp. MSI and HSI)
%     P1,P2,Pm: spatial and spectral degratation matrices
%     ranks: specified multilinear rank
%     opts: options structure 
% OUTPUT ARGUMENTS:
%     SRI_hat: estimated SRI
%     info: informative structure
%
% Copyright (c) 2018 Clemence Prevost, Konstantin Usevich, Pierre Comon, David Brie
% https://github.com/cprevost4/HSR_Tucker
% Contact: clemence.prevost@univ-lorraine.fr

if nargin==6
    lambda = 1; Display = 'false'; alpha = 0;
elseif nargin==7
    lambda = opts.lambda; Display = opts.Display; alpha = opts.alpha;
end

[U, ~, ~] = svds(tens2mat(MSI,1,[]),R(1));
[V, ~, ~] = svds(tens2mat(MSI,2,[]),R(2));
[W, ~, ~] = svds(tens2mat(HSI,3,[]), R(3));

%%%%%%%%%%%%%%%%%%%
% ALWAYS WORKS BUT SLOWER  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A = kron(eye(R(3)),kron(V'*(P2'*P2)*V, U'*(P1'*P1)*U)) + lam * kron(W'*(Pm'*Pm)*W, eye(R(1)*R(2)));
% b = tmprod(HSI,{U'*P1', V'*P2', W'},[1,2,3]) + lam * tmprod(MSI,{U', V', W'*Pm'},[1,2,3]);
% S = reshape(A\ b(:),R);

%%%%%%%%%%%%%%%%%%
% ONLY WORKS UNDER IDENTIFIABILITY CONDITIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = kron(V'*(P2'*P2)*V, U'*(P1'*P1)*U);
B = lambda* W'*(Pm'*Pm)*W;
b_old = tmprod(HSI,{U'*P1', V'*P2', W'},[1,2,3]) + lambda * tmprod(MSI,{U', V', W'*Pm'},[1,2,3]);
C = reshape(b_old, R(1)*R(2), R(3));
S = reshape(sylvester((A+alpha*norm(A,2)^2*ones(size(A,2))),B,C),R);

SRI_hat = lmlragen({U,V,W},S);
info.factors = {U,V,W};
info.core = {S};
info.rank = {ranks};

end

