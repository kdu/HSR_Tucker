function [SRI_hat, info] = run_hosvd(HSI, MSI, P1, P2, Pm, ranks, opts)

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
    lambda = 1; Display = 'false'; alpha = 0;
elseif nargin==7
    lambda = opts.lambda; Display = opts.Display; alpha = opts.alpha;
end

[U, ~, ~] = svds(tens2mat(MSI,1,[]),ranks(1));
[V, ~, ~] = svds(tens2mat(MSI,2,[]),ranks(2));
[W, ~, ~] = svds(tens2mat(HSI,3,[]), ranks(3));

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
C = reshape(b_old, ranks(1)*ranks(2), ranks(3));
S = reshape(sylvester((A+alpha*norm(A,2)^2*ones(size(A,2))),B,C),ranks);

SRI_hat = lmlragen({U,V,W},S);
info.factors = {'U','V','W'};
info.core = {'S'};
info.rank = {'ranks'};

end

