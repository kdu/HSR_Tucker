function [SRI_hat,cost, err] = run_hosvd(SRI,MSI,HSI,R,P1,P2,Pm, alpha)

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

tic 

[U, ~, ~] = svds(tens2mat(MSI,1,[]),R(1));
[V, ~, ~] = svds(tens2mat(MSI,2,[]),R(2));
[W, ~, ~] = svds(tens2mat(HSI,3,[]), R(3));

lam = 1;
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
B = lam* W'*(Pm'*Pm)*W;
b_old = tmprod(HSI,{U'*P1', V'*P2', W'},[1,2,3]) + lam * tmprod(MSI,{U', V', W'*Pm'},[1,2,3]);
C = reshape(b_old, R(1)*R(2), R(3));
S = reshape(sylvester(A,B,C),R);

%%%%%%%%%%%
% IF REGULARIZATION %
%%%%%%%%%%%%%%%%%%%%%
%S = reshape((A + alpha*norm(A,2)^2*ones(size(A,2)))\ b(:),R); %alpha = regularization weight

SRI_hat = lmlragen({U,V,W},S);
time = toc 
cost = frob(HSI - lmlragen({P1*U,P2*V,W}, S),'squared') + frob(MSI - lmlragen({U,V,Pm*W},S),'squared');
err = {cost nmse(SRI,SRI_hat), cc(SRI,SRI_hat), sam(SRI,SRI_hat), ergas(SRI,SRI_hat,1/4), r_snr(SRI,SRI_hat),time};
%err = r_snr(SRI,SRI_hat);

end

