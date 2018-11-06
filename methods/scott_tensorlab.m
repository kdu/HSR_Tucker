function [SRI_hat,info] = scott_tensorlab(HSI, MSI, P1, P2, Pm, R, opts)

% SCOTT_TENSORLAB runs the SCOTT algorithm for specified rank R followed by 
% TensorLab optimization
% [SRI_hat,info] = SCOTT_TENSORLAB(HSI, MSI, P1, P2, Pm, R, opts) returns 
% estimation of SRI and informative structure
% 
% INPUT ARGUMENTS:
%     MSI, HSI: input datasets (resp. MSI and HSI)
%     R: specified multilinear rank
%     P1,P2,Pm: spatial and spectral degratation matrices
%     opts: options structure
% OUTPUT ARGUMENTS:
%     SRI_hat: estimated SRI
%     info: informative structure
%     
% Copyright (c) 2018 Clemence Prevost, Konstantin Usevich, Pierre Comon, David Brie
% https://github.com/cprevost4/HSR_Tucker
% Contact: clemence.prevost@univ-lorraine.fr

if nargin==6
    lambda = 1; Display = 'false'; alpha = 0; opti_Niter = 30;
elseif nargin==7
    lambda = opts.lambda; Display = opts.Display; alpha = opts.alpha;
    opti_Niter = opts.opti_Niter;
end

[U, ~, ~] = svds(tens2mat(MSI,1,[]),R(1));
[V, ~, ~] = svds(tens2mat(MSI,2,[]),R(2));
[W, ~, ~] = svds(tens2mat(HSI,3,[]), R(3));

A = kron(V'*(P2'*P2)*V, U'*(P1'*P1)*U);
B = lambda* W'*(Pm'*Pm)*W;
b_old = tmprod(HSI,{U'*P1', V'*P2', W'},[1,2,3]) + lambda * tmprod(MSI,{U', V', W'*Pm'},[1,2,3]);
C = reshape(b_old, R(1)*R(2), R(3));
S = reshape(sylvester((A+alpha*norm(A,2)^2*ones(size(A,2))),B,C),R);

model = struct;
model.variables = {U,V,W,S};
model.factors.A = 1; model.factors.B = 2; model.factors.C = 3; model.factors.S = 4;
model.factors.D = {1, @(U,task) struct_matvec(U, task, P1, eye(R(1)))};
model.factors.E = {2, @(V,task) struct_matvec(V, task, P2, eye(R(2)))};
model.factors.F = {3, @(W,task) struct_matvec(W, task, Pm, eye(R(3)))};
model.factorizations{1}.data = HSI;
model.factorizations{1}.weight = 2;
model.factorizations{1}.lmlra  = {'D', 'E', 'C','S'};
model.factorizations{2}.data = MSI;
model.factorizations{2}.weight = 2;
model.factorizations{2}.lmlra  = {'A', 'B', 'F', 'S'};
[sol, output] = sdf_nls(model, 'MaxIter', opti_Niter);

U = sol.variables{1}; V = sol.variables{2}; W = sol.variables{3}; S = sol.variables{4};
SRI_hat = lmlragen({U,V,W},S); 
info.factors = {U,V,W};
info.core = {S};
info.rank = {ranks};
info.opti = {output};

end

