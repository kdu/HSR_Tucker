function [SRI_hat,info] = run_sdf(HSI, MSI, P1, P2, Pm, ranks, opts)

% RUN_SDF runs the HOSVD algorithm for specified rank R followed by 
% TensorLab optimization
% [SRI_hat,cost, err] = RUN_SDF(MSI, HSI, SRI, ranks, options, P1,P2,Pm) returns 
% estimation of SRI, value of cost function and metrics in the cell array err
% 
% INPUT ARGUMENTS:
%     SRI, MSI, HSI: input datasets (resp. groundtruth SRI, MSI and HSI
%     ranks: specified multilinear rank
%     P1,P2,Pm: spatial and spectral degratation matrices
%     options: number of iterations, verbosity
% OUTPUT ARGUMENTS:
%     SRI_hat: estimated SRI
%     cost: value of the cost function
%     err: cell array of metrics
% Copyright (c) 2018 Clemence Prevost, Konstantin Usevich, Pierre Comon, David Brie
% https://github.com/cprevost4/HSR_Tucker
% Contact: clemence.prevost@univ-lorraine.fr

if nargin==6
    lambda = 1; Display = 'false'; alpha = 0; opti_Niter = 30;
elseif nargin==7
    lambda = opts.lambda; Display = opts.Display; alpha = opts.alpha;
    opti_Niter = opts.opti_Niter;
end

[U, ~, ~] = svds(tens2mat(MSI,1,[]),ranks(1));
[V, ~, ~] = svds(tens2mat(MSI,2,[]),ranks(2));
[W, ~, ~] = svds(tens2mat(HSI,3,[]), ranks(3));

A = kron(V'*(P2'*P2)*V, U'*(P1'*P1)*U);
B = lambda* W'*(Pm'*Pm)*W;
b_old = tmprod(HSI,{U'*P1', V'*P2', W'},[1,2,3]) + lambda * tmprod(MSI,{U', V', W'*Pm'},[1,2,3]);
C = reshape(b_old, ranks(1)*ranks(2), ranks(3));
S = reshape(sylvester((A+alpha*norm(A,2)^2*ones(size(A,2))),B,C),ranks);

model = struct;
model.variables = {U,V,W,S};
model.factors.A = 1; model.factors.B = 2; model.factors.C = 3; model.factors.S = 4;
model.factors.D = {1, @(U,task) struct_matvec(U, task, P1, eye(ranks(1)))};
model.factors.E = {2, @(V,task) struct_matvec(V, task, P2, eye(ranks(2)))};
model.factors.F = {3, @(W,task) struct_matvec(W, task, Pm, eye(ranks(3)))};
model.factorizations{1}.data = HSI;
model.factorizations{1}.weight = 2;
model.factorizations{1}.lmlra  = {'D', 'E', 'C','S'};
model.factorizations{2}.data = MSI;
model.factorizations{2}.weight = 2;
model.factorizations{2}.lmlra  = {'A', 'B', 'F', 'S'};
[sol, output] = sdf_nls(model, 'MaxIter', opti_Niter);

U = sol.variables{1}; V = sol.variables{2}; W = sol.variables{3}; S = sol.variables{4};
SRI_hat = lmlragen({U,V,W},S); 
info.factors = {'U','V','W'};
info.core = {'S'};
info.rank = {'ranks'};
info.opti = {'output'};

end

