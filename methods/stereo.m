function [SRI_hat, info] = stereo(HSI, MSI, P1, P2, Pm, ranks, opts)

% STEREO provides estimation of SRI with an AO algorithm
% [SRI_hat, info] = STEREO(HSI, MSI, P1, P2, Pm, ranks, opts) returns
% SRI_hat = [A,B,C] from HSI and MSI
% 
% INPUT ARGUMENTS:
%     F: tensor rank of estimation
%     HSI, MSI: lower-resolution images
%     P1,P2,Pm: degradation matrices
%     opts: structure containing parameters
% OUTPUT ARGUMENTS:
%     SRI_hat: estimation of SRI such that SRI_hat = [A,B,C]
%     info: informative structure about the method
% 
% SEE ALSO: TENREC
% Copyright (c) 2018 Clemence Prevost, Konstantin Usevich, Pierre Comon, David Brie
% https://github.com/cprevost4/HSR_Tucker
% Contact: clemence.prevost@univ-lorraine.fr

if nargin==6
    lambda = 1; Niter = 10; opts.CPD_Niter = 25; Display = 'false';
    [A, B, ~,~, C] = tenRec(MSI, HSI, P1, P2, Pm, ranks, opts);
elseif nargin==7
    lambda = opts.lambda; Niter = opts.Niter; Display = opts.Display;
    A = opts.factors{1}; B = opts.factors{2}; C = opts.factors{3};
end

opts2.POSDEF = true; opts2.SYM = true;

Yh1 = tens2mat(HSI,[],1); Yh2 = tens2mat(HSI,[],2); Yh3 = tens2mat(HSI,[],3);
Ym1 = tens2mat(MSI,[],1); Ym2 = tens2mat(MSI,[],2); Ym3 = tens2mat(MSI,[],3);

for n = 1:Niter
    
    disp('A...')
    R = lambda*Ym1'*kr(Pm*C,B) + P1'*Yh1'*kr(C,P2*B);
    Q = lambda*kron(((Pm*C)'*(Pm*C)).*(B'*B), eye(size(P1,2))) + kron((C'*C).*((P2*B)'*(P2*B)), P1'*P1);
    A = reshape(linsolve(Q, R(:), opts2),[size(P1,2) ranks]);
    %A = reshape(Q\R(:),[size(P1,2) ranks]);
    
    disp('B...')
    R = lambda*Ym2'*kr(Pm*C,A) + P2'*Yh2'*kr(C,P1*A);
    Q = lambda*kron(((Pm*C)'*(Pm*C)).*(A'*A), eye(size(P2,2))) + kron((C'*C).*((P1*A)'*(P1*A)), P2'*P2);
    B = reshape(linsolve(Q, R(:), opts2),[size(P2,2) ranks]);
    %B = reshape(Q\R(:),[size(P2,2) ranks]);
    
    disp('C...')
    R = lambda*Pm'*Ym3'*kr(B,A) + Yh3'*kr(P2*B,P1*A);
    Q = lambda*kron((B'*B).*(A'*A), Pm'*Pm) + kron(((P2*B)'*(P2*B)).*((P1*A)'*(P1*A)), eye(size(Pm,2)));
    C = reshape(linsolve(Q, R(:), opts2),[size(Pm,2) ranks]);
    %C = reshape(Q\R(:),[size(Pm,2) ranks]);
    
end

SRI_hat = cpdgen({A,B,C});
info.factors = {'A','B','C'};
info.rank = {'ranks'};
info.Niter = {'Niter'};

end

