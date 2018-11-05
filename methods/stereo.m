function [A,B,C, SRI_hat, err, cost] = stereo(lambda, F, A0,B0,C0,SRI, HSI, MSI, P1,P2,Pm, Niter)

% STEREO provides estimation of SRI with an AO algorithm
% [A,B,C, SRI_hat] = STEREO(lambda, F, A0,B0,C0, HSI, MSI, P1,P2,Pm,Niter) returns
% SRI_hat = [A,B,C] from HSI and MSI
% 
% INPUT ARGUMENTS:
%     lambda: weighting argument between HSI and MSI
%     F: tensor rank of estimation
%     A0,B0,C0: factor matrices obtained from initialization TENREC
%     HSI, MSI: lower-resolution images
%     P1,P2,Pm: degradation matrices
%     Niter: number of iterations (optional)
% OUTPUT ARGUMENTS:
%     A,B,C: factor matrices for CPD of estimation
%     SRI_hat: estimation of SRI such that SRI_hat = [A,B,C]
%     err: structure containing NMSE, SAM, ERGAS, R-SNR, CC, CPU TIME
% 
% SEE ALSO: TENREC
% Copyright (c) 2018 Clemence Prevost, Konstantin Usevich, Pierre Comon, David Brie
% https://github.com/cprevost4/HSR_Tucker
% Contact: clemence.prevost@univ-lorraine.fr


A = A0; B = B0; C = C0; opts.POSDEF = true; opts.SYM = true;

Yh1 = tens2mat(HSI,[],1); Yh2 = tens2mat(HSI,[],2); Yh3 = tens2mat(HSI,[],3);
Ym1 = tens2mat(MSI,[],1); Ym2 = tens2mat(MSI,[],2); Ym3 = tens2mat(MSI,[],3);

cost_A = []; cost_B = []; cost_C = [];

for n = 1:Niter
    
    disp('A...')
    R = lambda*Ym1'*kr(Pm*C,B) + P1'*Yh1'*kr(C,P2*B);
    Q = lambda*kron(((Pm*C)'*(Pm*C)).*(B'*B), eye(size(P1,2))) + kron((C'*C).*((P2*B)'*(P2*B)), P1'*P1);
    A = reshape(linsolve(Q, R(:), opts),[size(P1,2) F]);
    %A = reshape(Q\R(:),[size(P1,2) F]);
    
    disp('B...')
    
    R = lambda*Ym2'*kr(Pm*C,A) + P2'*Yh2'*kr(C,P1*A);
    Q = lambda*kron(((Pm*C)'*(Pm*C)).*(A'*A), eye(size(P2,2))) + kron((C'*C).*((P1*A)'*(P1*A)), P2'*P2);
    B = reshape(linsolve(Q, R(:), opts),[size(P2,2) F]);
    %B = reshape(Q\R(:),[size(P2,2) F]);
    
    disp('C...')
    R = lambda*Pm'*Ym3'*kr(B,A) + Yh3'*kr(P2*B,P1*A);
    Q = lambda*kron((B'*B).*(A'*A), Pm'*Pm) + kron(((P2*B)'*(P2*B)).*((P1*A)'*(P1*A)), eye(size(Pm,2)));
    C = reshape(linsolve(Q, R(:), opts),[size(Pm,2) F]);
    %C = reshape(Q\R(:),[size(Pm,2) F]);
    
    SRI_hat = cpdgen({A,B,C});
    mse = nmse(SRI, SRI_hat);
s = sam(SRI, SRI_hat);
glob_err = ergas(SRI, SRI_hat,1/8);
snr = r_snr(SRI, SRI_hat);
crossco = cc(SRI, SRI_hat);
    
end


err = {mse, s, glob_err, snr, crossco}; cost = {cost_A, cost_B, cost_C};


end

