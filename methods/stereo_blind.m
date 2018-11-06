function [SRI_hat, info] = stereo_blind(HSI, MSI, Pm, ranks, opts)

% STEREO_BLIND provides estimation of SRI with an AO algorithm
% [A,B,C, SRI_hat] = STEREO_blind(lambda, F, A0,B0,C0,A_tilde0, B_tilde0, HSI, MSI,Pm, Niter) returns
% SRI_hat = [A,B,C] from HSI and MSI
% 
% INPUT ARGUMENTS:
%     lambda: weighting argument between HSI and MSI
%     F: tensor rank of estimation
%     A0,B0,C0,A_tilde0, B_tilde0: factor matrices obtained from initialization TENREC_BLIND
%     HSI, MSI: lower-resolution images
%     SRI: super-resolution image
%     Pm: spectral degradation matrix
%     Niter: number of iterations (optional)
% OUTPUT ARGUMENTS:
%     A,B,C: factor matrices for CPD of estimation
%     SRI_hat: estimation of SRI such that SRI_hat = [A,B,C]
%     err: structure containing NMSE, SAM, ERGAS, R-SNR, CPU TIME
% 
% SEE ALSO: TENREC_BLIND, STEREO

F = ranks;

if nargin==4
    lambda = 1; Niter = 10; opts.CPD_Niter = 25; Display = 'false';
    [A, B, A_tilde, B_tilde, C] = stereo_init_blind(MSI,HSI, ranks, Pm);
elseif nargin==5
    lambda = opts.lambda; Niter = opts.Niter; Display = opts.Display;
    A = opts.factors{1}; B = opts.factors{2}; C = opts.factors{3}; 
    A_tilde = opts.factors{4}; B_tilde = opts.factors{5};
end

options.POSDEF = true; options.SYM = true;

Yh1 = tens2mat(HSI,[],1); Yh2 = tens2mat(HSI,[],2); Yh3 = tens2mat(HSI,[],3);
Ym1 = tens2mat(MSI,[],1); Ym2 = tens2mat(MSI,[],2); Ym3 = tens2mat(MSI,[],3);

for n = 1:Niter

    disp('C...')
    mat1 = kr(B,A); mat2 = kr(B_tilde,A_tilde);
    R = lambda*Pm'*Ym3'*mat1 + Yh3'*mat2;
    Q = lambda*kron((B'*B).*(A'*A), Pm'*Pm) + kron((B_tilde'*B_tilde).*(A_tilde'*A_tilde), eye(size(Pm,2)));
    C = reshape(linsolve(Q, R(:), options),[size(Pm,2) F]);
    
    disp('A tilde...')
    R = Yh1'*kr(C,B_tilde);
    Q = kron((C'*C).*(B_tilde'*B_tilde), eye(size(A_tilde,1)));
    A_tilde = reshape(linsolve(Q, R(:), options),[size(HSI,1) F]); 
    
    disp('B tilde...')
    R = Yh2'*kr(C,A_tilde);
    Q = kron((C'*C).*(A_tilde'*A_tilde), eye(size(B_tilde,1)));
    B_tilde = reshape(linsolve(Q, R(:), options),[size(HSI,2) F]);
    
    disp('A...')
    mat1 = Pm*C;
    R = Ym1'*kr(mat1,B);
    Q = kron((mat1'*mat1).*(B'*B), eye(size(MSI,1))); 
    A = reshape(linsolve(Q, R(:), options),[size(MSI,1) F]);
    
    disp('B...')
    mat1 = Pm*C;
    R = Ym2'*kr(mat1,A);
    Q = kron((mat1'*mat1).*(A'*A), eye(size(B,1)));
    B = reshape(linsolve(Q, R(:), options),[size(MSI,2) F]);
     

end

SRI_hat = cpdgen({A,B,C});
info.factors = {A,B,C,A_tilde, B_tilde};
info.rank = {ranks};
info.Niter = {Niter};



end

