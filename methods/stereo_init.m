function [A, B, A_tilde, B_tilde, C] = stereo_init(MSI, HSI, P1, P2, Pm, ranks)

% STEREO_INIT initialization of STEREO
% [A, B, A_tilde, B_tilde, C] = STEREO_INIT((MSI, HSI, P1, P2, Pm, ranks, opts)
% performs CPD of DATA 
% in the case where spatial degradation factors are known
% 
% INPUT ARGUMENTS:
%     MSI: multispectral image
%     HSI: hyperspectral image
%     ranks: Tensor rank of MSI rank decomposition
%     P1,P2,Pm: spectral and spatial degradation matrices 
%     opts: options structure
% OUTPUT ARGUMENTS:
%     A: matrix of size ImxF
%     B: matrix of size JmxF
%     A_tilde: matrix of size IhxF
%     B_tilde: matrix of size JhxF
%     C: matrix of size KhxF
% such that HSI = [A_tilde, B_tilde, C] and MSI = [A,B,Pm*C]
% 
% Copyright (c) 2018 Clemence Prevost, Konstantin Usevich, Pierre Comon, David Brie
% https://github.com/cprevost4/HSR_Tucker
% Contact: clemence.prevost@univ-lorraine.fr


options.MaxIter = 25; %options.Display = true;
U = cpd(MSI,ranks,options);
A = cell2mat(U(1)); B = cell2mat(U(2)); C_tilde = cell2mat(U(3));

A_tilde = P1*A; B_tilde = P2*B; 

%C = solveC_qr(B_tilde, A_tilde, C_tilde, 200, HSI, F, Pm);
C = solveC_normal(A_tilde, B_tilde, C_tilde, Pm, HSI);

end

