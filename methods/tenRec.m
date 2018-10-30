function [A, B, A_tilde, B_tilde, C] = tenRec(DATA, HSI, F, P1,P2,Pm)

% TENREC initialization of STEREO
% [A_tilde, B_tilde, C] = TENREC(DATA, F, P1,P2,Pm) performs CPD of DATA 
% in the case where spatial degradation factors are known
% 
% INPUT ARGUMENTS:
%     DATA: dataset on which CPD is performed
%     HSI: hyperspectral image
%     F: Tensor rank of MSI rank decomposition
%     P1,P2,Pm: spectral and spatial degradation matrices 
% OUTPUT ARGUMENTS:
%     A: matrix of size ImxF
%     B: matrix of size JmxF
%     A_tilde: matrix of size IhxF
%     B_tilde: matrix of size JhxF
%     C: matrix of size KhxF
% such that HSI = [A_tilde, B_tilde, C] and MSI = [A,B,Pm*C]
% 
% Copyright (c) 2018 Clémence Prévost, Konstantin Usevich, Pierre Comon, David Brie
% https://github.com/cprevost4/HSR_Tucker


options.MaxIter = 25; %options.Display = true;
U = cpd(DATA,F,options);
A = cell2mat(U(1)); B = cell2mat(U(2)); C_tilde = cell2mat(U(3));

A_tilde = P1*A; B_tilde = P2*B; 

%C = solveC_qr(B_tilde, A_tilde, C_tilde, 200, HSI, F, Pm);
C = solveC_normal(A_tilde, B_tilde, C_tilde, Pm, HSI);

end

