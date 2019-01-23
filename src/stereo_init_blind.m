function [A, B, A_tilde, B_tilde, C] = stereo_init_blind(MSI,HSI, ranks, Pm)

% TENREC initialization of STEREO_BLIND
% [A_tilde, B_tilde, C] = TENREC_BLIND(DATA, F, Pm) performs CPD of DATA
% in the case where spatial degradation factors are unknown
% 
% INPUT ARGUMENTS:
%     DATA: dataset on which CPD is performed
%     HSI: Hyperspectral image
%     F: Tensor rank of MSI rank decomposition
%     Pm: spectral degradation matrix
% OUTPUT ARGUMENTS:
%     A: matrix of size ImxF
%     B: matrix of size JmxF
%     A_tilde: matrix of size IhxF
%     B_tilde: matrix of size JhxF
%     C: matrix of size KhxF
% such that HSI = [A_tilde, B_tilde, C] and MSI = [A,B,Pm*C]
% 
% SEE ALSO: TENREC
F = ranks;

options.MaxIter = 25;
U = cpd(MSI,F,options);
A = cell2mat(U(1)); B = cell2mat(U(2)); C_tilde = cell2mat(U(3));

d = 4; Im = size(MSI,1); Ih = Im/d; 
Jm = size(MSI,2); Jh = Jm/d;
A_tilde = zeros(Ih,F); B_tilde = zeros(Jh,F);

for i=1:Ih
    for k=d*(i-1)+1:d*i
        A_tilde(i,:) = A_tilde(i,:) + A(k,:);
    end
end
for j=1:Jh
    for k=d*(j-1)+1:d*j
        B_tilde(j,:) = B_tilde(j,:) + B(k,:);
    end
end

[Q,R] = qr(kr(B_tilde,A_tilde),0);
b = [reshape(tens2mat(HSI,[],3)'*Q,[],1); reshape(C_tilde,[],1)];
X = [kron(R,eye(size(HSI,3))); kron(eye(F),Pm)];
C = reshape(b\X,[size(HSI,3),F]);


end

