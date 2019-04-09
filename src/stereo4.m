function [SRI_hat, err] = stereo4( HSI, MSI, P1,P2,Pm, ranks, opts)

if ~exist('opts','var')
    opts = struct();
end
if ~isfield(opts,'factors') || isempty(opts.factors)
    [A, B, ~,~, C] = stereo_init(MSI, HSI, P1, P2, Pm, ranks);
    %[A,B,C,~,~,~]=TenRec(MSI,tens2mat(HSI,[],3),opts.CPD_Niter,ranks,P1,P2);
else
    A = opts.factors{1}; B = opts.factors{2}; C = opts.factors{3};
end

lamda=1;
s_iter=10; % same as us
[ A1,B1,C1,~ ] = STEREO( HSI,MSI,P1,P2,Pm,s_iter,lamda,A,B,C,Pm*C);

SRI_hat = cpdgen({A1,B1,C1});
err = {};


end

