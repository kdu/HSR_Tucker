function [SRI_hat, err] = stereo3( HSI, MSI, P1,P2,Pm, ranks, opts)

%STEREO using Bartels-Stewart

if ~exist('opts','var')
    opts = struct();
end
if ~isfield(opts,'lambda') || isempty(opts.lambda)
    opts.lambda = 1;
end
if ~isfield(opts,'Niter') || isempty(opts.Niter)
    opts.Niter = 10;
end
if ~isfield(opts,'CPD_Niter') || isempty(opts.CPD_Niter)
    opts.CPD_Niter = 25;
end
if ~isfield(opts,'factors') || isempty(opts.factors)
    [A, B, ~,~, C] = stereo_init(MSI, HSI, P1, P2, Pm, ranks);
    %[A,B,C,~,~,~]=TenRec(MSI,tens2mat(HSI,[],3),opts.CPD_Niter,ranks,P1,P2);
else
    A = opts.factors{1}; B = opts.factors{2}; C = opts.factors{3};
end

opts2.POSDEF = true; opts2.SYM = true;

Yh1 = tens2mat(HSI,[],1); Yh2 = tens2mat(HSI,[],2); Yh3 = tens2mat(HSI,[],3);
Ym1 = tens2mat(MSI,[],1); Ym2 = tens2mat(MSI,[],2); Ym3 = tens2mat(MSI,[],3);

for n=1:opts.Niter
    %n 
    
    %disp('A...')
    mat1 = (C'*C).*((P2*B)'*(P2*B));
    mat2 = opts.lambda* (((Pm*C)'*(Pm*C)).*(B'*B));
    Z = opts.lambda*Ym1'*kr(Pm*C,B) + P1'*Yh1'*kr(C,P2*B);
    %A = bartelsStewart(P1'*P1, mat1', [], mat2',Z);
    A = bartelsStewart(mat2, [], mat1, P1'*P1,Z'); A = A';
    
    %disp('B...')
     mat1 = (C'*C).*((P1*A)'*(P1*A));
     mat2 = opts.lambda* (((Pm*C)'*(Pm*C)).*(A'*A));
     Z = opts.lambda*Ym2'*kr(Pm*C,A) + P2'*Yh2'*kr(C,P1*A);
     %B = bartelsStewart(P2'*P2, mat1', [], mat2',Z);
     B = bartelsStewart(mat2, [], mat1, P2'*P2,Z'); B = B';
     
    %disp('C...')
    mat1 = (B'*B).*(A'*A);
    mat2 = ((P2*B)'*(P2*B)).*((P1*A)'*(P1*A));
    Z = opts.lambda*Pm'*Ym3'*kr(B,A) + Yh3'*kr(P2*B,P1*A);
    %C = bartelsStewart(lambda*(Pm'*Pm), mat1', [], mat2',Z);
    C = bartelsStewart(mat2, [], mat1, opts.lambda*Pm'*Pm,Z'); C = C';
    
     

end

       SRI_hat = cpdgen({A,B,C});
%err = {mse, s, glob_err, snr, crossco, cost};
err = {};

end

