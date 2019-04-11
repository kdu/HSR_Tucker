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

%P1 = sparse(P1); P2 = sparse(P2); Pm = sparse(Pm);

for n=1:opts.Niter
    %n 
    
    %disp('A...')
    %tic;
    C_tilde = Pm*C; B_tilde = P2*B;
    mat1 = (C'*C).*(B_tilde'*B_tilde);
    mat2 = opts.lambda* ((C_tilde'*C_tilde).*(B'*B));
    Z = opts.lambda*Ym1'*kr(C_tilde,B) + P1'*Yh1'*kr(C,B_tilde);
    %A = bartelsStewart(P1'*P1, mat1', [], mat2',Z);
    A = bartelsStewart(mat2, [], mat1, P1'*P1,Z'); A = A'; %toc;
    
    %disp('B...')
    %tic;
    A_tilde = P1*A;
     mat1 = (C'*C).*(A_tilde'*A_tilde);
     mat2 = opts.lambda* ((C_tilde'*C_tilde).*(A'*A));
     Z = opts.lambda*Ym2'*kr(C_tilde,A) + P2'*Yh2'*kr(C,A_tilde);
     %B = bartelsStewart(P2'*P2, mat1', [], mat2',Z);
     B = bartelsStewart(mat2, [], mat1, P2'*P2,Z'); B = B'; %toc;
     
    %disp('C...')
    B_tilde = P2*B;
    %disp('Computing KR products'); tic;
    mat1 = (B'*B).*(A'*A);
    mat2 = (B_tilde'*B_tilde).*(A_tilde'*A_tilde); %toc
    %disp('Computing Z'); tic;
    Z = opts.lambda*Pm'*Ym3'*kr(B,A) + Yh3'*kr(B_tilde,A_tilde); %toc
    %C = bartelsStewart(lambda*(Pm'*Pm), mat1', [], mat2',Z);
    %disp('Solving equation'); tic;
    C = bartelsStewart(mat2, [], mat1, opts.lambda*Pm'*Pm,Z'); C = C'; %toc
    %C = bartelsStewart(eye(size(HSI,3)), mat2, opts.lambda*Pm'*Pm, mat1, Z); %toc
     

end

       SRI_hat = cpdgen({A,B,C});
%err = {mse, s, glob_err, snr, crossco, cost};
err = {};

end

