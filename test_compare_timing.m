%% compare time between our and Kana's implementation

%% Load Indian Pines

SRI = cell2mat(struct2cell(load('Indian_pines.mat')));
SRI(:,:,[104:108 150:163 220]) = []; %Regions of water absorption
SRI(1,:,:) = []; SRI(:,1,:) = [];
Pm = spectral_deg(SRI,"LANDSAT");
MSI = tmprod(SRI,Pm,3);
d1 = 4; d2 = 4; q = 9;
[P1,P2] = spatial_deg(SRI, q, d1, d2);
HSI = tmprod(tmprod(SRI,P1,1),P2,2);

%% Initialization

F = 50:50:400;  t1 = []; t2 = [];
for k=1:length(F)
    tic; [A, B, ~,~, C] = stereo_init(MSI, HSI, P1, P2, Pm, F(k)); t1(k) = toc;
    tic; [A,B,C, A_tilde,B_tilde,C_tilde]=TenRec(MSI,tens2mat(HSI,[],3),cpd_iter,F(k),P1,P2); t2(k) = toc;
end

figure
plot(F,t1); hold on; plot(F,t2); legend('Our','Kana')

%% Full algo

F = 50:50:400; cpd_iter = 25;  MAXIT = 1; t1 = []; t2 = []; lamda = 1;
for k=1:length(F)
    [A,B,C, ~,~,C_tilde]=TenRec(MSI,tens2mat(HSI,[],3),cpd_iter,F(k),P1,P2); 
    opts.factors = {A,B,C}; opts.Niter = 1;
    tic; [ A,B,C,cost ] = STEREO( HSI,MSI,P1,P2,Pm,MAXIT,lamda,A,B,C,C_tilde ); t1(k) = toc;
    tic; [SRI_hat, err] = stereo3( HSI, MSI, P1,P2,Pm, F(k), opts); t2(k) = toc;
end

figure
plot(F,t1); hold on; plot(F,t2); legend('Kana','Our')

%% Kana's iterations

F = 50; cpd_iter = 25; lamda = 1;

[Ih,Jh,K]=size(HSI);
[I,J,Km]=size(MSI);
H1=reshape(HSI,[Ih,Jh*K])';
H2=reshape(permute(HSI,[2 1 3]),[Jh,K*Ih])';
H3=reshape(HSI,[Ih*Jh,K]);
nH=norm(H3,'fro');
M1=reshape(MSI,[I,J*Km])';
M2=reshape(permute(MSI,[2 1 3]),[J,Km*I])';
M3=reshape(MSI,[I*J,Km]);
nM=norm(M3,'fro');

[A,B,C, A_tilde,B_tilde,C_tilde]=TenRec(MSI,tens2mat(HSI,[],3),cpd_iter,F,P1,P2);

    B_tilde=P2*B;

    tic;
    temp1h=khatri_rao(C,B_tilde);
    temp1m=khatri_rao(C_tilde,B);
    Kp=(C_tilde'*C_tilde).*(B'*B);
    K= (C'*C).*(B_tilde'*B_tilde);
    inv_K=pinv(K);
    As=lamda*(P1'*P1);
    Bs=Kp*inv_K;
    Cs=(lamda*P1'*H1'*temp1h+M1'*temp1m)*inv_K;
    A=sylvester(full(As),Bs,Cs);
    A_tilde=P1*A; tA1 = toc
    
%     tic;
%     temp2h=khatri_rao(C,A_tilde);
%     temp2m=khatri_rao(C_tilde,A);
%     Kp=(C_tilde'*C_tilde).*(A'*A);
%     K= (C'*C).*(A_tilde'*A_tilde);
%     inv_K=pinv(K);
%     As=lamda*(P2'*P2);
%     Bs=Kp*inv_K;
%     Cs=(lamda*P2'*H2'*temp2h+M2'*temp2m)*inv_K;
%     B=sylvester(full(As),Bs,Cs);
%     B_tilde=P2*B; tB1 = toc
%      
%     tic;
%     temp3h=khatri_rao(B_tilde,A_tilde);
%     temp3m=khatri_rao(B,A);
%     Kp=lamda*(B_tilde'*B_tilde).*(A_tilde'*A_tilde);
%     K= (B'*B).*(A'*A);
%     inv_K=pinv(K);
%     As=Pm'*Pm;
%     Bs=Kp*inv_K;
%     Cs=(lamda*H3'*temp3h+Pm'*M3'*temp3m)*inv_K;
%     C=sylvester(full(As),Bs,Cs);
%     C_tilde=Pm*C; tC1 = toc
    
    %% Our iterations
    
    F = 50; cpd_iter = 25; lambda = 1;
    [A,B,C, ~,~,~]=TenRec(MSI,tens2mat(HSI,[],3),cpd_iter,F,P1,P2);
    
    Yh1 = tens2mat(HSI,[],1); Yh2 = tens2mat(HSI,[],2); Yh3 = tens2mat(HSI,[],3);
    Ym1 = tens2mat(MSI,[],1); Ym2 = tens2mat(MSI,[],2); Ym3 = tens2mat(MSI,[],3);

    tic;
    mat1 = (C'*C).*((P2*B)'*(P2*B));
    mat2 = lambda* (((Pm*C)'*(Pm*C)).*(B'*B));
    Z = lambda*Ym1'*kr(Pm*C,B) + P1'*Yh1'*kr(C,P2*B);
    A = bartelsStewart(mat2, [], mat1, P1'*P1,Z'); A = A'; tA2 = toc
    
%     tic;
%      mat1 = (C'*C).*((P1*A)'*(P1*A));
%      mat2 = lambda* (((Pm*C)'*(Pm*C)).*(A'*A));
%      Z = lambda*Ym2'*kr(Pm*C,A) + P2'*Yh2'*kr(C,P1*A);
%      B = bartelsStewart(mat2, [], mat1, P2'*P2,Z'); B = B'; tB2 = toc
%      
%      tic;
%     mat1 = (B'*B).*(A'*A);
%     mat2 = ((P2*B)'*(P2*B)).*((P1*A)'*(P1*A));
%     Z = lambda*Pm'*Ym3'*kr(B,A) + Yh3'*kr(P2*B,P1*A);
%     C = bartelsStewart(mat2, [], mat1, lambda*Pm'*Pm,Z'); C = C'; tC2 = toc
    
    

