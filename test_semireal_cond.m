%%  Load data
SRI = cell2mat(struct2cell(load('SalinasA.mat')));
SRI(:,:,[108:112 154:167 224]) = []; %Regions of water absorption (Salinas)
SRI = crop(SRI,[80,84,size(SRI,3)]);
% 2. degradation
Pm = spectral_deg(SRI,"LANDSAT");
d1 = 4; d2 = 4; q = 9;
[P1,P2] = spatial_deg(SRI, q, d1, d2);
MSI = tmprod(SRI,Pm,3); HSI = tmprod(tmprod(SRI,P1,1),P2,2);

%% Condition numbers

figure
conds1 = zeros(size(P1,1),1);
for i=1:size(conds1,1) 
  matp = P1 * U(:,1:i);
  conds1(i) = cond(matp);
end
plot(conds1)

conds2 = zeros(size(P2,1),1);
for i=1:size(conds2,1)
  matp = P2 * V(:,1:i);
  conds2(i) = cond(matp);
end
hold on
plot(conds2); legend('Cond of P_1U','Cond of P_2V')

figure % total condition of Kronecker product
plot(conds1(1:20).*conds2(1:20)); legend('Cond of Kron. product')



conds = zeros(size(Pm,1),1);
for i=1:size(conds,1), 
  matp = Pm * W(:,1:i);
  conds(i) = cond(matp);
end
figure
plot(conds); legend('Total cond')


%% FIGURE 14 FOR SALINAS-A

R3 = size(MSI,3); cond_XtX = [];
for R1=1:size(HSI,1)
    [U, ~, ~] = svds(tens2mat(MSI,1,[]),R1);
    [V, ~, ~] = svds(tens2mat(MSI,2,[]),R1);
    [W, ~, ~] = svds(tens2mat(HSI,3,[]), R3);
    

    cond_XtX(R1) = (svds(Pm*W,1,'largest')^2 + svds(P1*U,1,'largest')^2*svds(P2*V,1,'largest')^2)....
        /(svds(P1*U,1,'smallest')^2*svds(P2*V,1,'smallest')^2);
    
end

figure
plot(cond_XtX); xlim([1 25]); xlabel('R_1=R_2');  

