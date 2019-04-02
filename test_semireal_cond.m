%%  Load data
SRI = cell2mat(struct2cell(load('SalinasA.mat')));
SRI(:,:,[108:112 154:167 224]) = []; %Regions of water absorption (Salinas)
SRI = crop(SRI,[80,84,size(SRI,3)]);
% 2. degradation
Pm = spectral_deg(SRI,"LANDSAT");
d1 = 4; d2 = 4; q = 9;
[P1,P2] = spatial_deg(SRI, q, d1, d2);
MSI = tmprod(SRI,Pm,3); HSI = tmprod(tmprod(SRI,P1,1),P2,2);

[U, S1, ~] = svd(tens2mat(MSI,1,[]),0);
[V, S2, ~] = svd(tens2mat(MSI,2,[]),0);

[W, S3, ~] = svd(tens2mat(HSI,3,[]),0);

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



%% Some other stuff
figure 

[Q1, R1] = qr(P1*U);
plot(abs(diag(R1)))
hold on
[Q2, R2] = qr(P2*V);
plot(abs(diag(R2)))

%%

minmax1 = zeros(20, 2)

figure
for i=1:20
  S1 = svd(P1 * U(:,1:i));
  minmax1(i,:) = [S1(1) S1(end)];
  plot(S1);
  hold on;
end
hold off;

minmax2 = zeros(21, 2)


figure
for i=1:21 
  S2 = svd(P2 * V(:,1:i));
  minmax2(i,:) = [S2(1) S2(end)];
  plot(S2);
  hold on;
end
hold off;


figure
for j=1:204 
  S3 = svd(Pm * W(:,1:j));
  plot(S3);
  minmax3(j,:) = [S3(1) S3(end)];
  hold on;
end
hold off;

conds = zeros(20,204);
%% 
for i=1:20
  for j=1:204
    conds(i,j) = (minmax3(j,1)^2 + (minmax1(i,1)*minmax2(i,1))^2) / ...
                 (minmax3(j,2)^2 + (minmax1(i,2)*minmax2(i,2))^2);
  end
end
    
figure
surfc(20:204,1:20,conds(:,20:204));
