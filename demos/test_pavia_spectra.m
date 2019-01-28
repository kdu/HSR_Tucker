% EXPERIMENT 1- SPECTRAL SIGNATURE FOR PAVIA %
% Copyright (c) 2018 Clemence Prevost, Konstantin Usevich, Pierre Comon, David Brie
% https://github.com/cprevost4/HSR_Tucker
% Contact: clemence.prevost@univ-lorraine.fr

%% Load Pavia + groundtruth

load PaviaU.mat
SRI = paviaU; clear paviaU; SRI(1:2,:,:) = []; SRI(:,1:4,:) = [];
d1 = 4; d2 = 4; %Spatial downsampling ratios
[P1,P2] = spatial_deg(SRI,7,d1,d2); Pm = spectral_deg(SRI,"Quickbird");
%Generate HSI and MSI
HSI = tmprod(SRI,{P1,P2},[1,2]); MSI = tmprod(SRI,Pm,3);
for k=1:size(HSI,3)
    HSI(:,:,k) = awgn(HSI(:,:,k),15);
end
for k=1:size(MSI,3)
    MSI(:,:,k) = awgn(MSI(:,:,k),25);
end
load paviaU_gt.mat
G = double(paviaU_gt); clear paviaU_gt; G(1:2,:,:) = []; G(:,1:4,:) = [];

%% Run algorithms

methods2 = {'SCUBA'  'scuba_adaptor' '[120,3]' [4,4]; ...
           'B-SCOTT' 'bscott_adaptor' '[size(HSI,1),size(HSI,2),3]' [4,4]; ...                 
           'SCUBA'  'scuba_adaptor' '[300,9]' [1,1]; ... 
           'B-SCOTT' 'bscott_adaptor' '[size(HSI,1),size(HSI,2),size(MSI,3)]' [1,1]; ...
           'SCUBA'  'scuba_adaptor' '[120,3]' [8,8]; ... 
           'B-SCOTT' 'bscott_adaptor' '[size(MSI,1)/8,size(MSI,2)/8,3]' [8,8];};      
[~,est] = compare_methods(SRI, HSI, MSI, struct('Pm', Pm), [d1 d2], methods2);

%% Reconstruction of spectral signatures - material 2

X = tens2mat(SRI,3,[]);
X2 = tens2mat(est{1},3,[]); X3 = tens2mat(est{2},3,[]);
X4 = tens2mat(est{3},3,[]); X5 = tens2mat(est{4},3,[]);
X6 = tens2mat(est{5},3,[]); X7 = tens2mat(est{6},3,[]);

ind = find(G(:) == 2);
spec = mean(X(:,ind),2); 
spec2 = mean(X2(:,ind),2);  spec3 = mean(X3(:,ind),2); 
spec4 = mean(X4(:,ind),2);  spec5 = mean(X5(:,ind),2); 
spec6 = mean(X6(:,ind),2);  spec7 = mean(X7(:,ind),2); 

figure
subplot(1,3,1); plot(spec,'Linewidth',1); hold on; plot(spec2,'Linewidth',1); hold on
plot(spec3,'Linewidth',1); title('(4,4) blocks'); legend('Groundtruth','SCUBA','BSCOTT');
xlabel('Spectral bands'); ylabel('Reflectance')
subplot(1,3,2); plot(spec,'Linewidth',1); hold on; plot(spec4,'Linewidth',1); hold on
plot(spec5,'Linewidth',1); title('(1,1) blocks'); legend('Groundtruth','SCUBA','BSCOTT')
xlabel('Spectral bands'); ylabel('Reflectance')
subplot(1,3,3); plot(spec,'Linewidth',1); hold on; plot(spec6,'Linewidth',1); hold on
plot(spec7,'Linewidth',1); title('(8,8) blocks'); legend('Groundtruth','SCUBA','BSCOTT')
xlabel('Spectral bands'); ylabel('Reflectance')

%% Reconstruction of spectral signatures - material 9

X = tens2mat(SRI,3,[]);
X2 = tens2mat(est{1},3,[]); X3 = tens2mat(est{2},3,[]);
X4 = tens2mat(est{3},3,[]); X5 = tens2mat(est{4},3,[]);
X6 = tens2mat(est{5},3,[]); X7 = tens2mat(est{6},3,[]);

ind = find(G(:) == 9);
spec = mean(X(:,ind),2); 
spec2 = mean(X2(:,ind),2);  spec3 = mean(X3(:,ind),2); 
spec4 = mean(X4(:,ind),2);  spec5 = mean(X5(:,ind),2); 
spec6 = mean(X6(:,ind),2);  spec7 = mean(X7(:,ind),2); 

figure
subplot(1,3,1); plot(spec,'Linewidth',1); hold on; plot(spec2,'Linewidth',1); hold on
plot(spec3,'Linewidth',1); title('(4,4) blocks'); legend('Groundtruth','SCUBA','BSCOTT');
xlabel('Spectral bands'); ylabel('Reflectance')
subplot(1,3,2); plot(spec,'Linewidth',1); hold on; plot(spec4,'Linewidth',1); hold on
plot(spec5,'Linewidth',1); title('(1,1) blocks'); legend('Groundtruth','SCUBA','BSCOTT')
xlabel('Spectral bands'); ylabel('Reflectance')
subplot(1,3,3); plot(spec,'Linewidth',1); hold on; plot(spec6,'Linewidth',1); hold on
plot(spec7,'Linewidth',1); title('(8,8) blocks'); legend('Groundtruth','SCUBA','BSCOTT')
xlabel('Spectral bands'); ylabel('Reflectance')