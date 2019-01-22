%% LOAD CUPRITE DATASET

load Cuprite.mat
X = double(X); %Convert groundtruth to double

d1 = 4; d2 = 4; %Spatial downsampling ratios
[P1,P2] = spatial_deg(X,7,d1,d2); Pm = spectral_deg(X,"LANDSAT");

%Generate HSI and MSI
HSI = tmprod(X,{P1,P2},[1,2]); MSI = tmprod(X,Pm,3);
for k=1:size(HSI,3)
    HSI(:,:,k) = awgn(HSI(:,:,k),15);
end
for k=1:size(MSI,3)
    MSI(:,:,k) = awgn(MSI(:,:,k),25);
end

%% LOAD PAVIA DATASET

load PaviaU.mat
X = paviaU; clear paviaU; X(1:2,:,:) = []; X(:,1:4,:) = [];

d1 = 4; d2 = 4; %Spatial downsampling ratios
[P1,P2] = spatial_deg(X,7,d1,d2); Pm = spectral_deg(X,"Quickbird");

%Generate HSI and MSI
HSI = tmprod(X,{P1,P2},[1,2]); MSI = tmprod(X,Pm,3);
for k=1:size(HSI,3)
    HSI(:,:,k) = awgn(HSI(:,:,k),15);
end
for k=1:size(MSI,3)
    MSI(:,:,k) = awgn(MSI(:,:,k),25);
end

%% RUN ALGORITHMS

opts.Nblocks = [4,4]; opts2.Nblocks = [1,1]; opts3.Nblocks = [8,8];
%tic; [Xhat1, ~] = stereo_blind3(HSI, MSI, Pm, 300); t1 = toc;
tic; [Xhat2,~] = scuba(MSI,HSI,Pm, [120,3], opts); t2 = toc;
tic; [Xhat3,~] = bscott(MSI,HSI,Pm, [size(HSI,1),size(HSI,2),3], opts); t3 = toc;

tic; [Xhat4,~] = scuba(MSI,HSI,Pm, [300,9], opts2); t4 = toc;
tic; [Xhat5,~] = bscott(MSI,HSI,Pm, [size(HSI,1),size(HSI,2),4], opts2); t5 = toc;

tic; [Xhat6,~] = scuba(MSI,HSI,Pm, [120,3], opts3); t6 = toc;
tic; [Xhat7,~] = bscott(MSI,HSI,Pm, [76,42,3], opts3); t7 = toc;

% %err1 = compute_metrics(X,Xhat1,d1,d2);
% err2 = compute_metrics(X,Xhat2,d1,d2); err2 = [cell2mat(err2) t2];
% err3 = compute_metrics(X,Xhat3,d1,d2); err3 = [cell2mat(err3) t3];
% err4 = compute_metrics(X,Xhat4,d1,d2); err4 = [cell2mat(err4) t4];
% err5 = compute_metrics(X,Xhat5,d1,d2); err5 = [cell2mat(err5) t5];
% err6 = compute_metrics(X,Xhat6,d1,d2); err6 = [cell2mat(err6) t6];
% err7 = compute_metrics(X,Xhat7,d1,d2); err7 = [cell2mat(err7) t7];
% T = [err2; err3; err4; err5; err6; err7]; save('table_BCuprite.txt','T','-ASCII')

% figure
% subplot(2,3,1); imagesc(Xhat2(:,:,44)); title('SCUBA (4,4)');
% subplot(2,3,2); imagesc(Xhat4(:,:,44)); title('SCUBA single');
% subplot(2,3,3); imagesc(Xhat6(:,:,44)); title('SCUBA (8,8)');
% subplot(2,3,4); imagesc(Xhat3(:,:,44)); title('BSCOTT (4,4)');
% subplot(2,3,5); imagesc(Xhat5(:,:,44)); title('BSCOTT single');
% subplot(2,3,6); imagesc(Xhat7(:,:,44)); title('BSCOTT (8,8)');

%% SVD FOR GUESS ON RANK

s1 = svd(tens2mat(MSI,1,[])); s2 = svd(tens2mat(MSI,2,[]));  s3 = svd(tens2mat(HSI,3,[])); 
figure
% subplot(1,3,1); plot(s1); 
% subplot(1,3,2); plot(s2);  subplot(1,3,3); plot(s3); 
%subplot(1,3,1); plot(log(s1)); subplot(1,3,2); plot(log(s2));  subplot(1,3,3); plot(log(s3)); 
subplot(1,3,1); semilogy(s1); 
subplot(1,3,2); semilogy(s2);  subplot(1,3,3); semilogy(s3);

%% RECONSTRUCTION OF SPECTRAL SIGNATURES

load paviaU_gt.mat
G = double(paviaU_gt); clear paviaU_gt; G(1:2,:,:) = []; G(:,1:4,:) = [];

mat = tens2mat(X,3,[]);
mat2 = tens2mat(Xhat2,3,[]); mat3 = tens2mat(Xhat3,3,[]);
mat4 = tens2mat(Xhat4,3,[]); mat5 = tens2mat(Xhat5,3,[]);
mat6 = tens2mat(Xhat6,3,[]); mat7 = tens2mat(Xhat7,3,[]);

for n=9:9
    ind = find(G(:) == n);
    spec = mean(mat(:,ind),2); 
    spec2 = mean(mat2(:,ind),2);  spec3 = mean(mat3(:,ind),2); 
    spec4 = mean(mat4(:,ind),2);  spec5 = mean(mat5(:,ind),2); 
    spec6 = mean(mat6(:,ind),2);  spec7 = mean(mat7(:,ind),2); 
    
    figure
    subplot(1,3,1); plot(spec,'Linewidth',1); hold on; plot(spec2,'Linewidth',1); 
    plot(spec3,'Linewidth',1); title('(4,4) blocks'); legend('GT','SCUBA','BSCOTT');
    xlabel('Spectral bands'); ylabel('Reflectance')
    subplot(1,3,2); plot(spec,'Linewidth',1); hold on; plot(spec4,'Linewidth',1);
    plot(spec5,'Linewidth',1); title('(1,1) blocks'); legend('GT','SCUBA','BSCOTT')
    xlabel('Spectral bands'); ylabel('Reflectance')
    subplot(1,3,3); plot(spec,'Linewidth',1); hold on; plot(spec6,'Linewidth',1);
    plot(spec7,'Linewidth',1); title('(8,8) blocks'); legend('GT','SCUBA','BSCOTT')
    xlabel('Spectral bands'); ylabel('Reflectance')
end






