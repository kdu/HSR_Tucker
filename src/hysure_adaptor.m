function [SRI_hat,info] = hysure_adaptor(HSI, MSI, Pm, R, opts)
%Changes the order of parameters 
  %[SRI_hat,info] = hysure(HSI,MSI, opts);
%   SRI = cell2mat(struct2cell(load('PaviaU.mat')));
%   SRI(1:2,:,:) = []; SRI(:,1:4,:) = [];
load Cuprite.mat
SRI = double(X); clear X; %Convert groundtruth to double
SRI(:,:,[1:2 104:113 148:167 221:224]) = [];
  [P1,P2] = spatial_deg(SRI, 9, 4,4);
  [SRI_hat,info] = hysure2(HSI,MSI,P1,P2,Pm, opts);
end

