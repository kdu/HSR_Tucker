% EXPERIMENTS ON WOOD DATA %
% Copyright (c) 2018 Clemence Prevost, Konstantin Usevich, Pierre Comon, David Brie
% https://github.com/cprevost4/HSR_Tucker
% Contact: clemence.prevost@univ-lorraine.fr

%% Load wood data

% 1. load data
SRI = cell2mat(struct2cell(load('SRI1.mat')));
SRI = permute(SRI,[2,3,1]);
SRI = SRI(:,3:end,:);
HSI = cell2mat(struct2cell(load('HSID2.mat'))); %Hyperspectral image, d=2
HSI = permute(HSI,[2,3,1]);
MSI = cell2mat(struct2cell(load('MSID8.mat'))); %Multispectral image, d=2
MSI = permute(MSI,[2,3,1]); 
% 2. degradation
Pm = zeros(size(MSI,3),size(SRI,3)); d = size(SRI,3)/size(MSI,3);
for i=1:size(MSI,3)
    Pm(i, d*(i-1)+1:d*i) = 1/d;
end
P1 = zeros(size(HSI,1),size(SRI,1)); d1 = size(SRI,1)/size(HSI,1);
P2 = zeros(size(HSI,2),size(SRI,2)); d2 = size(SRI,2)/size(HSI,2);
for i=1:size(P1,1)
        P1(i,1+d1*(i-1)) = 1;
end
for j=1:size(P2,1)
        P2(j,1+d2*(j-1)) = 1;
end

%% Make table 

methods = {'STEREO' 'stereo3' '20' []; ...
           'STEREO' 'stereo3' '30' []; ...
           'STEREO' 'stereo3' '50' []; ...
           'STEREO' 'stereo3' '100' []; ...
           'SCOTT' 'scott' '[50,50,8]' [];...
           'SCOTT' 'scott' '[30,30,15]' [];...
           'HySure' 'hysure_b1_adaptor' '[]' 8;};      
DegMat = struct('Pm', Pm, 'P1', P1, 'P2', P2);         
res = compare_methods(SRI, HSI, MSI, DegMat, [d1 d2], methods); 
T = cell2mat(res(:,2:end)); save('exp3_table2_wood.txt','T','-ASCII')

%% Make figure

methods = {'STEREO' 'stereo3' '50' []; ...
           'SCOTT' 'scott' '[50,50,8]' [];...
           'HySure' 'hysure_b1_adaptor' '[]' 8;};      
DegMat = struct('Pm', Pm, 'P1', P1, 'P2', P2);         
[~,est] = compare_methods(SRI, HSI, MSI, DegMat, [d1 d2], methods); 

figure
subplot(2,2,1); imagesc(SRI(:,:,160)); title('Groundtruth SRI')
subplot(2,2,2); imagesc(real(est{1}(:,:,160))); title(sprintf('%s, %s',methods{1,1},methods{1,3}))
subplot(2,2,3); imagesc(real(est{2}(:,:,160))); title(sprintf('%s, %s',methods{2,1},methods{2,3}))
subplot(2,2,4); imagesc(real(est{3}(:,:,160))); title(sprintf('%s',methods{3,1}))




