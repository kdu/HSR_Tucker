% EXPERIMENT 3: COMPARISON OF METRICS WITH STEREO %
% Copyright (c) 2018 Clemence Prevost, Konstantin Usevich
% https://github.com/cprevost4/HSR_Tucker
% Contact: clemence.prevost@univ-lorraine.fr
addpath ../utils
addpath ../methods
addpath ../metrics


%% Indian Pines (removed regions of water absorption):
SRI = cell2mat(struct2cell(load('Indian_pines.mat')));
SRI(:,:,[104:108 150:163 220]) = []; %Regions of water absorption
SRI(1,:,:) = []; SRI(:,1,:) = [];
Pm = spectral_deg(SRI,"LANDSAT");
MSI = tmprod(SRI,Pm,3);
d1 = 4; d2 = 4; q = 9;
[P1,P2] = spatial_deg(SRI, q, d1, d2);
HSI = tmprod(tmprod(SRI,P1,1),P2,2);

% 3. run metrics
methods = {'STEREO' 'stereo3' '50'; ...
           'STEREO' 'stereo3' '100'; ...
           'SCOTT' 'scott' '[40,40,6]'; ...
           'SCOTT' 'scott' '[30,30,16]'; ...
           'SCOTT' 'scott' '[24,24,25]';...
           'B-SCOTT' 'bscott_b1_adaptor' '[40,40,6]'; ...                 
           'B-SCOTT' 'bscott_b1_adaptor' '[60,60,6]'; ...                 
           'B-SCOTT' 'bscott_b1_adaptor' '[100,100,6]'};      
DegMat = struct('Pm', Pm, 'P1', P1, 'P2', P2);         
res = compare_methods(SRI, HSI, MSI, DegMat, [d1 d2], methods);           

%% Salinas-A scene
% 1. load data
SRI = cell2mat(struct2cell(load('Salinas.mat')));
SRI = crop(SRI,[80,84,size(SRI,3)]);
SRI(:,:,[108:112 154:167 224]) = []; %Regions of water absorption (Salinas)
% 2. degradation
Pm = spectral_deg(SRI,"LANDSAT");
d1 = 4; d2 = 4; q = 9;
[P1,P2] = spatial_deg(SRI, q, d1, d2);
MSI = tmprod(SRI,Pm,3); HSI = tmprod(tmprod(SRI,P1,1),P2,2);

% 3. run metrics
methods = {'STEREO' 'stereo3' '50'; ...
           'STEREO' 'stereo3' '100'; ...
           'SCOTT' 'scott' '[40,40,6]'; ...
           'SCOTT' 'scott' '[14,14,15]'; ...
           'SCOTT' 'scott' '[10,15,25]';...
           'SCOTT' 'scott' '[30,30,6]';...                
           'SCOTT' 'scott' '[58,58,6]'};      
DegMat = struct('Pm', Pm, 'P1', P1, 'P2', P2);         
res = compare_methods(SRI, HSI, MSI, DegMat, [d1 d2], methods); 

%% Pavia University
% 1. load data
SRI = cell2mat(struct2cell(load('PaviaU.mat')));
SRI(1:2,:,:) = []; SRI(:,1:4,:) = [];
% 2. degradation
Pm = spectral_deg(SRI,"Quickbird");
MSI = tmprod(SRI,Pm,3);
d1 = 4; d2 = 4; q = 9;
[P1,P2] = spatial_deg(SRI, q, d1, d2);
HSI = tmprod(tmprod(SRI,P1,1),P2,2);

% 3. run metrics
methods2 = {'SCUBA'  'scuba_adaptor' '[120,3]' [4,4]; ...
           'B-SCOTT' 'bscott_adaptor' '[100,100,4]' [1,1]; ...                 
           'B-SCOTT' 'bscott_adaptor' '[60,60,3]' [4,4]; ...                 
           'B-SCOTT' 'bscott_adaptor' '[120,60,4]' [4,4]};      
res = compare_methods(SRI, HSI, MSI, struct('Pm', Pm), [d1 d2], methods2);           

%% Wood
% 1. load data
SRI = cell2mat(struct2cell(load('SRI1.mat')));
SRI = permute(SRI,[2,3,1]);
SRI = SRI(:,3:end,:);
HSI = cell2mat(struct2cell(load('HSID2.mat'))); %Hyperspectral image, d=2
HSI = permute(HSI,[2,3,1]);
MSI = cell2mat(struct2cell(load('MSID8.mat'))); %Multispectral image, d=2
MSI = permute(MSI,[2,3,1]); w_msi = cell2mat(struct2cell(load('W_MSID8.mat')));
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

% 3. run metrics
methods = {'STEREO' 'stereo3' '20'; ...
           'STEREO' 'stereo3' '30'; ...
           'STEREO' 'stereo3' '50'; ...
           'STEREO' 'stereo3' '100'; ...
           'SCOTT' 'scott' '[50,50,8]';...
           'SCOTT' 'scott' '[30,30,15]';};      
DegMat = struct('Pm', Pm, 'P1', P1, 'P2', P2);         
res = compare_methods(SRI, HSI, MSI, DegMat, [d1 d2], methods); 

%% Indian Pines - pansharpening
% 1. load data
SRI = cell2mat(struct2cell(load('Indian_pines.mat')));
SRI(:,:,[104:108 150:163 220]) = []; %Regions of water absorption (Indian pines)
SRI(1,:,:) = []; SRI(:,1,:) = [];
% 2. degradation
d1 = 4; d2 = 4; q = 9;
[P1,P2] = spatial_deg(SRI, q, d1, d2);
HSI = tmprod(tmprod(SRI,P1,1),P2,2);
Km = 1; Pm = zeros(Km,size(SRI,3));
spec_range = linspace(400,2500,size(SRI,3));
Pm(1,:) = 1/length(spec_range);
MSI = tmprod(SRI,Pm,3); 

% 3. run metrics
methods = {'SCOTT' 'scott' '[24,24,25]';...
           'SCOTT' 'scott' '[30,30,16]';...
           'SCOTT' 'scott' '[35,35,6]';};      
DegMat = struct('Pm', Pm, 'P1', P1, 'P2', P2);         
res = compare_methods(SRI, HSI, MSI, DegMat, [d1 d2], methods); 

%% Pavia University - blind (v.2)
% 1. load data
SRI = cell2mat(struct2cell(load('PaviaU.mat')));
SRI(1:2,:,:) = []; SRI(:,1:4,:) = [];
% 2. degradation
Pm = spectral_deg(SRI,"Quickbird");
MSI = tmprod(SRI,Pm,3);
d1 = 4; d2 = 4; q = 7;
[P1,P2] = spatial_deg(SRI, q, d1, d2);
HSI = tmprod(tmprod(SRI,P1,1),P2,2);
% 3. add noise
for k=1:size(HSI,3)
    HSI(:,:,k) = awgn(HSI(:,:,k),15);
end
for k=1:size(MSI,3)
    MSI(:,:,k) = awgn(MSI(:,:,k),25);
end

% 4. run metrics
methods2 = {'SCUBA'  'scuba_adaptor' '[120,3]' [4,4]; ...
           'B-SCOTT' 'bscott_adaptor' '[size(HSI,1),size(HSI,2),3]' [4,4]; ...                 
           'SCUBA'  'scuba_adaptor' '[300,9]' [1,1]; ... 
           'B-SCOTT' 'bscott_adaptor' '[size(HSI,1),size(HSI,2),size(MSI,3)]' [1,1]; ...
           'SCUBA'  'scuba_adaptor' '[120,3]' [8,8]; ... 
           'B-SCOTT' 'bscott_adaptor' '[size(MSI,1)/8,size(MSI,2)/8,3]' [8,8];};      
res = compare_methods(SRI, HSI, MSI, struct('Pm', Pm), [d1 d2], methods2);

%% Cuprite - blind
% 1. load data
load Cuprite.mat
SRI = double(X); clear X %Convert groundtruth to double
% 2. degradation
d1 = 4; d2 = 4; %Spatial downsampling ratios
[P1,P2] = spatial_deg(SRI,7,d1,d2); Pm = spectral_deg(SRI,"LANDSAT");
HSI = tmprod(SRI,{P1,P2},[1,2]); MSI = tmprod(SRI,Pm,3);
% 3. add noise
for k=1:size(HSI,3)
    HSI(:,:,k) = awgn(HSI(:,:,k),15);
end
for k=1:size(MSI,3)
    MSI(:,:,k) = awgn(MSI(:,:,k),25);
end

% 4. run metrics
methods2 = {'SCUBA'  'scuba_adaptor' '[45,3]' [4,4]; ...
           'B-SCOTT' 'bscott_adaptor' '[45,45,3]' [4,4]; ...                 
           'SCUBA'  'scuba_adaptor' '[150,10]' [1,1]; ... 
           'B-SCOTT' 'bscott_adaptor' '[150,150,6]' [1,1]; ...
           'SCUBA'  'scuba_adaptor' '[45,3]' [8,8]; ... 
           'B-SCOTT' 'bscott_adaptor' '[45,45,3]' [8,8];};      
res = compare_methods(SRI, HSI, MSI, struct('Pm', Pm), [d1 d2], methods2);



