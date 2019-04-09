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

%% 
F = 50;

tic; [SRI_hat1, ~] = stereo3( HSI, MSI, P1,P2,Pm, F); toc
tic; [SRI_hat2, ~] = stereo4( HSI, MSI, P1,P2,Pm, F); toc


    
    

