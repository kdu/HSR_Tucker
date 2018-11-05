% EXPERIMENT 3: COMPARISON OF METRICS WITH STEREO %
% Copyright (c) 2018 ClÃ©mence PrÃ©vost, Konstantin Usevich, Pierre Comon, David Brie
% https://github.com/cprevost4/HSR_Tucker

addpath ../utils
addpath ../metrics
addpath ../methods

% Indian Pines
SRI = cell2mat(struct2cell(load('Indian_pines.mat')));
SRI(:,:,[104:108 150:163 220]) = []; %Regions of water absorption (Indian pines)
SRI(1,:,:) = []; SRI(:,1,:) = [];
% 
Pm = spectral_deg(SRI);
MSI = tmprod(SRI,Pm,3);
 
d1 = 4; d2 = 4; q = 9;
[P1,P2] = spatial_deg(SRI, q, d1, d2);
HSI = tmprod(tmprod(SRI,P1,1),P2,2);


[S7,~, err7] = run_hosvd_blind(SRI,MSI,HSI,[40,40,6],P1,P2,Pm);
[S8,~, err8] = run_hosvd_blind(SRI,MSI,HSI,[60,60,6],P1,P2,Pm);
[S9,~, err9] = run_hosvd_blind(SRI,MSI,HSI,[100,100,6],P1,P2,Pm);

%% MAKE TABLE FROM RESULTS

 T2 = [ err7{6} err7{2} err7{3} err7{4} err7{5} err7{7}; ...
    err8{6} err8{2} err8{3} err8{4} err8{5} err8{7};...
    err9{6} err9{2} err9{3} err9{4} err9{5} err9{7};...
];
T2 = table(T2);
writetable(T2, 'exp3_table2_ip_ku')
