% EXPERIMENT 3: COMPARISON OF METRICS WITH STEREO %
% Copyright (c) 2018 Clemence Prevost, Konstantin Usevich, Pierre Comon, David Brie
% https://github.com/cprevost4/HSR_Tucker
% Contact: clemence.prevost@univ-lorraine.fr

addpath ../utils
addpath ../metrics
addpath ../methods

% %Salinas
% SRI = cell2mat(struct2cell(load('Salinas.mat')));
% SRI = crop(SRI,[80,84,size(SRI,3)]);
% SRI(:,:,[108:112 154:167 224]) = []; %Regions of water absorption (Salinas)

% Pavia dataset
SRI = cell2mat(struct2cell(load('PaviaU.mat')));
SRI(1:2,:,:) = []; SRI(:,1:4,:) = [];
%SRI = SRI(:,1:320,:);
%SRI(1,:,:) = []; SRI(:,1,:) = [];
% 
Pm = spectral_deg(SRI,"Quickbird");
MSI = tmprod(SRI,Pm,3);
 
d1 = 4; d2 = 4; q = 9;
[P1,P2] = spatial_deg(SRI, q, d1, d2);
HSI = tmprod(tmprod(SRI,P1,1),P2,2);


% R1 = [40,40,6]; R2 = [14,14,15]; R3 = [10,15,25];
% R4 = [30,30,6]; R5 = [58,58,6];

% F = 20;
% [A0, B0, ~,~, C0] = tenRec(MSI, HSI, F, P1,P2,Pm);
% [~,~,~,Sa] = stereo(1, F, A0,B0,C0, HSI, MSI, P1,P2,Pm, 10);
% erra = {nmse(SRI,Sa), SAM(SRI,Sa), ergas(SRI,Sa), r_snr(SRI,Sa), cc(SRI,Sa)};
% F = 30;
% [A0, B0, ~,~, C0] = tenRec(MSI, HSI, F, P1,P2,Pm);
% [~,~,~,Sb] = stereo(1, F, A0,B0,C0, HSI, MSI, P1,P2,Pm, 10);
% errb = {nmse(SRI,Sb), SAM(SRI,Sb), ergas(SRI,Sb), r_snr(SRI,Sb), cc(SRI,Sb)};
% F = 50;
% [A0, B0, ~,~, C0] = tenRec(MSI, HSI, F, P1,P2,Pm);
% [~,~,~,S0] = stereo(1, F, A0,B0,C0, SRI, HSI, MSI, P1,P2,Pm, 10);
% err0 = {nmse(SRI,S0), sam(SRI,S0), ergas(SRI,S0,1/d1), r_snr(SRI,S0), cc(SRI,S0)};
% F = 100;
% [A0, B0, ~,~, C0] = tenRec(MSI, HSI, F, P1,P2,Pm);
% [~,~,~,S1] = stereo(1, F, A0,B0,C0, SRI, HSI, MSI, P1,P2,Pm, 10);
% err1 = {nmse(SRI,S1), sam(SRI,S1), ergas(SRI,S1,1/d1), r_snr(SRI,S1), cc(SRI,S1)};
opts.Nblocks = [4,4];
[S1,~] = scuba(MSI,HSI,120,3,Pm,opts);

[S2,~] = scott(HSI, MSI, P1, P2, Pm,[60,60,4]);

[S3,~] = bscott(MSI,HSI,[100,100,4],Pm);

[S4,~] = bscott(MSI,HSI,[60,60,3],Pm,opts);

[S5,~] = bscott(MSI,HSI,[120,60,4],Pm,opts);

figure(1)
subplot(2,2,1)
imagesc(SRI(:,:,44))
title('Groundtruth SRI')
subplot(2,2,2)
imagesc(S1(:,:,44))
title('SCUBA')
subplot(2,2,3)
imagesc(S4(:,:,44))
title('block-HOSVD [60,60,3]')
subplot(2,2,4)
imagesc(S5(:,:,44))
title('block-HOSVD [120,60,4]')




%% MAKE TABLE FROM RESULTS

%errH = load('metrics_hs_ip'); %This is obtained from Hysure and it must contain the same metrics

 T2 = [% real([erra{4}(end) erra{1}(end) erra{5}(end) erra{2}(end) erra{3}(end)]);
 %    real([errb{4}(end) errb{1}(end) errb{5}(end) errb{2}(end) errb{3}(end)]);
     err1{6} err1{2} err1{3} err1{4} err1{5} err1{7};
     %err2{6} err2{2} err2{3} err2{4} err2{5} err2{7};
    err3{6} err3{2} err3{3} err3{4} err3{5} err3{7};
    err4{6} err4{2} err4{3} err4{4} err4{5} err4{7};
    err5{6} err5{2} err5{3} err5{4} err5{5} err5{7};
%     err5{6} err5{2} err5{3} err5{4} err5{5};
%     err6{6} err6{2} err6{3} err6{4} err6{5};
%    cell2mat(errH.err_hysure_ip);
%     cell2mat(errH.err);
%     cell2mat(errC.err)]; 
]
T2 = table(T2);
writetable(T2, 'exp3_table2_pavia');

