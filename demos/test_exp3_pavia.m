% EXPERIMENT 3: COMPARISON OF METRICS WITH STEREO %
% Copyright (c) 2018 Clemence Prevost, Konstantin Usevich, Pierre Comon, David Brie
% https://github.com/cprevost4/HSR_Tucker
% Contact: clemence.prevost@univ-lorraine.fr

addpath ../utils
addpath ../metrics
addpath ../methods

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


%%

methods_nb = {'stereo','scott'};
names_nb = {'STEREO F=20','STEREO F=30','STEREO F=50','STEREO F=100',...
    'SCOTT = [60,60,4]','NaN','NaN','NaN'};
ranks_nb = {'20','[60,60,4]';
            '30','0';
            '50','0';
            '100','0';
            };       
c = 0;
        
 for i=1:size(methods_nb,2)
     for j=1:size(ranks_nb,1)
         c = c+1;
         if ranks_nb{j,i}~="0"
            tic,
            eval(sprintf('[Y_hat,info] = %s(HSI, MSI, P1, P2, Pm, %s);', methods_nb{i},ranks_nb{j,i}));
            time = toc;
            eval('err = [compute_metrics(SRI,Y_hat,d1,d2), time];')
            res{1,c} = cell2mat(err)';
         else
            res{1,c} = NaN*ones(1,5)';
         end
     end
 end       
 save_text_table('res_pavia_nonblind.txt', names_nb, cell2mat(res));
 
 methods_b = {'scuba','bscott'};
 names_b = {'SCUBA [120,3]','NaN','NaN', 'BSCOTT R=[100,100,4]',...
     'BSCOTT R=[60,60,3]', 'BSCOTT R=[120,60,4]'};
 ranks_b = {'[120,3]','[100,100,4]';
            '0','[60,60,3]';
            '0','[120,60,4]';
            };
opts.Nblocks = [4,4]; c = 0;

 for i=1:size(methods_b,2)
     for j=1:size(ranks_b,1)
         c = c+1;
         if ranks_b{j,i} == "[100,100,4]"
             opts.Nblocks = [1,1];
         end
         if ranks_b{j,i}~="0"
            tic,
            eval(sprintf('[Y_hat,info] = %s(MSI, HSI, Pm, %s, opts);', methods_b{i},ranks_b{j,i}));
            time = toc;
            eval('err = [compute_metrics(SRI,Y_hat,d1,d2), time];')
            res{1,c} = cell2mat(err)';
         else
            res{1,c} = NaN*ones(1,5)';
         end
     end
 end  
 save_text_table('res_pavia_blind.txt', names_b, cell2mat(res));


%%%%%%%%%%%%
% SHOW RECONSTRUCTION 
%%%%%%%%%%%%%%%%%%%%%

% figure(1)
% subplot(2,2,1)
% imagesc(SRI(:,:,44))
% title('Groundtruth SRI')
% subplot(2,2,2)
% imagesc(S1(:,:,44))
% title('SCUBA')
% subplot(2,2,3)
% imagesc(S4(:,:,44))
% title('block-HOSVD [60,60,3]')
% subplot(2,2,4)
% imagesc(S5(:,:,44))
% title('block-HOSVD [120,60,4]')

